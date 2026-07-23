# -*- coding: utf-8 -*-
"""
core.py
Core execution loop for the LIQLEV boundary layer simulation.

This module calculates the transient liquid level rise and two-phase fluid
behavior during zero-gravity venting. It uses a one-dimensional piston-like
displacement assumption where vapor generated along the wetted walls
(boundary layer) and bulk liquid flashing causes volumetric swelling.

Performance notes:
  - Inner solver loop is Numba JIT-compiled (~10-30x faster than pure Python)
  - Pre-computed property tables eliminate per-step CoolProp calls
  - Vectorized coefficient arrays for boundary layer sums
  - Pre-allocated result arrays avoid per-step dict/list overhead
"""

import numpy as np
import pandas as pd
from numba import njit
from thermo_utils import sli, Tsat, Psat, DensitySat, Cpsat, LHoV, dPdTsat


# ═══════════════════════════════════════════════════════════════════════════════
#  NUMBA JIT-COMPILED SOLVER CORE
# ═══════════════════════════════════════════════════════════════════════════════
@njit(cache=True)
def _interp(x, xp, fp):
    """Numba-compatible linear interpolation (equivalent to np.interp)."""
    if x <= xp[0]:
        return fp[0]
    if x >= xp[-1]:
        return fp[-1]
    # Binary search
    lo, hi = 0, len(xp) - 1
    while lo < hi - 1:
        mid = (lo + hi) >> 1
        if xp[mid] <= x:
            lo = mid
        else:
            hi = mid
    t = (x - xp[lo]) / (xp[hi] - xp[lo]) if xp[hi] != xp[lo] else 0.0
    return fp[lo] + t * (fp[hi] - fp[lo])


@njit(cache=True)
def _hydrogen_props(t1):
    """Hydrogen empirical polynomials (legacy AS-203 correlations)."""
    rhol = (0.1709 + 0.7454 * t1 - 0.04421 * t1**2
            + 0.001248 * t1**3 - 1.738e-5 * t1**4 + 9.424e-8 * t1**5)
    rhov = (-0.2511 + 0.04294 * t1 - 0.00286 * t1**2
            + 9.159e-5 * t1**3 - 1.422e-6 * t1**4 + 1.001e-8 * t1**5)
    cs = 0.078 * (t1 - 34.0) + 2.12
    hfg = -2.0 * (t1 - 34.0) + 194.5
    dpdts = 2.49 - 0.22 * t1 + 0.00407 * t1**2 + 5.22e-5 * t1**3
    return rhol, rhov, cs, hfg, dpdts


@njit(cache=True)
def _solver_loop(
    delta, pfinal, tinit, htzero, xmlzro, perim, ac, volt,
    tvmdot, xvmdot,
    neps, teps, xeps,
    tspal, xspacl, tspav, xspacv,
    tggo, xggo,
    is_hydrogen,
    pt_t, pt_rhol, pt_rhov, pt_cs, pt_hfg, pt_dpdts,
    use_prop_table,
    coeff_s1, exp_s1, coeff_sm, exp_sm, coeff_vbl, exp_vbl,
    dtank, pinit,
    grav_samples_t, grav_samples_g, use_grav_samples,
):
    """JIT-compiled inner solver loop.

    Returns
    -------
    result_arr : ndarray (N, 29)
    step_count : int
    """
    # Pre-allocate
    est_steps = int((tvmdot[-1] / delta) * 1.2) + 100
    n_cols = 29
    res = np.empty((est_steps, n_cols))

    # State
    dhdt = 0.0
    vbl1 = 0.0
    hldak3 = 0.0
    p1 = pinit
    t1 = tinit
    h1 = htzero
    theta1 = 0.0
    xml1 = xmlzro
    zht1 = htzero
    xmvap1 = 0.0
    is_first = True
    cumulative_vap_kg = 0.0
    step = 0
    LBM_TO_KG = 0.453592

    while p1 > pfinal:
        # --- Thermodynamic properties ---
        if is_hydrogen:
            rhol, rhov, cs, hfg, dpdts = _hydrogen_props(t1)
        elif use_prop_table:
            rhol = _interp(t1, pt_t, pt_rhol)
            rhov = _interp(t1, pt_t, pt_rhov)
            cs = _interp(t1, pt_t, pt_cs)
            hfg = _interp(t1, pt_t, pt_hfg)
            dpdts = _interp(t1, pt_t, pt_dpdts)
        else:
            # Fallback — should not reach here in JIT path
            rhol = rhov = cs = hfg = dpdts = 1.0

        # Volumes and masses
        volliq = xml1 / rhol + vbl1 if rhol != 0 else vbl1
        volgas = volt - volliq
        xmvap3 = volgas * rhov
        if is_first:
            xmvap1 = xmvap3

        # Time advance & interpolation
        dtdps = 1.0 / dpdts if dpdts != 0 else 0.0
        theta2 = theta1 + delta
        thetav = 0.5 * (theta1 + theta2)
        vmdot = _interp(thetav, tvmdot, xvmdot)

        # Depressurization physics
        if hfg * p1 * t1 != 0:
            denom = xml1 * cs * dtdps / hfg + xmvap1 * (1.0 / p1 - dtdps / t1)
        else:
            denom = 1.0
        dpdtha = -vmdot / denom if denom != 0 else 0.0
        delp = dpdtha * delta
        p2 = p1 + delp
        t2 = t1 + dtdps * delp

        # Bulk flashing
        delme = xml1 * cs * (t2 - t1) / hfg if hfg != 0 else 0.0
        xml2 = xml1 + delme
        delmv = vmdot * delta

        # Epsilon
        if neps > 0:
            eps = _interp(thetav, teps, xeps)
        else:
            eps_denom = perim * h1 + ac
            eps = (perim * h1) / eps_denom if eps_denom != 0 else 0.0

        # Bubble spacing & gravity
        spacv = _interp(thetav, tspav, xspacv)
        spacl = _interp(thetav, tspal, xspacl)
        if use_grav_samples:
            ggo_ft_s2 = _interp(thetav, grav_samples_t, grav_samples_g)
        else:
            ggo_ft_s2 = _interp(thetav, tggo, xggo)

        # Boundary layer coefficients
        if rhol != 0:
            ak1_term = 10.8 * (1 + spacl) * (1 + spacv) * ggo_ft_s2 * (rhol - rhov) / rhol
        else:
            ak1_term = 0.0
        ak1 = 1.089 * (ak1_term ** 0.5) if ak1_term > 0 else 0.0
        ak2 = -eps * cs * rhol * dtdps * dpdtha / rhov / hfg if rhov * hfg != 0 else 0.0
        ak3 = ak2 / ak1 if ak1 != 0 else 0.0

        # ── Boundary layer volume solver (Secant + Newton-Raphson) ──
        nconv = 0
        ak4 = 0.0
        fvbl4 = 0.0
        savak3 = 0.0
        svfvbl = 0.0
        solver_active = True

        while solver_active and nconv < 80:
            zht2_s = zht1 + dhdt * delta
            if ak3 < 0:
                ak3 = hldak3

            # Initial guess
            arg = 0.375 * dtank * ak3 * zht2_s
            delblz = arg ** (2.0 / 3.0) if arg > 0 else 0.0

            # Newton-Raphson inner loop (vectorized sums via dot products)
            for _ in range(20):
                # sum1 = coeff_s1 . delblz^exp_s1
                s1 = 0.0
                for k in range(10):
                    s1 += coeff_s1[k] * (delblz ** exp_s1[k])
                fdelt = 8.0 * s1 / ak3 - zht2_s if ak3 != 0 else 1e30
                if abs(fdelt) <= 1e-5 * zht2_s:
                    break
                sm = 0.0
                for k in range(10):
                    sm += coeff_sm[k] * (delblz ** exp_sm[k])
                fpdelt = 4.0 * sm / ak3 if ak3 != 0 else 1e30
                delblz -= fdelt / (fpdelt if fpdelt != 0 else 1e-9)

            # BL volume
            svbl = 0.0
            for k in range(10):
                svbl += coeff_vbl[k] * (delblz ** exp_vbl[k])
            vbl2 = svbl * np.pi / ak3 if ak3 != 0 else 0.0

            # Error function
            if rhol != 0:
                fvbl = (vbl2 - ak2 * xml1 * delta / rhol
                        + 2.1 * ak1 * dtank * (delblz ** 1.5) * delta - vbl1)
            else:
                fvbl = 1e30

            # Convergence check
            if abs(fvbl) <= 0.001 * vbl2:
                solver_active = False
            else:
                if nconv <= 1:
                    if nconv == 0:
                        nconv = 1
                    savak3 = ak3
                    svfvbl = fvbl
                    if fvbl > 0:
                        ak3 *= 0.1
                    else:
                        ak3 *= 2.0
                        nconv = 0
                else:
                    if not ((fvbl > 0 and svfvbl < 0) or (fvbl < 0 and svfvbl > 0)):
                        savak3 = ak4
                        svfvbl = fvbl4
                    ak4 = ak3
                    fvbl4 = fvbl
                    fvbl_diff = fvbl - svfvbl
                    if fvbl_diff != 0:
                        secant = ak3 - fvbl * (ak3 - savak3) / fvbl_diff
                    else:
                        secant = ak3
                    ak3 = 0.5 * (secant + savak3)
                    nconv += 1

        # ── Post-solver kinematics ──
        if rhol * ac * delta != 0:
            dhdt = ((vbl2 - vbl1) + (xml2 - xml1) / rhol) / ac / delta
        else:
            dhdt = 0.0

        zht2 = zht1 + dhdt * delta
        h2 = zht2
        delh = h2 - h1
        hratio = (h2 - htzero) / htzero if htzero != 0 else 0.0

        # Beta
        if hfg != 0:
            beta_exp = -cs * (tinit - t2) / hfg
            # Use safe exp to avoid overflow
            if beta_exp > 500:
                beta_denom = -1.0  # exp is huge, 1-exp ≈ -exp
            elif beta_exp < -500:
                beta_denom = 1.0   # exp ≈ 0
            else:
                beta_denom = 1.0 - np.exp(beta_exp)
        else:
            beta_denom = 0.0
        beta = vbl2 * rhov / xmlzro / beta_denom if xmlzro * beta_denom != 0 else 0.0

        xmvbl2 = vbl2 * rhov
        xmdtbl = 2.1 * ak1 * dtank * (delblz ** 1.5) * delta * rhov

        # Vapor generation
        mass_gen_lbm = xmdtbl - (1.0 - eps) * delme
        vap_gen_rate = (mass_gen_lbm / delta) * LBM_TO_KG

        if not is_first:
            cumulative_vap_kg += mass_gen_lbm * LBM_TO_KG

        # Ullage mass
        xmvap2 = xmvap1 - (1.0 - eps) * delme - delmv + xmdtbl

        # Superheat (approximate via slope when using prop table)
        if use_prop_table and not is_hydrogen:
            t_sat_r = t2 - delp * dtdps if dtdps != 0 else t2
        else:
            # For hydrogen or no-table, superheat ≈ 0 in this model
            t_sat_r = t2
        superheat_k = (t2 - t_sat_r) / 1.8

        # Write row
        if step >= len(res):
            # Grow array (rare — only if estimate was too small)
            new_arr = np.empty((len(res) + 500, n_cols))
            new_arr[:len(res)] = res
            res = new_arr

        res[step, 0] = theta2
        res[step, 1] = p2
        res[step, 2] = t2
        res[step, 3] = xml2
        res[step, 4] = xmvap2
        res[step, 5] = h2
        res[step, 6] = dhdt
        res[step, 7] = delh
        res[step, 8] = hratio
        res[step, 9] = dpdtha
        res[step, 10] = delp
        res[step, 11] = eps
        res[step, 12] = beta
        res[step, 13] = vbl2
        res[step, 14] = delblz
        res[step, 15] = ak3
        res[step, 16] = vmdot
        res[step, 17] = ak1
        res[step, 18] = ak2
        res[step, 19] = ak2 / ak1 if ak1 != 0 else 0.0
        res[step, 20] = xmvbl2
        res[step, 21] = xmdtbl
        res[step, 22] = xmvap3
        res[step, 23] = float(nconv)
        res[step, 24] = ggo_ft_s2 / 32.174
        res[step, 25] = superheat_k
        res[step, 26] = vap_gen_rate
        res[step, 27] = cumulative_vap_kg
        res[step, 28] = 0.0 if not solver_active else 1.0  # 0=converged, 1=failed
        step += 1

        # Advance state
        p1 = p2
        t1 = t2
        h1 = h2
        theta1 = theta2
        xml1 = xml2
        xmvap1 = xmvap2
        vbl1 = vbl2
        zht1 = zht2
        hldak3 = ak3
        is_first = False

        if theta2 >= tvmdot[-1] - delta:
            break

    return res, step


# ═══════════════════════════════════════════════════════════════════════════════
#  PUBLIC API (Python wrapper)
# ═══════════════════════════════════════════════════════════════════════════════
_COL_NAMES = [
    'Time', 'Press', 'Temp', 'Liq Mass', 'Ullage Mass',
    'Height', 'dh/dt', 'delh', 'Hratio', 'dP/dtha',
    'delp', 'eps', 'beta', 'VBL vol', 'BL thick',
    'AK3', 'Vent Rate', 'AK1', 'AK2', 'AK2/AK1',
    'Vapor in BL', 'BL Vap Out', 'Ullage from Calc', 'Conv Iterations',
    'Gravity_g', 'Superheat', 'Vap Gen Rate (kg/s)', 'Total Vap Gen (kg)',
    'Conv Failed',
]


def liqlev_simulation(inputs, verbose=True, prop_table=None, progress_cb=None):
    """Main simulation logic converted from the LIQLEV VBA subroutine.

    Parameters
    ----------
    inputs : dict       Simulation inputs (see config.py / gui.py for structure).
    verbose : bool      If False, suppress per-timestep print output.
    prop_table : tuple  Optional pre-computed property arrays from
                        thermo_utils.build_property_table().
    progress_cb : callable or None
        Called periodically with a dict of live solver stats.
    """
    # ---------------------------------------------------------
    # 1. UNPACK INPUTS & INITIALIZE GEOMETRY
    # ---------------------------------------------------------
    inputs = {k.lower(): v for k, v in inputs.items()}
    delta = inputs['delta']
    tvmdot, xvmdot = inputs['tvmdot'], inputs['xvmdot']
    neps = int(inputs['neps'])
    teps, xeps = inputs['teps'], inputs['xeps']
    tspal, xspacl = inputs['tspal'], inputs['xspacl']
    tspav, xspacv = inputs['tspav'], inputs['xspacv']
    tggo, xggo = inputs['tggo'], inputs['xggo']
    gravity_function = inputs.get('gravity_function', None)
    fluid = inputs['liquid']
    units = inputs['units']

    if units == "SI":
        dtank = inputs['dtank'] * 3.28
        htzero = inputs['htzero'] * 3.28
        volt = inputs['volt'] * (3.28 ** 3)
        xmlzro = inputs['xmlzro'] * 2.20462
        pinit = inputs['pinit'] * 0.145
        pfinal = inputs['pfinal'] * 0.145
        tinit = inputs['tinit'] * 1.8
    else:
        dtank = inputs['dtank']
        htzero = inputs['htzero']
        volt = inputs['volt']
        xmlzro = inputs['xmlzro']
        pinit = inputs['pinit']
        pfinal = inputs['pfinal']
        tinit = inputs['tinit']

    perim = np.pi * dtank
    ac = 0.7854 * (dtank ** 2)
    htank = volt / ac
    fill = htzero * ac / volt
    print(f"Initial % fill volume: {fill * 100:.2f}%")

    # Pre-compute vectorized coefficients for boundary layer sums
    _L = np.arange(1, 11, dtype=np.float64)
    coeff_s1 = (4.0 ** (_L - 1)) / (dtank ** _L) / (2.0 ** _L + 1)
    exp_s1 = _L + 0.5
    coeff_sm = (4.0 ** (_L - 1)) / (dtank ** _L)
    exp_sm = _L - 0.5
    coeff_vbl = (2.0 * _L + 1) / (_L + 1.5) / (dtank ** (_L - 1))
    exp_vbl = _L + 1.5

    # Ensure arrays are float64 contiguous for Numba
    tvmdot = np.ascontiguousarray(tvmdot, dtype=np.float64)
    xvmdot = np.ascontiguousarray(xvmdot, dtype=np.float64)
    tspal = np.ascontiguousarray(tspal, dtype=np.float64)
    xspacl = np.ascontiguousarray(xspacl, dtype=np.float64)
    tspav = np.ascontiguousarray(tspav, dtype=np.float64)
    xspacv = np.ascontiguousarray(xspacv, dtype=np.float64)
    tggo = np.ascontiguousarray(tggo, dtype=np.float64)
    xggo = np.ascontiguousarray(xggo, dtype=np.float64)

    # Handle neps=0 case — provide dummy arrays for Numba
    if neps == 0 or teps is None:
        teps_arr = np.array([0.0, 1.0], dtype=np.float64)
        xeps_arr = np.array([0.0, 0.0], dtype=np.float64)
        neps_jit = 0
    else:
        teps_arr = np.ascontiguousarray(teps, dtype=np.float64)
        xeps_arr = np.ascontiguousarray(xeps, dtype=np.float64)
        neps_jit = neps

    # Property table arrays (or empty dummies)
    is_hydrogen = (fluid == "Hydrogen")
    use_prop_table = (prop_table is not None) and not is_hydrogen
    if use_prop_table:
        pt_t, pt_rhol, pt_rhov, pt_cs, pt_hfg, pt_dpdts = prop_table
        pt_t = np.ascontiguousarray(pt_t, dtype=np.float64)
        pt_rhol = np.ascontiguousarray(pt_rhol, dtype=np.float64)
        pt_rhov = np.ascontiguousarray(pt_rhov, dtype=np.float64)
        pt_cs = np.ascontiguousarray(pt_cs, dtype=np.float64)
        pt_hfg = np.ascontiguousarray(pt_hfg, dtype=np.float64)
        pt_dpdts = np.ascontiguousarray(pt_dpdts, dtype=np.float64)
    else:
        # Dummy arrays (Numba needs typed arrays even if unused)
        _d = np.array([0.0, 1.0], dtype=np.float64)
        pt_t = pt_rhol = pt_rhov = pt_cs = pt_hfg = pt_dpdts = _d

    # Handle gravity_function by pre-sampling into an array
    use_grav_samples = False
    if gravity_function is not None:
        n_grav = int(tvmdot[-1] / delta) + 10
        grav_t = np.linspace(0, tvmdot[-1], n_grav)
        grav_g = np.array([gravity_function(t) for t in grav_t],
                          dtype=np.float64)
        use_grav_samples = True
    else:
        grav_t = np.array([0.0, 1.0], dtype=np.float64)
        grav_g = np.array([0.0, 0.0], dtype=np.float64)

    # ---------------------------------------------------------
    # 2. PROGRESS CALLBACK WRAPPER
    # ---------------------------------------------------------
    # For verbose/progress_cb, we run the JIT solver and then optionally
    # report progress. The JIT function itself doesn't do I/O.
    import time as _time
    _wall_start = _time.perf_counter()

    if verbose:
        print("  [JIT solver running...]")

    # ---------------------------------------------------------
    # 3. RUN JIT SOLVER
    # ---------------------------------------------------------
    res_arr, n_steps = _solver_loop(
        delta, pfinal, tinit, htzero, xmlzro, perim, ac, volt,
        tvmdot, xvmdot,
        neps_jit, teps_arr, xeps_arr,
        tspal, xspacl, tspav, xspacv,
        tggo, xggo,
        is_hydrogen,
        pt_t, pt_rhol, pt_rhov, pt_cs, pt_hfg, pt_dpdts,
        use_prop_table,
        coeff_s1, exp_s1, coeff_sm, exp_sm, coeff_vbl, exp_vbl,
        dtank, pinit,
        grav_t, grav_g, use_grav_samples,
    )

    elapsed = _time.perf_counter() - _wall_start
    if verbose:
        print(f"  [JIT solver complete: {n_steps} steps in {elapsed:.2f}s]")

    # Final progress callback
    if progress_cb is not None and n_steps > 0:
        progress_cb({
            "sim_time": res_arr[n_steps - 1, 0],
            "pressure": res_arr[n_steps - 1, 1],
            "height": res_arr[n_steps - 1, 5],
            "dh_h0": res_arr[n_steps - 1, 8],
            "wall_time": elapsed,
            "pct_complete": 100.0,
        })

    return pd.DataFrame(res_arr[:n_steps], columns=_COL_NAMES)
