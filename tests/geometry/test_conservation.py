"""F6 multi-step conservation and BL inventory identity tests.

Plan 4.3 / design §10.2 instrument for the rise-dt investigation (2026-07-25).

Global conservation
-------------------
For a vent transient the tank+vent inventory is

    m_liq(t) + m_ullage(t) + m_BL_vap(t) + m_vented(t)

where ``m_BL_vap = Vapor in BL`` (V_BL · ρ_v) is tank inventory not yet handed
to the free ullage, and ``m_vented`` is the cumulative ``Vent Rate · Δt``.
The three-term form without BL vapour misses that hold-up; both are reported.

Volume identity (custom geometry)
---------------------------------
    V(h) = m_liq / ρ_ℓ + V_BL

and the two ullage estimates

    Ullage Mass  (ODE)    vs    (V_tank − V_occ) · ρ_v    (volumetric)

Measured-vs-set tolerances (stock solver @ 74a0d7d, pinned env; just above noise):

| Case                         | Metric                         | Measured max | Set tol   |
|------------------------------|--------------------------------|-------------:|----------:|
| H2 custom cylinder, 0.25 s, 8 s | |Δm_3| (liq+ull+vent)      | 2.42e-4     | 5.0e-4    |
| H2 custom cylinder, 0.25 s, 8 s | |Δm_4| (+BL vap)            | 1.70e-6     | 5.0e-6    |
| H2 custom cylinder               | |V(h)−(m/ρ+V_BL)|          | 1.65e-6     | 5.0e-6    |
| H2 custom cylinder               | rel ullage ODE vs volumetric | 8.3e-7      | 5.0e-6    |
| LOX G3, Δt=0.02, 12 s            | |Δm_4| (+BL vap)            | 2.06e-2     | 3.0e-2    |
| LOX G3, Δt=0.02, 12 s            | |V(h)−(m/ρ+V_BL)|          | 4.0e-6      | 1.0e-5    |
| LOX G3, Δt=0.02, 12 s            | rel ullage ODE vs volumetric | 1.84e-2     | 2.5e-2    |

BL inventory identity (rate-scaled)
----------------------------------
Per step, reconstructing S·Δt and E·Δt from CSV columns the same way as
``evidence/10x_i1_residual_check.txt`` (i1-mechanism):

    S·Δt = AK2 · LiqMass_{prev} · Δt / ρ_ℓ(T)
    E·Δt = (BL Vap Out) / ρ_v(T)
    f_VBL = ΔV_BL − (S − E)·Δt

Rate-scaled acceptance (correct nondimensionalisation of the residual gate):

    |f_VBL| ≤ 0.001 · (|S| + |E|) · Δt + ABS_FLOOR

where the coefficient 0.001 matches the recommended solver gate and ABS_FLOOR
absorbs post-hoc CSV reconstruction lag (measured max |f_VBL| ≈ 8.8e-6 at
true quasi-steady Δt=0.02, t∈[25, 40] s → floor 2e-5).

- Δt = 0.02, QS window: GREEN (true E≈S root).
- Δt = 0.005: GREEN after the custom-mode rate-scaled residual gate fix
  (2026-07-25 s5a-gate-fix). Was a strict-xfail tripwire under the absolute
  0.001·V_BL gate, which accepted a false quasi-steady BL (E ≠ S, V_BL frozen).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from core import liqlev_simulation
from liqlev.geometry.fixtures import cylinder_kernel
from liqlev.geometry.jit import eval_ppoly, invert_monotone_volume
from liqlev.geometry.package import load_geometry_package
from thermo_utils import build_property_table
from validation.lox_vent_cases import (
    G3,
    TANK_TOTAL_VOLUME_FT3,
    run_lox_vent_case,
)


ROOT = Path(__file__).resolve().parents[2]
GEOMETRY_DIAMETER_FT = 4.0
TANK_HEIGHT_FT = 8.0
TINIT_R = 38.3

# --- Measured-vs-set tolerances (see module docstring) ---
CYL_MASS3_TOL = 5.0e-4
CYL_MASS4_TOL = 5.0e-6
CYL_VOLUME_TOL = 5.0e-6
CYL_ULLAGE_REL_TOL = 5.0e-6

LOX_MASS4_TOL = 3.0e-2
LOX_VOLUME_TOL = 1.0e-5
LOX_ULLAGE_REL_TOL = 2.5e-2

# Rate-scaled BL residual: 0.001 · (|S|+|E|)·Δt + absolute recon floor.
RATE_SCALED_COEFF = 0.001
BL_RECON_ABS_FLOOR = 2.0e-5
BL_SEED_SKIP_STEPS = 10


def _hydrogen_density(t_rankine: float) -> tuple[float, float]:
    liquid = (
        0.1709
        + 0.7454 * t_rankine
        - 0.04421 * t_rankine**2
        + 0.001248 * t_rankine**3
        - 1.738e-5 * t_rankine**4
        + 9.424e-8 * t_rankine**5
    )
    vapor = (
        -0.2511
        + 0.04294 * t_rankine
        - 0.00286 * t_rankine**2
        + 9.159e-5 * t_rankine**3
        - 1.422e-6 * t_rankine**4
        + 1.001e-8 * t_rankine**5
    )
    return liquid, vapor


def _custom_cylinder_multistep_inputs(
    *,
    fill: float = 0.5,
    delta_s: float = 0.25,
    duration_s: float = 8.0,
    vent_rate_lbm_s: float = 5.0e-5,
    gravity_g: float = 0.001,
) -> tuple[dict[str, object], object]:
    """Analytic custom-cylinder vent transient (multi-step, plateau-regime Δt)."""

    kernel = cylinder_kernel(
        GEOMETRY_DIAMETER_FT,
        TANK_HEIGHT_FT,
        node_count=1025,
    )
    liquid_density, _ = _hydrogen_density(TINIT_R)
    liquid_mass = liquid_density * fill * kernel.total_volume_ft3
    h0 = float(
        invert_monotone_volume(
            fill * kernel.total_volume_ft3,
            kernel.height_ft,
            kernel.volume_ft3,
            kernel.volume_coefficients,
        )
    )
    time = np.array([0.0, duration_s], dtype=np.float64)
    ones = np.ones(2, dtype=np.float64)
    inputs: dict[str, object] = {
        "Liquid": "Hydrogen",
        "Units": "British",
        "Delta": delta_s,
        "Dtank": 1.25,
        "Htzero": h0,
        "Volt": kernel.total_volume_ft3,
        "Xmlzro": liquid_mass,
        "Pinit": 19.5,
        "Pfinal": 10.0,
        "Tinit": TINIT_R,
        "Neps": 0,
        "Tvmdot": time,
        "Xvmdot": np.full(2, vent_rate_lbm_s, dtype=np.float64),
        "Teps": None,
        "Xeps": None,
        "Tspal": time,
        "Xspacl": ones,
        "Tspav": time,
        "Xspacv": ones,
        "Tggo": time,
        "Xggo": np.full(2, gravity_g * 32.174, dtype=np.float64),
        "GeometryMode": 1,
        "FillFraction": fill,
        "GeomHeight": kernel.height_ft,
        "GeomVolume": kernel.volume_ft3,
        "GeomVolumeCoefficients": kernel.volume_coefficients,
        "GeomAreaSamples": kernel.section_area_ft2,
        "GeomPerimeter": kernel.perimeter_ft,
        "GeomSidewallArea": kernel.sidewall_area_ft2,
        "GeomSidewallCoefficients": kernel.sidewall_coefficients,
    }
    return inputs, kernel


def _cumulative_vented(dataframe, delta_s: float) -> np.ndarray:
    return np.cumsum(dataframe["Vent Rate"].to_numpy(dtype=float) * delta_s)


def _mass_inventories(dataframe, delta_s: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (m3, m4) series: without / with BL vapour hold-up."""

    liq = dataframe["Liq Mass"].to_numpy(dtype=float)
    ull = dataframe["Ullage Mass"].to_numpy(dtype=float)
    bl_vap = dataframe["Vapor in BL"].to_numpy(dtype=float)
    vented = _cumulative_vented(dataframe, delta_s)
    m3 = liq + ull + vented
    m4 = liq + ull + bl_vap + vented
    return m3, m4


def _max_abs_drift(series: np.ndarray) -> float:
    return float(np.max(np.abs(series - series[0])))


def test_analytic_cylinder_multistep_global_mass_and_volume_conservation() -> None:
    """Custom-cylinder multi-step vent: mass + volume + ullage identities.

    Tolerances measured at Δt=0.25 s, duration 8 s (see module docstring).
    """

    delta_s = 0.25
    duration_s = 8.0
    inputs, kernel = _custom_cylinder_multistep_inputs(
        delta_s=delta_s, duration_s=duration_s
    )
    dataframe = liqlev_simulation(inputs, verbose=False)

    assert len(dataframe) >= 8
    assert int(dataframe["Conv Failed"].sum()) == 0
    assert np.isfinite(dataframe.to_numpy(dtype=float)).all()

    m3, m4 = _mass_inventories(dataframe, delta_s)
    m3_drift = _max_abs_drift(m3)
    m4_drift = _max_abs_drift(m4)
    # Measured: m3≈2.42e-4, m4≈1.70e-6 → set 5e-4 / 5e-6
    assert m3_drift <= CYL_MASS3_TOL
    assert m4_drift <= CYL_MASS4_TOL

    temps = dataframe["Temp"].to_numpy(dtype=float)
    rhol = np.array([_hydrogen_density(t)[0] for t in temps])
    rhov = np.array([_hydrogen_density(t)[1] for t in temps])
    liq = dataframe["Liq Mass"].to_numpy(dtype=float)
    vbl = dataframe["VBL vol"].to_numpy(dtype=float)
    heights = dataframe["Height"].to_numpy(dtype=float)
    v_occ = liq / rhol + vbl
    v_from_h = np.array(
        [
            eval_ppoly(h, kernel.height_ft, kernel.volume_coefficients)
            for h in heights
        ]
    )
    volume_err = float(np.max(np.abs(v_from_h - v_occ)))
    # Measured ≈ 1.65e-6 → set 5e-6
    assert volume_err <= CYL_VOLUME_TOL

    # Two ullage estimates: ODE mass vs volumetric (V_tank - V_occ)*rhov
    ull = dataframe["Ullage Mass"].to_numpy(dtype=float)
    ull_vol = (float(kernel.total_volume_ft3) - v_occ) * rhov
    rel_ull = np.abs(ull - ull_vol) / np.maximum(np.abs(ull_vol), 1e-30)
    # Measured ≈ 8.3e-7 → set 5e-6
    assert float(np.max(rel_ull)) <= CYL_ULLAGE_REL_TOL


def test_lox_global_mass_and_volume_conservation_at_plateau_dt() -> None:
    """LOX NASA geometry global conservation at plateau Δt=0.02 s (short run).

    Duration 12 s keeps the suite fast while covering BL growth. Tolerances
    measured on stock solver (see module docstring).
    """

    delta_s = 0.02
    duration_s = 12.0
    dataframe, summary = run_lox_vent_case(
        G3, timestep_s=delta_s, duration_s=duration_s
    )

    assert summary.rows >= 100
    assert summary.conv_failed_total == 0
    assert summary.finite
    assert summary.timestep_s == pytest.approx(delta_s)

    m3, m4 = _mass_inventories(dataframe, delta_s)
    m4_drift = _max_abs_drift(m4)
    # Complete inventory (incl. BL vapour). Measured ≈ 0.0206 → set 0.03
    assert m4_drift <= LOX_MASS4_TOL

    geometry = load_geometry_package(
        ROOT / "geometry" / "tables" / "nhq01-m21a-0201_LIQLEV_GEOMETRY.npz"
    )
    temps = dataframe["Temp"].to_numpy(dtype=float)
    table = build_property_table("Oxygen", 30.0, 45.0)
    t_grid, rhol_grid, rhov_grid = table[0], table[1], table[2]
    rhol = np.interp(temps, t_grid, rhol_grid)
    rhov = np.interp(temps, t_grid, rhov_grid)
    liq = dataframe["Liq Mass"].to_numpy(dtype=float)
    vbl = dataframe["VBL vol"].to_numpy(dtype=float)
    heights = dataframe["Height"].to_numpy(dtype=float)
    v_occ = liq / rhol + vbl
    v_from_h = np.array(
        [
            eval_ppoly(h, geometry.height_ft, geometry.volume_coefficients)
            for h in heights
        ]
    )
    volume_err = float(np.max(np.abs(v_from_h - v_occ)))
    # Measured ≈ 4e-6 → set 1e-5
    assert volume_err <= LOX_VOLUME_TOL

    ull = dataframe["Ullage Mass"].to_numpy(dtype=float)
    ull_vol = (TANK_TOTAL_VOLUME_FT3 - v_occ) * rhov
    rel_ull = np.abs(ull - ull_vol) / np.maximum(np.abs(ull_vol), 1e-30)
    # Measured ≈ 0.0184 → set 0.025 (known ullage-closure gap; not relaxed further)
    assert float(np.max(rel_ull)) <= LOX_ULLAGE_REL_TOL


def _reconstruct_bl_residual(
    dataframe, delta_s: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Mirror i1 residual reconstruction (10x_i1_residual_check.txt).

    Returns ``(t_end_of_step, f_VBL, scale_SE_dt, rate_imbalance)`` for steps
    1..N-1, using Liq Mass at previous row as xml1.
    """

    table = build_property_table("Oxygen", 30.0, 45.0)
    t_grid, rhol_grid, rhov_grid = table[0], table[1], table[2]
    temps = dataframe["Temp"].to_numpy(dtype=float)
    rhol = np.interp(temps, t_grid, rhol_grid)
    rhov = np.interp(temps, t_grid, rhov_grid)
    ak2 = dataframe["AK2"].to_numpy(dtype=float)
    liq = dataframe["Liq Mass"].to_numpy(dtype=float)
    vbl = dataframe["VBL vol"].to_numpy(dtype=float)
    bl_out = dataframe["BL Vap Out"].to_numpy(dtype=float)
    t = dataframe["Time"].to_numpy(dtype=float)

    # S*dt from (AK2_i, LiqMass_{i-1}, rhol_i); E*dt from BL Vap Out_i / rhov_i
    source_dt = ak2[1:] * liq[:-1] * delta_s / rhol[1:]
    exit_dt = bl_out[1:] / rhov[1:]
    d_vbl = np.diff(vbl)
    fvbl = d_vbl - source_dt + exit_dt
    scale = np.abs(source_dt) + np.abs(exit_dt)
    imbalance = np.abs(source_dt - exit_dt) / np.maximum(scale, 1e-30)
    return t[1:], fvbl, scale, imbalance


def _assert_rate_scaled_bl_identity(
    dataframe,
    delta_s: float,
    *,
    t_min: float,
    t_max: float,
    seed_skip: int = BL_SEED_SKIP_STEPS,
) -> None:
    t, fvbl, scale, imbalance = _reconstruct_bl_residual(dataframe, delta_s)
    mask = (t >= t_min) & (t <= t_max)
    # also skip early seed steps by index
    idx = np.arange(len(t))
    mask &= idx >= seed_skip
    assert np.any(mask), f"empty BL window t∈[{t_min},{t_max}] seed_skip={seed_skip}"

    tol = RATE_SCALED_COEFF * scale[mask] + BL_RECON_ABS_FLOOR
    residual = np.abs(fvbl[mask])
    worst = float(np.max(residual - tol))
    assert np.all(residual <= tol), (
        f"rate-scaled BL identity violated: max(residual-tol)={worst:.3e}, "
        f"max|f_VBL|={float(np.max(residual)):.3e}, "
        f"max rate-imbalance |S-E|/(|S|+|E|)={float(np.max(imbalance[mask])):.4f}, "
        f"window t∈[{t_min},{t_max}] Δt={delta_s}"
    )


def test_lox_bl_inventory_identity_rate_scaled_at_plateau_dt() -> None:
    """BL inventory identity at Δt=0.02 (true quasi-steady regime) — GREEN.

    Bounded window t∈[25, 40] s after seed skip: E≈S, residual within
    0.001·(|S|+|E|)·Δt + recon floor (measured max |f_VBL|≈8.8e-6).
    """

    delta_s = 0.02
    dataframe, summary = run_lox_vent_case(
        G3, timestep_s=delta_s, duration_s=40.0
    )
    assert summary.conv_failed_total == 0
    assert summary.finite
    # Sanity: plateau-regime inventory is not the fine-dt freeze
    assert summary.final_vbl_vol_ft3 > 0.25
    assert summary.max_ak3 > 0.08

    _assert_rate_scaled_bl_identity(
        dataframe, delta_s, t_min=25.0, t_max=40.0
    )


def test_lox_bl_inventory_identity_rate_scaled_fine_dt_tripwire() -> None:
    """Rate-scaled BL inventory identity at Δt=0.005 (post gate-fix).

    Was a strict-xfail tripwire under the absolute 0.001·V_BL gate (freeze by
    ~5 s, |S−E|/(|S|+|E|)≈0.60). With the rate-scaled custom-mode gate the
    window t∈[5, 12] s satisfies the identity and BL inventory grows.
    """

    delta_s = 0.005
    dataframe, summary = run_lox_vent_case(
        G3, timestep_s=delta_s, duration_s=12.0
    )
    assert summary.conv_failed_total == 0
    assert summary.finite
    # Post-fix sanity: fine-dt no longer freezes the BL (stock had V_BL<0.15,
    # max_ak3<0.04 with |S-E|/(|S|+|E|)≈0.60). Gate fix recovers growth.
    assert summary.final_vbl_vol_ft3 > 0.15
    assert summary.max_ak3 > 0.04

    _assert_rate_scaled_bl_identity(
        dataframe, delta_s, t_min=5.0, t_max=12.0
    )
