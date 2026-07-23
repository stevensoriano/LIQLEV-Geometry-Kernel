# -*- coding: utf-8 -*-
"""
thermo_utils.py
Helper and thermodynamic functions for the LIQLEV simulation.
"""

import numpy as np
from CoolProp.CoolProp import PropsSI

def dPdTsat(fluid, press_kpa):
    """Evaluates slope of the saturated P-T curve (kPa/K)."""
    t_plus = Tsat(fluid, press_kpa + 0.1)
    t_minus = Tsat(fluid, press_kpa - 0.1)
    if t_plus == t_minus: return 0
    return 0.2 / (t_plus - t_minus)

def Cpsat(fluid, phase, press_kpa):
    """Evaluates saturation specific heat at constant pressure (kJ/kg-K)."""
    quality = 1 if phase == "vapor" else 0
    return PropsSI("C", "P", press_kpa * 1000, "Q", quality, fluid) / 1000

def DensitySat(fluid, phase, press_kpa):
    """Evaluates saturation density (kg/m^3)."""
    quality = 1 if phase == "vapor" else 0
    return PropsSI("Dmass", "P", press_kpa * 1000, "Q", quality, fluid)

def EnthalpySat(fluid, phase, press_kpa):
    """Evaluates saturation enthalpy (kJ/kg)."""
    quality = 1 if phase == "vapor" else 0
    return PropsSI("H", "P", press_kpa * 1000, "Q", quality, fluid) / 1000

def LHoV(fluid, press_kpa):
    """Evaluates latent heat of vaporization (kJ/kg)."""
    h_vap = EnthalpySat(fluid, "vapor", press_kpa)
    h_liq = EnthalpySat(fluid, "liquid", press_kpa)
    return h_vap - h_liq

def Psat(fluid, temp_sat_k):
    """Evaluates saturation pressure (kPa)."""
    return PropsSI("P", "T", temp_sat_k, "Q", 1, fluid) / 1000

def Tsat(fluid, press_kpa):
    """Evaluates saturation temperature (K)."""
    return PropsSI("T", "P", press_kpa * 1000, "Q", 1, fluid)

def sli(arg, x, y):
    """Simple Linear Interpolation using numpy."""
    return np.interp(arg, x, y)


# ── Pre-computed Property Table for Fast Simulation ──────────────────────────
def build_property_table(fluid, p_min_psia, p_max_psia, n_points=400):
    """
    Pre-compute saturation properties over the simulation pressure range.

    Instead of calling CoolProp every timestep (~10ms each), this builds
    interpolation arrays once (~1s), then np.interp runs in ~1μs per call.
    Gives ~10,000x speedup on property evaluation.

    Returns: (temps_R, rhol, rhov, cs, hfg, dpdts) — all numpy arrays
             indexed by temperature in Rankine.
    """
    PSI_TO_KPA = 6.89475729

    # Get temperature bounds (with safety margin)
    t_min_K = Tsat(fluid, p_min_psia * PSI_TO_KPA)
    t_max_K = Tsat(fluid, p_max_psia * PSI_TO_KPA)

    t_min_R = (t_min_K - 2.0) * 1.8
    t_max_R = (t_max_K + 2.0) * 1.8

    temps_R = np.linspace(t_min_R, t_max_R, n_points)
    rhol_arr = np.empty(n_points)
    rhov_arr = np.empty(n_points)
    cs_arr = np.empty(n_points)
    hfg_arr = np.empty(n_points)
    dpdts_arr = np.empty(n_points)

    for i, t_R in enumerate(temps_R):
        t_k = t_R / 1.8
        try:
            # Single Psat call (was called 4x before via dPdTsat + LHoV)
            ps_pa = PropsSI("P", "T", t_k, "Q", 1, fluid)
            ps_kpa = ps_pa / 1000.0

            # Batch property calls — direct PropsSI instead of wrapper overhead
            rhol_arr[i] = PropsSI("Dmass", "P", ps_pa, "Q", 0, fluid) * 0.0624279606
            rhov_arr[i] = PropsSI("Dmass", "P", ps_pa, "Q", 1, fluid) * 0.0624279606
            cs_arr[i]   = PropsSI("C", "P", ps_pa, "Q", 0, fluid) / 1000 * 0.2388458966
            h_liq = PropsSI("H", "P", ps_pa, "Q", 0, fluid) / 1000
            h_vap = PropsSI("H", "P", ps_pa, "Q", 1, fluid) / 1000
            hfg_arr[i]  = (h_vap - h_liq) * 0.4299226

            # dP/dT via finite difference (2 calls instead of 4)
            t_p = PropsSI("T", "P", (ps_kpa + 0.1) * 1000, "Q", 1, fluid)
            t_m = PropsSI("T", "P", (ps_kpa - 0.1) * 1000, "Q", 1, fluid)
            dpdts_arr[i] = (0.2 / (t_p - t_m) if t_p != t_m else 0) * 0.08057652094
        except Exception:
            # Near phase boundaries CoolProp can fail — use neighbor value
            if i > 0:
                rhol_arr[i] = rhol_arr[i - 1]
                rhov_arr[i] = rhov_arr[i - 1]
                cs_arr[i] = cs_arr[i - 1]
                hfg_arr[i] = hfg_arr[i - 1]
                dpdts_arr[i] = dpdts_arr[i - 1]
            else:
                rhol_arr[i] = rhov_arr[i] = cs_arr[i] = hfg_arr[i] = dpdts_arr[i] = 0.0

    return (temps_R, rhol_arr, rhov_arr, cs_arr, hfg_arr, dpdts_arr)