"""Canonical physics cases for LIQLEV result-preservation checks.

These helpers intentionally avoid importing the GUI. They build the same input
shape expected by core.liqlev_simulation so baseline checks can run headlessly.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Tuple

import numpy as np
import pandas as pd

from core import liqlev_simulation
from thermo_utils import DensitySat, Psat, Tsat, build_property_table


PSI_TO_KPA = 6.89475729
G_TO_FT_S2 = 32.174


@dataclass(frozen=True)
class PhysicsCase:
    name: str
    description: str
    fluid: str
    pinit_psia: float
    pfinal_psia: float
    dtank_ft: float
    htank_ft: float
    fill_fraction: float
    duration_s: float
    delta_t_s: float
    vent_rate_lbm_s: float
    epsilon_mode: str
    gravity_g: float
    ramp_duration_s: Optional[float] = None
    ramp_target_factor: float = 1.0
    epsilon_value: Optional[float] = None
    xmlzro_override_lbm: Optional[float] = None
    tinit_override_r: Optional[float] = None


def hydrogen_density_lbm_ft3(t_rankine: float) -> float:
    """Legacy hydrogen saturated-liquid density polynomial."""
    return (
        0.1709
        + 0.7454 * t_rankine
        - 0.04421 * t_rankine**2
        + 0.001248 * t_rankine**3
        - 1.738e-5 * t_rankine**4
        + 9.424e-8 * t_rankine**5
    )


def epsilon_schedule(
    mode: str,
    duration_s: float,
    epsilon_value: Optional[float] = None,
) -> Tuple[int, Optional[np.ndarray], Optional[np.ndarray]]:
    """Return the Neps/Teps/Xeps triple expected by the solver input dict."""
    if mode == "height_dep":
        return 0, None, None

    if mode == "bulk_fake":
        return (
            2,
            np.array([0.0, duration_s], dtype=np.float64),
            np.array([50.0, 50.0], dtype=np.float64),
        )

    if mode == "AS-203 Schedule":
        return (
            11,
            np.array(
                [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, duration_s],
                dtype=np.float64,
            ),
            np.array(
                [0.0000, 0.0513, 0.1780, 0.2800, 0.3620, 0.4220, 0.4700, 0.5200, 0.5600, 0.6000, 0.6000],
                dtype=np.float64,
            ),
        )

    if mode == "custom":
        if epsilon_value is None:
            raise ValueError("custom epsilon mode requires epsilon_value")
        value = float(epsilon_value)
        return (
            2,
            np.array([0.0, duration_s], dtype=np.float64),
            np.array([value, value], dtype=np.float64),
        )

    raise ValueError(f"Unknown epsilon mode: {mode}")


def build_case_inputs(case: PhysicsCase) -> Dict[str, object]:
    """Build a solver input dictionary for a canonical physics case."""
    volt = (np.pi / 4.0) * case.dtank_ft**2 * case.htank_ft
    ac = 0.7854 * case.dtank_ft**2
    htzero = case.fill_fraction * case.htank_ft

    press_kpa = case.pinit_psia * PSI_TO_KPA
    tinit = Tsat(case.fluid, press_kpa) * 1.8
    if case.tinit_override_r is not None:
        tinit = case.tinit_override_r

    if case.fluid == "Hydrogen":
        rhol = hydrogen_density_lbm_ft3(tinit)
    else:
        t_k = tinit / 1.8
        ps_kpa = Psat(case.fluid, t_k)
        rhol = DensitySat(case.fluid, "liquid", ps_kpa) * 0.0624279606

    xmlzro = rhol * (htzero * ac)
    if case.xmlzro_override_lbm is not None:
        xmlzro = case.xmlzro_override_lbm

    ramp_duration = case.ramp_duration_s
    if ramp_duration is None:
        ramp_duration = case.duration_s

    if case.duration_s > ramp_duration:
        tvmdot = np.array([0.0, ramp_duration, case.duration_s], dtype=np.float64)
        xvmdot = np.array(
            [
                case.vent_rate_lbm_s,
                case.vent_rate_lbm_s * case.ramp_target_factor,
                case.vent_rate_lbm_s * case.ramp_target_factor,
            ],
            dtype=np.float64,
        )
    else:
        slope = (case.vent_rate_lbm_s * case.ramp_target_factor - case.vent_rate_lbm_s) / ramp_duration
        end_rate = case.vent_rate_lbm_s + slope * case.duration_s
        tvmdot = np.array([0.0, case.duration_s], dtype=np.float64)
        xvmdot = np.array([case.vent_rate_lbm_s, end_rate], dtype=np.float64)

    neps, teps, xeps = epsilon_schedule(case.epsilon_mode, case.duration_s, case.epsilon_value)

    g_ft_s2 = case.gravity_g * G_TO_FT_S2
    tggo = np.array([0.0, case.duration_s], dtype=np.float64)
    xggo = np.array([g_ft_s2, g_ft_s2], dtype=np.float64)

    return {
        "Title": f"Baseline case: {case.name}",
        "Liquid": case.fluid,
        "Units": "British",
        "Delta": case.delta_t_s,
        "Dtank": case.dtank_ft,
        "Htzero": htzero,
        "Volt": volt,
        "Xmlzro": xmlzro,
        "Pinit": case.pinit_psia,
        "Pfinal": case.pfinal_psia,
        "Tinit": tinit,
        "Thetin": 0.0,
        "Nvmd": len(tvmdot),
        "Neps": neps,
        "Nlattm": 2,
        "Nvertm": 2,
        "Nggo": 2,
        "Tvmdot": tvmdot,
        "Xvmdot": xvmdot,
        "Teps": teps,
        "Xeps": xeps,
        "Tspal": np.array([0.0, case.duration_s], dtype=np.float64),
        "Xspacl": np.array([1.0, 1.0], dtype=np.float64),
        "Tspav": np.array([0.0, case.duration_s], dtype=np.float64),
        "Xspacv": np.array([1.0, 1.0], dtype=np.float64),
        "Tggo": tggo,
        "Xggo": xggo,
        "gravity_function": None,
    }


def iter_cases() -> Iterable[PhysicsCase]:
    """Return canonical cases covering main solver paths."""
    return (
        PhysicsCase(
            name="as203_default_high_vent",
            description="AS-203 hydrogen default with measured mass/temperature and schedule epsilon.",
            fluid="Hydrogen",
            pinit_psia=19.5,
            pfinal_psia=13.8,
            dtank_ft=21.670,
            htank_ft=28.18,
            fill_fraction=0.5116,
            duration_s=400.0,
            delta_t_s=10.0,
            vent_rate_lbm_s=3.3069,
            epsilon_mode="AS-203 Schedule",
            gravity_g=0.00000963,
            ramp_duration_s=400.0,
            xmlzro_override_lbm=16300.0,
            tinit_override_r=38.3,
        ),
        PhysicsCase(
            name="hydrogen_height_dep_mid_fill",
            description="Hydrogen height-dependent epsilon path at mid fill.",
            fluid="Hydrogen",
            pinit_psia=19.5,
            pfinal_psia=13.8,
            dtank_ft=21.670,
            htank_ft=28.18,
            fill_fraction=0.50,
            duration_s=300.0,
            delta_t_s=10.0,
            vent_rate_lbm_s=0.0015,
            epsilon_mode="height_dep",
            gravity_g=0.001,
            ramp_duration_s=300.0,
        ),
        PhysicsCase(
            name="nitrogen_custom_epsilon",
            description="Non-hydrogen property table path with custom epsilon.",
            fluid="Nitrogen",
            pinit_psia=30.0,
            pfinal_psia=20.0,
            dtank_ft=5.0,
            htank_ft=10.0,
            fill_fraction=0.45,
            duration_s=160.0,
            delta_t_s=5.0,
            vent_rate_lbm_s=0.02,
            epsilon_mode="custom",
            epsilon_value=0.4,
            gravity_g=0.002,
            ramp_duration_s=160.0,
        ),
    )


def get_case(name: str) -> PhysicsCase:
    for case in iter_cases():
        if case.name == name:
            return case
    available = ", ".join(case.name for case in iter_cases())
    raise KeyError(f"Unknown case {name!r}. Available cases: {available}")


def run_case(case: PhysicsCase) -> pd.DataFrame:
    """Run one canonical case and return the solver DataFrame."""
    inputs = build_case_inputs(case)
    prop_table = None
    if case.fluid != "Hydrogen":
        prop_table = build_property_table(case.fluid, case.pfinal_psia, case.pinit_psia)
    return liqlev_simulation(inputs, verbose=False, prop_table=prop_table)


def dataframe_summary(df: pd.DataFrame) -> Dict[str, float]:
    """Return compact engineering metrics for diagnostics and logs."""
    if df.empty:
        return {
            "rows": 0,
            "max_dh_h0": 0.0,
            "time_to_peak_s": 0.0,
            "final_pressure_psia": 0.0,
            "final_height_ft": 0.0,
            "convergence_failures": 0,
        }

    peak_idx = df["Hratio"].idxmax()
    return {
        "rows": int(len(df)),
        "max_dh_h0": float(df.loc[peak_idx, "Hratio"]),
        "time_to_peak_s": float(df.loc[peak_idx, "Time"]),
        "final_pressure_psia": float(df["Press"].iloc[-1]),
        "final_height_ft": float(df["Height"].iloc[-1]),
        "convergence_failures": int(df["Conv Failed"].sum()) if "Conv Failed" in df.columns else 0,
    }

