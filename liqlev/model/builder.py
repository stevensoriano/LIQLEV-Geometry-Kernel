"""Solver input construction and small physics-adjacent helpers.

These functions preserve the legacy GUI input shape expected by
``core.liqlev_simulation``. They do not change solver formulas.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import Any

import numpy as np

from thermo_utils import DensitySat, Psat, Tsat


PSI_TO_KPA = 6.89475729
G_TO_FT_S2 = 32.174


def build_inputs(
    fluid: str,
    pinit_psia: float,
    pfinal_psia: float,
    dtank: float,
    htank: float,
    fill_fraction: float,
    duration: float,
    delta_t: float,
    vent_rate: float,
    neps: int | None,
    teps: np.ndarray | None,
    xeps: np.ndarray | None,
    ramp_duration: float,
    ramp_target_factor: float,
    nggo: int,
    tggo: np.ndarray,
    xggo: np.ndarray,
    gravity_function: Callable[[float], float] | None = None,
    xmlzro_override: float | None = None,
    tinit_override: float | None = None,
) -> dict[str, Any]:
    """Build the input dictionary consumed by ``core.liqlev_simulation``."""
    volt = (np.pi / 4) * (dtank**2) * htank
    ac = 0.7854 * (dtank**2)
    htzero = fill_fraction * htank

    press_kpa = pinit_psia * PSI_TO_KPA
    tinit = Tsat(fluid, press_kpa) * 1.8

    if tinit_override is not None:
        tinit = tinit_override

    if fluid == "Hydrogen":
        rhol = (
            0.1709
            + 0.7454 * tinit
            - 0.04421 * tinit**2
            + 0.001248 * tinit**3
            - 1.738e-5 * tinit**4
            + 9.424e-8 * tinit**5
        )
    else:
        t_k = tinit / 1.8
        ps_kpa = Psat(fluid, t_k)
        rhol = DensitySat(fluid, "liquid", ps_kpa) * 0.0624279606

    xmlzro = rhol * (htzero * ac)
    if xmlzro_override is not None:
        xmlzro = xmlzro_override

    if duration > ramp_duration:
        tvmdot = np.array([0.0, ramp_duration, duration])
        xvmdot = np.array(
            [
                vent_rate,
                vent_rate * ramp_target_factor,
                vent_rate * ramp_target_factor,
            ]
        )
    else:
        slope = (vent_rate * ramp_target_factor - vent_rate) / ramp_duration
        end_rate = vent_rate + slope * duration
        tvmdot = np.array([0.0, duration])
        xvmdot = np.array([vent_rate, end_rate])

    if neps is None:
        neps = 11
        teps = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, duration])
        xeps = np.array(
            [0.0, 0.0513, 0.178, 0.28, 0.362, 0.422, 0.47, 0.52, 0.56, 0.6, 0.6]
        )

    return {
        "Title": "LIQLEV GUI Simulation",
        "Liquid": fluid,
        "Units": "British",
        "Delta": delta_t,
        "Dtank": dtank,
        "Htzero": htzero,
        "Volt": volt,
        "Xmlzro": xmlzro,
        "Pinit": pinit_psia,
        "Pfinal": pfinal_psia,
        "Tinit": tinit,
        "Thetin": 0.0,
        "Nvmd": len(tvmdot),
        "Neps": neps,
        "Nlattm": 2,
        "Nvertm": 2,
        "Nggo": nggo,
        "Tvmdot": tvmdot,
        "Xvmdot": xvmdot,
        "Teps": teps,
        "Xeps": xeps,
        "Tspal": np.array([0.0, duration]),
        "Xspacl": np.array([1.0, 1.0]),
        "Tspav": np.array([0.0, duration]),
        "Xspacv": np.array([1.0, 1.0]),
        "Tggo": tggo,
        "Xggo": xggo,
        "gravity_function": gravity_function,
    }


def calculate_epsilon(h: float, dtank: float) -> float:
    """Return geometry-based wall-contact epsilon for the tank fill height."""
    perim = np.pi * dtank
    a_wall = perim * h
    ac = 0.7854 * (dtank**2)
    return a_wall / (a_wall + ac) if (a_wall + ac) != 0 else 0


def epsilon_schedule(
    eps_spec: str | float, duration: float
) -> tuple[int, np.ndarray | None, np.ndarray | None, str]:
    """Return the legacy Neps/Teps/Xeps schedule tuple for an epsilon selection."""
    if eps_spec == "height_dep":
        return 0, None, None, "height_dep"
    if eps_spec == "bulk_fake":
        return 2, np.array([0.0, duration]), np.array([50.0, 50.0]), "bulk_fake"
    if eps_spec == "AS-203 Schedule":
        return (
            11,
            np.array(
                [
                    0.0,
                    20.0,
                    40.0,
                    60.0,
                    80.0,
                    100.0,
                    120.0,
                    140.0,
                    160.0,
                    180.0,
                    duration,
                ]
            ),
            np.array(
                [
                    0.0000,
                    0.0513,
                    0.1780,
                    0.2800,
                    0.3620,
                    0.4220,
                    0.4700,
                    0.5200,
                    0.5600,
                    0.6000,
                    0.6000,
                ]
            ),
            "AS-203 Schedule",
        )

    value = float(eps_spec)
    return 2, np.array([0.0, duration]), np.array([value, value]), f"{value:.4f}"


_GRAV_SAFE_DICT = {
    "sin": math.sin,
    "cos": math.cos,
    "tan": math.tan,
    "asin": math.asin,
    "acos": math.acos,
    "atan": math.atan,
    "exp": math.exp,
    "log": math.log,
    "log10": math.log10,
    "sqrt": math.sqrt,
    "abs": abs,
    "pi": math.pi,
    "e": math.e,
    "min": min,
    "max": max,
    "pow": pow,
}


def safe_eval_gravity(expr: str, t: float) -> float:
    """Evaluate a gravity expression in g units with the legacy safe namespace."""
    safe = dict(_GRAV_SAFE_DICT)
    safe["t"] = t
    safe["__builtins__"] = {}
    return float(eval(expr, safe))


def make_gravity_function(expr: str) -> Callable[[float], float]:
    """Create a callable returning gravity in ft/s^2 from an expression in g units."""

    def gravity_func(t_sec: float) -> float:
        g_level = safe_eval_gravity(expr, t_sec)
        return g_level * G_TO_FT_S2

    return gravity_func
