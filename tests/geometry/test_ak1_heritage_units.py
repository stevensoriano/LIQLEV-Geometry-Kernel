"""Heritage-comparison regression for AK1 gravity units (finding F1).

The calibrated LIQLEV lineage evaluates the AK1 buoyancy correlation with gravity
in dimensionless standard-g units (VBA ``Ggo``). The kernel currently multiplies
by ``ggo_ft_s2`` at ``core.py:197`` without dividing by ``STD_GRAVITY_FT_S2``,
so AK1 is exactly ``sqrt(32.174) ≈ 5.672213`` times too large.

This test is TDD step 1 (RED on current code). It does not apply the fix.

Reference-input reconstruction
------------------------------
``core._solver_loop`` evaluates properties from state temperature ``t1`` at the
start of each step (``core.py`` hydrogen path: ``_hydrogen_props(t1)``), then
stores end-of-step ``t2`` in the result column ``'Temp'``. Therefore:

* row 0: ``t1 = Tinit`` (input)
* row i > 0: ``t1 = result['Temp'].iloc[i - 1]`` (previous step's ``t2``)

Liquid/vapor densities for the heritage reference are obtained with the same
AS-203 hydrogen polynomials as ``core._hydrogen_props`` (no property table for
``Liquid='Hydrogen'``). Gravity and spacings are the constant schedule values
passed into the legacy-cylinder case.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from core import liqlev_simulation

# Standard gravity (ft/s^2), same constant as heritage STD_GRAVITY_FT_S2.
STD_GRAVITY_FT_S2 = 32.174
EXPECTED_AK1_RATIO_BUG = math.sqrt(STD_GRAVITY_FT_S2)  # ≈ 5.672213

# Relative tolerance required after the core.py:197 fix lands.
POST_FIX_RTOL = 1e-9


# ---------------------------------------------------------------------------
# Vendored heritage reference (NO liqlev_heritage import)
# Provenance: LIQLEV-Heritage-Python @ f84de57
#   src/liqlev_heritage/solver.py:29-50
#   file sha256: ebb55df9b14222abf7de040de7d944100717d4ac77b3acc8e6af9a3c77359ec5
# ---------------------------------------------------------------------------
def boundary_layer_ak1(
    *,
    lateral_spacing: float,
    vertical_spacing: float,
    gravity_ft_s2: float,
    liquid_density_lbm_ft3: float,
    vapor_density_lbm_ft3: float,
) -> float:
    """Evaluate the VBA AK1 correlation with gravity expressed in standard g."""
    if liquid_density_lbm_ft3 == 0.0:
        return 0.0
    gravity_g = gravity_ft_s2 / STD_GRAVITY_FT_S2
    term = (
        10.8
        * (1.0 + lateral_spacing)
        * (1.0 + vertical_spacing)
        * gravity_g
        * (liquid_density_lbm_ft3 - vapor_density_lbm_ft3)
        / liquid_density_lbm_ft3
    )
    return 1.089 * math.sqrt(term) if term > 0.0 else 0.0


def _hydrogen_props(t1: float) -> tuple[float, float]:
    """Same AS-203 polynomials as ``core._hydrogen_props`` (rhol, rhov only)."""
    rhol = (
        0.1709
        + 0.7454 * t1
        - 0.04421 * t1**2
        + 0.001248 * t1**3
        - 1.738e-5 * t1**4
        + 9.424e-8 * t1**5
    )
    rhov = (
        -0.2511
        + 0.04294 * t1
        - 0.00286 * t1**2
        + 9.159e-5 * t1**3
        - 1.422e-6 * t1**4
        + 1.001e-8 * t1**5
    )
    return rhol, rhov


def _t1_for_row(result, row: int, tinit: float) -> float:
    """Reconstruct the solver's start-of-step temperature used for properties."""
    if row == 0:
        return float(tinit)
    return float(result["Temp"].iloc[row - 1])


def _legacy_cylinder_inputs(
    *,
    gravity_g: float = 0.001,
    duration_s: float = 40.0,
    delta_s: float = 10.0,
    spacl: float = 1.0,
    spacv: float = 1.0,
    tinit: float = 38.459461035098116,
) -> dict:
    """Small deterministic legacy-cylinder case (adapted from probe_ak1_units)."""
    time2 = np.array([0.0, duration_s], dtype=np.float64)
    g_ft = gravity_g * STD_GRAVITY_FT_S2
    return {
        "Title": "ak1 heritage units probe",
        "Delta": delta_s,
        "Thetin": 0.0,
        "Units": "British",
        "Liquid": "Hydrogen",
        "Dtank": 21.67,
        "Htzero": 14.09,
        "Volt": (np.pi / 4.0) * 21.67**2 * 28.18,
        "Xmlzro": 3000.0,
        "Pinit": 19.5,
        "Pfinal": 13.8,
        "Tinit": tinit,
        "Nvmd": 2,
        "Tvmdot": time2,
        "Xvmdot": np.full(2, 0.0015, dtype=np.float64),
        "Neps": 0,
        "Teps": None,
        "Xeps": None,
        "Nlattm": 2,
        "Tspal": time2,
        "Xspacl": np.full(2, spacl, dtype=np.float64),
        "Nvertm": 2,
        "Tspav": time2,
        "Xspacv": np.full(2, spacv, dtype=np.float64),
        "Nggo": 2,
        "Tggo": time2,
        "Xggo": np.full(2, g_ft, dtype=np.float64),
    }


def test_ak1_matches_heritage_standard_g_units() -> None:
    """Kernel AK1 must match heritage boundary_layer_ak1 (standard-g gravity).

    On current code (F1 unfixed) this fails with kernel/ref ≈ sqrt(32.174).
    """
    gravity_g = 0.001
    spacl = 1.0
    spacv = 1.0
    tinit = 38.459461035098116
    g_ft = gravity_g * STD_GRAVITY_FT_S2
    inputs = _legacy_cylinder_inputs(
        gravity_g=gravity_g,
        spacl=spacl,
        spacv=spacv,
        tinit=tinit,
    )
    result = liqlev_simulation(inputs, verbose=False)
    assert len(result) >= 2, "need at least two result rows for the comparison"

    for row in (0, 1):
        t1 = _t1_for_row(result, row, tinit)
        rhol, rhov = _hydrogen_props(t1)
        reference = boundary_layer_ak1(
            lateral_spacing=spacl,
            vertical_spacing=spacv,
            gravity_ft_s2=g_ft,
            liquid_density_lbm_ft3=rhol,
            vapor_density_lbm_ft3=rhov,
        )
        kernel_ak1 = float(result["AK1"].iloc[row])
        assert reference > 0.0
        ratio = kernel_ak1 / reference
        rel = abs(kernel_ak1 - reference) / abs(reference)
        assert kernel_ak1 == pytest.approx(reference, rel=POST_FIX_RTOL), (
            f"row {row}: kernel AK1={kernel_ak1:.15g} reference={reference:.15g} "
            f"ratio(kernel/ref)={ratio:.15g} expected_bug_ratio={EXPECTED_AK1_RATIO_BUG:.15g} "
            f"rel_err={rel:.3e} (post-fix rtol={POST_FIX_RTOL})"
        )


# NOTE (review, 2026-07-25): a zero-gravity parity control was drafted here and removed —
# the legacy solver at exactly g = 0 is practically non-terminating on this case (finding
# F3 territory; >90 s for a 4-step schedule). Zero-g behaviour is exercised by the custom
# geometry path once plan Phase 2 lands; legacy g = 0 stays baseline-locked by decision.
