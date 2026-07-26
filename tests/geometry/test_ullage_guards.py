"""Phase 4.1 / F4 ullage-mass acceptance guards — four-criteria coverage.

Each guard test documents which of the lead-mandated criteria it covers:

1. FIRES ON A BAD CASE (synthetic, or real where one exists).
2. ASSERTS THE RECOVERED STATE, never the defect state.
3. DRIVES THE REAL CODE PATH — no reimplementation of the residual logic.
4. PROVES REACHABILITY ON THE LIVE PATH — production entry points users hit.

Real call chains exercised here
--------------------------------
LOX production path (G0 fire / G3 pass)::

    run_lox_vent_case
      -> run_single_case
      -> summarize_lox_vent_result
           -> ullage_closure_metric          # single residual implementation
           -> assess_lox_dataframe_physicality
                -> ullage_mass_is_acceptable # positivity + 5% threshold

NASA physicality path (negative-ullage synthetic)::

    _dataframe_is_physical
      -> ullage_mass_is_acceptable           # same predicate as LOX

Threshold is the approved 5% (``ULLAGE_CLOSURE_RELATIVE_TOLERANCE``). G0 is
never special-cased; the 5% band is never relaxed.
"""

from __future__ import annotations

from types import SimpleNamespace

import pandas as pd
import pytest

from validation.custom_geometry_cases import _dataframe_is_physical
from validation.lox_vent_cases import (
    G0,
    G3,
    TANK_TOTAL_HEIGHT_FT,
    TIMESTEP_S,
    ULLAGE_CLOSURE_RELATIVE_TOLERANCE,
    ULLAGE_MASS_CLOSURE_CLASSIFICATION,
    assess_lox_dataframe_physicality,
    run_lox_vent_case,
    ullage_closure_metric,
    ullage_mass_is_acceptable,
)


# Committed G0 anchor at 7d6a744 (rate-scaled BL gate matrix): 5.2468...%.
G0_COMMITTED_CLOSURE = 0.05246846816058793


def test_ullage_guard_fires_on_real_g0_via_lox_entry_point() -> None:
    """REAL bad case: G0 full schedule through the live LOX validation entry.

    Criteria covered
    ----------------
    1. FIRES ON BAD CASE — G0 ullage closure ~5.247% > 5% (committed 7d6a744).
    2. RECOVERY-FORM — asserts classification ``ullage_mass_closure`` and
       ``physical is False`` / ``ullage_closure_within_tolerance is False``
       (the recovered guard verdict), not that a broken residual equals a
       hard-coded bad number as a pass condition.
    3. REAL CODE PATH — calls ``run_lox_vent_case`` only; residual comes from
       ``ullage_closure_metric`` inside ``summarize_lox_vent_result``.
    4. REACHABILITY — production call chain documented in the module docstring;
       this test is the user-facing LOX validation entry point.

    Schedule matches the committed G0 matrix row (dt=0.02, duration 60 s).
    """

    assert ULLAGE_CLOSURE_RELATIVE_TOLERANCE == 0.05  # never relaxed

    dataframe, summary = run_lox_vent_case(
        G0, timestep_s=TIMESTEP_S, duration_s=60.0
    )

    # Evidence the case is still the known-bad G0 (band check, not a pin).
    assert summary.ullage_closure_max_relative > ULLAGE_CLOSURE_RELATIVE_TOLERANCE
    assert summary.ullage_closure_max_relative == pytest.approx(
        G0_COMMITTED_CLOSURE, rel=0.05, abs=5e-4
    )
    # Same residual the summary stores (one implementation).
    assert summary.ullage_closure_max_relative == pytest.approx(
        ullage_closure_metric(dataframe)
    )

    # Guard FIRED on the recovered classification / verdict fields.
    assert not summary.ullage_closure_within_tolerance
    assert not summary.physical
    assert ULLAGE_MASS_CLOSURE_CLASSIFICATION in summary.failure_classifications
    assert summary.failure_classifications == (
        ULLAGE_MASS_CLOSURE_CLASSIFICATION,
    )


def test_ullage_guard_passes_on_g3_via_lox_entry_point() -> None:
    """Recovery-form pass case: G3 (~1.24% closure) on the same live path.

    Criteria covered
    ----------------
    2. RECOVERY-FORM — asserts recovered-good state (physical, no
       ``ullage_mass_closure``, within-tolerance True), not a defect pin.
    3. REAL CODE PATH — ``run_lox_vent_case`` + summary threshold wiring.
    4. REACHABILITY — same production entry as the G0 fire test.
    """

    dataframe, summary = run_lox_vent_case(
        G3, timestep_s=TIMESTEP_S, duration_s=60.0
    )

    assert summary.ullage_closure_max_relative <= ULLAGE_CLOSURE_RELATIVE_TOLERANCE
    assert summary.ullage_closure_max_relative == pytest.approx(
        ullage_closure_metric(dataframe)
    )
    assert summary.ullage_closure_within_tolerance is True
    assert summary.physical is True
    assert ULLAGE_MASS_CLOSURE_CLASSIFICATION not in summary.failure_classifications
    assert summary.failure_classifications == ()


def test_ullage_guard_fires_on_synthetic_negative_ullage_via_real_predicates() -> None:
    """Synthetic negative-ullage case for the ``Ullage Mass > 0`` clause.

    Criteria covered
    ----------------
    1. FIRES ON BAD CASE — crafted frame with negative ullage mass.
    2. RECOVERY-FORM — asserts predicate False + classification present.
    3. REAL CODE PATH — drives ``ullage_mass_is_acceptable``,
       ``assess_lox_dataframe_physicality``, and NASA ``_dataframe_is_physical``
       (no reimplemented residual math).
    4. REACHABILITY — both live predicates used by production validation.
    """

    # Minimal two-row frame with otherwise finite physical columns.
    frame = pd.DataFrame(
        {
            "Height": [0.8, 0.81],
            "VBL vol": [0.01, 0.012],
            "BL thick": [0.001, 0.0012],
            "Ullage Mass": [-0.1, -0.05],
            "Ullage from Calc": [0.5, 0.5],
        }
    )
    assert not ullage_mass_is_acceptable(frame)

    physical, classifications = assess_lox_dataframe_physicality(frame)
    assert physical is False
    assert ULLAGE_MASS_CLOSURE_CLASSIFICATION in classifications

    # NASA live predicate accepts any carrier with .dataframe (SingleCaseResult).
    carrier = SimpleNamespace(dataframe=frame)
    assert _dataframe_is_physical(carrier, TANK_TOTAL_HEIGHT_FT) is False  # type: ignore[arg-type]
