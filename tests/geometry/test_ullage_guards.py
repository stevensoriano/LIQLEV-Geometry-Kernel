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

from functools import lru_cache
from types import SimpleNamespace

import pandas as pd
import pytest

from liqlev.runner.single import run_single_case
from validation.custom_geometry_cases import (
    _dataframe_is_physical,
    build_nasa_tank_config,
)
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

# Reported (NOT gated) wall-film mass-drift anchors for the same G0 run, from
# the closed mass-ledger audit run of record
# LIQLEV_vent_study_spyder/results_vapor_balance/2026-07-29_204715/report.txt
# ("A SECOND UNRECONCILED RESERVOIR — the wall film (surfaced, NOT gated)"):
#   row   retention    EOS inv      drift  drift/inv %
#   G0      0.99838    0.94055    0.05783         6.15
G0_AUDIT_FILM_RETENTION_LBM = 0.99838
G0_AUDIT_FILM_EOS_INVENTORY_LBM = 0.94055
G0_AUDIT_FILM_DRIFT_LBM = 0.05783
G0_AUDIT_FILM_DRIFT_OVER_INVENTORY = 0.0615


@lru_cache(maxsize=1)
def _g0_production_run():
    """The ONE full-schedule G0 production run shared by the G0 guards.

    Both the ullage-closure fire test and the film-drift report test assert
    against the same committed G0 matrix row (dt = 0.02, duration 60 s); running
    it once keeps the suite from paying for a second 60 s solve.
    """

    return run_lox_vent_case(G0, timestep_s=TIMESTEP_S, duration_s=60.0)


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

    dataframe, summary = _g0_production_run()

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


def test_g0_film_drift_is_reported_not_gated_and_matches_the_audit() -> None:
    """Production anchor for the REPORTED wall-film mass-drift metric at G0.

    Reference values are the audit run of record cited beside the constants
    above (closed mass ledger to 6e-13 lbm), reconstructed here through the
    production summary path — same shared G0 run as the closure guard.

    The point of this test is twofold:
    1. the magnitude is surfaced per run (assumption 6 / F4: "must be reported
       per run, not assumed small"), and
    2. surfacing it changes NO acceptance verdict — G0's only classification is
       still ``ullage_mass_closure``, exactly as before the metric existed.
    """

    _dataframe, summary = _g0_production_run()

    assert summary.film_retention_lbm == pytest.approx(
        G0_AUDIT_FILM_RETENTION_LBM, rel=0.05
    )
    assert summary.film_eos_inventory_lbm == pytest.approx(
        G0_AUDIT_FILM_EOS_INVENTORY_LBM, rel=0.05
    )
    assert summary.film_drift_lbm == pytest.approx(
        G0_AUDIT_FILM_DRIFT_LBM, rel=0.05
    )
    assert summary.film_drift_lbm == pytest.approx(
        summary.film_retention_lbm - summary.film_eos_inventory_lbm
    )
    assert summary.film_drift_over_inventory == pytest.approx(
        summary.film_drift_lbm / summary.film_eos_inventory_lbm
    )
    assert summary.film_drift_over_inventory == pytest.approx(
        G0_AUDIT_FILM_DRIFT_OVER_INVENTORY, rel=0.05
    )
    # The audit's headline comparison — drift = 3.8x the ullage gap D at G0 — is
    # mass against mass, and D is not a summary field (only its normalized form
    # ullage_closure_max_relative is), so it is recorded in docs/vent-study.md
    # rather than asserted here against a dimensionless ratio.
    #
    # What IS asserted: surfacing the drift changed no verdict. G0 still fails on
    # ullage closure alone, and the 5% band is untouched.
    assert summary.failure_classifications == (
        ULLAGE_MASS_CLOSURE_CLASSIFICATION,
    )
    assert not summary.ullage_closure_within_tolerance
    assert not summary.physical
    assert ULLAGE_CLOSURE_RELATIVE_TOLERANCE == 0.05  # never relaxed


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


@pytest.mark.xfail(
    strict=True,
    reason=(
        "F4: hydrogen stand-in (3000x volume-mismatched acceptance case) fails "
        "strict ullage closure at high fill — measured 14.6% @ 0.75 / ~121% + "
        "negative ullage @ 0.90 (pre-campaign 11.5%/34.4%; deepened by the "
        "approved F1 physics correction). Fires when the bookkeeping is "
        "reconciled or the stand-in is retired (path-forward: LOX case replaces "
        "it)."
    ),
)
def test_nasa_hydrogen_standin_ullage_passes_strict_guard() -> None:
    """F4 quarantine tripwire: hydrogen stand-in fill 0.90 must pass strict ullage.

    Lead ruling (supervisor attribution confirmed): the guard stays STRICT on
    LOX/custom; this 3000x volume-mismatched NASA hydrogen stand-in is
    quarantined as xfail(strict). Pre-existing refinement tests are untouched.

    Criteria covered
    ----------------
    1. FIRES ON BAD CASE — live production fill 0.90 currently fails strict
       ullage (closure ~121%, negative ullage mass); xfail(strict) records that
       pathology until reconcile/retire.
    2. RECOVERY-FORM — asserts the recovered-good state (predicate True, no
       ``ullage_mass_closure``, all Ullage Mass > 0). Fails today → xfail
       satisfied; XPASS later forces marker removal (tripwire).
    3. REAL CODE PATH — ``build_nasa_tank_config`` (same builder as
       ``run_nasa_tank_validation``) + ``run_single_case`` at fill 0.90 on the
       production geometry package; strict predicate is the s6a LOX function
       ``assess_lox_dataframe_physicality`` / ``ullage_mass_is_acceptable``
       (no residual reimplementation).
    4. REACHABILITY — production NASA config builder + single-case runner +
       LOX strict ullage gate used on the live LOX path.
    """

    assert ULLAGE_CLOSURE_RELATIVE_TOLERANCE == 0.05  # never relaxed

    # Same builder run_nasa_tank_validation uses; production grid only (one fill).
    config = build_nasa_tank_config()
    result = run_single_case(config, fill_fraction=0.90)
    dataframe = result.dataframe

    # s6a LOX strict ullage predicate + classification (drive, do not copy).
    assert ullage_mass_is_acceptable(dataframe) is True
    physical, classifications = assess_lox_dataframe_physicality(dataframe)
    assert physical is True
    assert ULLAGE_MASS_CLOSURE_CLASSIFICATION not in classifications
    assert classifications == ()
    assert bool((dataframe["Ullage Mass"] > 0.0).all())
