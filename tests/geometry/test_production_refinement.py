"""Phase 4.2 / F5 production(92)-vs-reference(1025) refinement metric.

Does NOT touch the frozen NASA result manifest or the D1 test. The existing
``maximum_refinement_difference`` (coarse 513 vs reference 1025) stays
bit-identical; this suite only asserts the companion production-vs-reference
metric exposed on ``NasaTankValidation``.

Four-criteria coverage
----------------------
1. FIRES ON BAD CASE — not applicable as a failure-mode guard; the metric is a
   pass gate vs 2e-3. A synthetic over-tolerance would reimplement comparison
   math, so this file pins the recovery-form PASS on the live path instead.
2. RECOVERY-FORM — asserts the live metric is finite, positive, and <= 2e-3
   (recovered-good production grid), expected order ~5e-5 per the review probe.
3. REAL CODE PATH — calls ``run_nasa_tank_validation()``; reads
   ``maximum_production_refinement_difference`` computed inside that function
   via the same ``_refinement_metrics`` helper used for coarse-vs-reference.
4. REACHABILITY — ``run_nasa_tank_validation`` is the production NASA
   validation entry (also used by the D1 fixture and CLI ``main()``).
"""

from __future__ import annotations

import numpy as np
import pytest

from validation.custom_geometry_cases import (
    REFINEMENT_RELATIVE_TOLERANCE,
    run_nasa_tank_validation,
)


def test_production_vs_reference_refinement_within_approved_limit() -> None:
    """Live NASA path: production(92) vs reference(1025) worst relative diff.

    Expected ~5.3e-05 (review probe) against REFINEMENT_RELATIVE_TOLERANCE
    (2.0e-3). Coarse-vs-reference ``maximum_refinement_difference`` is left
    untouched for the frozen D1 manifest comparison.
    """

    validation = run_nasa_tank_validation()

    metric = validation.maximum_production_refinement_difference
    assert np.isfinite(metric)
    assert metric >= 0.0
    # Recovery-form pass vs the approved grid-refinement band (2.0e-3).
    assert metric <= REFINEMENT_RELATIVE_TOLERANCE
    # Live order is ~1e-6..1e-4 (review probe quoted ~5.3e-5; post-F1/F2
    # integration sits near ~6e-6). Gate is the 2e-3 tolerance, not the probe.
    assert metric < 1.0e-3
    assert metric > 0.0

    # Per-fill companion map is populated and consistent with the aggregate.
    assert set(validation.production_refinement_by_fill) == set(
        validation.refinement_by_fill
    )
    per_fill_max = max(
        m.maximum_relative_difference
        for m in validation.production_refinement_by_fill.values()
    )
    assert metric == pytest.approx(per_fill_max, rel=0.0, abs=0.0)

    # D1 boundary: coarse-vs-reference field remains the one the frozen
    # manifest pins; this test never writes or regenerates that manifest.
    assert np.isfinite(validation.maximum_refinement_difference)
    assert "production_evaluation_grid_refinement" not in (
        validation.failure_classifications
    )
