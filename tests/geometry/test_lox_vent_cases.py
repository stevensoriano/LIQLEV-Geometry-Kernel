"""Bounded tests for the LOX 43 L vent acceptance-case machinery."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from liqlev.geometry.package import load_geometry_package
from liqlev.geometry.jit import invert_monotone_volume
from liqlev.runner.single import SingleCaseResult, run_single_case
from validation.lox_vent_cases import (
    DURATION_S,
    FILL_FRACTION,
    FLUID,
    FLUID_STEP_PATH,
    G0,
    G3,
    GEOMETRY_NPZ_PATH,
    INITIAL_HEIGHT_FT,
    LOX_VENT_MODULE_PATH,
    PFINAL_PSIA,
    PINIT_PSIA,
    TANK_TOTAL_HEIGHT_FT,
    TIMESTEP_S,
    VENT_RATE_LBM_S,
    _summary_manifest,
    build_lox_vent_config,
    film_drift_metric,
    run_lox_vent_case,
    sha256_file,
    summarize_lox_vent_result,
    ullage_closure_metric,
    write_lox_vent_manifest,
)


# Manifest keys as committed BEFORE the reported film-drift metric was added.
# Asserted as a subset so the feature is proven additive-only.
PRE_FILM_DRIFT_MANIFEST_KEYS = frozenset(
    {
        "gravity_g",
        "timestep_s",
        "rows",
        "t_end_s",
        "p_end_psia",
        "h_initial_ft",
        "h_final_ft",
        "dh_ft",
        "dh_over_h0",
        "final_vbl_vol_ft3",
        "final_bl_thick_ft",
        "max_ak3",
        "conv_failed_total",
        "ullage_closure_max_relative",
        "ullage_closure_within_tolerance",
        "ullage_closure_relative_tolerance",
        "finite",
        "physical",
        "failure_classifications",
    }
)

FILM_DRIFT_FIELDS = (
    "film_retention_lbm",
    "film_eos_inventory_lbm",
    "film_drift_lbm",
    "film_drift_over_inventory",
)


@lru_cache(maxsize=1)
def _short_g3_result() -> SingleCaseResult:
    """One shared 1 s G3 run — the fast pattern used by the smoke tests.

    Shared (not re-run per test) so the film-drift tests cost one solver run and
    one property-table build; the no-mutation test below is what proves sharing
    is safe.
    """

    return run_single_case(
        build_lox_vent_config(G3, timestep_s=0.02, duration_s=1.0)
    )


def test_lox_vent_config_matches_test_definition() -> None:
    config = build_lox_vent_config(G3)

    assert config.fluid.name == FLUID
    assert config.fluid.initial_pressure_psia == PINIT_PSIA
    assert config.fluid.final_pressure_psia == PFINAL_PSIA
    assert config.tank.fill_fractions == (FILL_FRACTION,)
    assert config.tank.geometry_path == str(GEOMETRY_NPZ_PATH)
    assert config.vent.rates_lbm_s == (VENT_RATE_LBM_S,)
    assert config.vent.ramp_target_factor == 1.0
    assert config.gravity.mode == "Constant"
    assert config.gravity.constant_g == pytest.approx(G3)
    assert config.epsilon.mode == "height_dep"
    assert config.run.duration_s == DURATION_S
    assert config.run.timestep_s == TIMESTEP_S

    geometry = load_geometry_package(GEOMETRY_NPZ_PATH)
    h0 = float(
        invert_monotone_volume(
            FILL_FRACTION * geometry.total_volume_ft3,
            geometry.height_ft,
            geometry.volume_ft3,
            geometry.volume_coefficients,
        )
    )
    assert h0 == pytest.approx(INITIAL_HEIGHT_FT, rel=0.0, abs=5e-6)


def test_lox_vent_short_smoke_g3_summary_and_ullage_metric() -> None:
    """Few-second smoke at G3: physical rows + summary fields compute."""

    dataframe, summary = run_lox_vent_case(G3, timestep_s=0.02, duration_s=1.0)

    assert summary.rows >= 2
    assert summary.finite
    assert summary.conv_failed_total == 0
    assert np.isfinite(dataframe.to_numpy(dtype=float)).all()
    assert summary.h_initial_ft == pytest.approx(INITIAL_HEIGHT_FT, abs=5e-4)
    assert 0.0 < summary.h_final_ft < 1.806770
    assert summary.t_end_s == pytest.approx(1.0, abs=summary.timestep_s * 1.5)
    assert PFINAL_PSIA <= summary.p_end_psia <= PINIT_PSIA + 0.5
    assert np.isfinite(summary.final_vbl_vol_ft3)
    assert np.isfinite(summary.final_bl_thick_ft)
    assert np.isfinite(summary.max_ak3)
    assert summary.max_ak3 >= 0.0
    assert np.isfinite(summary.ullage_closure_max_relative)
    assert summary.ullage_closure_max_relative >= 0.0
    # Metric is measured only — not an acceptance threshold (Phase 4 owns that).
    recomputed = ullage_closure_metric(dataframe)
    assert summary.ullage_closure_max_relative == pytest.approx(recomputed)


def test_lox_vent_g0_accepted_single_short_step() -> None:
    """G0 is legal input after Step 3 custom-mode zero-g support."""

    # duration must exceed timestep (validator); one solver step uses 2*dt window.
    dataframe, summary = run_lox_vent_case(G0, timestep_s=0.02, duration_s=0.04)

    assert summary.rows >= 1
    assert summary.finite
    assert summary.conv_failed_total == 0
    assert np.isfinite(dataframe.to_numpy(dtype=float)).all()
    assert summary.gravity_g == 0.0


def test_film_drift_metric_on_short_g3_run() -> None:
    """Reported wall-film mass-drift metric: shape and sign on a live run.

    The film is carried across steps as a VOLUME (``core.py:567`` ``vbl1 =
    vbl2``) while saturated vapour density falls, so the cumulative mass the
    ullage was debited to fill it (retention) can only exceed the film's
    equation-of-state inventory ``vbl * rho_v``. Reported, never gated.
    """

    result = _short_g3_result()
    metric = film_drift_metric(result.dataframe, result.inputs)

    assert np.isfinite(metric.film_retention_lbm)
    assert np.isfinite(metric.film_eos_inventory_lbm)
    assert np.isfinite(metric.film_drift_lbm)
    assert np.isfinite(metric.film_drift_over_inventory)
    # retention >= inventory >= 0 while rho_v falls through the blowdown.
    assert metric.film_eos_inventory_lbm > 0.0
    assert metric.film_retention_lbm >= metric.film_eos_inventory_lbm
    assert metric.film_drift_lbm >= 0.0
    assert metric.film_drift_lbm == pytest.approx(
        metric.film_retention_lbm - metric.film_eos_inventory_lbm
    )
    assert metric.film_drift_over_inventory == pytest.approx(
        metric.film_drift_lbm / metric.film_eos_inventory_lbm
    )


def test_film_drift_metric_is_nan_on_degenerate_input() -> None:
    """Degenerate frames / missing run inputs yield NaN, mirroring closure."""

    inputs = {
        "Liquid": FLUID,
        "Tinit": 175.0,
        "Pinit": PINIT_PSIA,
        "Pfinal": PFINAL_PSIA,
    }
    one_row = pd.DataFrame({"Temp": [175.0], "VBL vol": [0.01]})
    degenerate = (
        film_drift_metric(pd.DataFrame(), inputs),          # empty frame
        film_drift_metric(one_row, None),                   # no run inputs
        film_drift_metric(pd.DataFrame({"Press": [40.0]}), inputs),  # no columns
        film_drift_metric(
            pd.DataFrame({"Temp": [np.nan], "VBL vol": [0.01]}), inputs
        ),                                                  # non-finite state
    )
    for metric in degenerate:
        assert not np.isfinite(metric.film_retention_lbm)
        assert not np.isfinite(metric.film_eos_inventory_lbm)
        assert not np.isfinite(metric.film_drift_lbm)
        assert not np.isfinite(metric.film_drift_over_inventory)


def test_film_drift_metric_does_not_mutate_the_dataframe() -> None:
    """Instrumentation only: the metric never writes to the run's dataframe."""

    result = _short_g3_result()
    before = result.dataframe.copy(deep=True)

    film_drift_metric(result.dataframe, result.inputs)

    pd.testing.assert_frame_equal(result.dataframe, before)


def test_summary_and_manifest_carry_the_film_drift_fields() -> None:
    """The 4 reported fields reach the summary and the manifest, additively."""

    result = _short_g3_result()
    summary = summarize_lox_vent_result(result, gravity_g=G3, timestep_s=0.02)
    metric = film_drift_metric(result.dataframe, result.inputs)

    # One measurement implementation: the summary stores exactly the metric.
    for field in FILM_DRIFT_FIELDS:
        assert getattr(summary, field) == pytest.approx(getattr(metric, field))

    payload = _summary_manifest(summary)
    for field in FILM_DRIFT_FIELDS:
        assert payload[field] == pytest.approx(getattr(summary, field))
    # Additive only: every pre-feature key is still present.
    assert PRE_FILM_DRIFT_MANIFEST_KEYS <= set(payload)
    # Reported, never gated: no new acceptance verdict rides along.
    assert summary.physical is True
    assert summary.failure_classifications == ()


def test_empty_dataframe_summary_reports_nan_film_drift() -> None:
    """Empty-frame branch mirrors ullage closure: NaN fields, None in manifest."""

    empty = SingleCaseResult(
        dataframe=pd.DataFrame(),
        inputs={},
        scenario_key="empty",
        vent_rate_lbm_s=VENT_RATE_LBM_S,
        fill_fraction=FILL_FRACTION,
        epsilon_label="height_dep",
        htank_ft=TANK_TOTAL_HEIGHT_FT,
        elapsed_s=0.0,
    )
    summary = summarize_lox_vent_result(empty, gravity_g=G0, timestep_s=0.02)

    for field in FILM_DRIFT_FIELDS:
        assert not np.isfinite(getattr(summary, field))
    payload = _summary_manifest(summary)
    for field in FILM_DRIFT_FIELDS:
        assert payload[field] is None
    assert payload["ullage_closure_max_relative"] is None


def test_lox_vent_manifest_refuses_dirty_tree(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(
        "validation.lox_vent_cases._git_worktree_is_dirty",
        lambda: True,
    )
    target = tmp_path / "lox_vent_manifest.json"
    with pytest.raises(RuntimeError, match="dirty git worktree"):
        write_lox_vent_manifest({}, target)


def test_lox_vent_manifest_hash_binding_and_clean_write(
    monkeypatch, tmp_path: Path
) -> None:
    monkeypatch.setattr(
        "validation.lox_vent_cases._git_worktree_is_dirty",
        lambda: False,
    )
    monkeypatch.setattr(
        "validation.lox_vent_cases._git_describe_dirty",
        lambda: "test-describe-dirty",
    )
    # Minimal synthetic summary so we do not need a solver run here.
    from validation.lox_vent_cases import LoxVentSummary

    summary = LoxVentSummary(
        gravity_g=G3,
        timestep_s=0.02,
        rows=2,
        t_end_s=0.02,
        p_end_psia=39.9,
        h_initial_ft=INITIAL_HEIGHT_FT,
        h_final_ft=INITIAL_HEIGHT_FT + 1e-4,
        dh_ft=1e-4,
        dh_over_h0=1e-4 / INITIAL_HEIGHT_FT,
        final_vbl_vol_ft3=0.01,
        final_bl_thick_ft=0.001,
        max_ak3=0.05,
        conv_failed_total=0,
        ullage_closure_max_relative=0.01,
        finite=True,
        ullage_closure_within_tolerance=True,
        physical=True,
        failure_classifications=(),
        film_retention_lbm=0.0125,
        film_eos_inventory_lbm=0.0120,
        film_drift_lbm=0.0005,
        film_drift_over_inventory=0.0005 / 0.0120,
    )
    target = tmp_path / "lox_vent_manifest.json"
    payload = write_lox_vent_manifest({"G3": summary}, target)

    assert target.is_file()
    on_disk = json.loads(target.read_text(encoding="utf-8"))
    assert on_disk == payload
    assert payload["schema"] == "liqlev.validation.lox_vent"
    assert payload["solver_describe"] == "test-describe-dirty"
    assert payload["geometry_npz_sha256"] == sha256_file(GEOMETRY_NPZ_PATH)
    assert payload["fluid_step_sha256"] == sha256_file(FLUID_STEP_PATH)
    assert payload["harness_module_sha256"] == sha256_file(LOX_VENT_MODULE_PATH)
    assert payload["harness_module"] == "validation/lox_vent_cases.py"
    assert "python" in payload["versions"]
    assert "numpy" in payload["versions"]
    assert "numba" in payload["versions"]
    assert payload["case_definition"]["fill_fraction"] == FILL_FRACTION
    assert payload["case_definition"]["vent_rate_lbm_s"] == VENT_RATE_LBM_S
    assert payload["results"]["G3"]["conv_failed_total"] == 0
    assert payload["results"]["G3"]["ullage_closure_max_relative"] == pytest.approx(
        0.01
    )
