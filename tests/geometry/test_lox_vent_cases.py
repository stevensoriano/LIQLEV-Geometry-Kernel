"""Bounded tests for the LOX 43 L vent acceptance-case machinery."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from liqlev.geometry.package import load_geometry_package
from liqlev.geometry.jit import invert_monotone_volume
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
    TIMESTEP_S,
    VENT_RATE_LBM_S,
    build_lox_vent_config,
    run_lox_vent_case,
    sha256_file,
    ullage_closure_metric,
    write_lox_vent_manifest,
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
