"""Bounded tests for the config-driven vent-study re-run interface.

Solver-running tests use 1 s durations (the pattern established by
``tests/geometry/test_lox_vent_cases.py``); everything else is pure config,
formatting or manifest work with monkeypatched provenance.
"""

from __future__ import annotations

import json
import platform
from pathlib import Path

import numpy as np
import pytest

from liqlev.analysis import vent_study
from liqlev.analysis.vent_study import (
    SpaceflightSpec,
    VentStudyConfigError,
    derive_gravity_levels,
    format_report,
    load_study_config,
    read_assumptions_block,
    run_study,
    select_rows,
    simulation_config_for,
    study_config_from_mapping,
    write_manifest,
)
from liqlev.geometry.package import load_geometry_package
from validation import lox_vent_cases as lox


REPO_ROOT = vent_study.ROOT
SHIPPED_CONFIG = REPO_ROOT / "configs" / "lox_43L_40to35_ssm3.json"
ASSUMPTIONS_DOC = REPO_ROOT / "docs" / "lox-vent-test-definition.md"

SSM3_SPACEFLIGHT = {
    "thrust_N": 1.60,
    "vehicle_mass_kg": 280.0,
    "on_ms": 120.0,
    "period_ms": 3400.0,
}


def base_payload(**overrides: object) -> dict[str, object]:
    """A minimal valid config; keyword arguments replace top-level keys."""

    payload: dict[str, object] = {
        "schema": "liqlev.analysis.vent_study",
        "version": 1,
        "name": "unit_case",
        "fluid": "Oxygen",
        "geometry_package": "geometry/tables/nhq01-m21a-0201_LIQLEV_GEOMETRY.npz",
        "fill": {"fraction": lox.FILL_FRACTION},
        "pressure": {"initial_psia": 40.0, "final_psia": 35.0},
        "vent_rate": {"lbm_per_s": lox.VENT_RATE_LBM_S},
        "duration_s": 1.0,
        "timestep_s": 0.02,
        "gravity": {"levels": {"G3": lox.G3}},
        "dt_plateau": {"enabled": False},
    }
    payload.update(overrides)
    return payload


def synthetic_result(config: vent_study.VentStudyConfig | None = None):
    """A one-row result built without running the solver."""

    config = config or study_config_from_mapping(base_payload())
    summary = lox.LoxVentSummary(
        gravity_g=lox.G3,
        timestep_s=0.02,
        rows=2,
        t_end_s=0.04,
        p_end_psia=39.9,
        h_initial_ft=lox.INITIAL_HEIGHT_FT,
        h_final_ft=lox.INITIAL_HEIGHT_FT + 1e-4,
        dh_ft=1e-4,
        dh_over_h0=1e-4 / lox.INITIAL_HEIGHT_FT,
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
    row = vent_study.VentStudyRow(
        name="G3",
        gravity_g=lox.G3,
        rise_mm=1e-4 * vent_study.FT_TO_MM,
        summary=summary,
    )
    return vent_study.VentStudyResult(
        config=config,
        rows=(row,),
        dt_plateau=None,
        assumptions=read_assumptions_block(),
    )


# --------------------------------------------------------------------------
# Config validation
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("mutation", "match"),
    [
        ({"fill": {"liters": 43.0, "fraction": 0.44}}, "fill requires exactly one"),
        ({"fill": {}}, "fill requires exactly one"),
        (
            {"vent_rate": {"g_per_s": 11.89, "lbm_per_s": 0.026212963}},
            "vent_rate requires exactly one",
        ),
        ({"vent_rate": {}}, "vent_rate requires exactly one"),
        (
            {
                "gravity": {
                    "levels": {"G3": lox.G3},
                    "spaceflight": dict(SSM3_SPACEFLIGHT),
                }
            },
            "gravity requires exactly one",
        ),
        ({"gravity": {}}, "gravity requires exactly one"),
        ({"gravity": {"levels": {}}}, "at least one row"),
        ({"gravity": {"levels": {"G3": -1.0}}}, "must be non-negative"),
        (
            {"pressure": {"initial_psia": 35.0, "final_psia": 40.0}},
            "initial_psia must exceed",
        ),
        ({"timestep_s": 2.0}, "must be smaller than duration_s"),
        ({"duration_s": 0.0}, "duration_s must be positive"),
        ({"timestep_s": "fast"}, "timestep_s must be a number"),
        ({"fluid": "Helium"}, "fluid must be one of"),
        ({"fill": {"fraction": 1.5}}, "outside"),
        ({"dt_plateau": {"row": "G9"}}, "not a configured gravity level"),
        ({"dt_plateau": {"enabled": "yes"}}, "must be true or false"),
        ({"schema": "liqlev.analysis.other"}, "schema must be"),
        ({"version": 7}, "unsupported config version"),
        ({"geometry_package": "geometry/tables/missing.npz"}, "not found"),
        ({"bogus_key": 1}, "unknown key"),
        ({"fill": {"litres": 43.0}}, "unknown key"),
    ],
)
def test_invalid_configs_are_rejected(mutation: dict, match: str) -> None:
    with pytest.raises(VentStudyConfigError, match=match):
        study_config_from_mapping(base_payload(**mutation))


@pytest.mark.parametrize(
    ("block", "match"),
    [
        (
            {**SSM3_SPACEFLIGHT, "duty_cycle": 0.035294},
            "exactly one of 'duty_cycle'",
        ),
        ({"thrust_N": 1.6, "vehicle_mass_kg": 280.0}, "exactly one of 'duty_cycle'"),
        ({"thrust_N": 1.6, "on_ms": 120.0, "period_ms": 3400.0}, "vehicle_mass_kg"),
        (
            {"thrust_N": 1.6, "vehicle_mass_kg": 280.0, "on_ms": 120.0},
            "both 'on_ms' and 'period_ms'",
        ),
        (
            {
                "thrust_N": 1.6,
                "vehicle_mass_kg": 280.0,
                "on_ms": 3400.0,
                "period_ms": 120.0,
            },
            "must not exceed period_ms",
        ),
        (
            {"thrust_N": 1.6, "vehicle_mass_kg": 280.0, "duty_cycle": 1.5},
            "outside",
        ),
    ],
)
def test_invalid_spaceflight_blocks_are_rejected(block: dict, match: str) -> None:
    with pytest.raises(VentStudyConfigError, match=match):
        study_config_from_mapping(base_payload(gravity={"spaceflight": block}))


def test_unknown_gravity_row_selection_is_rejected() -> None:
    config = study_config_from_mapping(base_payload())
    with pytest.raises(VentStudyConfigError, match="unknown gravity row"):
        select_rows(config, ["G9"])


# --------------------------------------------------------------------------
# Derived quantities
# --------------------------------------------------------------------------


def test_ssm3_gravity_derivation_is_exact() -> None:
    """280 kg / 1.60 N / 120:3400 duty reproduces the pinned G1 and G3."""

    spec = SpaceflightSpec(
        thrust_n=1.60,
        vehicle_mass_kg=280.0,
        duty_cycle=120.0 / 3400.0,
        on_ms=120.0,
        period_ms=3400.0,
    )
    levels = derive_gravity_levels(spec)

    assert tuple(levels) == ("G0", "G1", "G3")
    assert levels["G0"] == 0.0
    # Exact at the precision the case definition pins (7 significant figures).
    assert float(f"{levels['G3']:.6e}") == lox.G3 == 5.826950e-04
    assert float(f"{levels['G1']:.6e}") == lox.G1 == 2.056571e-05
    assert levels["G3"] == pytest.approx(lox.G3, rel=1e-6)
    assert levels["G1"] == pytest.approx(lox.G1, rel=1e-6)
    assert spec.peak_g == 1.60 / (280.0 * 9.80665)
    assert spec.duty_average_g == spec.peak_g * (120.0 / 3400.0)


def test_spaceflight_block_derives_the_gravity_matrix() -> None:
    config = study_config_from_mapping(
        base_payload(gravity={"spaceflight": dict(SSM3_SPACEFLIGHT)})
    )

    assert config.gravity_source == "spaceflight"
    assert config.row_names() == ("G0", "G1", "G3")
    assert config.dt_plateau_row == "G3"
    assert config.spaceflight is not None
    assert config.spaceflight.duty_cycle == pytest.approx(0.03529411764705882)
    assert float(f"{config.gravity_matrix['G3']:.6e}") == lox.G3
    assert float(f"{config.gravity_matrix['G1']:.6e}") == lox.G1

    duty_form = study_config_from_mapping(
        base_payload(
            gravity={
                "spaceflight": {
                    "thrust_N": 1.60,
                    "vehicle_mass_kg": 280.0,
                    "duty_cycle": 120.0 / 3400.0,
                }
            }
        )
    )
    assert duty_form.gravity_matrix == config.gravity_matrix


def test_fill_liters_resolve_against_the_geometry_package() -> None:
    config = study_config_from_mapping(base_payload(fill={"liters": 43.0}))
    geometry = load_geometry_package(lox.GEOMETRY_NPZ_PATH)
    expected = 43.0 / (geometry.total_volume_ft3 * vent_study.LITERS_PER_FT3)

    assert config.fill_liters == 43.0
    assert config.fill_fraction == expected
    assert config.fill_fraction == pytest.approx(lox.FILL_FRACTION, rel=1e-5)


def test_vent_rate_in_grams_matches_the_pinned_lbm_rate() -> None:
    config = study_config_from_mapping(base_payload(vent_rate={"g_per_s": 11.89}))

    assert config.vent_rate_lbm_s == 11.89 / vent_study.GRAMS_PER_LBM
    assert config.vent_rate_lbm_s == pytest.approx(lox.VENT_RATE_LBM_S, rel=1e-8)
    assert config.vent_rate_g_s == pytest.approx(11.89, rel=1e-12)


# --------------------------------------------------------------------------
# Shipped preset
# --------------------------------------------------------------------------


def test_shipped_config_reproduces_the_production_case() -> None:
    config = load_study_config(SHIPPED_CONFIG)

    assert config.fluid == lox.FLUID
    assert config.fill_fraction == lox.FILL_FRACTION
    assert config.initial_pressure_psia == lox.PINIT_PSIA
    assert config.final_pressure_psia == lox.PFINAL_PSIA
    assert config.vent_rate_lbm_s == lox.VENT_RATE_LBM_S
    assert config.duration_s == lox.DURATION_S
    assert config.timestep_s == lox.TIMESTEP_S
    assert config.geometry_path == lox.GEOMETRY_NPZ_PATH
    assert config.gravity_matrix == lox.GRAVITY_MATRIX
    assert config.row_names() == ("G0", "G1", "G2", "G3", "G4")
    assert config.dt_plateau_enabled is True
    assert config.dt_plateau_row == "G3"

    # The solver config for every row must be byte-identical to the committed
    # case builder — the preset re-runs the production case, it does not
    # redefine it.
    for name, gravity_g in config.gravity_levels:
        assert simulation_config_for(config, gravity_g) == lox.build_lox_vent_config(
            gravity_g
        ), name


def test_overrides_narrow_rows_and_keep_the_plateau_row_valid() -> None:
    config = load_study_config(SHIPPED_CONFIG)

    only_g3 = select_rows(config, ["G3"])
    assert only_g3.row_names() == ("G3",)
    assert only_g3.dt_plateau_row == "G3"

    low_rows = select_rows(config, ["G0", "G1"])
    assert low_rows.row_names() == ("G0", "G1")
    assert low_rows.dt_plateau_row == "G1"  # falls back to the highest level


# --------------------------------------------------------------------------
# Study execution (bounded)
# --------------------------------------------------------------------------


def test_run_study_smoke_two_rows_short_duration() -> None:
    config = study_config_from_mapping(
        base_payload(gravity={"levels": {"G3": lox.G3, "G4": lox.G4}})
    )
    result = run_study(config)

    assert [row.name for row in result.rows] == ["G3", "G4"]
    assert result.dt_plateau is None
    for row in result.rows:
        summary = row.summary
        assert summary.rows >= 2
        assert summary.finite
        assert summary.conv_failed_total == 0
        assert summary.h_initial_ft == pytest.approx(lox.INITIAL_HEIGHT_FT, abs=5e-4)
        assert summary.t_end_s == pytest.approx(1.0, abs=summary.timestep_s * 1.5)
        assert lox.PFINAL_PSIA <= summary.p_end_psia <= lox.PINIT_PSIA + 0.5
        assert row.rise_mm == pytest.approx(summary.dh_ft * 304.8)
        # Reused ullage/physicality machinery is wired through, not re-derived.
        assert np.isfinite(summary.ullage_closure_max_relative)
        assert summary.ullage_closure_max_relative >= 0.0
        assert summary.ullage_closure_within_tolerance == (
            summary.ullage_closure_max_relative
            <= lox.ULLAGE_CLOSURE_RELATIVE_TOLERANCE
        )
        assert isinstance(summary.failure_classifications, tuple)
        assert summary.physical == (not summary.failure_classifications)
    # Higher gravity settles harder: less swelling.
    assert result.rows[0].rise_mm > result.rows[1].rise_mm
    assert "Steady-state gravity" in result.assumptions


def test_dt_plateau_check_reruns_at_half_and_quarter_timestep() -> None:
    config = study_config_from_mapping(
        base_payload(dt_plateau={"enabled": True, "row": "G3"})
    )
    result = run_study(config)
    plateau = result.dt_plateau

    assert plateau is not None
    assert plateau.row == "G3"
    assert plateau.gravity_g == lox.G3
    assert [point.timestep_s for point in plateau.points] == pytest.approx(
        [0.02, 0.01, 0.005]
    )
    assert plateau.points[0].delta_percent == 0.0
    assert plateau.base_rise_mm == pytest.approx(result.rows[0].rise_mm)
    assert all(np.isfinite(point.rise_mm) for point in plateau.points)
    assert all(np.isfinite(point.delta_percent) for point in plateau.points)
    assert plateau.max_abs_delta_percent == max(
        abs(point.delta_percent) for point in plateau.points
    )

    report = format_report(result)
    assert "dt-plateau check" in report
    assert "max |delta| vs base" in report


# --------------------------------------------------------------------------
# Report and assumptions block
# --------------------------------------------------------------------------


def test_report_quotes_the_eight_assumptions_verbatim() -> None:
    block = read_assumptions_block()
    document = ASSUMPTIONS_DOC.read_text(encoding="utf-8")

    assert block in document  # verbatim slice, not a paraphrase
    assert block.startswith("## 5. Assumptions and limitations")
    for index in range(1, 9):
        assert f"{index}. **" in block
    for phrase in (
        "Steady-state gravity",
        "No boundary-layer time constant",
        "g = 0 is a bound, not a prediction",
        "Saturated initial state",
        "Constant average vent rate",
        "Heritage ullage formulation",
        "Numerical, not experimental",
        "No internal hardware",
    ):
        assert phrase in block

    report = format_report(synthetic_result())
    assert block in report
    assert "docs/lox-vent-test-definition.md" in report


def test_report_table_carries_the_reported_columns() -> None:
    report = format_report(synthetic_result())

    assert "LIQLEV vent study — unit_case" in report
    assert "rise (mm)" in report
    assert "ullage" in report
    assert "Ullage-closure acceptance gate: 5%" in report
    assert "dt-plateau check: disabled for this run." in report


# --------------------------------------------------------------------------
# Manifest provenance
# --------------------------------------------------------------------------


def test_manifest_write_is_hash_bound(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(vent_study, "_git_worktree_is_dirty", lambda: False)
    monkeypatch.setattr(vent_study, "_git_describe", lambda: "test-describe")
    result = synthetic_result()
    target = tmp_path / "study_manifest.json"

    payload = write_manifest(result, target)
    on_disk = json.loads(target.read_text(encoding="utf-8"))

    assert on_disk == payload
    assert payload["schema"] == "liqlev.analysis.vent_study"
    assert payload["version"] == 1
    assert payload["solver_describe"] == "test-describe"
    assert payload["geometry_npz"] == (
        "geometry/tables/nhq01-m21a-0201_LIQLEV_GEOMETRY.npz"
    )
    assert payload["geometry_npz_sha256"] == lox.sha256_file(lox.GEOMETRY_NPZ_PATH)
    assert payload["harness_module"] == "validation/lox_vent_cases.py"
    assert payload["harness_module_sha256"] == lox.sha256_file(
        lox.LOX_VENT_MODULE_PATH
    )
    assert payload["analysis_module"] == "liqlev/analysis/vent_study.py"
    assert payload["analysis_module_sha256"] == lox.sha256_file(
        vent_study.ANALYSIS_MODULE_PATH
    )
    assert payload["assumptions_doc"] == "docs/lox-vent-test-definition.md"
    assert payload["assumptions_doc_sha256"] == lox.sha256_file(
        vent_study.ASSUMPTIONS_DOC_PATH
    )
    assert payload["fluid_step_sha256"] == lox.sha256_file(lox.FLUID_STEP_PATH)
    assert payload["versions"]["python"] == platform.python_version()
    assert payload["versions"]["numpy"] == np.__version__
    assert payload["case_definition"]["fill_fraction"] == lox.FILL_FRACTION
    assert payload["case_definition"]["vent_rate_lbm_s"] == lox.VENT_RATE_LBM_S
    assert payload["case_definition"]["gravity_matrix_g"] == {"G3": lox.G3}
    assert payload["case_definition"]["epsilon"] == "height_dep"
    assert payload["dt_plateau"] is None
    row = payload["results"]["G3"]
    assert row["gravity_row"] == "G3"
    assert row["conv_failed_total"] == 0
    assert row["ullage_closure_relative_tolerance"] == (
        lox.ULLAGE_CLOSURE_RELATIVE_TOLERANCE
    )
    assert row["rise_mm"] == pytest.approx(1e-4 * 304.8)


def test_manifest_refuses_a_dirty_worktree(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(vent_study, "_git_worktree_is_dirty", lambda: True)
    target = tmp_path / "study_manifest.json"

    with pytest.raises(RuntimeError, match="dirty git worktree"):
        write_manifest(synthetic_result(), target)
    assert not target.exists()


def test_dirty_probe_delegates_to_the_validation_guard(monkeypatch) -> None:
    """The F8 guard is reused, not reimplemented."""

    monkeypatch.setattr(lox, "_git_worktree_is_dirty", lambda: True)
    assert vent_study._git_worktree_is_dirty() is True
    monkeypatch.setattr(lox, "_git_describe_dirty", lambda: "delegated")
    assert vent_study._git_describe() == "delegated"


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def test_cli_runs_a_bounded_study_and_writes_outputs(
    monkeypatch, tmp_path: Path, capsys
) -> None:
    monkeypatch.setattr(vent_study, "_git_worktree_is_dirty", lambda: False)
    monkeypatch.setattr(vent_study, "_git_describe", lambda: "test-describe")

    exit_code = vent_study.main(
        [
            "--config",
            str(SHIPPED_CONFIG),
            "--output-dir",
            str(tmp_path),
            "--gravity",
            "G3",
            "--duration-s",
            "1.0",
            "--skip-dt-plateau",
        ]
    )
    captured = capsys.readouterr()

    assert exit_code == 0
    assert "LIQLEV vent study" in captured.out
    assert "## 5. Assumptions and limitations" in captured.out
    assert "Steady-state gravity" in captured.out

    manifest_path = tmp_path / "lox_43L_40to35_ssm3_manifest.json"
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert list(payload["results"]) == ["G3"]
    assert payload["solver_describe"] == "test-describe"
    assert payload["case_definition"]["duration_s"] == 1.0
    assert payload["case_definition"]["gravity_matrix_g"] == {"G3": lox.G3}
    assert payload["dt_plateau"] is None
    assert payload["config_source"] == "configs/lox_43L_40to35_ssm3.json"
    assert payload["config_sha256"] == lox.sha256_file(SHIPPED_CONFIG)
    assert (tmp_path / "lox_43L_40to35_ssm3_report.txt").is_file()


def test_cli_reports_config_errors_without_traceback(tmp_path: Path, capsys) -> None:
    bad = tmp_path / "bad.json"
    bad.write_text(
        json.dumps(base_payload(fill={"liters": 43.0, "fraction": 0.44})),
        encoding="utf-8",
    )

    exit_code = vent_study.main(["--config", str(bad), "--no-manifest"])
    captured = capsys.readouterr()

    assert exit_code == 2
    assert "config error" in captured.err
    assert "exactly one" in captured.err
