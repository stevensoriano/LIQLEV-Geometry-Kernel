"""Bounded tests for the pulsed-gravity (GP) row of the vent-study harness.

Same discipline as ``tests/analysis/test_vent_study.py``: the square-wave
synthesis and the schema guards are pure work, and every solver-running test
uses a short duration.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from liqlev.analysis import vent_study
from liqlev.analysis.vent_study import (
    PULSED_GRAVITY_MODE,
    PULSED_ROW,
    PULSED_ROW_CAVEAT,
    SpaceflightSpec,
    VentStudyConfigError,
    format_report,
    load_study_config,
    prepare_pulsed_gravity,
    pulse_table_sha256,
    pulse_timing_notes,
    pulsed_gravity_table,
    read_assumptions_block,
    run_dt_plateau,
    run_pulsed_row,
    run_row,
    run_study,
    select_rows,
    study_config_from_mapping,
    write_manifest,
)
from liqlev.model.builder import G_TO_FT_S2
from liqlev.runner.single import PreparedGravity
from validation import lox_vent_cases as lox


SSM3_PULSED = {
    "thrust_N": 1.60,
    "vehicle_mass_kg": 280.0,
    "on_ms": 120.0,
    "period_ms": 3400.0,
    "pulsed_row": True,
}

PEAK_G = 1.60 / (280.0 * 9.80665)
PEAK_FT_S2 = PEAK_G * G_TO_FT_S2


def base_payload(**overrides: object) -> dict[str, object]:
    """A minimal valid config; keyword arguments replace top-level keys."""

    payload: dict[str, object] = {
        "schema": "liqlev.analysis.vent_study",
        "version": 1,
        "name": "pulsed_unit_case",
        "fluid": "Oxygen",
        "geometry_package": "geometry/tables/nhq01-m21a-0201_LIQLEV_GEOMETRY.npz",
        "fill": {"fraction": lox.FILL_FRACTION},
        "pressure": {"initial_psia": 40.0, "final_psia": 35.0},
        "vent_rate": {"lbm_per_s": lox.VENT_RATE_LBM_S},
        "duration_s": 1.0,
        "timestep_s": 0.02,
        "gravity": {"spaceflight": dict(SSM3_PULSED)},
        "dt_plateau": {"enabled": False},
    }
    payload.update(overrides)
    return payload


def spaceflight_payload(**block: object) -> dict[str, object]:
    """Base payload whose gravity block is the given spaceflight mapping."""

    return base_payload(gravity={"spaceflight": block})


def synthetic_summary(gravity_g: float) -> lox.LoxVentSummary:
    """A plausible summary built without running the solver."""

    return lox.LoxVentSummary(
        gravity_g=gravity_g,
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


def synthetic_pulsed_result(**overrides: object):
    """A two-row result (constant G3 plus GP) with no solver run."""

    config = study_config_from_mapping(base_payload(**overrides))
    spec = config.spaceflight
    _, descriptor = prepare_pulsed_gravity(config)
    constant = vent_study.VentStudyRow(
        name="G3",
        gravity_g=spec.peak_g,
        rise_mm=1e-4 * vent_study.FT_TO_MM,
        summary=synthetic_summary(spec.peak_g),
    )
    pulsed = vent_study.VentStudyRow(
        name=PULSED_ROW,
        gravity_g=spec.duty_average_g,
        rise_mm=2e-4 * vent_study.FT_TO_MM,
        summary=synthetic_summary(spec.duty_average_g),
        gravity_mode=PULSED_GRAVITY_MODE,
        pulse=descriptor,
    )
    return vent_study.VentStudyResult(
        config=config,
        rows=(constant, pulsed),
        dt_plateau=None,
        assumptions=read_assumptions_block(),
    )


# --------------------------------------------------------------------------
# Square-wave synthesis
# --------------------------------------------------------------------------


def test_pulsed_table_integrates_to_the_duty_average_over_whole_periods() -> None:
    """Trapezoid over N whole periods returns N * on_s * peak."""

    on_s, period_s, periods = 0.12, 3.4, 5
    duration_s = periods * period_s
    tggo, xggo = pulsed_gravity_table(PEAK_G, on_s, period_s, 0.0, duration_s)

    integral = float(np.trapezoid(xggo, tggo))
    assert integral == pytest.approx(periods * on_s * PEAK_FT_S2, rel=1e-6)

    # The time-mean of the synthesized train is the duty-averaged level.
    assert integral / duration_s == pytest.approx(
        PEAK_FT_S2 * (on_s / period_s), rel=1e-6
    )


def test_pulsed_table_nodes_only_take_the_two_square_wave_values() -> None:
    tggo, xggo = pulsed_gravity_table(PEAK_G, 0.12, 3.4, 0.0, 10.0)

    assert set(np.unique(xggo)) == {0.0, PEAK_FT_S2}
    assert xggo.dtype == np.float64
    assert xggo.flags["C_CONTIGUOUS"]


def test_pulsed_table_times_are_strictly_increasing_and_cover_the_duration() -> None:
    duration_s = 10.0
    tggo, xggo = pulsed_gravity_table(PEAK_G, 0.12, 3.4, 0.0, duration_s)

    assert tggo.dtype == np.float64
    assert tggo.flags["C_CONTIGUOUS"]
    assert tggo[0] == 0.0
    assert tggo[-1] == duration_s
    assert np.all(np.diff(tggo) > 0.0)
    assert len(tggo) == len(xggo)
    # First point is g(0) and the train starts OFF only when the phase delays it.
    assert xggo[0] == PEAK_FT_S2


def test_pulsed_table_edges_use_the_double_breakpoint_riser() -> None:
    riser_s = 1e-9
    tggo, xggo = pulsed_gravity_table(
        PEAK_G, 0.12, 3.4, 0.0, 3.4, riser_s=riser_s
    )

    # (0, peak) then the OFF edge as a double breakpoint, then the ON edge that
    # lands exactly on the duration boundary. 4 * periods + 1 nodes.
    assert len(tggo) == 5
    assert tggo[1] == pytest.approx(0.12, abs=1e-15)
    assert tggo[2] == pytest.approx(0.12 + riser_s, abs=1e-15)
    assert xggo[1] == PEAK_FT_S2
    assert xggo[2] == 0.0
    assert tggo[-1] == 3.4
    assert xggo[-1] == PEAK_FT_S2

    long_t, _ = pulsed_gravity_table(PEAK_G, 0.12, 3.4, 0.0, 4 * 3.4)
    assert len(long_t) == 4 * 4 + 1


def test_pulsed_table_phase_shifts_the_train() -> None:
    on_s, period_s, duration_s = 0.12, 3.4, 10.0
    phase_s = 1.25
    tggo, xggo = pulsed_gravity_table(PEAK_G, on_s, period_s, phase_s, duration_s)

    # Nothing fires before the phase offset; the first rising edge is at phase.
    assert xggo[0] == 0.0
    rising = tggo[np.flatnonzero(xggo == PEAK_FT_S2)[0]]
    assert rising == pytest.approx(phase_s + 1e-9, abs=1e-12)

    # A phase in (period - on, period) leaves a pulse already ON at t = 0.
    late_phase = period_s - on_s / 2.0
    _, late_x = pulsed_gravity_table(
        PEAK_G, on_s, period_s, late_phase, duration_s
    )
    assert late_x[0] == PEAK_FT_S2


def test_pulsed_table_at_full_duty_is_the_constant_table() -> None:
    """on == period collapses to the two-point constant table, bit for bit."""

    tggo, xggo = pulsed_gravity_table(PEAK_G, 3.4, 3.4, 0.0, 12.0)

    assert tggo.tolist() == [0.0, 12.0]
    assert xggo.tolist() == [PEAK_FT_S2, PEAK_FT_S2]
    assert xggo[0] == PEAK_G * G_TO_FT_S2


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"period_s": 0.0}, "period_s must be positive"),
        ({"on_s": 0.0}, "on_s must be positive"),
        ({"on_s": 4.0}, "on_s must not exceed period_s"),
        ({"phase_s": -0.1}, "phase_s must satisfy"),
        ({"phase_s": 3.4}, "phase_s must satisfy"),
        ({"duration_s": 0.0}, "duration_s must be positive"),
        ({"riser_s": 0.0}, "riser_s must be positive"),
        ({"riser_s": 5.0}, "riser_s must be smaller"),
        ({"peak_g": -1.0}, "peak_g must be non-negative"),
    ],
)
def test_pulsed_table_rejects_invalid_arguments(kwargs: dict, match: str) -> None:
    call = {
        "peak_g": PEAK_G,
        "on_s": 0.12,
        "period_s": 3.4,
        "phase_s": 0.0,
        "duration_s": 10.0,
    }
    call.update(kwargs)
    with pytest.raises(ValueError, match=match):
        pulsed_gravity_table(**call)


def test_pulse_table_sha256_binds_both_arrays() -> None:
    tggo, xggo = pulsed_gravity_table(PEAK_G, 0.12, 3.4, 0.0, 10.0)
    other_t, other_x = pulsed_gravity_table(PEAK_G, 0.12, 3.4, 0.5, 10.0)

    digest = pulse_table_sha256(tggo, xggo)
    assert digest == pulse_table_sha256(tggo, xggo)
    assert digest != pulse_table_sha256(other_t, other_x)
    assert len(digest) == 64


# --------------------------------------------------------------------------
# Schema: pulsed_row / phase_ms
# --------------------------------------------------------------------------


def test_pulsed_keys_default_off_and_stay_backward_compatible() -> None:
    config = study_config_from_mapping(
        spaceflight_payload(
            thrust_N=1.60, vehicle_mass_kg=280.0, on_ms=120.0, period_ms=3400.0
        )
    )
    spec = config.spaceflight

    assert spec is not None
    assert spec.pulsed_row is False
    assert spec.phase_ms == 0.0
    assert config.has_pulsed_row is False
    assert config.row_names() == ("G0", "G1", "G3")
    assert PULSED_ROW not in config.row_names()


def test_pulsed_row_does_not_change_the_gravity_matrix() -> None:
    pulsed = study_config_from_mapping(base_payload())
    duty_form = study_config_from_mapping(
        spaceflight_payload(
            thrust_N=1.60, vehicle_mass_kg=280.0, duty_cycle=120.0 / 3400.0
        )
    )

    assert pulsed.has_pulsed_row is True
    assert pulsed.gravity_matrix == duty_form.gravity_matrix
    assert pulsed.row_names() == ("G0", "G1", "G3")


def test_pulsed_spec_round_trips_through_a_config_file(tmp_path: Path) -> None:
    target = tmp_path / "pulsed.json"
    target.write_text(
        json.dumps(
            base_payload(
                gravity={"spaceflight": {**SSM3_PULSED, "phase_ms": 500.0}}
            )
        ),
        encoding="utf-8",
    )

    config = load_study_config(target)
    spec = config.spaceflight

    assert spec is not None
    assert spec.pulsed_row is True
    assert spec.phase_ms == 500.0
    assert spec.on_ms == 120.0
    assert spec.period_ms == 3400.0
    assert spec.on_s == pytest.approx(0.12)
    assert spec.period_s == pytest.approx(3.4)
    assert spec.phase_s == pytest.approx(0.5)

    payload = spec.as_dict()
    assert payload["pulsed_row"] is True
    assert payload["phase_ms"] == 500.0
    assert payload["peak_g"] == spec.peak_g
    assert payload["duty_average_g"] == spec.duty_average_g


@pytest.mark.parametrize(
    ("block", "match"),
    [
        (
            {
                "thrust_N": 1.6,
                "vehicle_mass_kg": 280.0,
                "duty_cycle": 120.0 / 3400.0,
                "pulsed_row": True,
            },
            "pulsed_row requires the 'on_ms'/'period_ms' form",
        ),
        (
            {**SSM3_PULSED, "pulsed_row": "yes"},
            "pulsed_row must be true or false",
        ),
        (
            {
                "thrust_N": 1.6,
                "vehicle_mass_kg": 280.0,
                "on_ms": 120.0,
                "period_ms": 3400.0,
                "phase_ms": 100.0,
            },
            "phase_ms requires pulsed_row",
        ),
        (
            {**SSM3_PULSED, "phase_ms": 3400.0},
            r"phase_ms must satisfy 0 <= phase_ms < period_ms",
        ),
        (
            {**SSM3_PULSED, "phase_ms": 9000.0},
            r"phase_ms must satisfy 0 <= phase_ms < period_ms",
        ),
        (
            {**SSM3_PULSED, "phase_ms": -1.0},
            "phase_ms must be non-negative",
        ),
        (
            {**SSM3_PULSED, "phase_ms": "later"},
            "phase_ms must be a number",
        ),
        (
            {**SSM3_PULSED, "spin_ms": 1.0},
            "unknown key",
        ),
    ],
)
def test_invalid_pulsed_blocks_are_rejected(block: dict, match: str) -> None:
    with pytest.raises(VentStudyConfigError, match=match):
        study_config_from_mapping(spaceflight_payload(**block))


# --------------------------------------------------------------------------
# Run-time timing guards
# --------------------------------------------------------------------------


def test_under_resolved_pulse_is_rejected_before_any_run() -> None:
    with pytest.raises(VentStudyConfigError, match="under-resolved"):
        study_config_from_mapping(base_payload(timestep_s=0.08))

    # Two steps per ON interval is the accepted floor, not an error.
    ok = study_config_from_mapping(base_payload(timestep_s=0.06))
    assert ok.timestep_s == 0.06


def test_timestep_override_re_checks_the_pulse_resolution(tmp_path: Path) -> None:
    target = tmp_path / "pulsed.json"
    target.write_text(json.dumps(base_payload()), encoding="utf-8")

    exit_code = vent_study.main(
        ["--config", str(target), "--timestep-s", "0.08", "--no-manifest"]
    )
    assert exit_code == 2


def test_aliasing_note_fires_only_when_the_timestep_does_not_divide() -> None:
    spec = SpaceflightSpec(
        thrust_n=1.60,
        vehicle_mass_kg=280.0,
        duty_cycle=120.0 / 3400.0,
        on_ms=120.0,
        period_ms=3400.0,
        pulsed_row=True,
    )

    # Production timesteps divide both 120 ms and 3400 ms exactly.
    for timestep_s in (0.02, 0.01, 0.005):
        assert pulse_timing_notes(spec, timestep_s) == ()

    notes = pulse_timing_notes(spec, 0.007)
    assert len(notes) == 2
    assert any("on_s" in note for note in notes)
    assert any("period_s" in note for note in notes)
    for note in notes:
        assert "aliasing" in note
        assert note.isascii()


# --------------------------------------------------------------------------
# Row selection and plateau wiring
# --------------------------------------------------------------------------


def test_select_rows_treats_the_pulsed_row_as_selectable() -> None:
    config = study_config_from_mapping(base_payload())

    both = select_rows(config, ["G3", PULSED_ROW])
    assert both.row_names() == ("G3",)
    assert both.has_pulsed_row is True

    constants_only = select_rows(config, ["G0", "G3"])
    assert constants_only.row_names() == ("G0", "G3")
    assert constants_only.has_pulsed_row is False

    pulsed_only = select_rows(config, [PULSED_ROW])
    assert pulsed_only.row_names() == ()
    assert pulsed_only.has_pulsed_row is True
    assert pulsed_only.dt_plateau_row == PULSED_ROW


def test_pulsed_row_is_not_selectable_without_pulsed_row_true() -> None:
    config = study_config_from_mapping(
        spaceflight_payload(
            thrust_N=1.60, vehicle_mass_kg=280.0, on_ms=120.0, period_ms=3400.0
        )
    )
    with pytest.raises(VentStudyConfigError, match="unknown gravity row"):
        select_rows(config, [PULSED_ROW])


def test_dt_plateau_row_may_name_the_pulsed_row() -> None:
    config = study_config_from_mapping(
        base_payload(dt_plateau={"enabled": True, "row": PULSED_ROW})
    )
    assert config.dt_plateau_row == PULSED_ROW

    with pytest.raises(VentStudyConfigError, match="not a configured gravity level"):
        study_config_from_mapping(
            spaceflight_payload(
                thrust_N=1.60,
                vehicle_mass_kg=280.0,
                on_ms=120.0,
                period_ms=3400.0,
            )
            | {"dt_plateau": {"enabled": True, "row": PULSED_ROW}}
        )


# --------------------------------------------------------------------------
# GP execution (bounded)
# --------------------------------------------------------------------------


def test_gp_row_feeds_the_synthesized_table_to_the_solver(monkeypatch) -> None:
    """The solver sees the square wave itself; result.inputs proves it."""

    config = study_config_from_mapping(base_payload())
    spec = config.spaceflight
    expected_t, expected_x = pulsed_gravity_table(
        spec.peak_g, 0.12, 3.4, 0.0, config.duration_s
    )

    captured: dict[str, object] = {}
    real_run = vent_study.run_single_case

    def spy(simulation, **kwargs):
        captured["prepared"] = kwargs.get("prepared_gravity")
        result = real_run(simulation, **kwargs)
        captured["result"] = result
        captured["simulation"] = simulation
        return result

    monkeypatch.setattr(vent_study, "run_single_case", spy)
    row = run_pulsed_row(config)

    prepared = captured["prepared"]
    assert isinstance(prepared, PreparedGravity)
    assert prepared.nggo == len(expected_t)
    assert prepared.gravity_function is None
    assert prepared.messages and "pulsed" in prepared.messages[0]

    # The builder echoes Tggo/Xggo into result.inputs: that echo is the
    # provenance truth for the gravity the solver actually integrated.
    inputs = captured["result"].inputs
    assert np.array_equal(inputs["Tggo"], expected_t)
    assert np.array_equal(inputs["Xggo"], expected_x)
    assert inputs["Nggo"] == len(expected_t)

    # The constant-g bookkeeping value is the true time mean, not the peak.
    assert captured["simulation"].gravity.constant_g == spec.duty_average_g

    assert row.name == PULSED_ROW
    assert row.gravity_mode == PULSED_GRAVITY_MODE
    assert row.gravity_g == spec.duty_average_g
    assert row.notes == ()
    assert row.pulse == {
        "thrust_N": 1.60,
        "vehicle_mass_kg": 280.0,
        "on_ms": 120.0,
        "period_ms": 3400.0,
        "phase_ms": 0.0,
        "peak_g": spec.peak_g,
        "mean_g": spec.duty_average_g,
        "n_table_points": len(expected_t),
        "table_sha256": pulse_table_sha256(expected_t, expected_x),
    }
    assert row.summary.finite
    assert row.summary.conv_failed_total == 0


def test_run_study_appends_the_pulsed_row_after_the_constant_rows() -> None:
    config = select_rows(study_config_from_mapping(base_payload()), ["G3", PULSED_ROW])
    result = run_study(config)

    assert [row.name for row in result.rows] == ["G3", PULSED_ROW]
    assert [row.gravity_mode for row in result.rows] == ["constant", "pulsed"]
    assert result.rows[1].pulse is not None
    assert result.dt_plateau is None


def test_dt_plateau_accepts_the_pulsed_row() -> None:
    config = select_rows(
        study_config_from_mapping(
            base_payload(dt_plateau={"enabled": True, "row": PULSED_ROW})
        ),
        [PULSED_ROW],
    )
    plateau = run_dt_plateau(config)

    assert plateau.row == PULSED_ROW
    assert plateau.gravity_g == config.spaceflight.duty_average_g
    assert [point.timestep_s for point in plateau.points] == pytest.approx(
        [0.02, 0.01, 0.005]
    )
    assert all(np.isfinite(point.rise_mm) for point in plateau.points)
    assert plateau.points[0].delta_percent == 0.0


def test_gp_rows_are_deterministic() -> None:
    config = study_config_from_mapping(base_payload())

    first = run_pulsed_row(config)
    second = run_pulsed_row(config)

    assert first.rise_mm == second.rise_mm
    assert first.summary == second.summary
    assert first.pulse == second.pulse


def test_full_duty_pulsed_row_reproduces_the_constant_peak_row() -> None:
    """Physics tripwire: a 100% duty square wave is the constant-g case.

    on_ms == period_ms makes the synthesized table constant at the peak, so the
    profile path must return the constant-g result with no arithmetic drift.
    """

    config = study_config_from_mapping(
        base_payload(
            duration_s=2.0,
            gravity={
                "spaceflight": {
                    "thrust_N": 1.60,
                    "vehicle_mass_kg": 280.0,
                    "on_ms": 120.0,
                    "period_ms": 120.0,
                    "pulsed_row": True,
                }
            },
        )
    )
    spec = config.spaceflight
    assert spec.duty_cycle == 1.0
    assert spec.duty_average_g == spec.peak_g

    constant = run_row(config, "G3", spec.peak_g)
    pulsed = run_pulsed_row(config)

    assert pulsed.pulse["n_table_points"] == 2
    assert pulsed.gravity_g == constant.gravity_g
    assert pulsed.rise_mm == constant.rise_mm  # bitwise
    assert pulsed.summary.dh_ft == constant.summary.dh_ft
    assert pulsed.summary.h_final_ft == pytest.approx(
        constant.summary.h_final_ft, rel=1e-15
    )
    assert pulsed.summary.max_ak3 == constant.summary.max_ak3


# --------------------------------------------------------------------------
# Report caveat and manifest
# --------------------------------------------------------------------------


def test_report_marks_the_pulsed_row_and_prints_the_mandatory_caveat() -> None:
    result = synthetic_pulsed_result()
    report = format_report(result)
    descriptor = result.rows[1].pulse

    assert PULSED_ROW_CAVEAT in report  # verbatim, not a paraphrase
    assert PULSED_ROW_CAVEAT.isascii()
    assert "—" not in PULSED_ROW_CAVEAT
    for phrase in (
        "SSM-3 square wave",
        "model-consistency experiment",
        "docs/lox-vent-test-definition.md",
        "section 3",
        "quasi-steady",
        "no relaxation timescale",
        "3.4 s cycle",
        "G1-G3",
        "not a validated prediction",
    ):
        assert phrase in PULSED_ROW_CAVEAT

    # The matrix table marks the mode and both g levels are visible.
    assert "mode" in report
    assert f"{PULSED_ROW:<6} {'pulsed':<9}" in report
    assert f"{'G3':<6} {'constant':<9}" in report
    assert f"peak g = {descriptor['peak_g']:.6e}" in report
    assert f"mean g = {descriptor['mean_g']:.6e}" in report
    assert descriptor["table_sha256"] in report
    assert f"{descriptor['n_table_points']} points" in report


def test_constant_only_reports_are_unchanged_by_the_feature() -> None:
    result = synthetic_pulsed_result()
    constant_only = vent_study.VentStudyResult(
        config=result.config,
        rows=(result.rows[0],),
        dt_plateau=None,
        assumptions=result.assumptions,
    )
    report = format_report(constant_only)

    assert PULSED_ROW_CAVEAT not in report
    assert "pulsed profile" not in report
    # No mode column: the header keeps its pre-feature layout.
    assert f"{'Row':<6} {'g (std)':>13}" in report


def test_manifest_carries_the_pulse_descriptor(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(vent_study, "_git_worktree_is_dirty", lambda: False)
    monkeypatch.setattr(vent_study, "_git_describe", lambda: "test-describe")
    result = synthetic_pulsed_result()
    target = tmp_path / "pulsed_manifest.json"

    payload = write_manifest(result, target)
    on_disk = json.loads(target.read_text(encoding="utf-8"))
    assert on_disk == payload

    constant_row = payload["results"]["G3"]
    assert constant_row["gravity_mode"] == "constant"
    assert "pulse" not in constant_row

    pulsed_row = payload["results"][PULSED_ROW]
    assert pulsed_row["gravity_mode"] == "pulsed"
    assert pulsed_row["pulse"] == result.rows[1].pulse
    assert pulsed_row["pulse"]["table_sha256"] == result.rows[1].pulse["table_sha256"]
    assert pulsed_row["gravity_g"] == result.config.spaceflight.duty_average_g

    spaceflight = payload["case_definition"]["spaceflight"]
    assert spaceflight["pulsed_row"] is True
    assert spaceflight["phase_ms"] == 0.0
    # The pulsed row never joins the constant-g matrix.
    assert list(payload["case_definition"]["gravity_matrix_g"]) == ["G0", "G1", "G3"]


def test_docs_record_the_pulsed_capability() -> None:
    vent_doc = (vent_study.ROOT / "docs" / "vent-study.md").read_text(
        encoding="utf-8"
    )
    definition = (
        vent_study.ROOT / "docs" / "lox-vent-test-definition.md"
    ).read_text(encoding="utf-8")

    assert "pulsed_row" in vent_doc
    assert "phase_ms" in vent_doc
    assert "GP" in vent_doc
    # §6 carries the dated additive note; §3 stays the method of record.
    deferred = definition.split("## 6. Deferred / not yet in scope")[1].split("## 7.")[
        0
    ]
    assert "2026-07-29" in deferred
    assert "feature/pulsed-gravity" in deferred
