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
    PULSED_ROW,
    SpaceflightSpec,
    VentStudyConfigError,
    load_study_config,
    pulse_table_sha256,
    pulse_timing_notes,
    pulsed_gravity_table,
    study_config_from_mapping,
)
from liqlev.model.builder import G_TO_FT_S2
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
