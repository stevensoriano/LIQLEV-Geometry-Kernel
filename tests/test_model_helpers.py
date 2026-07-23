"""Tests for extracted non-GUI model helpers."""

from __future__ import annotations

import numpy as np
import pytest

from liqlev.model.builder import G_TO_FT_S2, build_inputs, safe_eval_gravity
from liqlev.model.config import FluidConfig, RunControls, SimulationConfig, TankConfig
from liqlev.model.parsing import parse_numeric_array
from liqlev.model.units import convert_from_british, convert_to_british, display_unit
from liqlev.model.validation import InputValidationError, validate_simulation_config
from validation.physics_cases import build_case_inputs, epsilon_schedule, get_case


def test_parse_numeric_array_preserves_legacy_forms() -> None:
    assert parse_numeric_array("") == []
    assert parse_numeric_array("0.1, 0.3, 0.5") == [0.1, 0.3, 0.5]
    assert parse_numeric_array("0.1:0.1:0.3") == pytest.approx([0.1, 0.2, 0.3])
    assert parse_numeric_array("linspace(0.1, 0.5, 5)") == pytest.approx(
        [0.1, 0.2, 0.3, 0.4, 0.5]
    )


def test_unit_helpers_preserve_solver_british_contract() -> None:
    assert display_unit("pressure", si_mode=False) == "psia"
    assert display_unit("pressure", si_mode=True) == "bar"
    assert convert_from_british(10.0, "mass_flow", si_mode=False) == 10.0
    assert convert_from_british(10.0, "mass_flow", si_mode=True) == pytest.approx(
        4.53592
    )
    assert convert_to_british(4.53592, "mass_flow", si_mode=True) == pytest.approx(10.0)


def test_build_inputs_matches_validation_case_setup_for_as203() -> None:
    case = get_case("as203_default_high_vent")
    expected = build_case_inputs(case)
    neps, teps, xeps = epsilon_schedule(
        case.epsilon_mode, case.duration_s, case.epsilon_value
    )
    actual = build_inputs(
        fluid=case.fluid,
        pinit_psia=case.pinit_psia,
        pfinal_psia=case.pfinal_psia,
        dtank=case.dtank_ft,
        htank=case.htank_ft,
        fill_fraction=case.fill_fraction,
        duration=case.duration_s,
        delta_t=case.delta_t_s,
        vent_rate=case.vent_rate_lbm_s,
        neps=neps,
        teps=teps,
        xeps=xeps,
        ramp_duration=case.ramp_duration_s or case.duration_s,
        ramp_target_factor=case.ramp_target_factor,
        nggo=2,
        tggo=np.array([0.0, case.duration_s]),
        xggo=np.array([case.gravity_g * G_TO_FT_S2, case.gravity_g * G_TO_FT_S2]),
        xmlzro_override=case.xmlzro_override_lbm,
        tinit_override=case.tinit_override_r,
    )

    for key, value in expected.items():
        if key == "Title":
            continue
        if isinstance(value, np.ndarray):
            assert np.asarray(actual[key]) == pytest.approx(value)
        elif isinstance(value, float):
            assert actual[key] == pytest.approx(value)
        else:
            assert actual[key] == value


def test_safe_gravity_eval_and_validation_errors_are_field_level() -> None:
    assert safe_eval_gravity("0.001 + sin(t)", 0.0) == pytest.approx(0.001)

    bad_config = SimulationConfig(
        fluid=FluidConfig(initial_pressure_psia=10.0, final_pressure_psia=12.0),
        tank=TankConfig(diameter_ft=-1.0, height_ft=28.18, fill_fractions=(1.2,)),
        run=RunControls(duration_s=10.0, timestep_s=10.0),
    )

    with pytest.raises(InputValidationError) as exc_info:
        validate_simulation_config(bad_config)

    fields = {issue.field for issue in exc_info.value.issues}
    assert "fluid.initial_pressure_psia" in fields
    assert "tank.diameter_ft" in fields
    assert "tank.fill_fractions[0]" in fields
    assert "run.timestep_s" in fields
