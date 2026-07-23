"""Tests for headless runner and config I/O paths."""

from __future__ import annotations

import pandas as pd

from liqlev.io.config_json import load_simulation_config, save_simulation_config
from liqlev.model.config import (
    EpsilonConfig,
    FluidConfig,
    GravityProfileConfig,
    RunControls,
    SimulationConfig,
    TankConfig,
    VentProfileConfig,
)
from liqlev.runner.sweep import run_sweep
from validation.physics_cases import get_case, run_case


def config_from_case(name: str) -> SimulationConfig:
    case = get_case(name)
    return SimulationConfig(
        fluid=FluidConfig(
            name=case.fluid,
            initial_pressure_psia=case.pinit_psia,
            final_pressure_psia=case.pfinal_psia,
            initial_mass_lbm=case.xmlzro_override_lbm,
            initial_temperature_r=case.tinit_override_r,
        ),
        tank=TankConfig(
            diameter_ft=case.dtank_ft,
            height_ft=case.htank_ft,
            fill_fractions=(case.fill_fraction,),
        ),
        vent=VentProfileConfig(
            rates_lbm_s=(case.vent_rate_lbm_s,),
            ramp_duration_s=case.ramp_duration_s or case.duration_s,
            ramp_target_factor=case.ramp_target_factor,
        ),
        gravity=GravityProfileConfig(mode="Constant", constant_g=case.gravity_g),
        epsilon=(
            EpsilonConfig(mode="Custom", values=(case.epsilon_value,))
            if case.epsilon_mode == "custom"
            else EpsilonConfig(mode=case.epsilon_mode)
        ),
        run=RunControls(duration_s=case.duration_s, timestep_s=case.delta_t_s),
    )


def test_run_sweep_matches_validation_case_dataframe() -> None:
    config = config_from_case("as203_default_high_vent")
    expected = run_case(get_case("as203_default_high_vent"))

    result = run_sweep(config)

    assert result.run_count == 1
    scenario = result.scenarios["Fill 51%, eps=AS-203 Schedule"]
    actual = scenario["dfs"][0]
    pd.testing.assert_frame_equal(actual, expected)


def test_run_sweep_returns_legacy_compatible_scenario_shape() -> None:
    config = config_from_case("hydrogen_height_dep_mid_fill")
    result = run_sweep(config)

    assert list(result.scenarios) == ["Fill 50%, eps=height_dep"]
    scenario = result.scenarios["Fill 50%, eps=height_dep"]
    assert set(scenario) == {"dfs", "vent_rates", "fill", "eps_label", "htank"}
    assert scenario["vent_rates"] == [0.0015]
    assert scenario["fill"] == 0.5
    assert scenario["eps_label"] == "height_dep"


def test_simulation_config_json_round_trip(tmp_path) -> None:
    config = config_from_case("as203_default_high_vent")
    path = tmp_path / "as203_config.json"

    save_simulation_config(config, path)
    loaded = load_simulation_config(path)

    assert loaded == config
