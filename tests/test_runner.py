"""Tests for headless runner and config I/O paths."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from liqlev.geometry.fixtures import cylinder_kernel
from liqlev.geometry.package import save_geometry_package
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
from liqlev.runner.single import run_single_case
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


def custom_geometry_config(geometry_path: str) -> SimulationConfig:
    return SimulationConfig(
        fluid=FluidConfig(
            name="Hydrogen",
            initial_pressure_psia=19.5,
            final_pressure_psia=10.0,
            initial_temperature_r=38.3,
        ),
        tank=TankConfig(
            diameter_ft=25.0,
            height_ft=99.0,
            fill_fractions=(0.5,),
            geometry_path=geometry_path,
        ),
        vent=VentProfileConfig(
            rates_lbm_s=(1.0e-7,),
            ramp_duration_s=10_000.0,
        ),
        gravity=GravityProfileConfig(mode="Constant", constant_g=0.001),
        epsilon=EpsilonConfig(mode="height_dep"),
        run=RunControls(duration_s=10_000.0, timestep_s=5_000.0),
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


def test_geometry_path_round_trips_in_schema_v2(tmp_path) -> None:
    config = SimulationConfig(
        tank=TankConfig(
            diameter_ft=4.0,
            height_ft=8.0,
            fill_fractions=(0.5,),
            geometry_path="tank.npz",
        )
    )
    path = tmp_path / "config.json"

    save_simulation_config(config, path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    loaded = load_simulation_config(path)

    assert payload["schema_version"] == 2
    assert payload["tank"]["geometry_path"] == "tank.npz"
    assert loaded.tank.geometry_path == "tank.npz"


def test_schema_v1_loads_with_empty_geometry_path(tmp_path) -> None:
    path = tmp_path / "legacy-config.json"
    path.write_text(json.dumps({"schema_version": 1}), encoding="utf-8")

    loaded = load_simulation_config(path)

    assert loaded.schema_version == 1
    assert loaded.tank.geometry_path == ""


def test_single_case_uses_exact_custom_geometry_inputs_and_height(tmp_path) -> None:
    kernel = cylinder_kernel(4.0, 8.0)
    package_path = tmp_path / "cylinder.npz"
    save_geometry_package(kernel, package_path)

    result = run_single_case(custom_geometry_config(str(package_path)))

    assert result.inputs["GeometryMode"] == 1
    assert result.inputs["FillFraction"] == pytest.approx(0.5)
    geometry_keys = {
        key
        for key in result.inputs
        if key == "GeometryMode" or key.startswith("Geom")
    }
    assert geometry_keys == {
        "GeometryMode",
        "GeomHeight",
        "GeomVolume",
        "GeomVolumeCoefficients",
        "GeomAreaSamples",
        "GeomPerimeter",
        "GeomSidewallArea",
        "GeomSidewallCoefficients",
    }
    expected_arrays = {
        "GeomHeight": kernel.height_ft,
        "GeomVolume": kernel.volume_ft3,
        "GeomVolumeCoefficients": kernel.volume_coefficients,
        "GeomAreaSamples": kernel.section_area_ft2,
        "GeomPerimeter": kernel.perimeter_ft,
        "GeomSidewallArea": kernel.sidewall_area_ft2,
        "GeomSidewallCoefficients": kernel.sidewall_coefficients,
    }
    for key, expected in expected_arrays.items():
        np.testing.assert_array_equal(result.inputs[key], expected)
    assert result.inputs["Htzero"] == pytest.approx(4.0)
    assert result.inputs["Volt"] == pytest.approx(kernel.total_volume_ft3)
    assert result.htank_ft == pytest.approx(8.0)


def test_sweep_uses_geometry_package_total_height(tmp_path) -> None:
    package_path = tmp_path / "cylinder.npz"
    save_geometry_package(cylinder_kernel(4.0, 8.0), package_path)

    result = run_sweep(custom_geometry_config(str(package_path)))

    scenario = result.scenarios["Fill 50%, eps=height_dep"]
    assert scenario["htank"] == pytest.approx(8.0)
