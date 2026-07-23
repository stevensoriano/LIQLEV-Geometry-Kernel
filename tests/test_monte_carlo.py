"""Tests for headless Monte Carlo geometry integration."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import liqlev.model.validation as validation_module
import liqlev.runner.monte_carlo as monte_carlo_module
from liqlev.geometry.fixtures import cylinder_kernel
from liqlev.geometry.package import save_geometry_package
from liqlev.model.config import (
    EpsilonConfig,
    FluidConfig,
    GravityProfileConfig,
    RunControls,
    SimulationConfig,
    TankConfig,
    VentProfileConfig,
)
from liqlev.runner.monte_carlo import MonteCarloRequest


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


def test_monte_carlo_reuses_one_custom_geometry_for_every_sample(
    tmp_path, monkeypatch
) -> None:
    kernel = cylinder_kernel(4.0, 8.0, node_count=17)
    package_path = tmp_path / "cylinder.npz"
    save_geometry_package(kernel, package_path)
    captured_inputs: list[dict[str, object]] = []
    load_count = 0
    real_load = validation_module.load_geometry_package

    def counting_load(path):
        nonlocal load_count
        load_count += 1
        return real_load(path)

    def capture_solver(inputs, **_kwargs):
        captured_inputs.append(inputs)
        return pd.DataFrame({"Hratio": [float(inputs["FillFraction"])]})

    monkeypatch.setattr(validation_module, "load_geometry_package", counting_load)
    monkeypatch.setattr(monte_carlo_module, "liqlev_simulation", capture_solver)
    request = MonteCarloRequest(
        n=3,
        vent_min_lbm_s=1.0e-7,
        vent_max_lbm_s=2.0e-7,
        fill_min=0.2,
        fill_max=0.8,
        gravity_min_g=0.0005,
        gravity_max_g=0.0015,
        seed=42,
    )

    result = monte_carlo_module.run_monte_carlo(
        custom_geometry_config(str(package_path)),
        request,
    )

    assert load_count == 1
    assert len(captured_inputs) == request.n
    expected_arrays = {
        "GeomHeight": kernel.height_ft,
        "GeomVolume": kernel.volume_ft3,
        "GeomVolumeCoefficients": kernel.volume_coefficients,
        "GeomAreaSamples": kernel.section_area_ft2,
        "GeomPerimeter": kernel.perimeter_ft,
        "GeomSidewallArea": kernel.sidewall_area_ft2,
        "GeomSidewallCoefficients": kernel.sidewall_coefficients,
    }
    first_inputs = captured_inputs[0]
    for inputs, params in zip(captured_inputs, result.all_params, strict=True):
        assert inputs["GeometryMode"] == 1
        assert inputs["Volt"] == pytest.approx(kernel.total_volume_ft3)
        assert inputs["FillFraction"] == pytest.approx(params["fill"])
        assert np.asarray(inputs["GeomHeight"])[-1] == pytest.approx(
            kernel.total_height_ft
        )
        for key, expected in expected_arrays.items():
            assert inputs[key] is first_inputs[key]
            np.testing.assert_array_equal(inputs[key], expected)
