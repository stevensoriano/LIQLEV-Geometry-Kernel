from __future__ import annotations

import json
import platform

import numba
import numpy as np
import pytest

from validation.custom_geometry_cases import (
    FILL_FRACTIONS,
    FLUID_STEP_PATH,
    GEOMETRY_NPZ_PATH,
    MANIFEST_PATH,
    REFINEMENT_RELATIVE_TOLERANCE,
    build_nasa_tank_config,
    run_nasa_tank_validation,
    sha256_file,
)
from validation.physics_cases import build_case_inputs, get_case


@pytest.fixture(scope="module")
def nasa_validation():
    return run_nasa_tank_validation()


def test_nasa_tank_config_uses_deterministic_legacy_hydrogen_case() -> None:
    legacy = get_case("hydrogen_height_dep_mid_fill")
    legacy_inputs = build_case_inputs(legacy)

    config = build_nasa_tank_config()

    assert config.tank.geometry_path == str(GEOMETRY_NPZ_PATH)
    assert config.tank.fill_fractions == FILL_FRACTIONS
    assert config.fluid.name == legacy.fluid
    assert config.fluid.initial_pressure_psia == legacy.pinit_psia
    assert config.fluid.final_pressure_psia == legacy.pfinal_psia
    assert config.fluid.initial_temperature_r == legacy_inputs["Tinit"]
    assert config.epsilon.mode == "height_dep"
    assert config.gravity.mode == "Constant"
    assert config.gravity.constant_g == 0.001
    assert config.run.timestep_s == 10.0
    assert config.run.duration_s == 300.0


@pytest.mark.parametrize("fill_fraction", FILL_FRACTIONS)
def test_nasa_tank_solver_is_finite_bounded_and_converged(
    nasa_validation,
    fill_fraction: float,
) -> None:
    dataframe = nasa_validation.production_results[fill_fraction].dataframe
    kernel = nasa_validation.production_kernel

    assert not dataframe.empty
    assert np.isfinite(dataframe.to_numpy(dtype=float)).all()
    assert (dataframe["Height"] >= 0.0).all()
    assert (dataframe["Height"] <= kernel.total_height_ft).all()
    assert (dataframe["VBL vol"] >= 0.0).all()
    assert (dataframe["BL thick"] >= 0.0).all()
    assert dataframe["Conv Failed"].sum() == 0.0


def test_solver_evaluation_grid_refinement_is_within_approved_limit(
    nasa_validation,
) -> None:
    assert nasa_validation.passed
    assert nasa_validation.maximum_convergence_count == 0
    assert nasa_validation.refinement_authority_node_count == 92
    assert nasa_validation.coarse_node_count == 513
    assert nasa_validation.reference_node_count == 1025

    for fill_metrics in nasa_validation.refinement_by_fill.values():
        assert (
            fill_metrics.max_height_relative_difference
            <= REFINEMENT_RELATIVE_TOLERANCE
        )
        assert (
            fill_metrics.max_boundary_layer_volume_relative_difference
            <= REFINEMENT_RELATIVE_TOLERANCE
        )


def test_nasa_tank_result_manifest_matches_current_evidence(
    nasa_validation,
) -> None:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))

    assert manifest["geometry_npz_sha256"] == sha256_file(GEOMETRY_NPZ_PATH)
    assert manifest["fluid_step_sha256"] == sha256_file(FLUID_STEP_PATH)
    assert manifest["versions"] == {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "numba": numba.__version__,
    }
    assert manifest["fill_cases"] == list(FILL_FRACTIONS)
    assert manifest["maximum_convergence_count"] == (
        nasa_validation.maximum_convergence_count
    )
    assert manifest["maximum_refinement_difference"] == pytest.approx(
        nasa_validation.maximum_refinement_difference,
        rel=0.0,
        abs=1e-12,
    )
    assert manifest["passed"] is nasa_validation.passed
    assert manifest["solver_evaluation_grid_refinement"][
        "continuous_authority"
    ] == "committed_92_node_cad_derived_kernel"
    assert manifest["solver_evaluation_grid_refinement"][
        "relative_tolerance"
    ] == REFINEMENT_RELATIVE_TOLERANCE
    recorded_by_fill = manifest["solver_evaluation_grid_refinement"][
        "per_fill_maxima"
    ]
    for fill_fraction, metrics in nasa_validation.refinement_by_fill.items():
        recorded = recorded_by_fill[f"{fill_fraction:.2f}"]
        assert recorded["maximum_height_relative_difference"] == pytest.approx(
            metrics.max_height_relative_difference,
            rel=0.0,
            abs=1e-12,
        )
        assert recorded[
            "maximum_boundary_layer_volume_relative_difference"
        ] == pytest.approx(
            metrics.max_boundary_layer_volume_relative_difference,
            rel=0.0,
            abs=1e-12,
        )
