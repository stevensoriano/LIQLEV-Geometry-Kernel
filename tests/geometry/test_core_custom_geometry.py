from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from core import _solver_loop, liqlev_simulation
from liqlev.geometry.fixtures import cylinder_kernel
from tests.geometry.test_boundary_layer import report_cylinder_profile


ROOT = Path(__file__).resolve().parents[2]
TINIT_R = 38.3
TIMESTEP_S = 5000.0
RUN_DURATION_S = 2.0 * TIMESTEP_S
GEOMETRY_DIAMETER_FT = 4.0
TANK_HEIGHT_FT = 8.0


def _hydrogen_density(t_rankine: float) -> tuple[float, float]:
    liquid = (
        0.1709
        + 0.7454 * t_rankine
        - 0.04421 * t_rankine**2
        + 0.001248 * t_rankine**3
        - 1.738e-5 * t_rankine**4
        + 9.424e-8 * t_rankine**5
    )
    vapor = (
        -0.2511
        + 0.04294 * t_rankine
        - 0.00286 * t_rankine**2
        + 9.159e-5 * t_rankine**3
        - 1.422e-6 * t_rankine**4
        + 1.001e-8 * t_rankine**5
    )
    return liquid, vapor


def _custom_cylinder_inputs(
    fill: float,
    *,
    legacy_dtank_ft: float = 1.25,
    gravity_g: float = 0.001,
    liquid_mass_factor: float = 1.0,
) -> dict[str, object]:
    kernel = cylinder_kernel(
        GEOMETRY_DIAMETER_FT,
        TANK_HEIGHT_FT,
        node_count=1025,
    )
    liquid_density, _ = _hydrogen_density(TINIT_R)
    liquid_mass = (
        liquid_density
        * fill
        * kernel.total_volume_ft3
        * liquid_mass_factor
    )
    time = np.array([0.0, RUN_DURATION_S], dtype=np.float64)
    ones = np.ones(2, dtype=np.float64)
    return {
        "Liquid": "Hydrogen",
        "Units": "British",
        "Delta": TIMESTEP_S,
        "Dtank": legacy_dtank_ft,
        "Htzero": 0.33,
        "Volt": kernel.total_volume_ft3,
        "Xmlzro": liquid_mass,
        "Pinit": 19.5,
        "Pfinal": 10.0,
        "Tinit": TINIT_R,
        "Neps": 0,
        "Tvmdot": time,
        "Xvmdot": np.full(2, 1.0e-7, dtype=np.float64),
        "Teps": None,
        "Xeps": None,
        "Tspal": time,
        "Xspacl": ones,
        "Tspav": time,
        "Xspacv": ones,
        "Tggo": time,
        "Xggo": np.full(2, gravity_g * 32.174, dtype=np.float64),
        "GeometryMode": 1,
        "FillFraction": fill,
        "GeomHeight": kernel.height_ft,
        "GeomVolume": kernel.volume_ft3,
        "GeomVolumeCoefficients": kernel.volume_coefficients,
        "GeomAreaSamples": kernel.section_area_ft2,
        "GeomPerimeter": kernel.perimeter_ft,
        "GeomSidewallArea": kernel.sidewall_area_ft2,
        "GeomSidewallCoefficients": kernel.sidewall_coefficients,
    }


def _run_custom_cylinder(fill: float, **overrides):
    inputs = _custom_cylinder_inputs(fill, **overrides)
    result = liqlev_simulation(inputs, verbose=False)
    assert len(result) == 1
    return result.iloc[0]


def test_legacy_physics_baseline_checker_passes_at_strict_tolerances() -> None:
    result = subprocess.run(
        [sys.executable, "scripts/check_physics_baseline.py"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "Physics baseline check passed." in result.stdout


def test_custom_cylinder_uses_volume_and_wetted_geometry_not_legacy_dtank() -> None:
    fill = 0.5
    narrow_legacy = _run_custom_cylinder(fill, legacy_dtank_ft=1.25)
    wide_legacy = _run_custom_cylinder(fill, legacy_dtank_ft=25.0)
    initial_height = fill * TANK_HEIGHT_FT
    area = np.pi * GEOMETRY_DIAMETER_FT**2 / 4.0
    perimeter = np.pi * GEOMETRY_DIAMETER_FT
    expected_eps = perimeter * initial_height / (
        perimeter * initial_height + area
    )
    liquid_density, _ = _hydrogen_density(TINIT_R)

    for row in (narrow_legacy, wide_legacy):
        expected_occupied_volume = (
            row["Liq Mass"] / liquid_density + row["VBL vol"]
        )
        expected_height = expected_occupied_volume / area
        assert row["Conv Failed"] == 0.0
        assert row["eps"] == pytest.approx(expected_eps, rel=1e-12)
        assert row["Height"] == pytest.approx(expected_height, rel=1e-11)
        assert row["dh/dt"] == pytest.approx(
            (expected_height - initial_height) / TIMESTEP_S,
            rel=1e-11,
            abs=1e-15,
        )

    geometry_columns = [
        "Height",
        "dh/dt",
        "eps",
        "VBL vol",
        "BL thick",
        "BL Vap Out",
    ]
    np.testing.assert_allclose(
        narrow_legacy[geometry_columns].to_numpy(dtype=float),
        wide_legacy[geometry_columns].to_numpy(dtype=float),
        rtol=1e-12,
        atol=1e-14,
    )


def test_custom_cylinder_does_not_require_a_nonzero_legacy_dtank() -> None:
    zero_legacy = _run_custom_cylinder(0.5, legacy_dtank_ft=0.0)
    positive_legacy = _run_custom_cylinder(0.5, legacy_dtank_ft=1.25)
    geometry_columns = [
        "Height",
        "dh/dt",
        "eps",
        "VBL vol",
        "BL thick",
        "BL Vap Out",
    ]

    assert zero_legacy["Conv Failed"] == 0.0
    np.testing.assert_allclose(
        zero_legacy[geometry_columns].to_numpy(dtype=float),
        positive_legacy[geometry_columns].to_numpy(dtype=float),
        rtol=1e-12,
        atol=1e-14,
    )


@pytest.mark.parametrize("fill", [0.1, 0.25, 0.5, 0.8, 0.95])
def test_custom_cylinder_boundary_layer_matches_published_report(
    fill: float,
) -> None:
    row = _run_custom_cylinder(fill)
    initial_height = fill * TANK_HEIGHT_FT
    expected_delta, expected_vbl, expected_normalized_exit = (
        report_cylinder_profile(
            GEOMETRY_DIAMETER_FT,
            float(row["AK3"]),
            initial_height,
        )
    )
    _, vapor_density = _hydrogen_density(TINIT_R)
    actual_normalized_exit = row["BL Vap Out"] / (
        row["AK1"] * TIMESTEP_S * vapor_density
    )

    assert row["Conv Failed"] == 0.0
    np.testing.assert_allclose(
        [row["BL thick"], row["VBL vol"], actual_normalized_exit],
        [expected_delta, expected_vbl, expected_normalized_exit],
        rtol=1e-3,
        atol=1e-10,
    )


def test_custom_geometry_requires_matching_total_volume() -> None:
    inputs = _custom_cylinder_inputs(0.5)
    inputs["Volt"] = float(inputs["Volt"]) * (1.0 + 2.0e-10)

    with pytest.raises(ValueError, match="GeomVolume"):
        liqlev_simulation(inputs, verbose=False)


def test_custom_geometry_rejects_non_float64_arrays_before_jit() -> None:
    inputs = _custom_cylinder_inputs(0.5)
    inputs["GeomPerimeter"] = np.asarray(
        inputs["GeomPerimeter"],
        dtype=np.float32,
    )

    with pytest.raises(ValueError, match="GeomPerimeter.*float64"):
        liqlev_simulation(inputs, verbose=False)


def test_boundary_layer_failure_writes_one_bounded_diagnostic_row() -> None:
    fill = 0.5
    result = liqlev_simulation(
        _custom_cylinder_inputs(fill, gravity_g=0.0),
        verbose=False,
    )

    assert len(result) == 1
    assert result.loc[0, "Conv Failed"] == 1.0
    assert result.loc[0, "Height"] == pytest.approx(fill * TANK_HEIGHT_FT)
    assert 0.0 <= result.loc[0, "Height"] <= TANK_HEIGHT_FT


def test_volume_inversion_failure_writes_one_bounded_diagnostic_row() -> None:
    fill = 0.5
    result = liqlev_simulation(
        _custom_cylinder_inputs(fill, liquid_mass_factor=2.0),
        verbose=False,
    )

    assert len(result) == 1
    assert result.loc[0, "Conv Failed"] == 1.0
    assert result.loc[0, "Height"] == pytest.approx(fill * TANK_HEIGHT_FT)
    assert 0.0 <= result.loc[0, "Height"] <= TANK_HEIGHT_FT


def test_custom_solver_path_compiles_in_nopython_mode() -> None:
    row = _run_custom_cylinder(0.25)

    assert row["Conv Failed"] == 0.0
    assert _solver_loop.nopython_signatures
    assert all(
        "pyobject" not in str(signature).lower()
        for signature in _solver_loop.nopython_signatures
    )
