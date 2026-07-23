from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from liqlev.geometry.fixtures import cylinder_kernel, sphere_kernel
from liqlev.geometry.jit import integrate_boundary_layer


ROOT = Path(__file__).resolve().parents[2]


def _report_height(delta: float, diameter_ft: float, ak3: float) -> float:
    term = delta**1.5 / (3.0 * diameter_ft)
    total = term
    ratio_base = 4.0 * delta / diameter_ft
    for n in range(2, 1000):
        term *= ratio_base * (2.0 * n - 1.0) / (2.0 * n + 1.0)
        total += term
        if abs(term) <= 1e-16 * abs(total):
            break
    return 8.0 * total / ak3


def report_cylinder_profile(
    diameter_ft: float, ak3: float, height_ft: float
) -> tuple[float, float, float]:
    delta = (3.0 * diameter_ft * ak3 * height_ft / 8.0) ** (2.0 / 3.0)
    for _ in range(40):
        residual = _report_height(delta, diameter_ft, ak3) - height_ft
        derivative = delta**0.5 / (ak3 * (diameter_ft / 4.0 - delta))
        correction = residual / derivative
        delta -= correction
        if abs(correction) <= 1e-14 * max(1.0, abs(delta)):
            break

    term = delta**2.5 / 5.0
    vbl_series = term
    ratio_base = 4.0 * delta / diameter_ft
    for n in range(2, 1000):
        term *= ratio_base * (2.0 * n + 1.0) / (2.0 * n + 3.0)
        vbl_series += term
        if abs(term) <= 1e-16 * abs(vbl_series):
            break
    vbl = 8.0 * np.pi * vbl_series / ak3
    normalized_exit = (2.0 / 3.0) * np.pi * diameter_ft * delta**1.5
    return float(delta), float(vbl), float(normalized_exit)


@pytest.mark.parametrize("fill", [0.1, 0.25, 0.5, 0.8, 0.95])
def test_numeric_cylinder_boundary_layer_matches_published_report(
    fill: float,
) -> None:
    diameter = 4.0
    tank_height = 8.0
    kernel = cylinder_kernel(diameter, tank_height, node_count=1025)
    top_height = fill * tank_height
    expected = report_cylinder_profile(diameter, 0.015, top_height)
    actual = integrate_boundary_layer(
        0.015,
        top_height,
        kernel.height_ft,
        kernel.volume_coefficients,
        kernel.perimeter_ft,
        4,
    )
    assert actual[3] == 0
    np.testing.assert_allclose(actual[:3], expected, rtol=1e-3, atol=1e-10)


def test_cylinder_fixture_has_exact_arrays_and_metadata() -> None:
    diameter = 4.0
    height = 2.0
    kernel = cylinder_kernel(diameter, height, node_count=3)
    h = np.array([0.0, 1.0, 2.0])
    area = np.pi * diameter**2 / 4.0
    perimeter = np.pi * diameter

    np.testing.assert_array_equal(kernel.height_ft, h)
    np.testing.assert_allclose(kernel.volume_ft3, area * h)
    np.testing.assert_allclose(kernel.section_area_ft2, area)
    np.testing.assert_allclose(kernel.perimeter_ft, perimeter)
    np.testing.assert_allclose(kernel.sidewall_area_ft2, perimeter * h)
    expected_total = area + perimeter * h
    expected_total[-1] += area
    np.testing.assert_allclose(kernel.total_wetted_area_ft2, expected_total)
    _assert_analytic_metadata(kernel.metadata, "analytic-cylinder", height)


def test_sphere_fixture_has_exact_arrays_and_metadata() -> None:
    radius = 2.0
    kernel = sphere_kernel(radius, node_count=5)
    h = np.linspace(0.0, 2.0 * radius, 5)
    radial_term = np.maximum(0.0, 2.0 * radius * h - h**2)
    area = np.pi * radial_term
    volume = np.pi * h**2 * (radius - h / 3.0)
    perimeter = 2.0 * np.pi * np.sqrt(radial_term)
    sidewall = 2.0 * np.pi * radius * h

    np.testing.assert_array_equal(kernel.height_ft, h)
    np.testing.assert_allclose(kernel.volume_ft3, volume)
    np.testing.assert_allclose(kernel.section_area_ft2, area)
    np.testing.assert_allclose(kernel.perimeter_ft, perimeter)
    np.testing.assert_allclose(kernel.sidewall_area_ft2, sidewall)
    np.testing.assert_allclose(kernel.total_wetted_area_ft2, sidewall)
    _assert_analytic_metadata(kernel.metadata, "analytic-sphere", 2.0 * radius)


def _assert_analytic_metadata(metadata, geometry_id: str, height_ft: float) -> None:
    assert metadata.schema_version == 1
    assert metadata.geometry_id == geometry_id
    assert metadata.source_step_sha256 == "0" * 64
    assert metadata.fluid_step_sha256 == "0" * 64
    assert metadata.axis == "+Y"
    assert metadata.gravity_direction == "-Y"
    assert metadata.length_unit == "ft"
    assert metadata.area_unit == "ft^2"
    assert metadata.volume_unit == "ft^3"
    assert metadata.y_min_mm == 0.0
    assert metadata.y_max_mm == pytest.approx(height_ft * 304.8)


@pytest.mark.parametrize("top_height", [-0.1, 8.1, np.nan])
def test_boundary_layer_rejects_invalid_height(top_height: float) -> None:
    kernel = cylinder_kernel(4.0, 8.0, node_count=17)
    result = integrate_boundary_layer(
        0.015,
        top_height,
        kernel.height_ft,
        kernel.volume_coefficients,
        kernel.perimeter_ft,
        4,
    )
    assert result[3] == 1
    assert np.isnan(result[:3]).all()


@pytest.mark.parametrize("ak3", [0.0, -0.015])
def test_boundary_layer_rejects_nonpositive_ak3(ak3: float) -> None:
    kernel = cylinder_kernel(4.0, 8.0, node_count=17)
    result = integrate_boundary_layer(
        ak3,
        4.0,
        kernel.height_ft,
        kernel.volume_coefficients,
        kernel.perimeter_ft,
        4,
    )
    assert result[3] == 1
    assert np.isnan(result[:3]).all()


@pytest.mark.parametrize("substeps", [0, -1])
def test_boundary_layer_rejects_invalid_substeps(substeps: int) -> None:
    kernel = cylinder_kernel(4.0, 8.0, node_count=17)
    result = integrate_boundary_layer(
        0.015,
        4.0,
        kernel.height_ft,
        kernel.volume_coefficients,
        kernel.perimeter_ft,
        substeps,
    )
    assert result[3] == 1
    assert np.isnan(result[:3]).all()


def test_boundary_layer_integrator_compiles_in_nopython_mode() -> None:
    kernel = cylinder_kernel(4.0, 8.0, node_count=17)
    result = integrate_boundary_layer(
        0.015,
        4.0,
        kernel.height_ft,
        kernel.volume_coefficients,
        kernel.perimeter_ft,
        4,
    )
    assert result[3] == 0
    assert integrate_boundary_layer.nopython_signatures


def test_legacy_physics_baseline_checker_passes_unchanged() -> None:
    result = subprocess.run(
        [sys.executable, "scripts/check_physics_baseline.py"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "Physics baseline check passed." in result.stdout
