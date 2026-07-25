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


@pytest.mark.parametrize("ak3", [-0.015, -1.0, -1e-12])
def test_boundary_layer_rejects_negative_ak3(ak3: float) -> None:
    """AUTHORIZED EDIT #1 (plan 2.7c): reject only ak3 < 0; ak3 == 0 is exact."""
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


def test_boundary_layer_ak3_zero_is_exact_zero_solution() -> None:
    """ak3 == 0 => dq/dh = 0 => q = delta = V_BL = 0 with status 0 (plan 2.3/2.7c)."""
    kernel = cylinder_kernel(4.0, 8.0, node_count=17)
    delta, vbl, q, status = integrate_boundary_layer(
        0.0,
        4.0,
        kernel.height_ft,
        kernel.volume_coefficients,
        kernel.perimeter_ft,
        4,
    )
    assert status == 0
    assert delta == 0.0
    assert vbl == 0.0
    assert q == 0.0


def _lox_probe_geometry():
    """Committed NASA tank at the 43 L fill used in the F2 findings table."""
    from pathlib import Path

    from liqlev.geometry.jit import (
        eval_ppoly,
        eval_ppoly_derivative,
        interp_linear_nonnegative,
        invert_monotone_volume,
    )
    from liqlev.geometry.package import load_geometry_package

    geom = load_geometry_package(
        Path(__file__).resolve().parents[2]
        / "geometry"
        / "tables"
        / "nhq01-m21a-0201_LIQLEV_GEOMETRY.npz"
    )
    fill_ft3 = 43.0 / 28.316846592
    top_height = float(
        invert_monotone_volume(
            fill_ft3,
            geom.height_ft,
            geom.volume_ft3,
            geom.volume_coefficients,
        )
    )
    area = float(
        eval_ppoly_derivative(
            top_height, geom.height_ft, geom.volume_coefficients
        )
    )
    perimeter = float(
        interp_linear_nonnegative(
            top_height, geom.height_ft, geom.perimeter_ft
        )
    )
    volume = float(
        eval_ppoly(top_height, geom.height_ft, geom.volume_coefficients)
    )
    return geom, top_height, area, perimeter, volume


def test_boundary_layer_saturation_lock_in_at_ak3_1e2() -> None:
    """Plan 2.7a: AK3 = 1e2 locks δ → A/P (0.999245) with V_BL ≤ V(top)."""
    from liqlev.geometry.jit import eval_ppoly

    geom, top_height, area, perimeter, volume = _lox_probe_geometry()
    delta, vbl, q, status = integrate_boundary_layer(
        1.0e2,
        top_height,
        geom.height_ft,
        geom.volume_coefficients,
        geom.perimeter_ft,
        4,
    )
    assert status == 0
    assert area > 0.0 and perimeter > 0.0
    assert delta / (area / perimeter) >= 0.999
    # Lock the measured pre-fix saturation ratio (findings table).
    np.testing.assert_allclose(
        delta / (area / perimeter), 0.999245, rtol=0.0, atol=5e-6
    )
    assert vbl <= volume * (1.0 + 1e-6)
    assert vbl / volume >= 0.99


@pytest.mark.parametrize("ak3", [1.0e3, 1.0e4, 1.0e6])
def test_boundary_layer_rejects_extreme_ak3_garbage_class(ak3: float) -> None:
    """Plan 2.7b: documented 3.2×/15×/323× silent-garbage class must not return success-with-garbage.

    Pre-fix: status 0 with V_BL ≫ V. Post-fix: status 1 (self-consistency)
    or status 0 with physical saturation bounds (δ ≤ A/P, V_BL ≤ V).
    """
    geom, top_height, area, perimeter, volume = _lox_probe_geometry()
    delta, vbl, q, status = integrate_boundary_layer(
        ak3,
        top_height,
        geom.height_ft,
        geom.volume_coefficients,
        geom.perimeter_ft,
        4,
    )
    if status == 1:
        assert np.isnan(delta) and np.isnan(vbl) and np.isnan(q)
        return
    assert status == 0
    assert np.isfinite(delta) and np.isfinite(vbl)
    assert delta <= (area / perimeter) * (1.0 + 1e-6)
    assert vbl <= volume * (1.0 + 1e-6)
    # Extreme AK3 must be at the saturation limit, not a partial/garbage state.
    assert delta / (area / perimeter) >= 0.99
    assert vbl / volume >= 0.99


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


def test_boundary_layer_rejects_nonfinite_rk4_candidate() -> None:
    kernel = cylinder_kernel(4.0, 8.0, node_count=17)
    volume_coefficients = kernel.volume_coefficients.copy()
    volume_coefficients[2, :] = np.nan

    result = integrate_boundary_layer(
        0.015,
        4.0,
        kernel.height_ft,
        volume_coefficients,
        kernel.perimeter_ft,
        4,
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
