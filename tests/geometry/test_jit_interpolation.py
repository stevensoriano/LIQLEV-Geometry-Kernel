from __future__ import annotations

import numpy as np
import pytest

from liqlev.geometry.coefficients import pchip_coefficients
from liqlev.geometry.jit import (
    eval_ppoly,
    eval_ppoly_derivative,
    interp_linear_nonnegative,
    invert_monotone_volume,
)


def test_linear_volume_is_exact_and_invertible() -> None:
    height = np.linspace(0.0, 4.0, 9)
    volume = 3.25 * height
    coefficients = pchip_coefficients(height, volume)
    for h in np.linspace(0.0, 4.0, 33):
        assert eval_ppoly(h, height, coefficients) == pytest.approx(3.25 * h)
        assert eval_ppoly_derivative(h, height, coefficients) == pytest.approx(3.25)
        assert invert_monotone_volume(
            3.25 * h, height, volume, coefficients
        ) == pytest.approx(h, abs=1e-11)


def test_inverse_rejects_out_of_domain_volume() -> None:
    height = np.array([0.0, 1.0, 2.0])
    volume = height**3
    coefficients = pchip_coefficients(height, volume)
    assert np.isnan(invert_monotone_volume(-1.0, height, volume, coefficients))
    assert np.isnan(invert_monotone_volume(9.0, height, volume, coefficients))


def test_perimeter_interpolation_is_nonnegative() -> None:
    nodes = np.array([0.0, 1.0, 2.0])
    values = np.array([0.0, 2.0, 0.0])
    assert interp_linear_nonnegative(0.5, nodes, values) == pytest.approx(1.0)
    assert interp_linear_nonnegative(3.0, nodes, values) == pytest.approx(0.0)


@pytest.mark.parametrize(
    ("height", "volume", "target", "expected_height"),
    [
        (
            np.array([0.0, 1.0, 2.0, 3.0]),
            np.array([0.0, 0.0, 1.0, 2.0]),
            0.0,
            0.0,
        ),
        (
            np.array([0.0, 1.0, 2.0, 3.0]),
            np.array([0.0, 1.0, 1.0, 2.0]),
            1.0,
            2.0,
        ),
        (
            np.array([0.0, 1.0, 2.0, 3.0]),
            np.array([0.0, 1.0, 2.0, 2.0]),
            2.0,
            3.0,
        ),
    ],
)
def test_inverse_returns_defined_height_for_duplicate_volume_target(
    height: np.ndarray,
    volume: np.ndarray,
    target: float,
    expected_height: float,
) -> None:
    coefficients = pchip_coefficients(height, volume)
    assert invert_monotone_volume(
        target, height, volume, coefficients
    ) == pytest.approx(expected_height)


def test_inverse_resolves_height_near_interior_plateau() -> None:
    height = np.array([0.0, 1.0, 2.0, 3.0])
    volume = np.array([0.0, 1.0, 1.0, 2.0])
    coefficients = pchip_coefficients(height, volume)
    expected_height = 2.000001
    target = eval_ppoly(expected_height, height, coefficients)

    assert invert_monotone_volume(
        target, height, volume, coefficients
    ) == pytest.approx(expected_height, abs=1e-12 * max(1.0, height[-1]))


def test_runtime_functions_compile_without_object_mode() -> None:
    height = np.array([0.0, 1.0, 2.0])
    volume = height**2
    coefficients = pchip_coefficients(height, volume)
    invert_monotone_volume(1.0, height, volume, coefficients)
    eval_ppoly_derivative(1.0, height, coefficients)
    interp_linear_nonnegative(1.0, height, volume)
    assert invert_monotone_volume.nopython_signatures
    assert eval_ppoly.nopython_signatures
    assert eval_ppoly_derivative.nopython_signatures
    assert interp_linear_nonnegative.nopython_signatures
