from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True)
def _find_interval(x: float, nodes: np.ndarray) -> int:
    if x <= nodes[0]:
        return 0
    if x >= nodes[-1]:
        return len(nodes) - 2
    low = 0
    high = len(nodes) - 1
    while low < high - 1:
        middle = (low + high) >> 1
        if nodes[middle] <= x:
            low = middle
        else:
            high = middle
    return low


@njit(cache=True)
def eval_ppoly(x: float, breaks: np.ndarray, coefficients: np.ndarray) -> float:
    index = _find_interval(x, breaks)
    bounded = min(max(x, breaks[0]), breaks[-1])
    dx = bounded - breaks[index]
    return (
        (coefficients[0, index] * dx + coefficients[1, index]) * dx
        + coefficients[2, index]
    ) * dx + coefficients[3, index]


@njit(cache=True)
def eval_ppoly_derivative(
    x: float, breaks: np.ndarray, coefficients: np.ndarray
) -> float:
    index = _find_interval(x, breaks)
    bounded = min(max(x, breaks[0]), breaks[-1])
    dx = bounded - breaks[index]
    return (
        3.0 * coefficients[0, index] * dx + 2.0 * coefficients[1, index]
    ) * dx + coefficients[2, index]


@njit(cache=True)
def interp_linear_nonnegative(
    x: float, nodes: np.ndarray, values: np.ndarray
) -> float:
    index = _find_interval(x, nodes)
    bounded = min(max(x, nodes[0]), nodes[-1])
    width = nodes[index + 1] - nodes[index]
    fraction = (bounded - nodes[index]) / width
    value = values[index] + fraction * (values[index + 1] - values[index])
    return max(0.0, value)


@njit(cache=True)
def invert_monotone_volume(
    target_volume: float,
    height: np.ndarray,
    volume: np.ndarray,
    volume_coefficients: np.ndarray,
) -> float:
    if target_volume < volume[0] or target_volume > volume[-1]:
        return np.nan
    if target_volume == volume[0]:
        return height[0]
    if target_volume == volume[-1]:
        return height[-1]

    low_index = 0
    high_index = len(volume) - 1
    while low_index < high_index - 1:
        middle = (low_index + high_index) >> 1
        if volume[middle] <= target_volume:
            low_index = middle
        else:
            high_index = middle

    if target_volume == volume[low_index]:
        return height[low_index]

    lower = height[low_index]
    upper = height[low_index + 1]
    volume_width = volume[low_index + 1] - volume[low_index]
    if volume_width <= 0.0:
        return np.nan
    guess = lower + (
        (target_volume - volume[low_index])
        / volume_width
        * (upper - lower)
    )
    tolerance = 1e-12 * max(1.0, height[-1])
    volume_ulp = np.spacing(max(abs(target_volume), 1.0))
    result = guess

    for _ in range(64):
        residual = eval_ppoly(guess, height, volume_coefficients) - target_volume
        if residual <= 0.0:
            lower = guess
        else:
            upper = guess
        derivative = eval_ppoly_derivative(guess, height, volume_coefficients)
        candidate = guess - residual / derivative if derivative > 0.0 else np.nan
        if not np.isfinite(candidate) or candidate <= lower or candidate >= upper:
            candidate = 0.5 * (lower + upper)
        guess = candidate
        if upper - lower <= tolerance:
            result = 0.5 * (lower + upper)
            result_residual = abs(
                eval_ppoly(result, height, volume_coefficients) - target_volume
            )
            lower_residual = abs(
                eval_ppoly(lower, height, volume_coefficients) - target_volume
            )
            upper_residual = abs(
                eval_ppoly(upper, height, volume_coefficients) - target_volume
            )
            if lower_residual < result_residual:
                result = lower
                result_residual = lower_residual
            if upper_residual < result_residual:
                result = upper
                result_residual = upper_residual
            local_derivative = abs(
                eval_ppoly_derivative(result, height, volume_coefficients)
            )
            height_ulp = np.spacing(max(abs(result), 1.0))
            allowed_volume_error = max(
                volume_ulp, 2.0 * local_derivative * height_ulp
            )
            if result_residual <= allowed_volume_error:
                return result
    return result


@njit(cache=True)
def _boundary_layer_derivatives(
    h: float,
    q: float,
    height: np.ndarray,
    volume_coefficients: np.ndarray,
    perimeter_values: np.ndarray,
    ak3: float,
) -> tuple[float, float]:
    perimeter = interp_linear_nonnegative(h, height, perimeter_values)
    if perimeter <= 0.0:
        if q <= 0.0:
            delta = 0.0
        else:
            return np.nan, np.nan
    else:
        delta = (1.5 * max(q, 0.0) / perimeter) ** (2.0 / 3.0)
    area = eval_ppoly_derivative(h, height, volume_coefficients)
    return ak3 * (area - perimeter * delta), perimeter * delta


@njit(cache=True)
def integrate_boundary_layer(
    ak3: float,
    top_height: float,
    height: np.ndarray,
    volume_coefficients: np.ndarray,
    perimeter_values: np.ndarray,
    substeps: int,
) -> tuple[float, float, float, int]:
    if (
        ak3 <= 0.0
        or not np.isfinite(ak3)
        or not np.isfinite(top_height)
        or top_height < height[0]
        or top_height > height[-1]
        or substeps <= 0
    ):
        return np.nan, np.nan, np.nan, 1

    q = 0.0
    vbl = 0.0
    interval_count = len(height) - 1
    for interval in range(interval_count):
        lower = height[interval]
        if lower >= top_height:
            break
        upper = min(height[interval + 1], top_height)
        step = (upper - lower) / substeps
        h = lower
        for _ in range(substeps):
            k1_q, k1_vbl = _boundary_layer_derivatives(
                h,
                q,
                height,
                volume_coefficients,
                perimeter_values,
                ak3,
            )
            k2_q, k2_vbl = _boundary_layer_derivatives(
                h + 0.5 * step,
                q + 0.5 * step * k1_q,
                height,
                volume_coefficients,
                perimeter_values,
                ak3,
            )
            k3_q, k3_vbl = _boundary_layer_derivatives(
                h + 0.5 * step,
                q + 0.5 * step * k2_q,
                height,
                volume_coefficients,
                perimeter_values,
                ak3,
            )
            k4_q, k4_vbl = _boundary_layer_derivatives(
                h + step,
                q + step * k3_q,
                height,
                volume_coefficients,
                perimeter_values,
                ak3,
            )
            next_q = (
                q
                + step * (k1_q + 2.0 * k2_q + 2.0 * k3_q + k4_q) / 6.0
            )
            next_vbl = (
                vbl
                + step
                * (k1_vbl + 2.0 * k2_vbl + 2.0 * k3_vbl + k4_vbl)
                / 6.0
            )
            if not np.isfinite(next_q) or not np.isfinite(next_vbl):
                return np.nan, np.nan, np.nan, 1
            q = max(0.0, next_q)
            vbl = max(0.0, next_vbl)
            h += step

    perimeter_top = interp_linear_nonnegative(
        top_height, height, perimeter_values
    )
    if perimeter_top <= 0.0:
        if q > 0.0:
            return np.nan, np.nan, np.nan, 1
        delta_top = 0.0
    else:
        delta_top = (1.5 * q / perimeter_top) ** (2.0 / 3.0)
    if not np.isfinite(delta_top):
        return np.nan, np.nan, np.nan, 1
    return delta_top, vbl, q, 0
