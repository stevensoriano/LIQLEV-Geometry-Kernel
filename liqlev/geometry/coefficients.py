from __future__ import annotations

import numpy as np


def pchip_coefficients(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y) or len(x) < 3:
        raise ValueError("x and y must be equal-length 1D arrays with N >= 3")
    h = np.diff(x)
    if np.any(h <= 0.0):
        raise ValueError("x must be strictly increasing")
    delta = np.diff(y) / h
    slopes = np.zeros_like(y)
    for index in range(1, len(y) - 1):
        if delta[index - 1] * delta[index] > 0.0:
            w1 = 2.0 * h[index] + h[index - 1]
            w2 = h[index] + 2.0 * h[index - 1]
            slopes[index] = (w1 + w2) / (
                w1 / delta[index - 1] + w2 / delta[index]
            )
    slopes[0] = ((2.0 * h[0] + h[1]) * delta[0] - h[0] * delta[1]) / (
        h[0] + h[1]
    )
    slopes[-1] = (
        (2.0 * h[-1] + h[-2]) * delta[-1] - h[-1] * delta[-2]
    ) / (h[-1] + h[-2])
    if slopes[0] * delta[0] <= 0.0:
        slopes[0] = 0.0
    elif abs(slopes[0]) > 3.0 * abs(delta[0]):
        slopes[0] = 3.0 * delta[0]
    if slopes[-1] * delta[-1] <= 0.0:
        slopes[-1] = 0.0
    elif abs(slopes[-1]) > 3.0 * abs(delta[-1]):
        slopes[-1] = 3.0 * delta[-1]
    coefficients = np.empty((4, len(x) - 1), dtype=np.float64)
    coefficients[0] = (slopes[:-1] + slopes[1:] - 2.0 * delta) / h**2
    coefficients[1] = (3.0 * delta - 2.0 * slopes[:-1] - slopes[1:]) / h
    coefficients[2] = slopes[:-1]
    coefficients[3] = y[:-1]
    return np.ascontiguousarray(coefficients)
