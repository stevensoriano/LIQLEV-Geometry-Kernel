"""Display-independent parsing helpers for LIQLEV inputs."""

from __future__ import annotations

import re

import numpy as np


def parse_numeric_array(text: str) -> list[float]:
    """Parse comma lists, inclusive ranges, or linspace expressions."""
    text = text.strip()
    if not text:
        return []

    match = re.match(r"linspace\(\s*([^,]+),\s*([^,]+),\s*([^)]+)\)", text, re.I)
    if match:
        return np.linspace(
            float(match.group(1)),
            float(match.group(2)),
            int(match.group(3)),
        ).tolist()

    parts = text.split(":")
    if len(parts) == 3:
        start, step, stop = float(parts[0]), float(parts[1]), float(parts[2])
        return np.arange(start, stop + step * 0.5, step).tolist()

    return [float(item.strip()) for item in text.split(",") if item.strip()]
