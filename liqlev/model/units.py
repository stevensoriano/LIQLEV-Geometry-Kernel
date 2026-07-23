"""Unit conversion helpers for display and input normalization."""

from __future__ import annotations

from collections.abc import Iterable


UNIT_CONV = {
    "length_ft": ("ft", "m", 0.3048),
    "length_in": ("in", "cm", 2.54),
    "pressure": ("psia", "bar", 0.0689476),
    "mass_flow": ("lbm/s", "kg/s", 0.453592),
    "temperature": ("R", "K", 1 / 1.8),
    "area": ("ft\u00b2", "m\u00b2", 0.092903),
    "volume": ("ft\u00b3", "m\u00b3", 0.028317),
    "density": ("lbm/ft\u00b3", "kg/m\u00b3", 16.0185),
}


def display_unit(key: str, si_mode: bool) -> str:
    """Return the unit label for the requested display mode."""
    british_unit, si_unit, _ = UNIT_CONV[key]
    return si_unit if si_mode else british_unit


def convert_from_british(value: float, key: str, si_mode: bool) -> float:
    """Convert a British-unit value to the requested display mode."""
    if not si_mode:
        return value
    _, _, factor = UNIT_CONV[key]
    return value * factor


def convert_to_british(value: float, key: str, si_mode: bool) -> float:
    """Convert a displayed value back to British solver units."""
    if not si_mode:
        return value
    _, _, factor = UNIT_CONV[key]
    return value / factor


def convert_array_to_british(
    values: Iterable[float], key: str, si_mode: bool
) -> list[float]:
    """Convert displayed values back to British solver units."""
    return [convert_to_british(value, key, si_mode) for value in values]


def convert_array_from_british(
    values: Iterable[float], key: str, si_mode: bool
) -> list[float]:
    """Convert British-unit values to the requested display mode."""
    return [convert_from_british(value, key, si_mode) for value in values]
