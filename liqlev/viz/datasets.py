"""Plot-ready dataset transforms derived from solver DataFrame columns."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Trace:
    name: str
    x: np.ndarray
    y: np.ndarray
    y_label: str
    unit: str = ""


def _array(df: pd.DataFrame, column: str) -> np.ndarray:
    if column not in df.columns:
        raise KeyError(f"Missing result column: {column}")
    return df[column].to_numpy(dtype=float)


def normalized(values: np.ndarray) -> np.ndarray:
    """Normalize values to 0..1, preserving flat series at zero."""
    values = np.asarray(values, dtype=float)
    min_value = float(np.nanmin(values))
    max_value = float(np.nanmax(values))
    span = max_value - min_value
    if span == 0:
        return np.zeros_like(values)
    return (values - min_value) / span


def time_trace(
    df: pd.DataFrame, column: str, name: str, y_label: str, unit: str = ""
) -> Trace:
    """Return a single time-history trace for a solver column."""
    return Trace(
        name=name,
        x=_array(df, "Time"),
        y=_array(df, column),
        y_label=y_label,
        unit=unit,
    )


def pressure_level_trace(df: pd.DataFrame) -> Trace:
    """Return pressure versus liquid level rise."""
    return Trace(
        name="Pressure vs dh/h0",
        x=_array(df, "Hratio"),
        y=_array(df, "Press"),
        y_label="Pressure",
        unit="psia",
    )


def event_evolution_traces(
    df: pd.DataFrame, threshold_dh_h0: float | None = None
) -> list[Trace]:
    """Return normalized tank-state traces for linked event-evolution views."""
    time_s = _array(df, "Time")
    traces = [
        Trace(
            "Liquid level", time_s, normalized(_array(df, "Height")), "Normalized State"
        ),
        Trace(
            "Ullage mass",
            time_s,
            normalized(_array(df, "Ullage Mass")),
            "Normalized State",
        ),
        Trace("Pressure", time_s, normalized(_array(df, "Press")), "Normalized State"),
        Trace(
            "Vent rate", time_s, normalized(_array(df, "Vent Rate")), "Normalized State"
        ),
    ]
    if threshold_dh_h0 is not None:
        traces.append(
            Trace(
                "Risk threshold",
                time_s,
                np.full_like(time_s, threshold_dh_h0, dtype=float),
                "dh/h0",
            )
        )
    return traces


def boundary_layer_traces(df: pd.DataFrame) -> list[Trace]:
    """Return boundary-layer diagnostic traces."""
    return [
        time_trace(df, "VBL vol", "Boundary-layer volume", "Volume", "ft^3"),
        time_trace(df, "BL thick", "Boundary-layer thickness", "Thickness", "ft"),
        time_trace(df, "Vapor in BL", "Vapor in BL", "Mass", "lbm"),
        time_trace(df, "BL Vap Out", "BL vapor out", "Mass Flow", "lbm/s"),
    ]


def convergence_traces(df: pd.DataFrame) -> list[Trace]:
    """Return convergence health traces."""
    return [
        time_trace(df, "Conv Iterations", "Iterations", "Iterations"),
        time_trace(df, "Conv Failed", "Failure flag", "Failure"),
    ]
