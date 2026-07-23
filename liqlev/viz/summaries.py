"""Engineering summaries and threshold detection for LIQLEV results."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import pandas as pd


@dataclass(frozen=True)
class RiskSummary:
    max_dh_h0: float
    time_to_peak_s: float
    final_pressure_psia: float
    threshold_crossed: bool
    threshold_time_s: float | None
    tank_exceeded: bool
    tank_exceeded_time_s: float | None
    convergence_failures: int


def first_crossing_time(
    df: pd.DataFrame, column: str, threshold: float
) -> float | None:
    """Return the first time where column >= threshold."""
    crossed = df[df[column] >= threshold]
    if crossed.empty:
        return None
    return float(crossed["Time"].iloc[0])


def risk_summary(
    df: pd.DataFrame,
    *,
    htank_ft: float,
    threshold_dh_h0: float | None = None,
) -> RiskSummary:
    """Summarize peak level, pressure, convergence, and risk threshold crossings."""
    if df.empty:
        return RiskSummary(
            max_dh_h0=0.0,
            time_to_peak_s=0.0,
            final_pressure_psia=0.0,
            threshold_crossed=False,
            threshold_time_s=None,
            tank_exceeded=False,
            tank_exceeded_time_s=None,
            convergence_failures=0,
        )

    peak = df.loc[df["Hratio"].idxmax()]
    threshold_time = None
    if threshold_dh_h0 is not None:
        threshold_time = first_crossing_time(df, "Hratio", threshold_dh_h0)

    tank_time = first_crossing_time(df, "Height", htank_ft)
    conv_fails = int(df["Conv Failed"].sum()) if "Conv Failed" in df.columns else 0

    return RiskSummary(
        max_dh_h0=float(peak["Hratio"]),
        time_to_peak_s=float(peak["Time"]),
        final_pressure_psia=float(df["Press"].iloc[-1]),
        threshold_crossed=threshold_time is not None,
        threshold_time_s=threshold_time,
        tank_exceeded=tank_time is not None,
        tank_exceeded_time_s=tank_time,
        convergence_failures=conv_fails,
    )


def summary_rows(
    scenarios: dict[str, dict[str, Any]], threshold_dh_h0: float | None = None
) -> list[dict[str, Any]]:
    """Return one engineering summary row per scenario and vent rate."""
    rows: list[dict[str, Any]] = []
    for scenario_key, scenario in scenarios.items():
        for dataframe, rate in zip(scenario["dfs"], scenario["vent_rates"]):
            summary = risk_summary(
                dataframe,
                htank_ft=float(scenario["htank"]),
                threshold_dh_h0=threshold_dh_h0,
            )
            rows.append(
                {
                    "Scenario": scenario_key,
                    "Vent Rate (lbm/s)": rate,
                    "Max dh/h0": round(summary.max_dh_h0, 6),
                    "Time to Peak (s)": round(summary.time_to_peak_s, 2),
                    "Final Pressure (psia)": round(summary.final_pressure_psia, 4),
                    "Threshold Crossed": "YES" if summary.threshold_crossed else "No",
                    "Threshold Time (s)": summary.threshold_time_s,
                    "Tank Exceeded": "YES" if summary.tank_exceeded else "No",
                    "Tank Exceeded Time (s)": summary.tank_exceeded_time_s,
                    "Convergence Failures": summary.convergence_failures,
                }
            )
    return rows
