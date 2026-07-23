"""Tests for plot-ready result transforms and risk summaries."""

from __future__ import annotations

import numpy as np

from liqlev.viz.datasets import (
    boundary_layer_traces,
    event_evolution_traces,
    pressure_level_trace,
)
from liqlev.viz.summaries import risk_summary, summary_rows
from validation.physics_cases import get_case, run_case


def test_visualization_traces_are_column_derived() -> None:
    dataframe = run_case(get_case("as203_default_high_vent"))

    event_traces = event_evolution_traces(dataframe)
    pressure_level = pressure_level_trace(dataframe)
    boundary_layer = boundary_layer_traces(dataframe)

    assert {trace.name for trace in event_traces} == {
        "Liquid level",
        "Ullage mass",
        "Pressure",
        "Vent rate",
    }
    assert np.array_equal(pressure_level.x, dataframe["Hratio"].to_numpy(dtype=float))
    assert np.array_equal(pressure_level.y, dataframe["Press"].to_numpy(dtype=float))
    assert {trace.name for trace in boundary_layer} == {
        "Boundary-layer volume",
        "Boundary-layer thickness",
        "Vapor in BL",
        "BL vapor out",
    }


def test_risk_summary_detects_threshold_and_tank_crossings() -> None:
    dataframe = run_case(get_case("as203_default_high_vent"))
    summary = risk_summary(dataframe, htank_ft=28.18, threshold_dh_h0=0.01)

    assert summary.max_dh_h0 > 0.01
    assert summary.threshold_crossed is True
    assert summary.threshold_time_s is not None
    assert summary.tank_exceeded is False
    assert summary.convergence_failures == int(dataframe["Conv Failed"].sum())


def test_summary_rows_match_export_contract() -> None:
    dataframe = run_case(get_case("hydrogen_height_dep_mid_fill"))
    scenarios = {
        "Fill 50%, eps=height_dep": {
            "dfs": [dataframe],
            "vent_rates": [0.0015],
            "fill": 0.5,
            "eps_label": "height_dep",
            "htank": 28.18,
        }
    }

    rows = summary_rows(scenarios, threshold_dh_h0=0.1)

    assert rows[0]["Scenario"] == "Fill 50%, eps=height_dep"
    assert rows[0]["Vent Rate (lbm/s)"] == 0.0015
    assert rows[0]["Threshold Crossed"] == "No"
    assert rows[0]["Convergence Failures"] == int(dataframe["Conv Failed"].sum())
