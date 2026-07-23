"""Headless LIQLEV sweep execution."""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any

from liqlev.model.config import SimulationConfig
from liqlev.model.validation import validate_simulation_config
from liqlev.runner.progress import ProgressCallback, ProgressEvent, emit_progress
from liqlev.runner.single import (
    build_property_table_for_config,
    load_vent_profile,
    prepare_gravity,
    run_single_case,
)


@dataclass(frozen=True)
class SweepResult:
    scenarios: dict[str, dict[str, Any]]
    run_count: int
    elapsed_s: float


def epsilon_specs(config: SimulationConfig) -> tuple[str | float, ...]:
    """Return the sweep epsilon specs in legacy order."""
    if config.epsilon.mode == "Custom":
        return tuple(config.epsilon.values)
    return (config.epsilon.mode,)


def vent_rates(config: SimulationConfig) -> tuple[float, ...]:
    """Return the sweep vent rates, preserving CSV placeholder behavior."""
    if config.vent.csv_path:
        return (0.0,)
    return tuple(config.vent.rates_lbm_s)


def run_sweep(
    config: SimulationConfig, progress_cb: ProgressCallback | None = None
) -> SweepResult:
    """Run fill x epsilon x vent sweeps without starting a GUI."""
    validate_simulation_config(config)

    fills = tuple(config.tank.fill_fractions)
    eps_values = epsilon_specs(config)
    rates = vent_rates(config)
    total_runs = len(fills) * len(eps_values) * len(rates)
    run_index = 0
    scenarios: dict[str, dict[str, Any]] = {}
    prepared_gravity = prepare_gravity(config)
    vent_profile = load_vent_profile(config)
    prop_table = build_property_table_for_config(config)

    start = time.perf_counter()
    for message in prepared_gravity.messages:
        emit_progress(progress_cb, ProgressEvent(kind="log", message=message))

    for fill_fraction in fills:
        for epsilon_spec in eps_values:
            scenario_key = None
            scenario_dfs = []
            scenario_rates = []
            epsilon_label = ""

            for rate in rates:
                run_index += 1
                result = run_single_case(
                    config,
                    fill_fraction=fill_fraction,
                    vent_rate_lbm_s=rate,
                    epsilon_spec=epsilon_spec,
                    prepared_gravity=prepared_gravity,
                    vent_profile=vent_profile,
                    prop_table=prop_table,
                    progress_cb=progress_cb,
                    run_index=run_index,
                    total_runs=total_runs,
                )
                scenario_key = result.scenario_key
                epsilon_label = result.epsilon_label
                scenario_dfs.append(result.dataframe)
                scenario_rates.append(result.vent_rate_lbm_s)

            if scenario_key is not None:
                scenarios[scenario_key] = {
                    "dfs": scenario_dfs,
                    "vent_rates": scenario_rates,
                    "fill": fill_fraction,
                    "eps_label": epsilon_label,
                    "htank": config.tank.height_ft,
                }

    elapsed_s = time.perf_counter() - start
    emit_progress(
        progress_cb,
        ProgressEvent(
            kind="complete",
            message=f"Completed {run_index} LIQLEV run(s) in {elapsed_s:.2f}s",
            fraction=1.0,
            run_index=run_index,
            total_runs=total_runs,
        ),
    )
    return SweepResult(scenarios=scenarios, run_count=run_index, elapsed_s=elapsed_s)
