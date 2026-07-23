"""Headless single-case LIQLEV execution."""

from __future__ import annotations

import time
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from core import liqlev_simulation
from liqlev.geometry.schema import GeometryKernel
from liqlev.io.profiles import (
    ProfileData,
    load_gravity_profile_csv,
    load_vent_rate_profile_csv,
)
from liqlev.model.builder import (
    G_TO_FT_S2,
    build_inputs,
    epsilon_schedule,
    make_gravity_function,
    safe_eval_gravity,
)
from liqlev.model.config import SimulationConfig
from liqlev.model.validation import validate_simulation_config
from liqlev.runner.progress import ProgressCallback, ProgressEvent, emit_progress
from thermo_utils import build_property_table


@dataclass(frozen=True)
class PreparedGravity:
    nggo: int
    tggo: np.ndarray
    xggo: np.ndarray
    gravity_function: Callable[[float], float] | None = None
    messages: tuple[str, ...] = ()


@dataclass(frozen=True)
class SingleCaseResult:
    dataframe: pd.DataFrame
    inputs: dict[str, Any]
    scenario_key: str
    vent_rate_lbm_s: float
    fill_fraction: float
    epsilon_label: str
    htank_ft: float
    elapsed_s: float
    warnings: tuple[str, ...] = field(default_factory=tuple)


def prepare_gravity(config: SimulationConfig) -> PreparedGravity:
    """Build legacy gravity arrays/function from typed config."""
    duration = config.run.duration_s
    gravity = config.gravity

    if gravity.mode == "Constant":
        const_g_ft = gravity.constant_g * G_TO_FT_S2
        return PreparedGravity(
            nggo=2,
            tggo=np.array([0.0, duration]),
            xggo=np.array([const_g_ft, const_g_ft]),
            messages=(f"Gravity: Constant {gravity.constant_g} g",),
        )

    if gravity.mode == "Function of Time":
        test_times = [0.0, duration / 4, duration / 2, duration]
        messages = [f"Gravity: g(t) = {gravity.expression}"]
        for test_time in test_times:
            value = safe_eval_gravity(gravity.expression, test_time)
            messages.append(f"g({test_time:.1f}s) = {value:.6f} g")
        return PreparedGravity(
            nggo=2,
            tggo=np.array([0.0, duration]),
            xggo=np.array([0.0, 0.0]),
            gravity_function=make_gravity_function(gravity.expression),
            messages=tuple(messages),
        )

    if gravity.mode == "CSV Profile":
        profile = load_gravity_profile_csv(gravity.csv_path, duration, gravity.hold_g)
        return PreparedGravity(
            nggo=profile.point_count,
            tggo=profile.time_s,
            xggo=profile.values,
            messages=(f"Gravity: CSV Profile from {profile.source}",),
        )

    raise ValueError(f"Unsupported gravity mode: {gravity.mode}")


def load_vent_profile(config: SimulationConfig) -> ProfileData | None:
    """Load a vent profile override when configured."""
    if not config.vent.csv_path:
        return None
    return load_vent_rate_profile_csv(config.vent.csv_path)


def build_property_table_for_config(config: SimulationConfig):
    """Build a non-hydrogen property table once per execution context."""
    if config.fluid.name == "Hydrogen":
        return None
    return build_property_table(
        config.fluid.name,
        config.fluid.final_pressure_psia,
        config.fluid.initial_pressure_psia,
    )


def run_single_case(
    config: SimulationConfig,
    *,
    fill_fraction: float | None = None,
    vent_rate_lbm_s: float | None = None,
    epsilon_spec: str | float | None = None,
    prepared_gravity: PreparedGravity | None = None,
    vent_profile: ProfileData | None = None,
    prop_table=None,
    progress_cb: ProgressCallback | None = None,
    run_index: int | None = None,
    total_runs: int | None = None,
    geometry: GeometryKernel | None = None,
) -> SingleCaseResult:
    """Run one LIQLEV case without constructing any GUI objects."""
    geometry = validate_simulation_config(config, geometry=geometry)

    fill = config.tank.fill_fractions[0] if fill_fraction is None else fill_fraction
    vent_rate = (
        config.vent.rates_lbm_s[0] if vent_rate_lbm_s is None else vent_rate_lbm_s
    )
    if config.vent.csv_path and vent_rate_lbm_s is None:
        vent_rate = 0.0

    if epsilon_spec is None:
        epsilon_spec = (
            config.epsilon.values[0]
            if config.epsilon.mode == "Custom"
            else config.epsilon.mode
        )

    gravity = prepared_gravity or prepare_gravity(config)
    loaded_vent_profile = vent_profile
    if loaded_vent_profile is None:
        loaded_vent_profile = load_vent_profile(config)

    if prop_table is None:
        prop_table = build_property_table_for_config(config)

    neps, teps, xeps, epsilon_label = epsilon_schedule(
        epsilon_spec, config.run.duration_s
    )
    scenario_key = f"Fill {fill * 100:.0f}%, eps={epsilon_label}"
    emit_progress(
        progress_cb,
        ProgressEvent(
            kind="run_start",
            message=scenario_key,
            scenario_key=scenario_key,
            run_index=run_index,
            total_runs=total_runs,
        ),
    )

    inputs = build_inputs(
        fluid=config.fluid.name,
        pinit_psia=config.fluid.initial_pressure_psia,
        pfinal_psia=config.fluid.final_pressure_psia,
        dtank=config.tank.diameter_ft,
        htank=config.tank.height_ft,
        fill_fraction=fill,
        duration=config.run.duration_s,
        delta_t=config.run.timestep_s,
        vent_rate=vent_rate,
        neps=neps,
        teps=teps,
        xeps=xeps,
        ramp_duration=config.vent.ramp_duration_s,
        ramp_target_factor=config.vent.ramp_target_factor,
        nggo=gravity.nggo,
        tggo=gravity.tggo,
        xggo=gravity.xggo,
        gravity_function=gravity.gravity_function,
        xmlzro_override=config.fluid.initial_mass_lbm,
        tinit_override=config.fluid.initial_temperature_r,
        geometry=geometry,
    )

    if loaded_vent_profile is not None:
        inputs["Nvmd"] = loaded_vent_profile.point_count
        inputs["Tvmdot"] = loaded_vent_profile.time_s
        inputs["Xvmdot"] = loaded_vent_profile.values

    def on_solver_progress(stats: dict[str, float]) -> None:
        emit_progress(
            progress_cb,
            ProgressEvent(
                kind="solver_progress",
                scenario_key=scenario_key,
                run_index=run_index,
                total_runs=total_runs,
                stats=stats,
            ),
        )

    start = time.perf_counter()
    dataframe = liqlev_simulation(
        inputs,
        verbose=False,
        prop_table=prop_table,
        progress_cb=on_solver_progress if progress_cb is not None else None,
    )
    elapsed_s = time.perf_counter() - start

    warnings: list[str] = []
    if not dataframe.empty and "Conv Failed" in dataframe.columns:
        conv_fails = int(dataframe["Conv Failed"].sum())
        if conv_fails > 0:
            warnings.append(f"{conv_fails} timesteps had solver convergence failures")

    for warning in warnings:
        emit_progress(
            progress_cb,
            ProgressEvent(kind="warning", message=warning, scenario_key=scenario_key),
        )

    emit_progress(
        progress_cb,
        ProgressEvent(
            kind="run_complete",
            message=scenario_key,
            scenario_key=scenario_key,
            run_index=run_index,
            total_runs=total_runs,
            fraction=(run_index / total_runs if run_index and total_runs else None),
        ),
    )

    return SingleCaseResult(
        dataframe=dataframe,
        inputs=inputs,
        scenario_key=scenario_key,
        vent_rate_lbm_s=vent_rate,
        fill_fraction=fill,
        epsilon_label=epsilon_label,
        htank_ft=(
            geometry.total_height_ft
            if geometry is not None
            else config.tank.height_ft
        ),
        elapsed_s=elapsed_s,
        warnings=tuple(warnings),
    )
