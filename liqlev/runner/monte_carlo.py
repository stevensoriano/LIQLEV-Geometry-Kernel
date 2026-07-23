"""Headless Monte Carlo execution for LIQLEV sensitivity analysis."""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np

from liqlev.model.builder import G_TO_FT_S2, build_inputs, epsilon_schedule
from liqlev.model.config import SimulationConfig
from liqlev.model.validation import validate_simulation_config
from liqlev.runner.progress import ProgressCallback, ProgressEvent, emit_progress
from liqlev.runner.single import build_property_table_for_config
from core import liqlev_simulation


@dataclass(frozen=True)
class MonteCarloRequest:
    n: int
    vent_min_lbm_s: float
    vent_max_lbm_s: float
    fill_min: float
    fill_max: float
    gravity_min_g: float
    gravity_max_g: float
    seed: int | None = None


@dataclass(frozen=True)
class MonteCarloResult:
    n: int
    all_dh: list[float]
    all_params: list[dict[str, float]]
    max_dh: float
    mean_dh: float
    std_dh: float
    p95: float
    p99: float
    worst: dict[str, float]
    elapsed_s: float


def validate_monte_carlo_request(request: MonteCarloRequest) -> None:
    """Validate Monte Carlo sample count and ranges."""
    if request.n < 2:
        raise ValueError("N Samples must be at least 2.")
    if request.vent_min_lbm_s >= request.vent_max_lbm_s:
        raise ValueError("Vent Rate Min must be less than Max.")
    if request.fill_min >= request.fill_max:
        raise ValueError("Fill Frac Min must be less than Max.")
    if request.gravity_min_g >= request.gravity_max_g:
        raise ValueError("Gravity Min must be less than Max.")


def run_monte_carlo(
    config: SimulationConfig,
    request: MonteCarloRequest,
    progress_cb: ProgressCallback | None = None,
) -> MonteCarloResult:
    """Run legacy-compatible Monte Carlo sampling without GUI objects."""
    geometry = validate_simulation_config(config)
    validate_monte_carlo_request(request)

    prop_table = build_property_table_for_config(config)
    rng = np.random.default_rng(request.seed)
    start = time.perf_counter()

    eps_spec = (
        config.epsilon.values[0]
        if config.epsilon.mode == "Custom"
        else config.epsilon.mode
    )
    neps, teps, xeps, _ = epsilon_schedule(eps_spec, config.run.duration_s)

    all_dh: list[float] = []
    all_params: list[dict[str, float]] = []
    worst = {"dh": 0.0, "vent": 0.0, "fill": 0.0, "grav": 0.0}

    for index in range(request.n):
        vent = float(rng.uniform(request.vent_min_lbm_s, request.vent_max_lbm_s))
        fill = float(rng.uniform(request.fill_min, request.fill_max))
        grav = float(rng.uniform(request.gravity_min_g, request.gravity_max_g))

        grav_ft = grav * G_TO_FT_S2
        tggo = np.array([0.0, config.run.duration_s])
        xggo = np.array([grav_ft, grav_ft])

        # Preserve legacy GUI Monte Carlo behavior: sampled constant gravity,
        # first epsilon only, and no AS-203 measured mass/temperature overrides.
        inputs = build_inputs(
            fluid=config.fluid.name,
            pinit_psia=config.fluid.initial_pressure_psia,
            pfinal_psia=config.fluid.final_pressure_psia,
            dtank=config.tank.diameter_ft,
            htank=config.tank.height_ft,
            fill_fraction=fill,
            duration=config.run.duration_s,
            delta_t=config.run.timestep_s,
            vent_rate=vent,
            neps=neps,
            teps=teps,
            xeps=xeps,
            ramp_duration=config.vent.ramp_duration_s,
            ramp_target_factor=config.vent.ramp_target_factor,
            nggo=2,
            tggo=tggo,
            xggo=xggo,
            geometry=geometry,
        )

        dataframe = liqlev_simulation(inputs, verbose=False, prop_table=prop_table)
        max_dh = float(dataframe["Hratio"].max()) if not dataframe.empty else 0.0
        all_dh.append(max_dh)
        all_params.append({"vent": vent, "fill": fill, "grav": grav})

        if max_dh > worst["dh"]:
            worst = {"dh": max_dh, "vent": vent, "fill": fill, "grav": grav}

        fraction = (index + 1) / request.n
        emit_progress(
            progress_cb,
            ProgressEvent(
                kind="solver_progress",
                message=f"Monte Carlo sample {index + 1}/{request.n}",
                fraction=fraction,
                run_index=index + 1,
                total_runs=request.n,
                stats={"max_dh": max_dh, "max_so_far": worst["dh"]},
            ),
        )

    arr = np.array(all_dh)
    elapsed_s = time.perf_counter() - start
    result = MonteCarloResult(
        n=request.n,
        all_dh=all_dh,
        all_params=all_params,
        max_dh=float(arr.max()),
        mean_dh=float(arr.mean()),
        std_dh=float(arr.std()),
        p95=float(np.percentile(arr, 95)),
        p99=float(np.percentile(arr, 99)),
        worst=worst,
        elapsed_s=elapsed_s,
    )
    emit_progress(
        progress_cb,
        ProgressEvent(
            kind="complete",
            message=f"Monte Carlo complete ({request.n} samples, {elapsed_s:.2f}s)",
            fraction=1.0,
            run_index=request.n,
            total_runs=request.n,
        ),
    )
    return result
