"""Headless Monte Carlo execution for LIQLEV sensitivity analysis.

F17 / Phase 4.11 — abort integrity
---------------------------------
Custom-geometry samples can abort (``Conv Failed``) and emit a single bounded
diagnostic row with ``Hratio ≈ 0``. Including those rows in the Monte Carlo
population biases ``max`` / ``mean`` / ``p95`` / ``p99`` **low** — the
non-conservative direction for level-rise analysis.

Contract (measured, not silently relaxed):
- Per-sample ``Conv Failed`` detection (any positive sum → aborted sample).
- Aborted samples are **excluded** from ``all_dh`` and from dh statistics.
- ``MonteCarloResult.aborted_count`` and ``aborted_params`` surface every abort.
- If ``aborted_count / n > ABORT_FRACTION_THRESHOLD`` (10%), the run fails
  loudly with :class:`MonteCarloAbortError`. The 10% threshold is small enough
  that a non-trivial abort cluster cannot hide inside the population, while
  still tolerating a rare isolated diagnostic without discarding an entire
  study. It is a fixed contract of this module.
"""

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


# F17 contract: fail the whole Monte Carlo run when more than this fraction of
# samples abort (Conv Failed). Chosen at 10% — small enough to keep statistics
# representative, large enough that a single rare diagnostic does not discard
# a study. Do not silently relax this constant.
ABORT_FRACTION_THRESHOLD = 0.10


class MonteCarloAbortError(RuntimeError):
    """Raised when the Monte Carlo abort fraction exceeds the contract threshold.

    Attributes mirror the partial run so callers can report which draws aborted
    without re-running the study.
    """

    def __init__(
        self,
        message: str,
        *,
        aborted_count: int,
        n: int,
        aborted_params: list[dict[str, float]],
        threshold: float = ABORT_FRACTION_THRESHOLD,
    ) -> None:
        super().__init__(message)
        self.aborted_count = aborted_count
        self.n = n
        self.aborted_params = list(aborted_params)
        self.threshold = threshold
        self.abort_fraction = aborted_count / n if n else 1.0


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
    aborted_count: int = 0
    aborted_params: list[dict[str, float]] | None = None

    def __post_init__(self) -> None:
        # frozen dataclass: normalize None → empty list for callers
        if self.aborted_params is None:
            object.__setattr__(self, "aborted_params", [])


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


def _sample_is_aborted(dataframe) -> bool:
    """Return True when a sample produced any solver convergence failure."""
    if dataframe is None or getattr(dataframe, "empty", True):
        return False
    if "Conv Failed" not in dataframe.columns:
        return False
    return int(dataframe["Conv Failed"].sum()) > 0


def run_monte_carlo(
    config: SimulationConfig,
    request: MonteCarloRequest,
    progress_cb: ProgressCallback | None = None,
) -> MonteCarloResult:
    """Run legacy-compatible Monte Carlo sampling without GUI objects.

    Aborted custom-geometry samples (``Conv Failed``) are counted, excluded
    from dh statistics, and surfaced on the result. When the abort fraction
    exceeds :data:`ABORT_FRACTION_THRESHOLD`, raises
    :class:`MonteCarloAbortError` (F17).
    """
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
    aborted_params: list[dict[str, float]] = []
    worst = {"dh": 0.0, "vent": 0.0, "fill": 0.0, "grav": 0.0}

    for index in range(request.n):
        vent = float(rng.uniform(request.vent_min_lbm_s, request.vent_max_lbm_s))
        fill = float(rng.uniform(request.fill_min, request.fill_max))
        grav = float(rng.uniform(request.gravity_min_g, request.gravity_max_g))
        draw = {"vent": vent, "fill": fill, "grav": grav}

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
            boundary_layer_substeps=config.run.boundary_layer_substeps,
        )

        dataframe = liqlev_simulation(inputs, verbose=False, prop_table=prop_table)

        if _sample_is_aborted(dataframe):
            # F17: do not let Hratio≈0 diagnostic rows enter the dh population.
            aborted_params.append(draw)
            fraction = (index + 1) / request.n
            emit_progress(
                progress_cb,
                ProgressEvent(
                    kind="solver_progress",
                    message=(
                        f"Monte Carlo sample {index + 1}/{request.n} aborted "
                        f"(Conv Failed)"
                    ),
                    fraction=fraction,
                    run_index=index + 1,
                    total_runs=request.n,
                    stats={
                        "aborted": 1.0,
                        "aborted_count": float(len(aborted_params)),
                        "max_so_far": worst["dh"],
                    },
                ),
            )
            continue

        max_dh = float(dataframe["Hratio"].max()) if not dataframe.empty else 0.0
        all_dh.append(max_dh)
        all_params.append(draw)

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

    aborted_count = len(aborted_params)
    abort_fraction = aborted_count / request.n
    if abort_fraction > ABORT_FRACTION_THRESHOLD:
        raise MonteCarloAbortError(
            (
                f"Monte Carlo abort fraction {abort_fraction:.1%} "
                f"({aborted_count}/{request.n}) exceeds the "
                f"{ABORT_FRACTION_THRESHOLD:.0%} contract threshold "
                f"(F17); aborted draws are excluded from dh statistics but "
                f"too many aborts invalidate the study."
            ),
            aborted_count=aborted_count,
            n=request.n,
            aborted_params=aborted_params,
            threshold=ABORT_FRACTION_THRESHOLD,
        )

    if not all_dh:
        # Defensive: only reachable if threshold >= 1.0; keep loud failure.
        raise MonteCarloAbortError(
            "Monte Carlo produced no successful samples; all draws aborted.",
            aborted_count=aborted_count,
            n=request.n,
            aborted_params=aborted_params,
            threshold=ABORT_FRACTION_THRESHOLD,
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
        aborted_count=aborted_count,
        aborted_params=aborted_params,
    )
    emit_progress(
        progress_cb,
        ProgressEvent(
            kind="complete",
            message=(
                f"Monte Carlo complete ({request.n} samples, "
                f"{aborted_count} aborted, {elapsed_s:.2f}s)"
            ),
            fraction=1.0,
            run_index=request.n,
            total_runs=request.n,
            stats={"aborted_count": float(aborted_count)},
        ),
    )
    return result
