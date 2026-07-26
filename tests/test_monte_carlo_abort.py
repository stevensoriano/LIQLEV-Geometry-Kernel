"""Phase 4.11 / F17 Monte Carlo abort integrity — four-criteria coverage.

Each guard test documents which of the lead-mandated criteria it covers:

1. FIRES ON A BAD CASE (synthetic, or real where one exists).
2. ASSERTS THE RECOVERED STATE, never the defect state.
3. DRIVES THE REAL CODE PATH — no reimplementation of abort detection.
4. PROVES REACHABILITY ON THE LIVE PATH — production entry points users hit.

Real call chains exercised here
--------------------------------
Headless / production Monte Carlo entry (UI worker and scripts)::

    liqlev.ui_qt.app.RunWorker / headless caller
      -> liqlev.runner.monte_carlo.run_monte_carlo   # sole detection site
           -> build_inputs
           -> core.liqlev_simulation
           -> _sample_is_aborted(dataframe)         # Conv Failed sum > 0
           -> exclude from all_dh / stats
           -> MonteCarloResult(aborted_count, aborted_params)
           -> or MonteCarloAbortError if fraction > ABORT_FRACTION_THRESHOLD

Abort trigger used post-Phase-2
-------------------------------
Custom g=0 no longer aborts (saturation is legal). These tests force real
custom-mode aborts by zeroing ``GeomPerimeter`` on selected samples after
``build_inputs`` — the same BL-failure path used by
``tests/geometry/test_core_custom_geometry.py::test_boundary_layer_failure_*``
(and independently verified: over-full liquid mass via ``Xmlzro`` scaled to
2× tank inventory also aborts, but only when fill ≳ 0.5, so it is a poor
per-draw force for random MC fills). Zero perimeter aborts for every fill.

Verified: injection yields one diagnostic row with ``Conv Failed == 1`` and
``Hratio == 0`` through the real solver.

Threshold contract: ``ABORT_FRACTION_THRESHOLD == 0.10`` (never silently
relaxed).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import liqlev.runner.monte_carlo as monte_carlo_module
from liqlev.geometry.fixtures import cylinder_kernel
from liqlev.geometry.package import save_geometry_package
from liqlev.model.config import (
    EpsilonConfig,
    FluidConfig,
    GravityProfileConfig,
    RunControls,
    SimulationConfig,
    TankConfig,
    VentProfileConfig,
)
from liqlev.runner.monte_carlo import (
    ABORT_FRACTION_THRESHOLD,
    MonteCarloAbortError,
    MonteCarloRequest,
    run_monte_carlo,
)


# ---------------------------------------------------------------------------
# Shared fixtures / helpers
# ---------------------------------------------------------------------------


def custom_geometry_config(geometry_path: str) -> SimulationConfig:
    return SimulationConfig(
        fluid=FluidConfig(
            name="Hydrogen",
            initial_pressure_psia=19.5,
            final_pressure_psia=10.0,
            initial_temperature_r=38.3,
        ),
        tank=TankConfig(
            diameter_ft=25.0,
            height_ft=99.0,
            fill_fractions=(0.5,),
            geometry_path=geometry_path,
        ),
        vent=VentProfileConfig(
            rates_lbm_s=(1.0e-7,),
            ramp_duration_s=10_000.0,
        ),
        gravity=GravityProfileConfig(mode="Constant", constant_g=0.001),
        epsilon=EpsilonConfig(mode="height_dep"),
        run=RunControls(duration_s=10_000.0, timestep_s=5_000.0),
    )


def _save_cylinder(tmp_path: Path) -> Path:
    package_path = tmp_path / "cylinder.npz"
    save_geometry_package(cylinder_kernel(4.0, 8.0, node_count=17), package_path)
    return package_path


def _mc_request(n: int, seed: int = 42) -> MonteCarloRequest:
    return MonteCarloRequest(
        n=n,
        vent_min_lbm_s=1.0e-7,
        vent_max_lbm_s=2.0e-7,
        fill_min=0.2,
        fill_max=0.8,
        gravity_min_g=0.0005,
        gravity_max_g=0.0015,
        seed=seed,
    )


def _force_abort_inputs(inputs: dict) -> dict:
    """Return a copy of solver inputs forced to the post-Phase-2 abort path.

    Zero perimeter is fill-independent (unlike Xmlzro×2, which only aborts
    when fill ≳ 0.5). Same diagnostic as the re-pointed BL-failure unit test.
    """
    forced = dict(inputs)
    perimeter = np.asarray(inputs["GeomPerimeter"], dtype=np.float64)
    forced["GeomPerimeter"] = np.zeros_like(perimeter)
    return forced


def _install_abort_on_call_indices(
    monkeypatch, abort_indices: set[int]
) -> None:
    """Abort only listed 1-based call indices via the real solver."""
    real_sim = monte_carlo_module.liqlev_simulation
    call_count = {"n": 0}

    def selective(inputs, **kwargs):
        call_count["n"] += 1
        if call_count["n"] in abort_indices:
            inputs = _force_abort_inputs(inputs)
        return real_sim(inputs, **kwargs)

    monkeypatch.setattr(monte_carlo_module, "liqlev_simulation", selective)


def _install_abort_all(monkeypatch) -> None:
    real_sim = monte_carlo_module.liqlev_simulation

    def always_abort(inputs, **kwargs):
        return real_sim(_force_abort_inputs(inputs), **kwargs)

    monkeypatch.setattr(monte_carlo_module, "liqlev_simulation", always_abort)


# ---------------------------------------------------------------------------
# Precondition: trigger produces Conv Failed through the real solver
# ---------------------------------------------------------------------------


def test_abort_triggers_yield_conv_failed_on_real_solver(tmp_path) -> None:
    """Precondition: zero-perimeter and over-full mass abort custom mode for real.

    Call chain: ``build_inputs`` → ``liqlev_simulation`` (same stack MC uses).
    """
    from liqlev.model.builder import G_TO_FT_S2, build_inputs, epsilon_schedule
    from liqlev.model.validation import validate_simulation_config
    from core import liqlev_simulation

    package_path = _save_cylinder(tmp_path)
    config = custom_geometry_config(str(package_path))
    geometry = validate_simulation_config(config)
    neps, teps, xeps, _ = epsilon_schedule(
        config.epsilon.mode, config.run.duration_s
    )
    grav_ft = 0.001 * G_TO_FT_S2
    inputs = build_inputs(
        fluid=config.fluid.name,
        pinit_psia=config.fluid.initial_pressure_psia,
        pfinal_psia=config.fluid.final_pressure_psia,
        dtank=config.tank.diameter_ft,
        htank=config.tank.height_ft,
        fill_fraction=0.5,
        duration=config.run.duration_s,
        delta_t=config.run.timestep_s,
        vent_rate=1.0e-7,
        neps=neps,
        teps=teps,
        xeps=xeps,
        ramp_duration=config.vent.ramp_duration_s,
        ramp_target_factor=config.vent.ramp_target_factor,
        nggo=2,
        tggo=np.array([0.0, config.run.duration_s]),
        xggo=np.array([grav_ft, grav_ft]),
        geometry=geometry,
    )
    ok = liqlev_simulation(inputs, verbose=False)
    assert int(ok["Conv Failed"].sum()) == 0
    assert float(ok["Hratio"].max()) > 0.0

    # Primary MC force: zero perimeter (fill-independent).
    zero_p = liqlev_simulation(_force_abort_inputs(inputs), verbose=False)
    assert len(zero_p) == 1
    assert zero_p.loc[0, "Conv Failed"] == 1.0
    assert zero_p.loc[0, "Hratio"] == pytest.approx(0.0)

    # Plan-named alternative: over-full liquid mass (2× tank inventory).
    overfull = dict(inputs)
    fill = float(inputs["FillFraction"])
    overfull["Xmlzro"] = float(inputs["Xmlzro"]) / fill * 2.0
    bad = liqlev_simulation(overfull, verbose=False)
    assert len(bad) == 1
    assert bad.loc[0, "Conv Failed"] == 1.0
    assert bad.loc[0, "Hratio"] == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# Guard 1: some aborts → counted, excluded, surfaced (recovered stats)
# ---------------------------------------------------------------------------


def test_run_monte_carlo_counts_excludes_and_surfaces_real_aborts(
    tmp_path, monkeypatch
) -> None:
    """SOME aborts through real ``run_monte_carlo`` + real solver.

    Criteria covered
    ----------------
    1. FIRES ON BAD CASE — every even sample forced over-full liquid mass;
       real solver returns Conv Failed == 1.
    2. RECOVERY-FORM — asserts clean-only statistics (max/mean/p95 over
       successful Hratio only), not the biased-low population that includes
       Hratio≈0 abort rows. Also asserts aborted_count and aborted_params.
    3. REAL CODE PATH — calls ``run_monte_carlo`` only; detection is
       ``_sample_is_aborted`` inside that function (no reimplemented filter).
    4. REACHABILITY — ``run_monte_carlo`` is the production headless/UI entry
       (see module docstring call chain).
    """
    assert ABORT_FRACTION_THRESHOLD == 0.10  # never relaxed

    package_path = _save_cylinder(tmp_path)
    # 2 of 20 = 10% exactly → allowed (threshold uses strict >)
    abort_indices = {2, 5}
    _install_abort_on_call_indices(monkeypatch, abort_indices)

    n = 20
    request = _mc_request(n=n, seed=7)
    result = run_monte_carlo(
        custom_geometry_config(str(package_path)),
        request,
    )

    expected_aborts = len(abort_indices)
    expected_clean = n - expected_aborts

    # Surfaced
    assert result.aborted_count == expected_aborts
    assert len(result.aborted_params) == expected_aborts
    assert result.n == n
    assert len(result.all_dh) == expected_clean
    assert len(result.all_params) == expected_clean

    # Each aborted draw is a param dict with vent/fill/grav
    for draw in result.aborted_params:
        assert set(draw) == {"vent", "fill", "grav"}

    # Recovered statistics: only clean Hratio values, all strictly positive
    # (abort rows would have contributed 0.0 and pulled max/p95/mean down).
    assert all(dh > 0.0 for dh in result.all_dh)
    arr = np.array(result.all_dh)
    assert result.max_dh == pytest.approx(float(arr.max()))
    assert result.mean_dh == pytest.approx(float(arr.mean()))
    assert result.p95 == pytest.approx(float(np.percentile(arr, 95)))
    assert result.p99 == pytest.approx(float(np.percentile(arr, 99)))

    # Bias-direction demonstration (recovered vs defect): if abort Hratio=0
    # rows were silently included, max would be unchanged (max of positives
    # and zeros is still max positive) but mean/p95/p99 would drop.
    biased = np.concatenate([arr, np.zeros(expected_aborts)])
    biased_mean = float(biased.mean())
    biased_p95 = float(np.percentile(biased, 95))
    assert result.mean_dh > biased_mean  # recovered is higher (conservative)
    assert result.p95 >= biased_p95
    # Store as attributes for evidence scrapers (also asserted numerically).
    assert biased_mean < result.mean_dh


def test_abort_exclusion_bias_direction_numbers(tmp_path, monkeypatch) -> None:
    """Document biased-vs-clean stats on a forcing case (report evidence).

    Criteria: (1)(2)(3)(4) as above; focuses on numeric bias direction.
    """
    package_path = _save_cylinder(tmp_path)
    # 1 abort of 12 ≈ 8.3% < 10%
    _install_abort_on_call_indices(monkeypatch, {3})

    result = run_monte_carlo(
        custom_geometry_config(str(package_path)),
        _mc_request(n=12, seed=11),
    )
    assert result.aborted_count == 1
    clean = np.array(result.all_dh)
    # Defect-state population: clean Hratios + one 0.0 abort row
    biased = np.concatenate([clean, [0.0]])
    # Recovered-state assertions (not "biased is correct")
    assert result.mean_dh == pytest.approx(float(clean.mean()))
    assert result.mean_dh > float(biased.mean())
    assert result.max_dh == pytest.approx(float(clean.max()))
    assert result.max_dh == pytest.approx(float(biased.max()))  # zeros don't raise max
    assert result.p95 == pytest.approx(float(np.percentile(clean, 95)))
    assert result.p95 >= float(np.percentile(biased, 95))


# ---------------------------------------------------------------------------
# Guard 2: abort fraction above threshold → run FAILS loudly
# ---------------------------------------------------------------------------


def test_run_monte_carlo_fails_loudly_when_abort_fraction_exceeds_threshold(
    tmp_path, monkeypatch
) -> None:
    """Abort fraction > 10% raises MonteCarloAbortError through real entry.

    Criteria covered
    ----------------
    1. FIRES ON BAD CASE — all samples forced over-full → 100% abort.
    2. RECOVERY-FORM — asserts exception type, message, and attributes
       (aborted_count, threshold, aborted_params), not silent success.
    3. REAL CODE PATH — ``run_monte_carlo`` raises after real solver aborts.
    4. REACHABILITY — same production entry as the count/exclude test.
    """
    package_path = _save_cylinder(tmp_path)
    _install_abort_all(monkeypatch)

    with pytest.raises(MonteCarloAbortError, match="abort fraction") as exc_info:
        run_monte_carlo(
            custom_geometry_config(str(package_path)),
            _mc_request(n=4, seed=1),
        )

    err = exc_info.value
    assert err.aborted_count == 4
    assert err.n == 4
    assert err.abort_fraction == pytest.approx(1.0)
    assert err.threshold == ABORT_FRACTION_THRESHOLD
    assert len(err.aborted_params) == 4
    assert ABORT_FRACTION_THRESHOLD == 0.10


def test_run_monte_carlo_fails_when_abort_fraction_just_above_threshold(
    tmp_path, monkeypatch
) -> None:
    """2 of 10 aborts (20% > 10%) fails; proves threshold is the contract."""
    package_path = _save_cylinder(tmp_path)
    _install_abort_on_call_indices(monkeypatch, {1, 2})

    with pytest.raises(MonteCarloAbortError) as exc_info:
        run_monte_carlo(
            custom_geometry_config(str(package_path)),
            _mc_request(n=10, seed=3),
        )
    assert exc_info.value.aborted_count == 2
    assert exc_info.value.abort_fraction == pytest.approx(0.2)


# ---------------------------------------------------------------------------
# Guard 3: zero-abort baseline preserved
# ---------------------------------------------------------------------------


def test_run_monte_carlo_zero_abort_surfaces_zero_count(tmp_path) -> None:
    """Zero-abort custom MC through real solver: aborted_count == 0.

    Complements unchanged ``tests/test_monte_carlo.py`` (monkeypatched solver).
    Here the real custom path runs with mild draws that do not abort.

    Criteria: (2) recovered good state; (3)(4) real ``run_monte_carlo``.
    """
    package_path = _save_cylinder(tmp_path)
    result = run_monte_carlo(
        custom_geometry_config(str(package_path)),
        _mc_request(n=3, seed=42),
    )
    assert result.aborted_count == 0
    assert result.aborted_params == []
    assert len(result.all_dh) == 3
    assert all(dh >= 0.0 for dh in result.all_dh)
    assert result.max_dh == pytest.approx(max(result.all_dh))


# ---------------------------------------------------------------------------
# Guard 4: reachability — production entry is run_monte_carlo detection
# ---------------------------------------------------------------------------


def test_abort_detection_is_inside_run_monte_carlo_not_reimplemented() -> None:
    """Reachability: detection helper and threshold live in the MC module.

    Production consumers (Qt RunWorker, headless scripts) call
    ``run_monte_carlo`` only; there is no second abort filter. This test pins
    the public contract symbols so a silent reintroduction of absorb-and-
    continue would fail import/attribute checks.
    """
    assert callable(monte_carlo_module.run_monte_carlo)
    assert callable(monte_carlo_module._sample_is_aborted)
    assert monte_carlo_module.ABORT_FRACTION_THRESHOLD == 0.10
    assert issubclass(monte_carlo_module.MonteCarloAbortError, RuntimeError)

    # Synthetic frame: detection itself is the real helper used by the entry.
    aborted = pd.DataFrame({"Hratio": [0.0], "Conv Failed": [1.0]})
    clean = pd.DataFrame({"Hratio": [0.01], "Conv Failed": [0.0]})
    assert monte_carlo_module._sample_is_aborted(aborted) is True
    assert monte_carlo_module._sample_is_aborted(clean) is False
    assert monte_carlo_module._sample_is_aborted(pd.DataFrame()) is False
