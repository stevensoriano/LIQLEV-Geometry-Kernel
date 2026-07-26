"""F7: configurable boundary-layer RK4 substeps (plan Phase 4.4).

Plumbing (chosen surface)
-------------------------
``RunControls.boundary_layer_substeps`` (default 4) mirrors other numeric run
knobs (``timestep_s``, ``duration_s``). The headless runner threads it through
``build_inputs(..., boundary_layer_substeps=...)`` → inputs key ``BLSubsteps``
→ ``liqlev_simulation`` → ``_solver_loop`` → ``integrate_boundary_layer``.
When the key is absent (legacy ``physics_cases`` path) the solver defaults to 4,
preserving byte-identical baseline behaviour.

Four-criteria map
-----------------
| Guard / invariant                         | C1 bad-case fires | C2 recovered state | C3 real path | C4 reachability |
|-------------------------------------------|-------------------|--------------------|--------------|-----------------|
| Default-4 config round-trip               | n/a (parity)      | default==4         | runner       | ``test_f7_default_substeps_via_runner_inputs`` |
| 4-vs-8 rise/VBL invariance G0 @ dt=0.02   | n/a (parity)      | |Δ| ≤ measured tol | LOX config   | ``test_f7_substeps_4_vs_8_invariance_g0_g3`` |
| 4-vs-8 rise/VBL invariance G3 @ dt=0.02   | n/a (parity)      | |Δ| ≤ measured tol | LOX config   | same (parametrized; both rows mandatory) |
| substeps <= 0 rejected (validation)       | InputValidationError | positive required | validate/runner | ``test_f7_substeps_nonpositive_rejected_via_validation`` |
| substeps <= 0 rejected (solver entry)     | ValueError        | positive required  | liqlev_simulation | ``test_f7_substeps_nonpositive_rejected_via_liqlev_simulation`` |

Measured 4-vs-8 differences (pinned env, 2026-07-25, full LOX 43 L dual-term
blowdown at plateau Δt=0.02; evidence ``s6b_f7_invariance_measure.txt``):

| Row | |Δ rise| mm     | |Δ VBL| ft³     | set abs tol rise_mm | set abs tol VBL |
|-----|----------------:|---------------:|--------------------:|----------------:|
| G0  | 2.21e-09        | 1.41e-11       | 1.0e-05             | 1.0e-07         |
| G3  | 4.39e-06        | 3.81e-08       | 1.0e-05             | 1.0e-07         |

Tolerances sit just above the larger (G3) measured noise.
"""

from __future__ import annotations

from dataclasses import replace

import pytest

from core import liqlev_simulation
from liqlev.model.config import RunControls, SimulationConfig
from liqlev.model.validation import InputValidationError, validate_simulation_config
from liqlev.runner.single import run_single_case
from validation.lox_vent_cases import (
    G0,
    G3,
    build_lox_vent_config,
)

# Just above measured G3 noise (see module docstring).
RISE_MM_ABS_TOL = 1.0e-5
VBL_ABS_TOL = 1.0e-7


def _rise_mm_and_vbl(result) -> tuple[float, float]:
    df = result.dataframe
    h0 = float(df["Height"].iloc[0])
    h1 = float(df["Height"].iloc[-1])
    rise_mm = (h1 - h0) * 304.8
    vbl = float(df["VBL vol"].iloc[-1])
    return rise_mm, vbl


def test_f7_default_substeps_via_runner_inputs() -> None:
    """C2/C3/C4: real runner path emits BLSubsteps=4 by default."""

    config = build_lox_vent_config(G3, timestep_s=0.02, duration_s=0.04)
    assert config.run.boundary_layer_substeps == 4
    result = run_single_case(config)
    assert result.inputs["BLSubsteps"] == 4
    assert result.dataframe is not None
    assert len(result.dataframe) >= 1


@pytest.mark.parametrize(
    "gravity_g,row_label",
    [
        (G0, "G0"),
        (G3, "G3"),
    ],
    ids=["G0", "G3"],
)
def test_f7_substeps_4_vs_8_invariance_g0_g3(
    gravity_g: float, row_label: str
) -> None:
    """C2/C3/C4: 4-vs-8 rise/VBL invariance at plateau dt via real config path.

    Both G0 and G3 are mandatory (single-row would be a false pass).
    """

    base = build_lox_vent_config(gravity_g, timestep_s=0.02, duration_s=60.0)
    cfg4 = replace(base, run=replace(base.run, boundary_layer_substeps=4))
    cfg8 = replace(base, run=replace(base.run, boundary_layer_substeps=8))

    r4 = run_single_case(cfg4)
    r8 = run_single_case(cfg8)

    assert r4.inputs["BLSubsteps"] == 4
    assert r8.inputs["BLSubsteps"] == 8

    rise4, vbl4 = _rise_mm_and_vbl(r4)
    rise8, vbl8 = _rise_mm_and_vbl(r8)

    d_rise = abs(rise8 - rise4)
    d_vbl = abs(vbl8 - vbl4)
    assert d_rise <= RISE_MM_ABS_TOL, (
        f"{row_label} |Δrise|={d_rise!r} mm exceeds tol {RISE_MM_ABS_TOL}; "
        f"rise4={rise4!r} rise8={rise8!r}"
    )
    assert d_vbl <= VBL_ABS_TOL, (
        f"{row_label} |ΔVBL|={d_vbl!r} ft3 exceeds tol {VBL_ABS_TOL}; "
        f"vbl4={vbl4!r} vbl8={vbl8!r}"
    )


def test_f7_substeps_nonpositive_rejected_via_validation() -> None:
    """C1/C2/C3/C4: substeps<=0 raises on the real validation/runner entry."""

    config = build_lox_vent_config(G3, timestep_s=0.02, duration_s=1.0)
    bad = replace(config, run=replace(config.run, boundary_layer_substeps=0))

    with pytest.raises(InputValidationError) as exc_info:
        validate_simulation_config(bad)
    assert "run.boundary_layer_substeps" in str(exc_info.value)

    with pytest.raises(InputValidationError):
        run_single_case(bad)


def test_f7_substeps_nonpositive_rejected_via_liqlev_simulation() -> None:
    """C1/C2/C3/C4: BLSubsteps<=0 raises at the live solver entry."""

    config = build_lox_vent_config(G3, timestep_s=0.02, duration_s=0.04)
    result = run_single_case(config)
    inputs = dict(result.inputs)
    inputs["BLSubsteps"] = 0

    with pytest.raises(ValueError, match="BLSubsteps must be a positive integer"):
        liqlev_simulation(inputs, verbose=False)

    inputs["BLSubsteps"] = -1
    with pytest.raises(ValueError, match="BLSubsteps must be a positive integer"):
        liqlev_simulation(inputs, verbose=False)
