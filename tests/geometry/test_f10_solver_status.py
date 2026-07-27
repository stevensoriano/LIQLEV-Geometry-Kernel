"""F10: opt-in Solver Status column (plan Phase 4.5 / guard 4.5).

Lead ruling (option a): default OFF, zero baseline impact; a test MUST exercise
the opt-in path so the feature cannot ship inert.

Plumbing (chosen surface)
-------------------------
``RunControls.include_solver_status`` (default False) mirrors other run knobs.
Headless runner threads it through
``build_inputs(..., include_solver_status=...)`` → inputs key
``IncludeSolverStatus`` → ``liqlev_simulation`` (default False if key absent).
``_solver_loop`` always writes an internal 30th column; the public DataFrame
slices to the legacy 29-column ``_COL_NAMES`` contract unless the flag is set.
When ON: 30 columns, name ``'Solver Status'``, dtype float.

Status codes (internal col 29)
------------------------------
0 ok; 1 AK3 non-convergence (bracket exhausted); 2 BL integration failure
(jit status != 0); 3 volume inversion out of domain; 4 BL saturated at A/P
(derived from δ_top vs A/P at the integration height — no change to
``integrate_boundary_layer`` return contract).

Four-criteria map (no-inert)
----------------------------
| Guard / invariant                         | C1 bad-case fires | C2 recovered state | C3 real path | C4 reachability |
|-------------------------------------------|-------------------|--------------------|--------------|-----------------|
| Opt-in ON clean custom short run          | n/a (parity)      | col present, all 0 | liqlev_simulation | ``test_f10_opt_in_clean_custom_all_zeros`` |
| Opt-in ON forced BL fail (zero perimeter) | status==2         | diagnostic row     | liqlev_simulation | ``test_f10_opt_in_zero_perimeter_status_code_2`` |
| Default OFF 29-col byte contract          | n/a (parity)      | names==_COL_NAMES  | liqlev_simulation | ``test_f10_default_off_29_column_contract`` |
| LOX/custom runner can opt in              | n/a (parity)      | key+column present | run_single_case   | ``test_f10_lox_runner_opt_in_reachable`` |

Legacy ``Conv Failed`` (index 28) is untouched. Baseline-checked legacy path
is never wired ON (physics_cases omits the key → default OFF).
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from core import _COL_NAMES, liqlev_simulation
from liqlev.geometry.fixtures import cylinder_kernel
from liqlev.runner.single import run_single_case
from validation.lox_vent_cases import G3, build_lox_vent_config


TINIT_R = 38.3
TIMESTEP_S = 5000.0
RUN_DURATION_S = 2.0 * TIMESTEP_S
GEOMETRY_DIAMETER_FT = 4.0
TANK_HEIGHT_FT = 8.0


def _hydrogen_density(t_rankine: float) -> float:
    return (
        0.1709
        + 0.7454 * t_rankine
        - 0.04421 * t_rankine**2
        + 0.001248 * t_rankine**3
        - 1.738e-5 * t_rankine**4
        + 9.424e-8 * t_rankine**5
    )


def _custom_cylinder_inputs(
    fill: float,
    *,
    gravity_g: float = 0.001,
    liquid_mass_factor: float = 1.0,
    include_solver_status: bool = False,
) -> dict[str, object]:
    """Minimal custom-geometry inputs (mirrors test_core_custom_geometry)."""
    kernel = cylinder_kernel(
        GEOMETRY_DIAMETER_FT,
        TANK_HEIGHT_FT,
        node_count=1025,
    )
    liquid_density = _hydrogen_density(TINIT_R)
    liquid_mass = (
        liquid_density
        * fill
        * kernel.total_volume_ft3
        * liquid_mass_factor
    )
    time = np.array([0.0, RUN_DURATION_S], dtype=np.float64)
    ones = np.ones(2, dtype=np.float64)
    inputs: dict[str, object] = {
        "Liquid": "Hydrogen",
        "Units": "British",
        "Delta": TIMESTEP_S,
        "Dtank": 1.25,
        "Htzero": 0.33,
        "Volt": kernel.total_volume_ft3,
        "Xmlzro": liquid_mass,
        "Pinit": 19.5,
        "Pfinal": 10.0,
        "Tinit": TINIT_R,
        "Neps": 0,
        "Tvmdot": time,
        "Xvmdot": np.full(2, 1.0e-7, dtype=np.float64),
        "Teps": None,
        "Xeps": None,
        "Tspal": time,
        "Xspacl": ones,
        "Tspav": time,
        "Xspacv": ones,
        "Tggo": time,
        "Xggo": np.full(2, gravity_g * 32.174, dtype=np.float64),
        "GeometryMode": 1,
        "FillFraction": fill,
        "GeomHeight": kernel.height_ft,
        "GeomVolume": kernel.volume_ft3,
        "GeomVolumeCoefficients": kernel.volume_coefficients,
        "GeomAreaSamples": kernel.section_area_ft2,
        "GeomPerimeter": kernel.perimeter_ft,
        "GeomSidewallArea": kernel.sidewall_area_ft2,
        "GeomSidewallCoefficients": kernel.sidewall_coefficients,
    }
    if include_solver_status:
        inputs["IncludeSolverStatus"] = True
    return inputs


def test_f10_opt_in_clean_custom_all_zeros() -> None:
    """C2/C3/C4: opt-in ON, clean short custom run → column present, all zeros."""

    result = liqlev_simulation(
        _custom_cylinder_inputs(0.5, include_solver_status=True),
        verbose=False,
    )

    assert list(result.columns) == list(_COL_NAMES) + ["Solver Status"]
    assert len(result.columns) == 30
    assert result["Solver Status"].dtype == float or np.issubdtype(
        result["Solver Status"].dtype, np.floating
    )
    assert (result["Solver Status"].to_numpy(dtype=float) == 0.0).all()
    # Legacy Conv Failed untouched and clean.
    assert (result["Conv Failed"].to_numpy(dtype=float) == 0.0).all()


def test_f10_opt_in_zero_perimeter_status_code_2() -> None:
    """C1/C2/C3/C4: zero-perimeter forces BL integration fail → status 2.

    Uses the same diagnostic trigger as
    ``test_boundary_layer_failure_writes_one_bounded_diagnostic_row``, through
    the real ``liqlev_simulation`` entry with IncludeSolverStatus ON.
    """

    inputs = _custom_cylinder_inputs(
        0.5, gravity_g=0.001, include_solver_status=True
    )
    inputs["GeomPerimeter"] = np.zeros_like(
        np.asarray(inputs["GeomPerimeter"], dtype=np.float64)
    )
    result = liqlev_simulation(inputs, verbose=False)

    assert len(result) == 1
    assert "Solver Status" in result.columns
    assert result.loc[0, "Conv Failed"] == 1.0
    assert result.loc[0, "Solver Status"] == 2.0  # BL integration failure


def test_f10_opt_in_volume_inversion_status_code_3() -> None:
    """C1/C2/C3: liquid-mass overflow forces volume inversion fail → status 3."""

    result = liqlev_simulation(
        _custom_cylinder_inputs(
            0.5,
            liquid_mass_factor=2.0,
            include_solver_status=True,
        ),
        verbose=False,
    )

    assert len(result) == 1
    assert result.loc[0, "Conv Failed"] == 1.0
    assert result.loc[0, "Solver Status"] == 3.0


def test_f10_default_off_29_column_contract() -> None:
    """C2/C3/C4: default OFF → exact 29 columns, names byte-identical to contract."""

    # No IncludeSolverStatus key (legacy path shape).
    result = liqlev_simulation(
        _custom_cylinder_inputs(0.5, include_solver_status=False),
        verbose=False,
    )

    assert list(result.columns) == list(_COL_NAMES)
    assert len(result.columns) == 29
    assert "Solver Status" not in result.columns
    # Explicit False via key also keeps 29 columns.
    inputs_false = _custom_cylinder_inputs(0.5)
    inputs_false["IncludeSolverStatus"] = False
    result_false = liqlev_simulation(inputs_false, verbose=False)
    assert list(result_false.columns) == list(_COL_NAMES)


def test_f10_lox_runner_opt_in_reachable() -> None:
    """C2/C3/C4: LOX/custom validation path can opt in via RunControls.

    Wires nothing into the baseline-checked legacy path: only this opt-in
    config sets IncludeSolverStatus. Default LOX config must remain 29-col.
    """

    base = build_lox_vent_config(G3, timestep_s=0.02, duration_s=0.04)
    assert base.run.include_solver_status is False

    # Default OFF through real runner → no status column, no inputs key.
    default_result = run_single_case(base)
    assert "IncludeSolverStatus" not in default_result.inputs
    assert list(default_result.dataframe.columns) == list(_COL_NAMES)
    assert "Solver Status" not in default_result.dataframe.columns

    # Opt-in ON via config → key present, 30th column, finite codes.
    opted = replace(
        base,
        run=replace(base.run, include_solver_status=True),
    )
    on_result = run_single_case(opted)
    assert on_result.inputs.get("IncludeSolverStatus") is True
    assert list(on_result.dataframe.columns) == list(_COL_NAMES) + [
        "Solver Status"
    ]
    status = on_result.dataframe["Solver Status"].to_numpy(dtype=float)
    assert status.shape[0] >= 1
    assert np.all(np.isfinite(status))
    # Short LOX smoke should be clean (codes 0; 4 only if fully saturated).
    assert np.all((status == 0.0) | (status == 4.0))
