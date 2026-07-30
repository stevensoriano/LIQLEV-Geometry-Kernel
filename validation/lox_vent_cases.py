"""LOX 43 L / 40→35 psia / SSM-3 acceptance-case machinery.

Case constants are pinned by
``docs/2026-07-24-lox-vent-test-definition.md`` §1–§3 (and the Step-4
reference). This module provides the config builder, a bounded runner that
returns a summary (including the ullage-closure metric and Phase 4.1
threshold verdict), and an F8-hardened result-manifest writer.

Phase 4.1 wires ``ULLAGE_CLOSURE_RELATIVE_TOLERANCE`` (5%) through
``ullage_mass_is_acceptable`` / ``assess_lox_dataframe_physicality`` into
``LoxVentSummary`` (``ullage_closure_within_tolerance``, ``physical``,
``failure_classifications``). The residual itself is computed only by
``ullage_closure_metric``.

``film_drift_metric`` adds a second F4-class measurement — the wall-film mass
drift — reported alongside the closure residual and deliberately wired to NO
acceptance predicate (``docs/lox-vent-test-definition.md`` §5 item 6: magnitude
reported per run, never assumed small).
"""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping

import numba
import numpy as np
import pandas as pd

from liqlev.model.config import (
    EpsilonConfig,
    FluidConfig,
    GravityProfileConfig,
    RunControls,
    SimulationConfig,
    TankConfig,
    VentProfileConfig,
)
from liqlev.runner.single import SingleCaseResult, run_single_case
from thermo_utils import build_property_table


ROOT = Path(__file__).resolve().parents[1]
GEOMETRY_NPZ_PATH = (
    ROOT / "geometry" / "tables" / "nhq01-m21a-0201_LIQLEV_GEOMETRY.npz"
).resolve()
FLUID_STEP_PATH = (
    ROOT / "geometry" / "output" / "nhq01-m21a-0201_LIQLEV_FLUID_DOMAIN.step"
).resolve()
LOX_VENT_MODULE_PATH = Path(__file__).resolve()
MANIFEST_PATH = (
    ROOT / "validation" / "results" / "lox_vent_manifest.json"
).resolve()

# §1 System
FLUID = "Oxygen"
FILL_FRACTION = 0.438286
INITIAL_HEIGHT_FT = 0.820444  # inverted from 43 L; geometry package authority
TANK_TOTAL_VOLUME_FT3 = 3.464700
TANK_TOTAL_HEIGHT_FT = 1.806770
# Equivalent diameter from max section area (~1.835); probe heritage used 1.84.
TANK_DIAMETER_FT = 1.84

# §2 Vent
PINIT_PSIA = 40.0
PFINAL_PSIA = 35.0
VENT_RATE_LBM_S = 0.026212963
DURATION_S = 60.0
TIMESTEP_S = 0.02

# §3 Gravity matrix (standard-g units; g0 = 9.80665 exactly)
G0 = 0.0
G1 = 2.056571e-05
G2 = 4.113141e-05
G3 = 5.826950e-04
G4 = 1.165390e-03
GRAVITY_MATRIX: dict[str, float] = {
    "G0": G0,
    "G1": G1,
    "G2": G2,
    "G3": G3,
    "G4": G4,
}

# Plan Phase 4.1 / finding F4 — approved 5% ullage-closure acceptance gate.
# Do not relax; G0 (~5.25%) is a real bad case that must fire.
ULLAGE_CLOSURE_RELATIVE_TOLERANCE = 0.05
ULLAGE_MASS_CLOSURE_CLASSIFICATION = "ullage_mass_closure"


@dataclass(frozen=True)
class LoxVentSummary:
    """Compact engineering summary for one LOX vent run."""

    gravity_g: float
    timestep_s: float
    rows: int
    t_end_s: float
    p_end_psia: float
    h_initial_ft: float
    h_final_ft: float
    dh_ft: float
    dh_over_h0: float
    final_vbl_vol_ft3: float
    final_bl_thick_ft: float
    max_ak3: float
    conv_failed_total: int
    ullage_closure_max_relative: float
    finite: bool
    # Phase 4.1 threshold verdict + physicality (wired to the one metric impl).
    ullage_closure_within_tolerance: bool
    physical: bool
    failure_classifications: tuple[str, ...]
    # REPORTED, NOT GATED — wall-film mass drift (:func:`film_drift_metric`).
    # Deliberately wired to no acceptance predicate; see that function.
    film_retention_lbm: float
    film_eos_inventory_lbm: float
    film_drift_lbm: float
    film_drift_over_inventory: float


def sha256_file(path: str | Path) -> str:
    """Return an uppercase SHA-256 digest for a validation artifact."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as file_obj:
        for block in iter(lambda: file_obj.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def build_lox_vent_config(
    gravity_g: float,
    *,
    timestep_s: float = TIMESTEP_S,
    duration_s: float = DURATION_S,
) -> SimulationConfig:
    """Build the §1–§3 LOX vent case at a given steady gravity level."""

    return SimulationConfig(
        fluid=FluidConfig(
            name=FLUID,
            initial_pressure_psia=PINIT_PSIA,
            final_pressure_psia=PFINAL_PSIA,
        ),
        tank=TankConfig(
            diameter_ft=TANK_DIAMETER_FT,
            height_ft=TANK_TOTAL_HEIGHT_FT,
            fill_fractions=(FILL_FRACTION,),
            geometry_path=str(GEOMETRY_NPZ_PATH),
        ),
        vent=VentProfileConfig(
            rates_lbm_s=(VENT_RATE_LBM_S,),
            ramp_duration_s=duration_s,
            ramp_target_factor=1.0,
        ),
        gravity=GravityProfileConfig(
            mode="Constant",
            constant_g=float(gravity_g),
        ),
        epsilon=EpsilonConfig(mode="height_dep"),
        run=RunControls(
            duration_s=float(duration_s),
            timestep_s=float(timestep_s),
        ),
    )


def ullage_closure_metric(dataframe: pd.DataFrame) -> float:
    """Max relative ullage closure residual (the single measurement implementation).

    Pairing is deliberate: ``Ullage Mass[k]`` vs ``Ullage from Calc[k+1]``
    because the volumetric estimate is computed from start-of-step state.
    Phase 4.1 wires the 5% threshold verdict via
    :func:`ullage_mass_is_acceptable` / :func:`assess_lox_dataframe_physicality`
    against this metric — do not reimplement the residual elsewhere.
    """

    if dataframe.empty or len(dataframe) < 2:
        return float("nan")
    if "Ullage Mass" not in dataframe.columns:
        return float("nan")
    if "Ullage from Calc" not in dataframe.columns:
        return float("nan")

    ullage_mass = dataframe["Ullage Mass"].to_numpy(dtype=float)
    ullage_from_calc = dataframe["Ullage from Calc"].to_numpy(dtype=float)
    numerator = np.abs(ullage_mass[:-1] - ullage_from_calc[1:])
    denominator = ullage_from_calc[1:]
    relative = np.full_like(numerator, np.nan, dtype=float)
    positive = np.isfinite(denominator) & (np.abs(denominator) > 0.0)
    relative[positive] = numerator[positive] / np.abs(denominator[positive])
    if not np.any(np.isfinite(relative)):
        return float("nan")
    return float(np.nanmax(relative))


@dataclass(frozen=True)
class FilmDriftMetric:
    """Endpoint wall-film mass drift for one run (lbm, plus the ratio).

    ``film_drift_lbm = film_retention_lbm - film_eos_inventory_lbm`` and
    ``film_drift_over_inventory`` is that drift divided by the inventory (a
    fraction, not a percent). All four are NaN together when the run state
    needed to reconstruct them is unavailable.
    """

    film_retention_lbm: float
    film_eos_inventory_lbm: float
    film_drift_lbm: float
    film_drift_over_inventory: float


_UNAVAILABLE_FILM_DRIFT = FilmDriftMetric(
    film_retention_lbm=float("nan"),
    film_eos_inventory_lbm=float("nan"),
    film_drift_lbm=float("nan"),
    film_drift_over_inventory=float("nan"),
)


@lru_cache(maxsize=8)
def _saturated_vapor_density_table(
    fluid: str, final_pressure_psia: float, initial_pressure_psia: float
) -> tuple[np.ndarray, np.ndarray]:
    """Memoized ``(temperature_R, rho_v)`` slice of the solver's property table.

    This is the identical :func:`thermo_utils.build_property_table` call
    :func:`liqlev.runner.single.build_property_table_for_config` makes, so
    interpolating it reproduces ``core.py:129-130`` instead of inventing a
    saturation curve. Memoized because the build costs seconds of CoolProp
    calls and is a pure function of its three arguments; the returned arrays
    are read-only so no caller can poison the cache.
    """

    temps, _rhol, rhov, *_ = build_property_table(
        fluid, final_pressure_psia, initial_pressure_psia
    )
    temps = np.array(temps, dtype=float)
    rhov = np.array(rhov, dtype=float)
    temps.setflags(write=False)
    rhov.setflags(write=False)
    return temps, rhov


def film_drift_metric(
    dataframe: pd.DataFrame, inputs: Mapping[str, Any] | None = None
) -> FilmDriftMetric:
    """Wall-film mass drift at ``t_end`` (the single measurement implementation).

    REPORTED, NEVER GATED. This mirrors ``ullage_closure_metric`` in discipline
    but is deliberately wired to no acceptance predicate: it is the assumption-6
    instruction of ``docs/lox-vent-test-definition.md`` §5 item 6 ("Its
    magnitude must be reported per run, not assumed small") applied to the wall
    film rather than to the ullage.

    The solver carries the film across steps as a VOLUME (``core.py:567``
    ``vbl1 = vbl2``) while saturated vapour density falls through the blowdown.
    So the cumulative mass the ullage was debited to fill the film,
    ``sum(d vbl * rho_v(t1))`` (the film volume balance ``core.py:306-311``
    evaluated at each step's own density), is NOT the mass the film holds under
    the equation of state, ``vbl(t) * rho_v(t)``. The difference is the drift.

    Reconstruction reference: the closed mass-ledger audit
    ``run_vapor_balance.py`` (run of record
    ``results_vapor_balance/2026-07-29_204715/report.txt``), which balances the
    full ledger to 6e-13 lbm. ``inputs`` is the solver input mapping carried by
    :class:`liqlev.runner.single.SingleCaseResult`; ``Tinit`` supplies the
    start-of-step temperature of step 0 and ``Liquid``/``Pinit``/``Pfinal``
    select the same property table the run interpolated. Without it (or for the
    hydrogen polynomial path, which has no table) every field is NaN.

    Does not mutate ``dataframe``.
    """

    if inputs is None or dataframe.empty:
        return _UNAVAILABLE_FILM_DRIFT
    if "Temp" not in dataframe.columns or "VBL vol" not in dataframe.columns:
        return _UNAVAILABLE_FILM_DRIFT
    fluid = str(inputs.get("Liquid", ""))
    if fluid == "Hydrogen":
        # core.py uses the _hydrogen_props polynomial fit, not a property table.
        return _UNAVAILABLE_FILM_DRIFT
    if not {"Tinit", "Pinit", "Pfinal"} <= set(inputs):
        return _UNAVAILABLE_FILM_DRIFT

    tinit = float(inputs["Tinit"])
    temp = dataframe["Temp"].to_numpy(dtype=float)
    vbl = dataframe["VBL vol"].to_numpy(dtype=float)
    if not (
        np.isfinite(tinit)
        and np.isfinite(temp).all()
        and np.isfinite(vbl).all()
    ):
        return _UNAVAILABLE_FILM_DRIFT

    pt_temps, pt_rhov = _saturated_vapor_density_table(
        fluid, float(inputs["Pfinal"]), float(inputs["Pinit"])
    )
    # Row k is the END of step k, so the start of step k is row k-1 — and for
    # step 0 it is the initial condition (core.py:110 vbl1 = 0.0, :113 t1).
    temp1 = np.concatenate([[tinit], temp[:-1]])
    vbl1 = np.concatenate([[0.0], vbl[:-1]])
    rhov1 = np.interp(temp1, pt_temps, pt_rhov)
    rhov2 = np.interp(temp, pt_temps, pt_rhov)

    retention = float(np.sum((vbl - vbl1) * rhov1))
    inventory = float(vbl[-1] * rhov2[-1])
    drift = retention - inventory
    over_inventory = drift / inventory if inventory != 0.0 else float("nan")
    if not (
        np.isfinite(retention) and np.isfinite(inventory) and np.isfinite(drift)
    ):
        return _UNAVAILABLE_FILM_DRIFT
    return FilmDriftMetric(
        film_retention_lbm=retention,
        film_eos_inventory_lbm=inventory,
        film_drift_lbm=drift,
        film_drift_over_inventory=float(over_inventory),
    )


def ullage_mass_is_acceptable(dataframe: pd.DataFrame) -> bool:
    """Acceptance predicate for ullage mass positivity and 5% closure.

    Requires ``(Ullage Mass > 0).all()`` and
    ``ullage_closure_metric(df) <= ULLAGE_CLOSURE_RELATIVE_TOLERANCE``.
    Shared by the LOX physicality path and NASA ``_dataframe_is_physical``.
    """

    if dataframe.empty:
        return False
    if "Ullage Mass" not in dataframe.columns:
        return False
    ullage = dataframe["Ullage Mass"].to_numpy(dtype=float)
    if not np.isfinite(ullage).all():
        return False
    if not bool((ullage > 0.0).all()):
        return False
    if len(dataframe) < 2:
        # Single-row frames cannot form a closure residual; positivity alone.
        return True
    closure = ullage_closure_metric(dataframe)
    return bool(
        np.isfinite(closure) and closure <= ULLAGE_CLOSURE_RELATIVE_TOLERANCE
    )


def assess_lox_dataframe_physicality(
    dataframe: pd.DataFrame,
    *,
    total_height_ft: float = TANK_TOTAL_HEIGHT_FT,
) -> tuple[bool, tuple[str, ...]]:
    """LOX live-path physicality + classifications (production validation entry).

    Returns ``(physical, failure_classifications)``. Ullage failures surface as
    the distinct classification ``ullage_mass_closure``; other bound/finiteness
    failures use ``physical_bounds``. Drives the same ullage predicate the NASA
    path imports — no reimplementation of the residual metric.
    """

    classifications: list[str] = []
    if dataframe.empty:
        return False, ("physical_bounds", ULLAGE_MASS_CLOSURE_CLASSIFICATION)

    numeric = dataframe.to_numpy(dtype=float)
    base_ok = bool(
        np.isfinite(numeric).all()
        and "Height" in dataframe.columns
        and "VBL vol" in dataframe.columns
        and "BL thick" in dataframe.columns
        and (dataframe["Height"] >= 0.0).all()
        and (dataframe["Height"] <= total_height_ft).all()
        and (dataframe["VBL vol"] >= 0.0).all()
        and (dataframe["BL thick"] >= 0.0).all()
    )
    if not base_ok:
        classifications.append("physical_bounds")

    if not ullage_mass_is_acceptable(dataframe):
        classifications.append(ULLAGE_MASS_CLOSURE_CLASSIFICATION)

    return (len(classifications) == 0, tuple(classifications))


def summarize_lox_vent_result(
    result: SingleCaseResult,
    *,
    gravity_g: float,
    timestep_s: float,
) -> LoxVentSummary:
    """Build the engineering summary from a single-case solver result."""

    dataframe = result.dataframe
    if dataframe.empty:
        physical, classifications = assess_lox_dataframe_physicality(dataframe)
        return LoxVentSummary(
            gravity_g=float(gravity_g),
            timestep_s=float(timestep_s),
            rows=0,
            t_end_s=float("nan"),
            p_end_psia=float("nan"),
            h_initial_ft=float("nan"),
            h_final_ft=float("nan"),
            dh_ft=float("nan"),
            dh_over_h0=float("nan"),
            final_vbl_vol_ft3=float("nan"),
            final_bl_thick_ft=float("nan"),
            max_ak3=float("nan"),
            conv_failed_total=0,
            ullage_closure_max_relative=float("nan"),
            finite=True,
            ullage_closure_within_tolerance=False,
            physical=physical,
            failure_classifications=classifications,
            film_retention_lbm=float("nan"),
            film_eos_inventory_lbm=float("nan"),
            film_drift_lbm=float("nan"),
            film_drift_over_inventory=float("nan"),
        )

    height = dataframe["Height"].to_numpy(dtype=float)
    h0 = float(height[0])
    h1 = float(height[-1])
    dh = h1 - h0
    conv_failed = (
        int(dataframe["Conv Failed"].sum())
        if "Conv Failed" in dataframe.columns
        else 0
    )
    values = dataframe.to_numpy(dtype=float)
    closure = ullage_closure_metric(dataframe)
    # Reported alongside closure and fed to NO acceptance predicate below.
    film = film_drift_metric(dataframe, result.inputs)
    physical, classifications = assess_lox_dataframe_physicality(dataframe)
    # Threshold verdict shares the same metric; do not re-derive residual.
    within = bool(
        np.isfinite(closure) and closure <= ULLAGE_CLOSURE_RELATIVE_TOLERANCE
    )
    return LoxVentSummary(
        gravity_g=float(gravity_g),
        timestep_s=float(timestep_s),
        rows=int(len(dataframe)),
        t_end_s=float(dataframe["Time"].iloc[-1]),
        p_end_psia=float(dataframe["Press"].iloc[-1]),
        h_initial_ft=h0,
        h_final_ft=h1,
        dh_ft=dh,
        dh_over_h0=(dh / h0) if h0 != 0.0 else float("nan"),
        final_vbl_vol_ft3=float(dataframe["VBL vol"].iloc[-1]),
        final_bl_thick_ft=float(dataframe["BL thick"].iloc[-1]),
        max_ak3=float(dataframe["AK3"].max()),
        conv_failed_total=conv_failed,
        ullage_closure_max_relative=closure,
        finite=bool(np.isfinite(values).all()),
        ullage_closure_within_tolerance=within,
        physical=physical,
        failure_classifications=classifications,
        film_retention_lbm=film.film_retention_lbm,
        film_eos_inventory_lbm=film.film_eos_inventory_lbm,
        film_drift_lbm=film.film_drift_lbm,
        film_drift_over_inventory=film.film_drift_over_inventory,
    )


def run_lox_vent_case(
    gravity_g: float,
    *,
    timestep_s: float = TIMESTEP_S,
    duration_s: float = DURATION_S,
) -> tuple[pd.DataFrame, LoxVentSummary]:
    """Run the LOX vent case and return ``(dataframe, summary)``."""

    config = build_lox_vent_config(
        gravity_g, timestep_s=timestep_s, duration_s=duration_s
    )
    result = run_single_case(config)
    summary = summarize_lox_vent_result(
        result, gravity_g=gravity_g, timestep_s=timestep_s
    )
    return result.dataframe, summary


def _git_describe_dirty() -> str:
    return subprocess.check_output(
        ["git", "describe", "--dirty", "--always"],
        cwd=ROOT,
        text=True,
    ).strip()


def _git_worktree_is_dirty() -> bool:
    """Return True when the repository has uncommitted changes."""

    porcelain = subprocess.check_output(
        ["git", "status", "--porcelain"],
        cwd=ROOT,
        text=True,
    )
    return bool(porcelain.strip())


def _finite_or_none(value: float) -> float | None:
    numeric = float(value)
    return numeric if np.isfinite(numeric) else None


def _summary_manifest(summary: LoxVentSummary) -> dict[str, float | int | bool | None]:
    return {
        "gravity_g": summary.gravity_g,
        "timestep_s": summary.timestep_s,
        "rows": summary.rows,
        "t_end_s": _finite_or_none(summary.t_end_s),
        "p_end_psia": _finite_or_none(summary.p_end_psia),
        "h_initial_ft": _finite_or_none(summary.h_initial_ft),
        "h_final_ft": _finite_or_none(summary.h_final_ft),
        "dh_ft": _finite_or_none(summary.dh_ft),
        "dh_over_h0": _finite_or_none(summary.dh_over_h0),
        "final_vbl_vol_ft3": _finite_or_none(summary.final_vbl_vol_ft3),
        "final_bl_thick_ft": _finite_or_none(summary.final_bl_thick_ft),
        "max_ak3": _finite_or_none(summary.max_ak3),
        "conv_failed_total": summary.conv_failed_total,
        "ullage_closure_max_relative": _finite_or_none(
            summary.ullage_closure_max_relative
        ),
        "ullage_closure_within_tolerance": summary.ullage_closure_within_tolerance,
        "ullage_closure_relative_tolerance": ULLAGE_CLOSURE_RELATIVE_TOLERANCE,
        "finite": summary.finite,
        "physical": summary.physical,
        "failure_classifications": list(summary.failure_classifications),
        # Reported wall-film mass drift — no gate, no tolerance key.
        "film_retention_lbm": _finite_or_none(summary.film_retention_lbm),
        "film_eos_inventory_lbm": _finite_or_none(summary.film_eos_inventory_lbm),
        "film_drift_lbm": _finite_or_none(summary.film_drift_lbm),
        "film_drift_over_inventory": _finite_or_none(
            summary.film_drift_over_inventory
        ),
    }


def build_lox_vent_manifest(
    results: Mapping[str, LoxVentSummary],
    *,
    solver_describe: str | None = None,
) -> dict[str, Any]:
    """Build the deterministic LOX vent result-manifest payload (F8 hardened)."""

    if solver_describe is None:
        solver_describe = _git_describe_dirty()
    return {
        "schema": "liqlev.validation.lox_vent",
        "version": 1,
        "geometry_npz": str(GEOMETRY_NPZ_PATH.relative_to(ROOT)).replace(
            "\\", "/"
        ),
        "geometry_npz_sha256": sha256_file(GEOMETRY_NPZ_PATH),
        "fluid_step": str(FLUID_STEP_PATH.relative_to(ROOT)).replace("\\", "/"),
        "fluid_step_sha256": sha256_file(FLUID_STEP_PATH),
        "harness_module": str(LOX_VENT_MODULE_PATH.relative_to(ROOT)).replace(
            "\\", "/"
        ),
        "harness_module_sha256": sha256_file(LOX_VENT_MODULE_PATH),
        "solver_describe": solver_describe,
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "numba": numba.__version__,
        },
        "case_definition": {
            "source": "2026-07-24-lox-vent-test-definition.md",
            "fluid": FLUID,
            "fill_fraction": FILL_FRACTION,
            "initial_height_ft": INITIAL_HEIGHT_FT,
            "initial_pressure_psia": PINIT_PSIA,
            "final_pressure_psia": PFINAL_PSIA,
            "vent_rate_lbm_s": VENT_RATE_LBM_S,
            "duration_s": DURATION_S,
            "timestep_s_primary": TIMESTEP_S,
            "epsilon": "height_dep",
            "gravity_matrix_g": dict(GRAVITY_MATRIX),
        },
        "results": {
            key: _summary_manifest(summary) for key, summary in results.items()
        },
    }


def write_lox_vent_manifest(
    results: Mapping[str, LoxVentSummary],
    path: str | Path | None = None,
    *,
    solver_describe: str | None = None,
) -> dict[str, Any]:
    """Write the LOX vent manifest, refusing a dirty git tree (F8).

    The production path is ``validation/results/lox_vent_manifest.json``; run
    engineers write it later. This function is the machinery — tests use
    ``tmp_path`` and must never commit a scratch manifest under results/.
    """

    if _git_worktree_is_dirty():
        raise RuntimeError(
            "Refusing to write LOX vent manifest from a dirty git worktree "
            "(F8 provenance guard). Commit or stash changes first."
        )

    target = Path(path) if path is not None else MANIFEST_PATH
    payload = build_lox_vent_manifest(results, solver_describe=solver_describe)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return payload
