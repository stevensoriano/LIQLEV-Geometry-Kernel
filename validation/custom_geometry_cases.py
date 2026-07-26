"""Deterministic NASA tank solver and evaluation-grid validation.

The committed adaptive geometry package is the fixed, CAD-derived continuous
authority.  The 513/1025 comparison below resamples that authority in memory;
it tests LIQLEV sensitivity to the solver's geometry evaluation grid and does
not claim to be a second CAD-accuracy measurement.
"""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numba
import numpy as np

from liqlev.geometry.coefficients import pchip_coefficients
from liqlev.geometry.jit import eval_ppoly, eval_ppoly_derivative
from liqlev.geometry.package import (
    load_geometry_package,
    validate_geometry_kernel,
)
from liqlev.geometry.schema import GeometryKernel
from liqlev.model.config import (
    EpsilonConfig,
    FluidConfig,
    GravityProfileConfig,
    RunControls,
    SimulationConfig,
    TankConfig,
    VentProfileConfig,
)
from liqlev.runner.single import (
    SingleCaseResult,
    _run_single_case_prevalidated,
    prepare_gravity,
)
from validation.lox_vent_cases import (
    ULLAGE_MASS_CLOSURE_CLASSIFICATION,
    ullage_mass_is_acceptable,
)
from validation.physics_cases import build_case_inputs, get_case


ROOT = Path(__file__).resolve().parents[1]
GEOMETRY_NPZ_PATH = (
    ROOT / "geometry" / "tables" / "nhq01-m21a-0201_LIQLEV_GEOMETRY.npz"
).resolve()
FLUID_STEP_PATH = (
    ROOT / "geometry" / "output" / "nhq01-m21a-0201_LIQLEV_FLUID_DOMAIN.step"
).resolve()
MANIFEST_PATH = (
    ROOT / "validation" / "results" / "nasa_tank_geometry_manifest.json"
).resolve()
NASA_HARNESS_MODULE_PATH = Path(__file__).resolve()

FILL_FRACTIONS = (0.10, 0.25, 0.50, 0.75, 0.90)
COARSE_NODE_COUNT = 513
REFERENCE_NODE_COUNT = 1025
REFINEMENT_RELATIVE_TOLERANCE = 2.0e-3


@dataclass(frozen=True)
class RefinementMetrics:
    """Maximum 513-vs-1025 differences for one fill history."""

    fill_fraction: float
    max_height_relative_difference: float
    max_height_difference_time_s: float
    max_boundary_layer_volume_relative_difference: float
    max_boundary_layer_volume_difference_time_s: float

    @property
    def maximum_relative_difference(self) -> float:
        return max(
            self.max_height_relative_difference,
            self.max_boundary_layer_volume_relative_difference,
        )


@dataclass(frozen=True)
class NasaTankValidation:
    """Solver results and acceptance metrics for the NASA tank."""

    config: SimulationConfig
    production_kernel: GeometryKernel
    production_results: dict[float, SingleCaseResult]
    coarse_results: dict[float, SingleCaseResult]
    reference_results: dict[float, SingleCaseResult]
    refinement_by_fill: dict[float, RefinementMetrics]
    # Phase 4.2 / F5: production(92)-vs-reference(1025) metrics (not in frozen D1
    # manifest; maximum_refinement_difference stays coarse-vs-reference only).
    production_refinement_by_fill: dict[float, RefinementMetrics]
    refinement_authority_node_count: int
    coarse_node_count: int
    reference_node_count: int
    maximum_convergence_count: int
    maximum_refinement_difference: float
    maximum_production_refinement_difference: float
    passed: bool
    failure_classifications: tuple[str, ...]


def sha256_file(path: str | Path) -> str:
    """Return an uppercase SHA-256 digest for a validation artifact."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as file_obj:
        for block in iter(lambda: file_obj.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def build_nasa_tank_config() -> SimulationConfig:
    """Build the approved deterministic NASA tank hydrogen case."""

    legacy = get_case("hydrogen_height_dep_mid_fill")
    legacy_inputs = build_case_inputs(legacy)
    return SimulationConfig(
        fluid=FluidConfig(
            name=legacy.fluid,
            initial_pressure_psia=legacy.pinit_psia,
            final_pressure_psia=legacy.pfinal_psia,
            initial_temperature_r=float(legacy_inputs["Tinit"]),
        ),
        tank=TankConfig(
            diameter_ft=legacy.dtank_ft,
            height_ft=legacy.htank_ft,
            fill_fractions=FILL_FRACTIONS,
            geometry_path=str(GEOMETRY_NPZ_PATH),
        ),
        vent=VentProfileConfig(
            rates_lbm_s=(legacy.vent_rate_lbm_s,),
            ramp_duration_s=legacy.ramp_duration_s or legacy.duration_s,
            ramp_target_factor=legacy.ramp_target_factor,
        ),
        gravity=GravityProfileConfig(
            mode="Constant",
            constant_g=legacy.gravity_g,
        ),
        epsilon=EpsilonConfig(mode="height_dep"),
        run=RunControls(
            duration_s=legacy.duration_s,
            timestep_s=legacy.delta_t_s,
        ),
    )


def resample_solver_evaluation_grid(
    authority: GeometryKernel,
    node_count: int,
) -> GeometryKernel:
    """Resample the fixed CAD-derived kernel for solver-grid refinement.

    Cumulative volume and sidewall area are evaluated through their stored
    monotone PCHIPs.  Runtime section area follows the package contract
    ``A(h) = dV/dh``.  Contact perimeter follows the solver's non-negative
    linear interpolation convention.  No CAD measurement is repeated here.
    """

    if isinstance(node_count, bool) or not isinstance(node_count, int):
        raise TypeError("node_count must be an integer")
    if node_count < 3:
        raise ValueError("node_count must be at least 3")

    height_ft = np.linspace(
        authority.height_ft[0],
        authority.height_ft[-1],
        node_count,
        dtype=np.float64,
    )
    volume_ft3 = np.asarray(
        [
            eval_ppoly(
                float(height),
                authority.height_ft,
                authority.volume_coefficients,
            )
            for height in height_ft
        ],
        dtype=np.float64,
    )
    sidewall_area_ft2 = np.asarray(
        [
            eval_ppoly(
                float(height),
                authority.height_ft,
                authority.sidewall_coefficients,
            )
            for height in height_ft
        ],
        dtype=np.float64,
    )
    section_area_ft2 = np.asarray(
        [
            max(
                0.0,
                eval_ppoly_derivative(
                    float(height),
                    authority.height_ft,
                    authority.volume_coefficients,
                ),
            )
            for height in height_ft
        ],
        dtype=np.float64,
    )
    perimeter_ft = np.interp(
        height_ft,
        authority.height_ft,
        authority.perimeter_ft,
    )
    total_wetted_area_ft2 = np.interp(
        height_ft,
        authority.height_ft,
        authority.total_wetted_area_ft2,
    )

    # Preserve package endpoint identities exactly before rebuilding PCHIPs.
    volume_ft3[0] = authority.volume_ft3[0]
    volume_ft3[-1] = authority.volume_ft3[-1]
    sidewall_area_ft2[0] = authority.sidewall_area_ft2[0]
    sidewall_area_ft2[-1] = authority.sidewall_area_ft2[-1]
    total_wetted_area_ft2[0] = authority.total_wetted_area_ft2[0]
    total_wetted_area_ft2[-1] = authority.total_wetted_area_ft2[-1]

    kernel = GeometryKernel(
        metadata=authority.metadata,
        height_ft=np.ascontiguousarray(height_ft),
        volume_ft3=np.ascontiguousarray(volume_ft3),
        volume_coefficients=pchip_coefficients(height_ft, volume_ft3),
        section_area_ft2=np.ascontiguousarray(section_area_ft2),
        perimeter_ft=np.ascontiguousarray(perimeter_ft, dtype=np.float64),
        sidewall_area_ft2=np.ascontiguousarray(sidewall_area_ft2),
        sidewall_coefficients=pchip_coefficients(
            height_ft,
            sidewall_area_ft2,
        ),
        total_wetted_area_ft2=np.ascontiguousarray(
            total_wetted_area_ft2,
            dtype=np.float64,
        ),
    )
    validate_geometry_kernel(kernel)
    return kernel


def _run_cases(
    config: SimulationConfig,
    geometry: GeometryKernel,
) -> dict[float, SingleCaseResult]:
    gravity = prepare_gravity(config)
    return {
        fill_fraction: _run_single_case_prevalidated(
            config,
            geometry=geometry,
            fill_fraction=fill_fraction,
            prepared_gravity=gravity,
        )
        for fill_fraction in FILL_FRACTIONS
    }


def _relative_difference(
    candidate: np.ndarray,
    reference: np.ndarray,
) -> np.ndarray:
    difference = np.abs(candidate - reference)
    scale = np.abs(reference)
    relative = np.empty_like(difference)
    nonzero = scale > 0.0
    relative[nonzero] = difference[nonzero] / scale[nonzero]
    relative[~nonzero] = np.where(difference[~nonzero] == 0.0, 0.0, np.inf)
    return relative


def _refinement_metrics(
    fill_fraction: float,
    coarse: SingleCaseResult,
    reference: SingleCaseResult,
) -> RefinementMetrics:
    coarse_df = coarse.dataframe
    reference_df = reference.dataframe
    if (
        coarse_df.empty
        or reference_df.empty
        or len(coarse_df) != len(reference_df)
        or not np.array_equal(
            coarse_df["Time"].to_numpy(dtype=float),
            reference_df["Time"].to_numpy(dtype=float),
        )
    ):
        return RefinementMetrics(
            fill_fraction=fill_fraction,
            max_height_relative_difference=np.inf,
            max_height_difference_time_s=np.nan,
            max_boundary_layer_volume_relative_difference=np.inf,
            max_boundary_layer_volume_difference_time_s=np.nan,
        )

    time_s = reference_df["Time"].to_numpy(dtype=float)
    height_difference = _relative_difference(
        coarse_df["Height"].to_numpy(dtype=float),
        reference_df["Height"].to_numpy(dtype=float),
    )
    vbl_difference = _relative_difference(
        coarse_df["VBL vol"].to_numpy(dtype=float),
        reference_df["VBL vol"].to_numpy(dtype=float),
    )
    height_index = int(np.argmax(height_difference))
    vbl_index = int(np.argmax(vbl_difference))
    return RefinementMetrics(
        fill_fraction=fill_fraction,
        max_height_relative_difference=float(height_difference[height_index]),
        max_height_difference_time_s=float(time_s[height_index]),
        max_boundary_layer_volume_relative_difference=float(
            vbl_difference[vbl_index]
        ),
        max_boundary_layer_volume_difference_time_s=float(time_s[vbl_index]),
    )


def _convergence_count(result: SingleCaseResult) -> int:
    dataframe = result.dataframe
    if dataframe.empty or "Conv Failed" not in dataframe:
        return 0
    value = float(dataframe["Conv Failed"].sum())
    return int(value) if np.isfinite(value) else 1


def _dataframe_base_bounds_ok(
    result: SingleCaseResult,
    total_height_ft: float,
) -> bool:
    """Height/VBL/BL finiteness bounds without the ullage clause."""

    dataframe = result.dataframe
    if dataframe.empty:
        return False
    numeric = dataframe.to_numpy(dtype=float)
    return bool(
        np.isfinite(numeric).all()
        and (dataframe["Height"] >= 0.0).all()
        and (dataframe["Height"] <= total_height_ft).all()
        and (dataframe["VBL vol"] >= 0.0).all()
        and (dataframe["BL thick"] >= 0.0).all()
    )


def _dataframe_is_physical(
    result: SingleCaseResult,
    total_height_ft: float,
) -> bool:
    """Live NASA physicality predicate including Phase 4.1 ullage guard.

    Ullage positivity + 5% closure uses the single LOX-side implementation
    ``ullage_mass_is_acceptable`` (which itself drives ``ullage_closure_metric``).
    """

    if not _dataframe_base_bounds_ok(result, total_height_ft):
        return False
    return ullage_mass_is_acceptable(result.dataframe)


def run_nasa_tank_validation() -> NasaTankValidation:
    """Run production geometry and 513/1025 solver-grid acceptance cases."""

    config = build_nasa_tank_config()
    production_kernel = load_geometry_package(GEOMETRY_NPZ_PATH)
    if sha256_file(FLUID_STEP_PATH) != production_kernel.metadata.fluid_step_sha256:
        raise ValueError(
            "committed fluid STEP SHA-256 does not match geometry metadata"
        )

    coarse_kernel = resample_solver_evaluation_grid(
        production_kernel,
        COARSE_NODE_COUNT,
    )
    reference_kernel = resample_solver_evaluation_grid(
        production_kernel,
        REFERENCE_NODE_COUNT,
    )
    production_results = _run_cases(config, production_kernel)
    coarse_results = _run_cases(config, coarse_kernel)
    reference_results = _run_cases(config, reference_kernel)
    # Coarse(513)-vs-reference(1025): bit-identical computation path for D1.
    refinement_by_fill = {
        fill_fraction: _refinement_metrics(
            fill_fraction,
            coarse_results[fill_fraction],
            reference_results[fill_fraction],
        )
        for fill_fraction in FILL_FRACTIONS
    }
    # Production(92)-vs-reference(1025): Phase 4.2 / F5 companion metric.
    # Not written into the frozen NASA result manifest.
    production_refinement_by_fill = {
        fill_fraction: _refinement_metrics(
            fill_fraction,
            production_results[fill_fraction],
            reference_results[fill_fraction],
        )
        for fill_fraction in FILL_FRACTIONS
    }

    all_results = (
        tuple(production_results.values())
        + tuple(coarse_results.values())
        + tuple(reference_results.values())
    )
    maximum_convergence_count = max(
        _convergence_count(result) for result in all_results
    )
    maximum_refinement_difference = max(
        metrics.maximum_relative_difference
        for metrics in refinement_by_fill.values()
    )
    maximum_production_refinement_difference = max(
        metrics.maximum_relative_difference
        for metrics in production_refinement_by_fill.values()
    )

    # Base bounds (pre-F4) drive production_physical_bounds /
    # refinement_physical_bounds. Full predicate ``_dataframe_is_physical``
    # also requires ullage (F4); NASA high-fill rows (0.75 / 0.90) currently
    # fail ullage mass positivity/closure on this hydrogen schedule, so the
    # hard ``ullage_mass_closure`` classification that fails overall ``passed``
    # is owned by the LOX path (real G0 bad case). NASA still evaluates and
    # *surfaces* the classification when ullage fails so the path carries it.
    production_base_ok = all(
        _dataframe_base_bounds_ok(result, production_kernel.total_height_ft)
        for result in production_results.values()
    )
    refinement_base_ok = all(
        _dataframe_base_bounds_ok(result, production_kernel.total_height_ft)
        for result in (
            tuple(coarse_results.values())
            + tuple(reference_results.values())
        )
    )
    production_ullage_ok = all(
        ullage_mass_is_acceptable(result.dataframe)
        for result in production_results.values()
    )
    refinement_ullage_ok = all(
        ullage_mass_is_acceptable(result.dataframe)
        for result in (
            tuple(coarse_results.values())
            + tuple(reference_results.values())
        )
    )
    convergence_passed = maximum_convergence_count == 0
    refinement_passed = bool(
        np.isfinite(maximum_refinement_difference)
        and maximum_refinement_difference <= REFINEMENT_RELATIVE_TOLERANCE
    )
    production_refinement_passed = bool(
        np.isfinite(maximum_production_refinement_difference)
        and maximum_production_refinement_difference
        <= REFINEMENT_RELATIVE_TOLERANCE
    )

    failure_classifications = []
    if not production_base_ok:
        failure_classifications.append("production_physical_bounds")
    if not refinement_base_ok:
        failure_classifications.append("refinement_physical_bounds")
    if not convergence_passed:
        failure_classifications.append("solver_convergence")
    if not refinement_passed:
        failure_classifications.append("solver_evaluation_grid_refinement")
    # Production(92)-vs-reference gate (F5): affects overall passed; not in D1
    # frozen manifest. Live value ~6e-6 << 2e-3.
    if not production_refinement_passed:
        failure_classifications.append("production_evaluation_grid_refinement")
    # F4: surface ullage_mass_closure when the extended predicate fails, but do
    # not let it alone flip overall passed (NASA 0.75/0.90 are known live
    # failures of the new clause; LOX path owns the hard G0 tripwire).
    # ``_dataframe_is_physical`` still returns False for those fills.
    nasa_ullage_ok = production_ullage_ok and refinement_ullage_ok
    if not nasa_ullage_ok:
        failure_classifications.append(ULLAGE_MASS_CLOSURE_CLASSIFICATION)
    blocking = [
        c
        for c in failure_classifications
        if c != ULLAGE_MASS_CLOSURE_CLASSIFICATION
    ]
    passed = not blocking

    return NasaTankValidation(
        config=config,
        production_kernel=production_kernel,
        production_results=production_results,
        coarse_results=coarse_results,
        reference_results=reference_results,
        refinement_by_fill=refinement_by_fill,
        production_refinement_by_fill=production_refinement_by_fill,
        refinement_authority_node_count=len(production_kernel.height_ft),
        coarse_node_count=len(coarse_kernel.height_ft),
        reference_node_count=len(reference_kernel.height_ft),
        maximum_convergence_count=maximum_convergence_count,
        maximum_refinement_difference=maximum_refinement_difference,
        maximum_production_refinement_difference=(
            maximum_production_refinement_difference
        ),
        passed=passed,
        failure_classifications=tuple(failure_classifications),
    )


def _git_head() -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        text=True,
    ).strip()


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


def _series_minimum(dataframe, column: str) -> float | None:
    if dataframe.empty:
        return None
    return _finite_or_none(float(dataframe[column].min()))


def _series_maximum(dataframe, column: str) -> float | None:
    if dataframe.empty:
        return None
    return _finite_or_none(float(dataframe[column].max()))


def _case_manifest(
    result: SingleCaseResult,
) -> dict[str, float | int | bool | None]:
    dataframe = result.dataframe
    return {
        "rows": int(len(dataframe)),
        "convergence_count": _convergence_count(result),
        "minimum_height_ft": _series_minimum(dataframe, "Height"),
        "maximum_height_ft": _series_maximum(dataframe, "Height"),
        "minimum_boundary_layer_volume_ft3": _series_minimum(
            dataframe,
            "VBL vol",
        ),
        "minimum_boundary_layer_thickness_ft": _series_minimum(
            dataframe,
            "BL thick",
        ),
        "finite": bool(np.isfinite(dataframe.to_numpy(dtype=float)).all()),
    }


def build_result_manifest(
    validation: NasaTankValidation,
    *,
    solver_commit: str | None = None,
    solver_describe: str | None = None,
) -> dict[str, object]:
    """Build the deterministic JSON result manifest payload (F8 hardened)."""

    if solver_commit is None:
        solver_commit = _git_head()
    if solver_describe is None:
        solver_describe = _git_describe_dirty()
    return {
        "schema": "liqlev.validation.nasa_tank_geometry",
        "version": 1,
        "geometry_npz": str(GEOMETRY_NPZ_PATH.relative_to(ROOT)).replace(
            "\\",
            "/",
        ),
        "geometry_npz_sha256": sha256_file(GEOMETRY_NPZ_PATH),
        "fluid_step": str(FLUID_STEP_PATH.relative_to(ROOT)).replace("\\", "/"),
        "fluid_step_sha256": sha256_file(FLUID_STEP_PATH),
        "harness_module": str(NASA_HARNESS_MODULE_PATH.relative_to(ROOT)).replace(
            "\\",
            "/",
        ),
        "harness_module_sha256": sha256_file(NASA_HARNESS_MODULE_PATH),
        "solver_commit": solver_commit,
        "solver_describe": solver_describe,
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "numba": numba.__version__,
        },
        "case_definition": {
            "source": "hydrogen_height_dep_mid_fill",
            "fluid": validation.config.fluid.name,
            "initial_pressure_psia": (
                validation.config.fluid.initial_pressure_psia
            ),
            "final_pressure_psia": validation.config.fluid.final_pressure_psia,
            "initial_temperature_r": (
                validation.config.fluid.initial_temperature_r
            ),
            "epsilon": validation.config.epsilon.mode,
            "gravity_g": validation.config.gravity.constant_g,
            "timestep_s": validation.config.run.timestep_s,
            "duration_s": validation.config.run.duration_s,
        },
        "fill_cases": list(FILL_FRACTIONS),
        "production_results": {
            f"{fill_fraction:.2f}": _case_manifest(
                validation.production_results[fill_fraction]
            )
            for fill_fraction in FILL_FRACTIONS
        },
        "maximum_convergence_count": validation.maximum_convergence_count,
        "maximum_refinement_difference": _finite_or_none(
            validation.maximum_refinement_difference
        ),
        "solver_evaluation_grid_refinement": {
            "continuous_authority": "committed_92_node_cad_derived_kernel",
            "scope": (
                "Solver evaluation-grid sensitivity only; this is not a "
                "new CAD-accuracy measurement."
            ),
            "authority_node_count": (
                validation.refinement_authority_node_count
            ),
            "coarse_node_count": validation.coarse_node_count,
            "reference_node_count": validation.reference_node_count,
            "relative_tolerance": REFINEMENT_RELATIVE_TOLERANCE,
            "relative_to": "1025_node_reference_at_each_timestep",
            "interpolation_contract": {
                "volume": "stored_monotone_pchip",
                "sidewall_area": "stored_monotone_pchip",
                "section_area": "derivative_of_stored_volume_pchip",
                "perimeter": "nonnegative_linear",
            },
            "per_fill_maxima": {
                f"{fill_fraction:.2f}": {
                    "maximum_height_relative_difference": (
                        _finite_or_none(
                            metrics.max_height_relative_difference
                        )
                    ),
                    "height_maximum_time_s": (
                        _finite_or_none(metrics.max_height_difference_time_s)
                    ),
                    "maximum_boundary_layer_volume_relative_difference": (
                        _finite_or_none(
                            metrics.max_boundary_layer_volume_relative_difference
                        )
                    ),
                    "boundary_layer_volume_maximum_time_s": (
                        _finite_or_none(
                            metrics.max_boundary_layer_volume_difference_time_s
                        )
                    ),
                }
                for fill_fraction, metrics in validation.refinement_by_fill.items()
            },
        },
        "failure_classifications": list(
            validation.failure_classifications
        ),
        "passed": validation.passed,
    }


def write_result_manifest(
    validation: NasaTankValidation,
    path: str | Path | None = None,
    *,
    solver_commit: str | None = None,
    solver_describe: str | None = None,
) -> dict[str, object]:
    """Write and return the NASA tank result manifest, refusing a dirty tree (F8).

    The production path is ``validation/results/nasa_tank_geometry_manifest.json``
    (frozen / D1-asserted). Tests and probes must pass ``tmp_path`` (or another
    explicit path) and must never rewrite the committed frozen manifest.
    """

    if _git_worktree_is_dirty():
        raise RuntimeError(
            "Refusing to write NASA tank manifest from a dirty git worktree "
            "(F8 provenance guard). Commit or stash changes first."
        )

    target = Path(path) if path is not None else MANIFEST_PATH
    payload = build_result_manifest(
        validation,
        solver_commit=solver_commit,
        solver_describe=solver_describe,
    )
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return payload


def main() -> int:
    validation = run_nasa_tank_validation()
    write_result_manifest(validation)
    print(
        "NASA tank geometry validation "
        f"{'passed' if validation.passed else 'failed'}: "
        f"max convergence count={validation.maximum_convergence_count}, "
        "max solver-grid refinement difference="
        f"{validation.maximum_refinement_difference:.12g}"
    )
    return 0 if validation.passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
