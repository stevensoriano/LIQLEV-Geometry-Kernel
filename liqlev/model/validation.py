"""Display-independent validation for LIQLEV simulation configuration."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from liqlev.geometry.package import GeometryPackageError, load_geometry_package
from liqlev.geometry.schema import GeometryKernel

from .builder import safe_eval_gravity
from .config import SimulationConfig


@dataclass(frozen=True)
class ValidationIssue:
    field: str
    message: str


class InputValidationError(ValueError):
    """Raised when a simulation config has one or more field-level errors."""

    def __init__(self, issues: list[ValidationIssue]):
        self.issues = issues
        super().__init__(
            "\n".join(f"{issue.field}: {issue.message}" for issue in issues)
        )


def validate_simulation_config(config: SimulationConfig) -> GeometryKernel | None:
    """Validate config values that affect run setup, without touching physics."""
    issues: list[ValidationIssue] = []
    geometry: GeometryKernel | None = None

    if config.fluid.initial_pressure_psia <= config.fluid.final_pressure_psia:
        issues.append(
            ValidationIssue(
                "fluid.initial_pressure_psia",
                "Initial pressure must exceed final pressure.",
            )
        )

    if config.tank.diameter_ft <= 0:
        issues.append(
            ValidationIssue("tank.diameter_ft", "Tank diameter must be positive.")
        )
    if config.tank.height_ft <= 0:
        issues.append(
            ValidationIssue("tank.height_ft", "Tank height must be positive.")
        )

    geometry_path_value = config.tank.geometry_path
    if not isinstance(geometry_path_value, str):
        issues.append(
            ValidationIssue(
                "tank.geometry_path", "Geometry package path must be a string."
            )
        )
    elif geometry_path_value:
        geometry_path = Path(geometry_path_value)
        metadata_path = geometry_path.with_suffix(".json")
        if geometry_path.suffix.lower() != ".npz":
            issues.append(
                ValidationIssue(
                    "tank.geometry_path",
                    f"Geometry package must be an .npz file: {geometry_path}",
                )
            )
        elif not geometry_path.is_file():
            issues.append(
                ValidationIssue(
                    "tank.geometry_path",
                    f"Geometry package not found: {geometry_path}",
                )
            )
        elif not metadata_path.is_file():
            issues.append(
                ValidationIssue(
                    "tank.geometry_path",
                    f"Geometry metadata not found: {metadata_path}",
                )
            )
        else:
            try:
                geometry = load_geometry_package(geometry_path)
            except GeometryPackageError as exc:
                issues.append(
                    ValidationIssue(
                        "tank.geometry_path", f"Invalid geometry package: {exc}"
                    )
                )

    if not config.tank.fill_fractions:
        issues.append(
            ValidationIssue(
                "tank.fill_fractions", "At least one fill fraction is required."
            )
        )
    else:
        for index, fill in enumerate(config.tank.fill_fractions):
            if fill <= 0 or fill > 1:
                issues.append(
                    ValidationIssue(
                        f"tank.fill_fractions[{index}]",
                        f"Fill fraction {fill} out of range (0, 1].",
                    )
                )

    if config.run.duration_s <= 0:
        issues.append(ValidationIssue("run.duration_s", "Duration must be positive."))
    if config.run.timestep_s <= 0:
        issues.append(ValidationIssue("run.timestep_s", "Time step must be positive."))
    if config.run.timestep_s >= config.run.duration_s:
        issues.append(
            ValidationIssue(
                "run.timestep_s",
                f"Time step ({config.run.timestep_s}s) must be less than duration ({config.run.duration_s}s).",
            )
        )

    if config.vent.csv_path:
        if not Path(config.vent.csv_path).exists():
            issues.append(
                ValidationIssue(
                    "vent.csv_path", f"Vent rate CSV not found: {config.vent.csv_path}"
                )
            )
    else:
        if not config.vent.rates_lbm_s:
            issues.append(
                ValidationIssue(
                    "vent.rates_lbm_s", "At least one vent rate is required."
                )
            )

    if config.vent.ramp_duration_s <= 0:
        issues.append(
            ValidationIssue("vent.ramp_duration_s", "Ramp duration must be positive.")
        )

    if config.epsilon.mode == "Custom" and not config.epsilon.values:
        issues.append(
            ValidationIssue(
                "epsilon.values", "At least one custom epsilon value is required."
            )
        )

    if config.gravity.mode == "Function of Time":
        if not config.gravity.expression.strip():
            issues.append(
                ValidationIssue(
                    "gravity.expression", "Gravity function expression cannot be empty."
                )
            )
        else:
            try:
                safe_eval_gravity(config.gravity.expression, 0.0)
            except Exception as exc:
                issues.append(
                    ValidationIssue(
                        "gravity.expression", f"Gravity function error at t=0: {exc}"
                    )
                )

    if config.gravity.mode == "CSV Profile" and not config.gravity.csv_path:
        issues.append(
            ValidationIssue("gravity.csv_path", "Please select a gravity CSV file.")
        )

    if config.run.threshold_dh_h0 is not None and config.run.threshold_dh_h0 <= 0:
        issues.append(
            ValidationIssue(
                "run.threshold_dh_h0", "Threshold must be positive when provided."
            )
        )

    # F7: BL RK4 substep count must be a positive integer (default 4).
    bl_substeps = config.run.boundary_layer_substeps
    if type(bl_substeps) is not int or isinstance(bl_substeps, bool):
        issues.append(
            ValidationIssue(
                "run.boundary_layer_substeps",
                "Boundary-layer RK4 substeps must be a positive integer.",
            )
        )
    elif bl_substeps <= 0:
        issues.append(
            ValidationIssue(
                "run.boundary_layer_substeps",
                "Boundary-layer RK4 substeps must be a positive integer.",
            )
        )

    if issues:
        raise InputValidationError(issues)
    return geometry
