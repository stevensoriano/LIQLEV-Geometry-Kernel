"""Display-independent validation for LIQLEV simulation configuration."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

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


def validate_simulation_config(config: SimulationConfig) -> None:
    """Validate config values that affect run setup, without touching physics."""
    issues: list[ValidationIssue] = []

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

    if issues:
        raise InputValidationError(issues)
