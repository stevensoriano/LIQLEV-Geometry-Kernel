"""Typed configuration objects for headless and GUI simulation setup."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal


FluidName = Literal["Nitrogen", "Oxygen", "Hydrogen", "Methane"]
EpsilonMode = Literal["height_dep", "bulk_fake", "AS-203 Schedule", "Custom"]
GravityMode = Literal["Constant", "Function of Time", "CSV Profile"]


@dataclass(frozen=True)
class FluidConfig:
    name: FluidName = "Hydrogen"
    initial_pressure_psia: float = 19.5
    final_pressure_psia: float = 13.8
    initial_mass_lbm: float | None = None
    initial_temperature_r: float | None = None


@dataclass(frozen=True)
class TankConfig:
    diameter_ft: float = 21.670
    height_ft: float = 28.18
    fill_fractions: tuple[float, ...] = (0.5116,)
    geometry_path: str = ""


@dataclass(frozen=True)
class VentProfileConfig:
    rates_lbm_s: tuple[float, ...] = (3.3069,)
    ramp_duration_s: float = 400.0
    ramp_target_factor: float = 1.0
    csv_path: str = ""


@dataclass(frozen=True)
class GravityProfileConfig:
    mode: GravityMode = "Constant"
    constant_g: float = 0.00000963
    expression: str = ""
    csv_path: str = ""
    hold_g: float = 0.0014


@dataclass(frozen=True)
class EpsilonConfig:
    mode: EpsilonMode = "AS-203 Schedule"
    values: tuple[float, ...] = (0.4,)


@dataclass(frozen=True)
class RunControls:
    duration_s: float = 400.0
    timestep_s: float = 10.0
    threshold_dh_h0: float | None = None
    # Fixed RK4 substeps per boundary-layer interval (custom geometry mode).
    # Default 4 matches the historical hard-coded core.py literal (finding F7).
    boundary_layer_substeps: int = 4
    # F10 / guard 4.5: when True, public DataFrame gains a 30th column
    # 'Solver Status' (float codes 0-4). Default False keeps the exact
    # 29-column baseline contract (check_physics_baseline column equality).
    include_solver_status: bool = False


@dataclass(frozen=True)
class SimulationConfig:
    fluid: FluidConfig = field(default_factory=FluidConfig)
    tank: TankConfig = field(default_factory=TankConfig)
    vent: VentProfileConfig = field(default_factory=VentProfileConfig)
    gravity: GravityProfileConfig = field(default_factory=GravityProfileConfig)
    epsilon: EpsilonConfig = field(default_factory=EpsilonConfig)
    run: RunControls = field(default_factory=RunControls)
    schema_version: int = 2
