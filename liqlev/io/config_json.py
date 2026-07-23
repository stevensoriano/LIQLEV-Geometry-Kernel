"""Versioned JSON load/save helpers for LIQLEV simulation configs."""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

from liqlev.model.config import (
    EpsilonConfig,
    FluidConfig,
    GravityProfileConfig,
    RunControls,
    SimulationConfig,
    TankConfig,
    VentProfileConfig,
)


CONFIG_SCHEMA_VERSION = 2


def simulation_config_to_dict(config: SimulationConfig) -> dict[str, Any]:
    """Return a JSON-serializable config dictionary with schema version."""
    payload = asdict(config)
    payload["schema_version"] = CONFIG_SCHEMA_VERSION
    return payload


def simulation_config_from_dict(payload: dict[str, Any]) -> SimulationConfig:
    """Build a typed config from a versioned dictionary."""
    version = payload.get("schema_version", 1)
    if type(version) is not int or version not in (1, CONFIG_SCHEMA_VERSION):
        raise ValueError(f"Unsupported config schema_version: {version}")

    fluid_data = dict(payload.get("fluid", {}))
    tank_data = dict(payload.get("tank", {}))
    vent_data = dict(payload.get("vent", {}))
    gravity_data = dict(payload.get("gravity", {}))
    epsilon_data = dict(payload.get("epsilon", {}))
    run_data = dict(payload.get("run", {}))

    if version == 1:
        tank_data["geometry_path"] = ""
    if "fill_fractions" in tank_data:
        tank_data["fill_fractions"] = tuple(tank_data["fill_fractions"])
    if "rates_lbm_s" in vent_data:
        vent_data["rates_lbm_s"] = tuple(vent_data["rates_lbm_s"])
    if "values" in epsilon_data:
        epsilon_data["values"] = tuple(epsilon_data["values"])

    return SimulationConfig(
        fluid=FluidConfig(**fluid_data),
        tank=TankConfig(**tank_data),
        vent=VentProfileConfig(**vent_data),
        gravity=GravityProfileConfig(**gravity_data),
        epsilon=EpsilonConfig(**epsilon_data),
        run=RunControls(**run_data),
        schema_version=version,
    )


def save_simulation_config(config: SimulationConfig, path: str | Path) -> None:
    """Write a versioned simulation config JSON file."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8") as file_obj:
        json.dump(simulation_config_to_dict(config), file_obj, indent=2)


def load_simulation_config(path: str | Path) -> SimulationConfig:
    """Read a versioned simulation config JSON file."""
    with Path(path).open("r", encoding="utf-8") as file_obj:
        payload = json.load(file_obj)
    return simulation_config_from_dict(payload)
