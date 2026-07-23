"""Smoke tests for solver-critical imports and validation case setup."""

from __future__ import annotations

import core
import thermo_utils
from validation.physics_cases import build_case_inputs, get_case


def test_solver_modules_and_validation_case_builder_import() -> None:
    assert callable(core.liqlev_simulation)
    assert callable(thermo_utils.build_property_table)

    case = get_case("as203_default_high_vent")
    inputs = build_case_inputs(case)

    assert inputs["Liquid"] == "Hydrogen"
    assert inputs["Pinit"] == case.pinit_psia
    assert inputs["Pfinal"] == case.pfinal_psia
    assert inputs["Nvmd"] == len(inputs["Tvmdot"])
