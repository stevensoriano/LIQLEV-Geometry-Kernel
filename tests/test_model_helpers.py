"""Tests for extracted non-GUI model helpers."""

from __future__ import annotations

import hashlib
import json

import numpy as np
import pytest

from liqlev.geometry.fixtures import cylinder_kernel
from liqlev.geometry.package import save_geometry_package
from liqlev.model.builder import G_TO_FT_S2, build_inputs, safe_eval_gravity
from liqlev.model.config import FluidConfig, RunControls, SimulationConfig, TankConfig
from liqlev.model.parsing import parse_numeric_array
from liqlev.model.units import convert_from_british, convert_to_british, display_unit
from liqlev.model.validation import InputValidationError, validate_simulation_config
from thermo_utils import build_property_table, sli
from validation.physics_cases import build_case_inputs, epsilon_schedule, get_case


def test_parse_numeric_array_preserves_legacy_forms() -> None:
    assert parse_numeric_array("") == []
    assert parse_numeric_array("0.1, 0.3, 0.5") == [0.1, 0.3, 0.5]
    assert parse_numeric_array("0.1:0.1:0.3") == pytest.approx([0.1, 0.2, 0.3])
    assert parse_numeric_array("linspace(0.1, 0.5, 5)") == pytest.approx(
        [0.1, 0.2, 0.3, 0.4, 0.5]
    )


def test_unit_helpers_preserve_solver_british_contract() -> None:
    assert display_unit("pressure", si_mode=False) == "psia"
    assert display_unit("pressure", si_mode=True) == "bar"
    assert convert_from_british(10.0, "mass_flow", si_mode=False) == 10.0
    assert convert_from_british(10.0, "mass_flow", si_mode=True) == pytest.approx(
        4.53592
    )
    assert convert_to_british(4.53592, "mass_flow", si_mode=True) == pytest.approx(10.0)


def test_build_inputs_matches_validation_case_setup_for_as203() -> None:
    case = get_case("as203_default_high_vent")
    expected = build_case_inputs(case)
    neps, teps, xeps = epsilon_schedule(
        case.epsilon_mode, case.duration_s, case.epsilon_value
    )
    actual = build_inputs(
        fluid=case.fluid,
        pinit_psia=case.pinit_psia,
        pfinal_psia=case.pfinal_psia,
        dtank=case.dtank_ft,
        htank=case.htank_ft,
        fill_fraction=case.fill_fraction,
        duration=case.duration_s,
        delta_t=case.delta_t_s,
        vent_rate=case.vent_rate_lbm_s,
        neps=neps,
        teps=teps,
        xeps=xeps,
        ramp_duration=case.ramp_duration_s or case.duration_s,
        ramp_target_factor=case.ramp_target_factor,
        nggo=2,
        tggo=np.array([0.0, case.duration_s]),
        xggo=np.array([case.gravity_g * G_TO_FT_S2, case.gravity_g * G_TO_FT_S2]),
        xmlzro_override=case.xmlzro_override_lbm,
        tinit_override=case.tinit_override_r,
    )

    assert set(actual) == set(expected)
    for key, value in expected.items():
        if key == "Title":
            continue
        if isinstance(value, np.ndarray):
            assert np.asarray(actual[key]) == pytest.approx(value)
        elif isinstance(value, float):
            assert actual[key] == pytest.approx(value)
        else:
            assert actual[key] == value


def test_lox_builder_density_matches_solver_property_table() -> None:
    """F9: non-hydrogen xmlzro density must match the solver property table.

    Pre-fix the builder used DensitySat while the solver interpolated a
    400-point table; the two sources disagree by the table interpolation
    residual (~9.5e-10 relative for LOX 40 psia). After the fix they share one
    source so initial liquid volume (xmlzro / rhol_table) equals fill·volt.
    """
    fluid = "Oxygen"
    pinit_psia = 40.0
    pfinal_psia = 35.0
    fill_fraction = 0.438286
    dtank = 1.84
    htank = 1.806770
    duration = 60.0

    inputs = build_inputs(
        fluid=fluid,
        pinit_psia=pinit_psia,
        pfinal_psia=pfinal_psia,
        dtank=dtank,
        htank=htank,
        fill_fraction=fill_fraction,
        duration=duration,
        delta_t=0.02,
        vent_rate=0.026212963,
        neps=0,
        teps=None,
        xeps=None,
        ramp_duration=duration,
        ramp_target_factor=1.0,
        nggo=2,
        tggo=np.array([0.0, duration]),
        xggo=np.array([0.0, 0.0]),
    )

    pt_t, pt_rhol, *_ = build_property_table(fluid, pfinal_psia, pinit_psia)
    rhol_table = float(sli(inputs["Tinit"], pt_t, pt_rhol))
    ac = 0.7854 * (dtank**2)
    htzero = fill_fraction * htank
    rhol_builder = float(inputs["Xmlzro"]) / (htzero * ac)

    # Recorded F9 finding number (builder DensitySat vs table, pre-unification).
    # After the fix this residual is exactly zero (same source).
    relative_divergence = abs(rhol_builder - rhol_table) / rhol_table
    assert rhol_builder == pytest.approx(rhol_table, rel=0.0, abs=1e-14), (
        f"F9 LOX density divergence: builder={rhol_builder!r} table={rhol_table!r} "
        f"rel={relative_divergence!r}"
    )
    assert float(inputs["Xmlzro"]) / rhol_table == pytest.approx(
        htzero * ac, rel=0.0, abs=1e-12
    )


def test_safe_gravity_eval_and_validation_errors_are_field_level() -> None:
    assert safe_eval_gravity("0.001 + sin(t)", 0.0) == pytest.approx(0.001)

    bad_config = SimulationConfig(
        fluid=FluidConfig(initial_pressure_psia=10.0, final_pressure_psia=12.0),
        tank=TankConfig(diameter_ft=-1.0, height_ft=28.18, fill_fractions=(1.2,)),
        run=RunControls(duration_s=10.0, timestep_s=10.0),
    )

    with pytest.raises(InputValidationError) as exc_info:
        validate_simulation_config(bad_config)

    fields = {issue.field for issue in exc_info.value.issues}
    assert "fluid.initial_pressure_psia" in fields
    assert "tank.diameter_ft" in fields
    assert "tank.fill_fractions[0]" in fields
    assert "run.timestep_s" in fields


def test_missing_geometry_npz_is_attributed_to_geometry_path(tmp_path) -> None:
    config = SimulationConfig(
        tank=TankConfig(geometry_path=str(tmp_path / "missing.npz"))
    )

    with pytest.raises(InputValidationError) as exc_info:
        validate_simulation_config(config)

    assert "tank.geometry_path" in {
        issue.field for issue in exc_info.value.issues
    }


def test_missing_geometry_metadata_is_attributed_to_geometry_path(tmp_path) -> None:
    package_path = tmp_path / "tank.npz"
    np.savez(package_path, placeholder=np.array([1.0]))
    config = SimulationConfig(tank=TankConfig(geometry_path=str(package_path)))

    with pytest.raises(InputValidationError) as exc_info:
        validate_simulation_config(config)

    assert "tank.geometry_path" in {
        issue.field for issue in exc_info.value.issues
    }


def test_geometry_package_hash_error_is_attributed_to_geometry_path(tmp_path) -> None:
    package_path = tmp_path / "tank.npz"
    save_geometry_package(cylinder_kernel(4.0, 8.0, node_count=17), package_path)
    metadata_path = package_path.with_suffix(".json")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["npz_sha256"] = "0" * 64
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
    config = SimulationConfig(tank=TankConfig(geometry_path=str(package_path)))

    with pytest.raises(InputValidationError) as exc_info:
        validate_simulation_config(config)

    issues = [
        issue
        for issue in exc_info.value.issues
        if issue.field == "tank.geometry_path"
    ]
    assert len(issues) == 1
    assert "SHA-256" in issues[0].message


def test_corrupt_geometry_archive_is_attributed_to_geometry_path(tmp_path) -> None:
    package_path = tmp_path / "tank.npz"
    save_geometry_package(cylinder_kernel(4.0, 8.0, node_count=17), package_path)
    package_path.write_bytes(b"not a valid numpy archive")
    metadata_path = package_path.with_suffix(".json")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["npz_sha256"] = hashlib.sha256(package_path.read_bytes()).hexdigest()
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
    config = SimulationConfig(tank=TankConfig(geometry_path=str(package_path)))

    with pytest.raises(InputValidationError) as exc_info:
        validate_simulation_config(config)

    assert "tank.geometry_path" in {
        issue.field for issue in exc_info.value.issues
    }


def test_non_string_geometry_path_is_attributed_to_geometry_path() -> None:
    config = SimulationConfig(tank=TankConfig(geometry_path=123))

    with pytest.raises(InputValidationError) as exc_info:
        validate_simulation_config(config)

    assert "tank.geometry_path" in {
        issue.field for issue in exc_info.value.issues
    }


@pytest.mark.parametrize(
    "metadata_case",
    [
        "invalid_json",
        "list_root",
        "scalar_root",
        "missing_hash",
        "null_hash",
        "numeric_hash",
        "missing_metadata_field",
    ],
)
def test_malformed_geometry_metadata_is_attributed_to_geometry_path(
    tmp_path,
    metadata_case: str,
) -> None:
    package_path = tmp_path / "tank.npz"
    save_geometry_package(cylinder_kernel(4.0, 8.0, node_count=17), package_path)
    metadata_path = package_path.with_suffix(".json")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))

    if metadata_case == "invalid_json":
        metadata_path.write_text("{", encoding="utf-8")
    elif metadata_case == "list_root":
        metadata_path.write_text(json.dumps([]), encoding="utf-8")
    elif metadata_case == "scalar_root":
        metadata_path.write_text(json.dumps(7), encoding="utf-8")
    elif metadata_case == "missing_hash":
        metadata.pop("npz_sha256")
        metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
    elif metadata_case == "null_hash":
        metadata["npz_sha256"] = None
        metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
    elif metadata_case == "numeric_hash":
        metadata["npz_sha256"] = 7
        metadata_path.write_text(json.dumps(metadata), encoding="utf-8")
    elif metadata_case == "missing_metadata_field":
        metadata.pop("axis")
        metadata_path.write_text(json.dumps(metadata), encoding="utf-8")

    config = SimulationConfig(tank=TankConfig(geometry_path=str(package_path)))

    with pytest.raises(InputValidationError) as exc_info:
        validate_simulation_config(config)

    assert "tank.geometry_path" in {
        issue.field for issue in exc_info.value.issues
    }
