from __future__ import annotations

import json

import numpy as np
import pytest

from liqlev.geometry.package import (
    GeometryPackageError,
    load_geometry_package,
    save_geometry_package,
    validate_geometry_kernel,
)
from liqlev.geometry.schema import GeometryKernel, GeometryMetadata


def valid_kernel() -> GeometryKernel:
    h = np.array([0.0, 1.0, 2.0], dtype=np.float64)
    volume = np.array([0.0, 3.0, 6.0], dtype=np.float64)
    area = np.array([3.0, 3.0, 3.0], dtype=np.float64)
    perimeter = np.array([2.0, 2.0, 2.0], dtype=np.float64)
    wall = np.array([0.0, 2.0, 4.0], dtype=np.float64)
    coeff = np.array(
        [[0.0, 0.0], [0.0, 0.0], [3.0, 3.0], [0.0, 3.0]],
        dtype=np.float64,
    )
    wall_coeff = np.array(
        [[0.0, 0.0], [0.0, 0.0], [2.0, 2.0], [0.0, 2.0]],
        dtype=np.float64,
    )
    return GeometryKernel(
        metadata=GeometryMetadata(
            schema_version=1,
            geometry_id="unit-cylinder",
            source_step_sha256="0" * 64,
            fluid_step_sha256="1" * 64,
            axis="+Y",
            gravity_direction="-Y",
            length_unit="ft",
            area_unit="ft^2",
            volume_unit="ft^3",
            y_min_mm=0.0,
            y_max_mm=609.6,
        ),
        height_ft=h,
        volume_ft3=volume,
        volume_coefficients=coeff,
        section_area_ft2=area,
        perimeter_ft=perimeter,
        sidewall_area_ft2=wall,
        sidewall_coefficients=wall_coeff,
        total_wetted_area_ft2=np.array([3.0, 5.0, 7.0]),
    )


def test_package_round_trip_preserves_arrays_and_metadata(tmp_path) -> None:
    target = tmp_path / "geometry.npz"
    source = valid_kernel()
    save_geometry_package(source, target)
    loaded = load_geometry_package(target)
    assert loaded.metadata == source.metadata
    for name in source.array_names():
        np.testing.assert_array_equal(getattr(loaded, name), getattr(source, name))
    metadata = json.loads(target.with_suffix(".json").read_text(encoding="utf-8"))
    assert metadata["npz_sha256"]


def test_validation_rejects_nonmonotone_volume() -> None:
    kernel = valid_kernel()
    bad = GeometryKernel(
        **{
            **kernel.as_mapping(),
            "volume_ft3": np.array([0.0, 4.0, 3.0], dtype=np.float64),
        }
    )
    with pytest.raises(GeometryPackageError, match="volume_ft3"):
        validate_geometry_kernel(bad)


def test_load_rejects_changed_binary_hash(tmp_path) -> None:
    target = tmp_path / "geometry.npz"
    save_geometry_package(valid_kernel(), target)
    with target.open("ab") as file_obj:
        file_obj.write(b"changed")
    with pytest.raises(GeometryPackageError, match="SHA-256"):
        load_geometry_package(target)
