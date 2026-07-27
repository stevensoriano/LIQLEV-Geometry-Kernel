from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from liqlev.geometry.package import (
    GeometryPackageError,
    load_geometry_package,
    save_geometry_package,
    validate_geometry_kernel,
)
from liqlev.geometry.schema import (
    V1_MIGRATED_AUDIT_STATUS,
    GeometryKernel,
    GeometryMetadata,
)

# Committed production package (on-disk schema_version=1; D1-asserted bytes).
COMMITTED_GEOMETRY_NPZ = (
    Path(__file__).resolve().parents[2]
    / "geometry"
    / "tables"
    / "nhq01-m21a-0201_LIQLEV_GEOMETRY.npz"
)


def valid_kernel(*, schema_version: int = 1, **metadata_overrides) -> GeometryKernel:
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
    meta_kwargs = {
        "schema_version": schema_version,
        "geometry_id": "unit-cylinder",
        "source_step_sha256": "0" * 64,
        "fluid_step_sha256": "1" * 64,
        "axis": "+Y",
        "gravity_direction": "-Y",
        "length_unit": "ft",
        "area_unit": "ft^2",
        "volume_unit": "ft^3",
        "y_min_mm": 0.0,
        "y_max_mm": 609.6,
    }
    meta_kwargs.update(metadata_overrides)
    return GeometryKernel(
        metadata=GeometryMetadata(**meta_kwargs),
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


@pytest.mark.parametrize("metadata_root", [[], 7])
def test_load_rejects_nonobject_metadata_root(tmp_path, metadata_root) -> None:
    target = tmp_path / "geometry.npz"
    save_geometry_package(valid_kernel(), target)
    target.with_suffix(".json").write_text(
        json.dumps(metadata_root),
        encoding="utf-8",
    )

    with pytest.raises(GeometryPackageError, match="root"):
        load_geometry_package(target)


def test_load_rejects_invalid_metadata_json(tmp_path) -> None:
    target = tmp_path / "geometry.npz"
    save_geometry_package(valid_kernel(), target)
    target.with_suffix(".json").write_text("{", encoding="utf-8")

    with pytest.raises(GeometryPackageError, match="JSON"):
        load_geometry_package(target)


def test_load_rejects_missing_metadata_hash(tmp_path) -> None:
    target = tmp_path / "geometry.npz"
    save_geometry_package(valid_kernel(), target)
    metadata_path = target.with_suffix(".json")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata.pop("npz_sha256")
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")

    with pytest.raises(GeometryPackageError, match="npz_sha256"):
        load_geometry_package(target)


@pytest.mark.parametrize("bad_hash", [None, 7])
def test_load_rejects_nonstring_metadata_hash(tmp_path, bad_hash) -> None:
    target = tmp_path / "geometry.npz"
    save_geometry_package(valid_kernel(), target)
    metadata_path = target.with_suffix(".json")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["npz_sha256"] = bad_hash
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")

    with pytest.raises(GeometryPackageError, match="npz_sha256"):
        load_geometry_package(target)


def test_load_rejects_malformed_metadata_fields(tmp_path) -> None:
    target = tmp_path / "geometry.npz"
    save_geometry_package(valid_kernel(), target)
    metadata_path = target.with_suffix(".json")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata.pop("axis")
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")

    with pytest.raises(GeometryPackageError, match="metadata fields"):
        load_geometry_package(target)


def test_load_rejects_corrupt_npz_with_matching_hash(tmp_path) -> None:
    target = tmp_path / "geometry.npz"
    save_geometry_package(valid_kernel(), target)
    target.write_bytes(b"not a valid numpy archive")
    metadata_path = target.with_suffix(".json")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["npz_sha256"] = hashlib.sha256(target.read_bytes()).hexdigest()
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")

    with pytest.raises(GeometryPackageError, match="NPZ archive"):
        load_geometry_package(target)


def test_load_rejects_non_path_value() -> None:
    with pytest.raises(GeometryPackageError, match="path"):
        load_geometry_package(123)


# ---------------------------------------------------------------------------
# F12 — GeometryMetadata schema v2 (four-criteria guard tests)
#
# Call chain (real entry points):
#   load_geometry_package
#     -> _metadata_from_payload  (v1 migrate / v2 accept / v3 reject)
#     -> GeometryKernel(...)
#     -> validate_geometry_kernel  (schema_version in {1, 2})
#   save_geometry_package
#     -> validate_geometry_kernel
#     -> write NPZ + JSON (tmp_path only in tests)
# ---------------------------------------------------------------------------


def test_load_rejects_schema_version_3_via_real_loader(tmp_path) -> None:
    """1 FIRES on bad (v3) · 2 recovered raise · 3/4 real load_geometry_package."""

    target = tmp_path / "geometry.npz"
    save_geometry_package(valid_kernel(schema_version=2), target)
    metadata_path = target.with_suffix(".json")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["schema_version"] = 3
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    with pytest.raises(
        GeometryPackageError,
        match=r"unsupported geometry schema_version 3",
    ):
        load_geometry_package(target)


def test_committed_v1_package_loads_with_migrated_v2_defaults() -> None:
    """1 v1 on-disk · 2 recovered sentinel defaults · 3/4 real load of committed package.

    Does not rewrite committed files; proves in-memory migration only.
    """

    assert COMMITTED_GEOMETRY_NPZ.is_file()
    kernel = load_geometry_package(COMMITTED_GEOMETRY_NPZ)

    on_disk = json.loads(
        COMMITTED_GEOMETRY_NPZ.with_suffix(".json").read_text(encoding="utf-8")
    )
    assert on_disk["schema_version"] == 1
    assert "tolerances" not in on_disk
    assert "audit_status" not in on_disk

    # Recovered in-memory v2-shaped metadata (sentinels; version stays 1 on-disk).
    assert kernel.metadata.schema_version == 1
    assert kernel.metadata.tolerances is None
    assert kernel.metadata.audit_status == V1_MIGRATED_AUDIT_STATUS
    assert kernel.metadata.geometry_id == on_disk["geometry_id"]
    assert len(kernel.height_ft) >= 3


def test_schema_v2_roundtrip_via_tmp_package(tmp_path) -> None:
    """1 v2 authoring · 2 recovered fields equal · 3/4 save+load real package path."""

    source = valid_kernel(
        schema_version=2,
        tolerances={"volume_relative": 5.0e-4, "sidewall_relative": 5.0e-4},
        audit_status="cad-audit-passed",
    )
    target = tmp_path / "geometry_v2.npz"
    save_geometry_package(source, target)
    loaded = load_geometry_package(target)

    assert loaded.metadata == source.metadata
    assert loaded.metadata.schema_version == 2
    assert loaded.metadata.tolerances == {
        "volume_relative": 5.0e-4,
        "sidewall_relative": 5.0e-4,
    }
    assert loaded.metadata.audit_status == "cad-audit-passed"
    for name in source.array_names():
        np.testing.assert_array_equal(getattr(loaded, name), getattr(source, name))

    disk_meta = json.loads(target.with_suffix(".json").read_text(encoding="utf-8"))
    assert disk_meta["schema_version"] == 2
    assert disk_meta["tolerances"] == source.metadata.tolerances
    assert disk_meta["audit_status"] == "cad-audit-passed"
    assert "npz_sha256" in disk_meta


def test_validate_rejects_schema_version_above_supported() -> None:
    """Direct validate path: unsupported version fires (supports reachability of gate)."""

    bad = valid_kernel(schema_version=3)
    with pytest.raises(GeometryPackageError, match=r"schema_version must be 1 or 2"):
        validate_geometry_kernel(bad)
