from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from zipfile import BadZipFile

import numpy as np

from .schema import GeometryKernel, GeometryMetadata


class GeometryPackageError(ValueError):
    pass


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file_obj:
        for block in iter(lambda: file_obj.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def _float64_contiguous(name: str, value: np.ndarray) -> np.ndarray:
    array = np.asarray(value)
    if array.dtype != np.float64:
        raise GeometryPackageError(f"{name} must use float64")
    if array.ndim not in (1, 2) or not array.flags.c_contiguous:
        raise GeometryPackageError(f"{name} must be a contiguous 1D or 2D array")
    if not np.all(np.isfinite(array)):
        raise GeometryPackageError(f"{name} contains non-finite values")
    return array


def validate_geometry_kernel(kernel: GeometryKernel) -> None:
    if kernel.metadata.schema_version != 1:
        raise GeometryPackageError("schema_version must equal 1")
    expected_units = ("+Y", "-Y", "ft", "ft^2", "ft^3")
    actual_units = (
        kernel.metadata.axis,
        kernel.metadata.gravity_direction,
        kernel.metadata.length_unit,
        kernel.metadata.area_unit,
        kernel.metadata.volume_unit,
    )
    if actual_units != expected_units:
        raise GeometryPackageError(f"axis/unit contract mismatch: {actual_units}")
    for name in kernel.array_names():
        _float64_contiguous(name, getattr(kernel, name))
    count = len(kernel.height_ft)
    for name in (
        "volume_ft3",
        "section_area_ft2",
        "perimeter_ft",
        "sidewall_area_ft2",
        "total_wetted_area_ft2",
    ):
        if len(getattr(kernel, name)) != count:
            raise GeometryPackageError(f"{name} length must match height_ft")
    if kernel.volume_coefficients.shape != (4, count - 1):
        raise GeometryPackageError("volume_coefficients must have shape (4, N-1)")
    if kernel.sidewall_coefficients.shape != (4, count - 1):
        raise GeometryPackageError("sidewall_coefficients must have shape (4, N-1)")
    if count < 3 or np.any(np.diff(kernel.height_ft) <= 0.0):
        raise GeometryPackageError("height_ft must contain at least 3 increasing nodes")
    for name in ("volume_ft3", "sidewall_area_ft2", "total_wetted_area_ft2"):
        if np.any(np.diff(getattr(kernel, name)) < 0.0):
            raise GeometryPackageError(f"{name} must be non-decreasing")
    for name in ("section_area_ft2", "perimeter_ft"):
        if np.any(getattr(kernel, name) < 0.0):
            raise GeometryPackageError(f"{name} must be non-negative")
    if abs(kernel.height_ft[0]) > 1e-12 or abs(kernel.volume_ft3[0]) > 1e-12:
        raise GeometryPackageError("height_ft and volume_ft3 must begin at zero")
    if kernel.total_height_ft <= 0.0 or kernel.total_volume_ft3 <= 0.0:
        raise GeometryPackageError("geometry totals must be positive")


def save_geometry_package(kernel: GeometryKernel, path: str | Path) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    validate_geometry_kernel(kernel)
    np.savez_compressed(
        target,
        **{name: getattr(kernel, name) for name in kernel.array_names()},
    )
    payload = kernel.metadata.as_dict()
    payload["npz_sha256"] = _sha256(target)
    target.with_suffix(".json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    csv_path = target.with_suffix(".csv")
    with csv_path.open("w", newline="", encoding="utf-8") as file_obj:
        writer = csv.writer(file_obj)
        writer.writerow(
            [
                "height_ft",
                "volume_ft3",
                "section_area_ft2",
                "perimeter_ft",
                "sidewall_area_ft2",
                "total_wetted_area_ft2",
            ]
        )
        writer.writerows(
            zip(
                kernel.height_ft,
                kernel.volume_ft3,
                kernel.section_area_ft2,
                kernel.perimeter_ft,
                kernel.sidewall_area_ft2,
                kernel.total_wetted_area_ft2,
                strict=True,
            )
        )


def load_geometry_package(path: str | Path) -> GeometryKernel:
    try:
        target = Path(path)
    except (TypeError, ValueError) as exc:
        raise GeometryPackageError(
            "geometry package path must be a string or path-like value"
        ) from exc

    metadata_path = target.with_suffix(".json")
    try:
        metadata_text = metadata_path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        raise GeometryPackageError(
            f"could not read geometry metadata JSON: {metadata_path}"
        ) from exc
    try:
        payload = json.loads(metadata_text)
    except json.JSONDecodeError as exc:
        raise GeometryPackageError("geometry metadata JSON is invalid") from exc
    if not isinstance(payload, dict):
        raise GeometryPackageError("geometry metadata JSON root must be an object")
    if "npz_sha256" not in payload:
        raise GeometryPackageError("geometry metadata is missing npz_sha256")
    expected_hash_value = payload.pop("npz_sha256")
    if not isinstance(expected_hash_value, str):
        raise GeometryPackageError("geometry metadata npz_sha256 must be a string")
    expected_hash = expected_hash_value.upper()
    try:
        actual_hash = _sha256(target)
    except OSError as exc:
        raise GeometryPackageError(f"could not read geometry NPZ: {target}") from exc
    if actual_hash != expected_hash:
        raise GeometryPackageError("NPZ SHA-256 does not match metadata")
    try:
        metadata = GeometryMetadata(**payload)
    except (TypeError, ValueError) as exc:
        raise GeometryPackageError(
            f"geometry metadata fields are invalid: {exc}"
        ) from exc
    try:
        with np.load(target, allow_pickle=False) as archive:
            arrays = {
                name: np.ascontiguousarray(archive[name], dtype=np.float64)
                for name in GeometryKernel.array_names()
            }
    except (
        OSError,
        EOFError,
        KeyError,
        TypeError,
        ValueError,
        OverflowError,
        BadZipFile,
    ) as exc:
        raise GeometryPackageError(f"invalid geometry NPZ archive: {exc}") from exc
    kernel = GeometryKernel(metadata=metadata, **arrays)
    validate_geometry_kernel(kernel)
    return kernel
