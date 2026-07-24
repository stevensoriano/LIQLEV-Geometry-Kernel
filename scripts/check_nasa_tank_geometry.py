"""Independently verify the published NASA tank geometry artifact set."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
import sys
import cadquery as cq
import numpy as np
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.GeomAbs import GeomAbs_Plane


REPOSITORY = Path(__file__).resolve().parents[1]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

from liqlev.geometry import load_geometry_package  # noqa: E402
from liqlev.cad.measure import (  # noqa: E402
    MM2_TO_FT2,
    MM3_TO_FT3,
    MM_TO_FT,
    measure_geometry,
)


DEFAULT_SOURCE_STEP = Path(
    r"C:\Users\sasorian\Documents\Eta_Space\geometry"
    r"\nhq01-m21a- 0201_TankAssy_NASA.STEP"
)
DEFAULT_SOURCE_SHA256 = (
    "0193EC296B5B754FFDC7B1652052F5B28BA99345A5BC89922D026DF011C837D5"
)
DEFAULT_GEOMETRY_ROOT = REPOSITORY / "geometry"
DEFAULT_Y_MIN_MM = -275.406791
DEFAULT_Y_MAX_MM = 275.296791
EXPECTED_FLUID_BOUNDS_MM = (
    -279.67000017112,
    279.67000017112,
    -275.4067913371,
    275.29679133713,
    297.81275947409995,
    857.15275961635,
)
EXPECTED_ENDPOINT_MOSAIC_FACE_COUNTS = (27, 3)
EXPECTED_FACE_COUNT = 40
STEP_RELATIVE = Path("output/nhq01-m21a-0201_LIQLEV_FLUID_DOMAIN.step")
PACKAGE_RELATIVE = Path("tables/nhq01-m21a-0201_LIQLEV_GEOMETRY.npz")
AUDIT_RELATIVE = Path("audit/nhq01-m21a-0201_LIQLEV_AUDIT.json")
PLANE_TOLERANCE_MM = 1e-5
ROUND_TRIP_RELATIVE_VOLUME_TOLERANCE = 1e-8
KERNEL_RELATIVE_TOLERANCE = 5e-4
EXPECTED_AXIS_CONTRACT = ("+Y", "-Y", "ft", "ft^2", "ft^3")


@dataclass(frozen=True)
class Check:
    name: str
    passed: bool
    detail: str


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file_obj:
        for block in iter(lambda: file_obj.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def _bounds(shape: cq.Shape) -> tuple[float, float, float, float, float, float]:
    box = shape.BoundingBox()
    return box.xmin, box.xmax, box.ymin, box.ymax, box.zmin, box.zmax


def _closure_plane_error(
    solid: cq.Solid,
    y_min_mm: float,
    y_max_mm: float,
) -> tuple[float, tuple[int, int]]:
    y_normal_planes: list[float] = []
    for face in solid.Faces():
        surface = BRepAdaptor_Surface(face.wrapped, True)
        if surface.GetType() != GeomAbs_Plane:
            continue
        plane = surface.Plane()
        direction = plane.Axis().Direction()
        if (
            math.hypot(direction.X(), direction.Z()) <= 1e-8
            and abs(abs(direction.Y()) - 1.0) <= 1e-8
        ):
            y_normal_planes.append(plane.Location().Y())
    deviations = [
        [
            abs(plane_y - closure_y)
            for plane_y in y_normal_planes
            if abs(plane_y - closure_y) <= PLANE_TOLERANCE_MM
        ]
        for closure_y in (y_min_mm, y_max_mm)
    ]
    counts = (len(deviations[0]), len(deviations[1]))
    if any(count == 0 for count in counts):
        return math.inf, counts
    return max(value for group in deviations for value in group), counts


def _numeric_dv_dh_volume(kernel) -> float:
    sample_height = np.linspace(
        kernel.height_ft[0],
        kernel.height_ft[-1],
        8193,
        dtype=np.float64,
    )
    interval = np.searchsorted(
        kernel.height_ft,
        sample_height,
        side="right",
    ) - 1
    interval = np.clip(interval, 0, len(kernel.height_ft) - 2)
    local = sample_height - kernel.height_ft[interval]
    coefficients = kernel.volume_coefficients
    derivative = (
        3.0 * coefficients[0, interval] * local**2
        + 2.0 * coefficients[1, interval] * local
        + coefficients[2, interval]
    )
    return float(np.trapezoid(derivative, sample_height))


def _eval_ppoly(
    x: np.ndarray,
    nodes: np.ndarray,
    coefficients: np.ndarray,
) -> np.ndarray:
    interval = np.searchsorted(nodes, x, side="right") - 1
    interval = np.clip(interval, 0, len(nodes) - 2)
    local = x - nodes[interval]
    return (
        (
            coefficients[0, interval] * local
            + coefficients[1, interval]
        )
        * local
        + coefficients[2, interval]
    ) * local + coefficients[3, interval]


def independent_cad_table_metrics(solid: cq.Solid, kernel) -> dict[str, float | int]:
    """Remeasure exact midpoint sections independently from the build pass."""

    midpoint_height_ft = 0.5 * (
        kernel.height_ft[:-1] + kernel.height_ft[1:]
    )
    absolute_midpoint_mm = (
        midpoint_height_ft / MM_TO_FT + kernel.metadata.y_min_mm
    )
    measurement_nodes_mm = np.concatenate(
        (
            [kernel.metadata.y_min_mm],
            absolute_midpoint_mm,
            [kernel.metadata.y_max_mm],
        )
    )
    samples = measure_geometry(
        solid,
        measurement_nodes_mm,
        y_min_mm=kernel.metadata.y_min_mm,
        y_max_mm=kernel.metadata.y_max_mm,
    )
    measured_volume_ft3 = samples.volume_mm3[1:-1] * MM3_TO_FT3
    measured_sidewall_ft2 = samples.sidewall_area_mm2[1:-1] * MM2_TO_FT2
    direct_area_ft2 = samples.section_area_mm2[1:-1] * MM2_TO_FT2
    predicted_volume_ft3 = _eval_ppoly(
        midpoint_height_ft,
        kernel.height_ft,
        kernel.volume_coefficients,
    )
    predicted_sidewall_ft2 = _eval_ppoly(
        midpoint_height_ft,
        kernel.height_ft,
        kernel.sidewall_coefficients,
    )
    interval = np.arange(len(midpoint_height_ft))
    local = midpoint_height_ft - kernel.height_ft[:-1]
    coefficients = kernel.volume_coefficients
    derivative_area_ft2 = (
        3.0 * coefficients[0, interval] * local**2
        + 2.0 * coefficients[1, interval] * local
        + coefficients[2, interval]
    )

    y_min_mm = kernel.metadata.y_min_mm
    y_max_mm = kernel.metadata.y_max_mm
    topology_tolerance_mm = 1e-8 * max(
        1.0,
        abs(y_min_mm),
        abs(y_max_mm),
    )
    topology_y_mm = np.asarray(
        [
            y_min_mm,
            *sorted(
                {
                    float(vertex.Center().y)
                    for vertex in solid.Vertices()
                    if y_min_mm + topology_tolerance_mm
                    < float(vertex.Center().y)
                    < y_max_mm - topology_tolerance_mm
                }
            ),
            y_max_mm,
        ],
        dtype=np.float64,
    )
    topology_coincidence = np.asarray(
        [
            np.min(np.abs(topology_y_mm - absolute_y))
            <= topology_tolerance_mm
            for absolute_y in absolute_midpoint_mm
        ]
    )
    endpoint_area_mm2 = np.asarray(
        [samples.section_area_mm2[0], samples.section_area_mm2[-1]]
    )
    degenerate_endpoint_y_mm = np.asarray(
        [y_min_mm, y_max_mm],
        dtype=np.float64,
    )[
        endpoint_area_mm2
        <= 1e-14 * max(1.0, float(np.max(samples.section_area_mm2)))
    ]
    endpoint_clearance_mm = max(
        topology_tolerance_mm,
        (y_max_mm - y_min_mm) / 32.0,
    )
    if len(degenerate_endpoint_y_mm):
        degenerate_neighborhood = np.asarray(
            [
                np.min(np.abs(degenerate_endpoint_y_mm - absolute_y))
                <= endpoint_clearance_mm
                for absolute_y in absolute_midpoint_mm
            ]
        )
    else:
        degenerate_neighborhood = np.zeros(
            len(absolute_midpoint_mm),
            dtype=bool,
        )
    positive_area = direct_area_ft2 > 1e-14 * max(
        1.0,
        float(np.max(direct_area_ft2)),
    )
    eligible = (
        ~topology_coincidence
        & ~degenerate_neighborhood
        & positive_area
    )
    if not np.any(eligible):
        raise ValueError("independent direct-area check has no eligible midpoints")

    cad_volume_ft3 = solid.Volume() * MM3_TO_FT3
    sidewall_scale = max(
        abs(float(kernel.sidewall_area_ft2[-1])),
        np.finfo(np.float64).tiny,
    )
    integrated_volume_ft3 = _numeric_dv_dh_volume(kernel)
    return {
        "cad_endpoint_volume_relative": float(
            abs(kernel.total_volume_ft3 - cad_volume_ft3) / cad_volume_ft3
        ),
        "integrated_dv_dh_relative": float(
            abs(integrated_volume_ft3 - kernel.total_volume_ft3)
            / kernel.total_volume_ft3
        ),
        "max_midpoint_volume_relative": float(
            np.max(np.abs(measured_volume_ft3 - predicted_volume_ft3))
            / cad_volume_ft3
        ),
        "max_midpoint_sidewall_relative": float(
            np.max(
                np.abs(measured_sidewall_ft2 - predicted_sidewall_ft2)
            )
            / sidewall_scale
        ),
        "max_direct_area_relative": float(
            np.max(
                np.abs(
                    derivative_area_ft2[eligible]
                    - direct_area_ft2[eligible]
                )
                / direct_area_ft2[eligible]
            )
        ),
        "eligible_direct_area_midpoints": int(np.count_nonzero(eligible)),
    }


def _read_csv_row_count(path: Path) -> tuple[list[str], int]:
    with path.open(newline="", encoding="utf-8") as file_obj:
        reader = csv.reader(file_obj)
        rows = list(reader)
    if not rows:
        return [], 0
    return rows[0], len(rows) - 1


def verify_geometry_root(
    *,
    source_step: Path,
    geometry_root: Path,
    expected_source_sha256: str,
    require_visual_review: bool = True,
    precomputed_measurement: dict[str, float | int] | None = None,
) -> list[Check]:
    """Return independent checks without invoking the build command."""

    checks: list[Check] = []

    def record(name: str, condition: bool, detail: str) -> None:
        checks.append(Check(name=name, passed=bool(condition), detail=detail))

    expected_source_sha256 = expected_source_sha256.upper()
    step_path = geometry_root / STEP_RELATIVE
    package_path = geometry_root / PACKAGE_RELATIVE
    metadata_path = package_path.with_suffix(".json")
    csv_path = package_path.with_suffix(".csv")
    audit_path = geometry_root / AUDIT_RELATIVE
    required = (source_step, step_path, package_path, metadata_path, csv_path, audit_path)
    missing = [str(path) for path in required if not path.is_file()]
    record("required files", not missing, f"missing={missing}")
    if missing:
        return checks

    source_hash = _sha256(source_step)
    record(
        "source SHA-256",
        source_hash == expected_source_sha256,
        f"actual={source_hash}",
    )
    try:
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        record("audit JSON", False, str(exc))
        return checks
    record("audit JSON", isinstance(audit, dict), f"path={audit_path}")
    if not isinstance(audit, dict):
        return checks

    record(
        "audit source preservation",
        audit.get("source_sha256") == source_hash
        and audit.get("source_preservation", {}).get("unchanged") is True
        and audit.get("source_preservation", {}).get("after_sha256") == source_hash
        and audit.get("source_preservation", {}).get("after_size_bytes")
        == source_step.stat().st_size,
        f"source_size_bytes={source_step.stat().st_size}",
    )

    expected_hashes = audit.get("artifact_sha256", {})
    artifact_paths = {
        "fluid_step": step_path,
        "geometry_npz": package_path,
        "geometry_json": metadata_path,
        "geometry_csv": csv_path,
    }
    artifact_hash_ok = isinstance(expected_hashes, dict)
    artifact_hash_details: list[str] = []
    for name, path in artifact_paths.items():
        actual = _sha256(path)
        expected = expected_hashes.get(name) if isinstance(expected_hashes, dict) else None
        artifact_hash_ok &= actual == expected
        artifact_hash_details.append(f"{name}={actual}")
    record(
        "core artifact SHA-256",
        artifact_hash_ok,
        ", ".join(artifact_hash_details),
    )

    try:
        imported = cq.importers.importStep(str(step_path))
        solids = imported.solids().vals()
    except Exception as exc:
        record("STEP re-import", False, repr(exc))
        return checks
    record("STEP re-import", len(solids) == 1, f"solid_count={len(solids)}")
    if len(solids) != 1 or not isinstance(solids[0], cq.Solid):
        return checks
    solid = solids[0]
    shells = solid.Shells()
    valid = solid.isValid() and BRepCheck_Analyzer(solid.wrapped).IsValid()
    closed = len(shells) == 1 and shells[0].wrapped.Closed()
    bounds = _bounds(solid)
    record(
        "BRep topology",
        len(solid.Solids()) == 1
        and len(shells) == 1
        and len(solid.Faces()) == EXPECTED_FACE_COUNT
        and valid
        and closed,
        (
            f"solids={len(solid.Solids())}, shells={len(shells)}, "
            f"faces={len(solid.Faces())}, valid={valid}, closed={closed}"
        ),
    )

    closure_planes = audit.get("closure_planes_mm", [])
    closure_contract = closure_planes == [
        DEFAULT_Y_MIN_MM,
        DEFAULT_Y_MAX_MM,
    ]
    cap_error, cap_counts = _closure_plane_error(
        solid,
        DEFAULT_Y_MIN_MM,
        DEFAULT_Y_MAX_MM,
    )
    record(
        "closure planes and endpoint mosaics",
        closure_contract
        and cap_error <= PLANE_TOLERANCE_MM
        and cap_counts == EXPECTED_ENDPOINT_MOSAIC_FACE_COUNTS
        and audit.get("endpoint_mosaic_face_counts")
        == list(EXPECTED_ENDPOINT_MOSAIC_FACE_COUNTS),
        (
            f"planes={closure_planes}, max_error_mm={cap_error}, "
            f"candidate_counts={cap_counts}"
        ),
    )

    recorded_bounds = audit.get("bounding_box_mm")
    bounds_error = (
        max(
            abs(float(after) - float(before))
            for before, after in zip(recorded_bounds, bounds, strict=True)
        )
        if isinstance(recorded_bounds, list) and len(recorded_bounds) == 6
        else math.inf
    )
    record(
        "recorded BRep metrics",
        audit.get("solid_count") == 1
        and audit.get("shell_count") == 1
        and audit.get("valid") is True
        and audit.get("closed_outer_shell") is True
        and bounds_error <= PLANE_TOLERANCE_MM
        and math.isclose(
            float(audit.get("volume_mm3", math.nan)),
            solid.Volume(),
            rel_tol=ROUND_TRIP_RELATIVE_VOLUME_TOLERANCE,
            abs_tol=1e-3,
        ),
        f"bounds_max_error_mm={bounds_error}, volume_mm3={solid.Volume()}",
    )
    expected_bounds_error = max(
        abs(actual - expected)
        for actual, expected in zip(
            bounds,
            EXPECTED_FLUID_BOUNDS_MM,
            strict=True,
        )
    )
    record(
        "fixed fluid bounds",
        expected_bounds_error <= PLANE_TOLERANCE_MM,
        f"max_error_mm={expected_bounds_error}, bounds={bounds}",
    )
    record(
        "round-trip volume tolerance",
        float(audit.get("round_trip_relative_volume_error", math.inf))
        <= ROUND_TRIP_RELATIVE_VOLUME_TOLERANCE,
        (
            "relative_error="
            f"{audit.get('round_trip_relative_volume_error', math.inf)}"
        ),
    )

    try:
        kernel = load_geometry_package(package_path)
    except Exception as exc:
        record("geometry package load", False, repr(exc))
        return checks
    record(
        "geometry package load",
        True,
        f"node_count={len(kernel.height_ft)}",
    )
    metadata_contract = (
        kernel.metadata.axis,
        kernel.metadata.gravity_direction,
        kernel.metadata.length_unit,
        kernel.metadata.area_unit,
        kernel.metadata.volume_unit,
    )
    record(
        "geometry metadata",
        metadata_contract == EXPECTED_AXIS_CONTRACT
        and kernel.metadata.source_step_sha256 == source_hash
        and kernel.metadata.fluid_step_sha256 == _sha256(step_path)
        and kernel.metadata.schema_version == 1
        and kernel.metadata.geometry_id
        == "nhq01-m21a-0201_LIQLEV_GEOMETRY"
        and kernel.metadata.y_min_mm == DEFAULT_Y_MIN_MM
        and kernel.metadata.y_max_mm == DEFAULT_Y_MAX_MM,
        (
            f"contract={metadata_contract}, "
            f"y=({kernel.metadata.y_min_mm}, {kernel.metadata.y_max_mm})"
        ),
    )
    expected_height_ft = (
        kernel.metadata.y_max_mm - kernel.metadata.y_min_mm
    ) / 304.8
    record(
        "height and node count",
        abs(kernel.total_height_ft - expected_height_ft) <= 1e-10
        and len(kernel.height_ft) == audit.get("node_count")
        and 3 <= len(kernel.height_ft) <= int(audit.get("max_nodes", 0)),
        (
            f"height_ft={kernel.total_height_ft}, "
            f"expected={expected_height_ft}, nodes={len(kernel.height_ft)}"
        ),
    )
    cad_volume_ft3 = solid.Volume() / 304.8**3
    cad_volume_error = abs(kernel.total_volume_ft3 - cad_volume_ft3) / cad_volume_ft3
    record(
        "CAD/table volume",
        cad_volume_error <= KERNEL_RELATIVE_TOLERANCE,
        f"relative_error={cad_volume_error}",
    )
    integrated_volume = _numeric_dv_dh_volume(kernel)
    dv_dh_error = (
        abs(integrated_volume - kernel.total_volume_ft3)
        / kernel.total_volume_ft3
    )
    record(
        "integrated dV/dh volume",
        dv_dh_error <= KERNEL_RELATIVE_TOLERANCE,
        (
            f"integrated_ft3={integrated_volume}, "
            f"relative_error={dv_dh_error}"
        ),
    )
    measurement = (
        precomputed_measurement
        if precomputed_measurement is not None
        else independent_cad_table_metrics(solid, kernel)
    )
    recorded_measurement = audit.get("independent_measurement", {})
    measurement_matches_audit = (
        isinstance(recorded_measurement, dict)
        and recorded_measurement.keys() == measurement.keys()
        and all(
            (
                int(recorded_measurement[name]) == int(value)
                if isinstance(value, int)
                else math.isclose(
                    float(recorded_measurement[name]),
                    float(value),
                    rel_tol=1e-12,
                    abs_tol=1e-15,
                )
            )
            for name, value in measurement.items()
        )
    )
    record(
        "independent CAD/table refinement",
        float(measurement.get("max_midpoint_volume_relative", math.inf))
        <= KERNEL_RELATIVE_TOLERANCE
        and float(measurement.get("max_midpoint_sidewall_relative", math.inf))
        <= KERNEL_RELATIVE_TOLERANCE
        and float(measurement.get("cad_endpoint_volume_relative", math.inf))
        <= KERNEL_RELATIVE_TOLERANCE
        and float(measurement.get("integrated_dv_dh_relative", math.inf))
        <= KERNEL_RELATIVE_TOLERANCE
        and float(measurement.get("max_direct_area_relative", math.inf))
        <= 2e-3
        and int(measurement.get("eligible_direct_area_midpoints", 0)) > 0
        and measurement_matches_audit,
        f"values={measurement}, matches_audit={measurement_matches_audit}",
    )

    csv_header, csv_rows = _read_csv_row_count(csv_path)
    with csv_path.open(newline="", encoding="utf-8") as file_obj:
        csv_payload = list(csv.reader(file_obj))
    csv_values = np.asarray(csv_payload[1:], dtype=np.float64)
    expected_csv_values = np.column_stack(
        (
            kernel.height_ft,
            kernel.volume_ft3,
            kernel.section_area_ft2,
            kernel.perimeter_ft,
            kernel.sidewall_area_ft2,
            kernel.total_wetted_area_ft2,
        )
    )
    record(
        "CSV table",
        csv_rows == len(kernel.height_ft)
        and csv_header
        == [
            "height_ft",
            "volume_ft3",
            "section_area_ft2",
            "perimeter_ft",
            "sidewall_area_ft2",
            "total_wetted_area_ft2",
        ]
        and np.array_equal(csv_values, expected_csv_values),
        f"rows={csv_rows}, values_match_npz={np.array_equal(csv_values, expected_csv_values)}",
    )

    image_paths = audit.get("section_images", [])
    image_hash_ok = isinstance(image_paths, list) and len(image_paths) == 10
    image_details: list[str] = []
    if isinstance(image_paths, list):
        for relative in image_paths:
            path = audit_path.parent / str(relative)
            exists = path.is_file() and path.stat().st_size > 0
            actual = _sha256(path) if exists else "MISSING"
            expected = (
                expected_hashes.get(str(relative))
                if isinstance(expected_hashes, dict)
                else None
            )
            image_hash_ok &= exists and actual == expected
            image_details.append(f"{relative}={actual}")
    record(
        "section image SHA-256",
        image_hash_ok,
        ", ".join(image_details),
    )
    render_modes = audit.get("section_render_modes", {})
    record(
        "section render modes",
        isinstance(render_modes, dict)
        and render_modes.get("sections/section_xz_h000.png")
        == "endpoint_closure_mosaic"
        and render_modes.get("sections/section_xz_h100.png")
        == "endpoint_closure_mosaic"
        and all(
            relative in render_modes
            for relative in image_paths
        ),
        f"modes={render_modes}",
    )
    tool_versions = audit.get("tool_versions", {})
    build_command = audit.get("build_command", {})
    record(
        "tool and command provenance",
        isinstance(tool_versions, dict)
        and all(
            isinstance(tool_versions.get(name), str)
            and bool(tool_versions[name])
            for name in ("python", "cadquery", "ocp", "numpy", "matplotlib")
        )
        and isinstance(build_command, dict)
        and isinstance(build_command.get("executable"), str)
        and bool(build_command["executable"])
        and isinstance(build_command.get("argv"), list)
        and bool(build_command["argv"]),
        f"tools={tool_versions}, command={build_command}",
    )
    visual_status = audit.get("visual_section_review")
    record(
        "visual section review",
        visual_status == "passed" if require_visual_review else visual_status in {"pending", "passed"},
        f"status={visual_status}",
    )
    record(
        "overall audit",
        audit.get("technical_validation_passed") is True
        and (
            audit.get("passed") is True
            if require_visual_review
            else audit.get("passed") in {False, True}
        ),
        (
            f"technical={audit.get('technical_validation_passed')}, "
            f"passed={audit.get('passed')}"
        ),
    )
    return checks


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-step", type=Path, default=DEFAULT_SOURCE_STEP)
    parser.add_argument(
        "--geometry-root",
        type=Path,
        default=DEFAULT_GEOMETRY_ROOT,
    )
    parser.add_argument(
        "--source-sha256",
        default=DEFAULT_SOURCE_SHA256,
    )
    parser.add_argument(
        "--allow-pending-visual-review",
        action="store_true",
        help="Allow a technically validated review set whose visual status is pending.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        checks = verify_geometry_root(
            source_step=args.source_step,
            geometry_root=args.geometry_root,
            expected_source_sha256=args.source_sha256,
            require_visual_review=not args.allow_pending_visual_review,
        )
    except Exception as exc:
        print(f"FAIL unexpected verification error: {exc}")
        return 1
    for check in checks:
        status = "PASS" if check.passed else "FAIL"
        print(f"{status} {check.name}: {check.detail}")
    if checks and all(check.passed for check in checks):
        print("NASA tank geometry audit passed.")
        return 0
    print("NASA tank geometry audit failed.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
