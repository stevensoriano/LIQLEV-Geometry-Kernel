"""Build the validated NASA tank CAD, numeric tables, audit, and sections."""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import sys
from tempfile import TemporaryDirectory
import time
from typing import Callable, TypeVar

import cadquery as cq
import matplotlib
import numpy as np
import OCP


matplotlib.use("Agg")
from matplotlib import pyplot as plt  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402


REPOSITORY = Path(__file__).resolve().parents[1]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

from check_nasa_tank_geometry import (  # noqa: E402
    AUDIT_RELATIVE,
    PACKAGE_RELATIVE,
    STEP_RELATIVE,
    independent_cad_table_metrics,
    verify_geometry_root,
)
from liqlev.cad import (  # noqa: E402
    build_fluid_domain,
    load_named_product,
    sha256_file,
    write_step_round_trip,
)
from liqlev.cad.measure import (  # noqa: E402
    build_geometry_kernel,
)
from liqlev.geometry import GeometryMetadata, save_geometry_package  # noqa: E402


DEFAULT_SOURCE_STEP = Path(
    r"C:\Users\sasorian\Documents\Eta_Space\geometry"
    r"\nhq01-m21a- 0201_TankAssy_NASA.STEP"
)
DEFAULT_OUTPUT_ROOT = REPOSITORY / "geometry"
DEFAULT_SOURCE_SHA256 = (
    "0193EC296B5B754FFDC7B1652052F5B28BA99345A5BC89922D026DF011C837D5"
)
DEFAULT_PRODUCT_NAME = "nhq01-m21a- 0202_short"
DEFAULT_Y_MIN_MM = -275.406791
DEFAULT_Y_MAX_MM = 275.296791
DEFAULT_MAX_NODES = 1025
PLANE_TOLERANCE_MM = 1e-5
REFINEMENT_RELATIVE_TOLERANCE = 5e-4
GEOMETRY_ID = "nhq01-m21a-0201_LIQLEV_GEOMETRY"
SECTION_SPECS = (
    ("orthogonal_xy.png", "XY"),
    ("orthogonal_yz.png", "YZ"),
    ("orthogonal_xz.png", "XZ"),
)
SECTION_HEIGHTS = (
    (0, 0.00),
    (10, 0.10),
    (25, 0.25),
    (50, 0.50),
    (75, 0.75),
    (90, 0.90),
    (100, 1.00),
)
T = TypeVar("T")


def _timed(
    timings: dict[str, float],
    name: str,
    operation: Callable[[], T],
) -> T:
    print(f"STAGE START {name}", file=sys.stderr, flush=True)
    started = time.perf_counter()
    try:
        result = operation()
    except Exception as exc:
        elapsed = round(time.perf_counter() - started, 6)
        print(
            f"STAGE FAIL {name} elapsed_seconds={elapsed}: {exc}",
            file=sys.stderr,
            flush=True,
        )
        raise
    timings[name] = round(time.perf_counter() - started, 6)
    print(
        f"STAGE PASS {name} elapsed_seconds={timings[name]}",
        file=sys.stderr,
        flush=True,
    )
    return result


def _sha256(path: Path) -> str:
    return sha256_file(path)


def _bounds(shape: cq.Shape) -> list[float]:
    box = shape.BoundingBox()
    return [box.xmin, box.xmax, box.ymin, box.ymax, box.zmin, box.zmax]


def _snapshot_brep_metrics(fluid: cq.Solid) -> dict[str, object]:
    shells = fluid.Shells()
    return {
        "solid_count": len(fluid.Solids()),
        "shell_count": len(shells),
        "face_count": len(fluid.Faces()),
        "valid": fluid.isValid(),
        "closed_outer_shell": (
            len(shells) == 1 and shells[0].wrapped.Closed()
        ),
        "volume_mm3": fluid.Volume(),
        "surface_area_mm2": fluid.Area(),
        "bounding_box_mm": _bounds(fluid),
    }


def _normalize_step_header_timestamp(step_path: Path) -> None:
    payload = step_path.read_bytes()
    normalized, replacements = re.subn(
        rb"(FILE_NAME\([^,]+,')"
        rb"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}"
        rb"(')",
        rb"\g<1>1970-01-01T00:00:00\g<2>",
        payload,
        count=1,
    )
    if replacements != 1:
        raise RuntimeError(
            f"STEP header has no unique FILE_NAME timestamp: {step_path}"
        )
    step_path.write_bytes(normalized)


def _reimport_single_step(step_path: Path) -> cq.Solid:
    imported = cq.importers.importStep(str(step_path))
    solids = imported.solids().vals()
    if len(solids) != 1 or not isinstance(solids[0], cq.Solid):
        raise RuntimeError(
            f"normalized STEP must reimport as one solid: {step_path}"
        )
    return solids[0]


def _configure_axes(ax, *, xlabel: str, ylabel: str, title: str) -> None:
    ax.set_facecolor("#07121c")
    ax.set_xlabel(xlabel, color="#d7e9f5")
    ax.set_ylabel(ylabel, color="#d7e9f5")
    ax.set_title(title, color="#f3f8fb", pad=16, fontsize=13)
    ax.tick_params(colors="#a8c6d9")
    for spine in ax.spines.values():
        spine.set_color("#36576b")
    ax.grid(True, color="#294554", linewidth=0.45, alpha=0.55)
    ax.set_aspect("equal", adjustable="box")


def _add_scale_bar(ax, x: np.ndarray, y: np.ndarray) -> None:
    scale_mm = 100.0
    x_min = float(np.min(x))
    x_max = float(np.max(x))
    y_min = float(np.min(y))
    y_max = float(np.max(y))
    x_start = x_min + 0.08 * (x_max - x_min)
    y_start = y_min + 0.08 * (y_max - y_min)
    ax.plot(
        [x_start, x_start + scale_mm],
        [y_start, y_start],
        color="#f5bd4f",
        linewidth=3.0,
        solid_capstyle="butt",
    )
    ax.text(
        x_start + 0.5 * scale_mm,
        y_start + 0.025 * (y_max - y_min),
        "100 mm",
        color="#f5bd4f",
        ha="center",
        va="bottom",
        fontsize=9,
    )


def _save_figure(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        path,
        dpi=160,
        facecolor="#050b10",
        bbox_inches="tight",
        metadata={"Software": "LIQLEV geometry publication"},
    )
    plt.close(fig)


def _render_projection(
    fluid: cq.Solid,
    path: Path,
    plane: str,
    *,
    y_min_mm: float,
    y_max_mm: float,
) -> None:
    vertices, triangles = fluid.tessellate(0.8, 0.08)
    points = np.asarray([vertex.toTuple() for vertex in vertices])
    triangle_indices = np.asarray(triangles, dtype=np.int64)
    mapping = {
        "XY": (0, 1, "Assembly X [mm]", "Absolute assembly Y [mm]"),
        "YZ": (2, 1, "Assembly Z [mm]", "Absolute assembly Y [mm]"),
        "XZ": (0, 2, "Assembly X [mm]", "Assembly Z [mm]"),
    }
    horizontal, vertical, xlabel, ylabel = mapping[plane]
    coordinates = points[:, [horizontal, vertical]]
    polygons = coordinates[triangle_indices]
    fig, ax = plt.subplots(figsize=(10, 8), facecolor="#050b10")
    collection = PolyCollection(
        polygons,
        facecolor="#146786",
        edgecolor="#61c3e6",
        linewidth=0.08,
        alpha=0.82,
    )
    ax.add_collection(collection)
    ax.autoscale_view()
    title = (
        f"NASA LIQLEV fluid domain — orthogonal {plane} view\n"
        f"Assembly axes | absolute Y = {y_min_mm:.6f}…{y_max_mm:.6f} mm | "
        "normalized h = 0…1"
    )
    _configure_axes(ax, xlabel=xlabel, ylabel=ylabel, title=title)
    _add_scale_bar(ax, coordinates[:, 0], coordinates[:, 1])
    fig.tight_layout()
    _save_figure(fig, path)


def _ordered_wire_points(
    wire: cq.Wire,
    *,
    samples_per_edge: int,
) -> np.ndarray:
    sampled_edges: list[np.ndarray] = []
    for edge in wire.Edges():
        points, _ = edge.sample(samples_per_edge)
        coordinates = np.asarray(
            [(point.x, point.z) for point in points],
            dtype=np.float64,
        )
        if len(coordinates) >= 2:
            sampled_edges.append(coordinates)
    if not sampled_edges:
        raise RuntimeError("wire has no drawable edges")

    ordered = [sampled_edges.pop(0)]
    while sampled_edges:
        end = ordered[-1][-1]
        distances = [
            (
                float(np.linalg.norm(edge[0] - end)),
                float(np.linalg.norm(edge[-1] - end)),
            )
            for edge in sampled_edges
        ]
        index, reverse = min(
            (
                (index, endpoint)
                for index, pair in enumerate(distances)
                for endpoint in range(2)
            ),
            key=lambda choice: distances[choice[0]][choice[1]],
        )
        selected = sampled_edges.pop(index)
        if reverse:
            selected = selected[::-1]
        ordered.append(selected[1:])
    result = np.vstack(ordered)
    if not np.allclose(result[0], result[-1], rtol=0.0, atol=1e-10):
        result = np.vstack((result, result[0]))
    return result


def _endpoint_closure_faces(
    fluid: cq.Solid,
    *,
    y_mm: float,
    normal_y: float,
) -> list[cq.Face]:
    faces: list[cq.Face] = []
    for face in fluid.Faces():
        if face.geomType() != "PLANE":
            continue
        center = face.Center()
        if abs(center.y - y_mm) > PLANE_TOLERANCE_MM:
            continue
        normal = face.normalAt(center)
        if (
            normal.y * normal_y <= 1.0 - 1e-9
            or abs(normal.x) > 1e-9
            or abs(normal.z) > 1e-9
        ):
            continue
        faces.append(face)
    if not faces:
        raise RuntimeError(
            f"no endpoint closure faces found at Y={y_mm:.12g} mm"
        )
    return faces


def _render_endpoint_closure_mosaic(
    fluid: cq.Solid,
    path: Path,
    *,
    normalized_height: float,
    y_mm: float,
) -> int:
    normal_y = -1.0 if normalized_height == 0.0 else 1.0
    faces = _endpoint_closure_faces(
        fluid,
        y_mm=y_mm,
        normal_y=normal_y,
    )
    fig, ax = plt.subplots(figsize=(9, 9), facecolor="#050b10")
    all_points: list[np.ndarray] = []
    for index, face in enumerate(faces):
        outer = _ordered_wire_points(
            face.outerWire(),
            samples_per_edge=160,
        )
        all_points.append(outer)
        ax.fill(
            outer[:, 0],
            outer[:, 1],
            color="#0c536f",
            alpha=0.55 + 0.08 * (index % 2),
        )
        for wire in face.Wires():
            boundary = _ordered_wire_points(
                wire,
                samples_per_edge=160,
            )
            all_points.append(boundary)
            ax.plot(
                boundary[:, 0],
                boundary[:, 1],
                color="#75daf6",
                linewidth=1.0,
            )
    combined = np.vstack(all_points)
    title = (
        "NASA LIQLEV fluid domain — exact X-Z endpoint closure mosaic\n"
        f"Assembly axes | absolute Y = {y_mm:.6f} mm | "
        f"normalized h = {normalized_height:.2f} | "
        f"closure faces = {len(faces)}"
    )
    _configure_axes(
        ax,
        xlabel="Assembly X [mm]",
        ylabel="Assembly Z [mm]",
        title=title,
    )
    _add_scale_bar(ax, combined[:, 0], combined[:, 1])
    fig.tight_layout()
    _save_figure(fig, path)
    return len(faces)


def _render_xz_section(
    fluid: cq.Solid,
    path: Path,
    *,
    normalized_height: float,
    y_min_mm: float,
    y_max_mm: float,
) -> str:
    design_y_mm = y_min_mm + normalized_height * (y_max_mm - y_min_mm)
    if normalized_height in (0.0, 1.0):
        _render_endpoint_closure_mosaic(
            fluid,
            path,
            normalized_height=normalized_height,
            y_mm=design_y_mm,
        )
        return "endpoint_closure_mosaic"
    sampled_y_mm = design_y_mm
    section = (
        cq.Workplane("XZ", origin=(0.0, sampled_y_mm, 0.0))
        .add(fluid)
        .section()
    )
    wires = section.wires().vals()
    if len(wires) != 1:
        raise RuntimeError(
            f"X-Z section at h={normalized_height:g} produced "
            f"{len(wires)} wires; exactly one outer fluid loop is required"
        )
    combined = _ordered_wire_points(
        wires[0],
        samples_per_edge=320,
    )
    fig, ax = plt.subplots(figsize=(9, 9), facecolor="#050b10")
    ax.fill(
        combined[:, 0],
        combined[:, 1],
        color="#0c536f",
        alpha=0.78,
    )
    ax.plot(
        combined[:, 0],
        combined[:, 1],
        color="#75daf6",
        linewidth=1.4,
    )
    title = (
        "NASA LIQLEV fluid domain — X-Z section\n"
        f"Assembly axes | absolute Y = {design_y_mm:.6f} mm | "
        f"normalized h = {normalized_height:.2f}"
    )
    _configure_axes(
        ax,
        xlabel="Assembly X [mm]",
        ylabel="Assembly Z [mm]",
        title=title,
    )
    _add_scale_bar(ax, combined[:, 0], combined[:, 1])
    fig.tight_layout()
    _save_figure(fig, path)
    return "interior_section"


def _render_sections(
    fluid: cq.Solid,
    section_directory: Path,
    *,
    y_min_mm: float,
    y_max_mm: float,
) -> tuple[list[str], dict[str, str], list[int]]:
    relative_paths: list[str] = []
    render_modes: dict[str, str] = {}
    for filename, plane in SECTION_SPECS:
        _render_projection(
            fluid,
            section_directory / filename,
            plane,
            y_min_mm=y_min_mm,
            y_max_mm=y_max_mm,
        )
        relative = f"sections/{filename}"
        relative_paths.append(relative)
        render_modes[relative] = "orthogonal_projection"
    for percentage, normalized_height in SECTION_HEIGHTS:
        filename = f"section_xz_h{percentage:03d}.png"
        mode = _render_xz_section(
            fluid,
            section_directory / filename,
            normalized_height=normalized_height,
            y_min_mm=y_min_mm,
            y_max_mm=y_max_mm,
        )
        relative = f"sections/{filename}"
        relative_paths.append(relative)
        render_modes[relative] = mode
    endpoint_counts = [
        len(
            _endpoint_closure_faces(
                fluid,
                y_mm=y_mm,
                normal_y=normal_y,
            )
        )
        for y_mm, normal_y in ((y_min_mm, -1.0), (y_max_mm, 1.0))
    ]
    return relative_paths, render_modes, endpoint_counts


def _managed_paths(section_images: list[str]) -> list[Path]:
    return [
        STEP_RELATIVE,
        PACKAGE_RELATIVE,
        PACKAGE_RELATIVE.with_suffix(".json"),
        PACKAGE_RELATIVE.with_suffix(".csv"),
        AUDIT_RELATIVE,
        *(Path("audit") / image for image in section_images),
    ]


def _manifest_relative_paths(values: object) -> list[Path]:
    if not isinstance(values, list):
        raise RuntimeError("promotion manifest paths must be a list")
    paths: list[Path] = []
    for value in values:
        if not isinstance(value, str):
            raise RuntimeError("promotion manifest path must be a string")
        path = Path(value)
        if path.is_absolute() or ".." in path.parts or not path.parts:
            raise RuntimeError(
                f"unsafe promotion manifest relative path: {value!r}"
            )
        paths.append(path)
    return paths


def _restore_promotion_backup(
    output_root: Path,
    backup_root: Path,
    *,
    managed_paths: list[Path],
    originally_existing: set[Path],
) -> None:
    for relative in reversed(managed_paths):
        backup = backup_root / relative
        target = output_root / relative
        if backup.is_file():
            if target.is_file():
                target.unlink()
            target.parent.mkdir(parents=True, exist_ok=True)
            os.replace(backup, target)
        elif relative not in originally_existing and target.is_file():
            target.unlink()


def _recover_interrupted_promotions(output_root: Path) -> int:
    output_root = output_root.resolve()
    parent = output_root.parent
    recovered = 0
    for backup_root in parent.glob(".nasa-tank-backup-*"):
        if not backup_root.is_dir():
            continue
        manifest_path = backup_root / "promotion-manifest.json"
        if not manifest_path.is_file():
            continue
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest_output = Path(manifest["output_root"]).resolve()
            managed_paths = _manifest_relative_paths(
                manifest["managed_paths"]
            )
            originally_existing = set(
                _manifest_relative_paths(
                    manifest["originally_existing"]
                )
            )
        except (OSError, KeyError, TypeError, ValueError) as exc:
            raise RuntimeError(
                f"invalid interrupted-promotion manifest: {manifest_path}"
            ) from exc
        if manifest_output != output_root:
            continue
        if backup_root.resolve().parent != parent.resolve():
            raise RuntimeError(
                f"unsafe interrupted-promotion backup path: {backup_root}"
            )
        _restore_promotion_backup(
            output_root,
            backup_root,
            managed_paths=managed_paths,
            originally_existing=originally_existing,
        )
        shutil.rmtree(backup_root)
        recovered += 1
    return recovered


def _promote_staged_outputs(
    staging_root: Path,
    output_root: Path,
    managed_paths: list[Path],
) -> None:
    output_root.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(
        prefix=".nasa-tank-backup-",
        dir=output_root.parent,
    ) as backup_text:
        backup_root = Path(backup_text)
        originally_existing = {
            relative
            for relative in managed_paths
            if (output_root / relative).is_file()
        }
        (backup_root / "promotion-manifest.json").write_text(
            json.dumps(
                {
                    "output_root": str(output_root.resolve()),
                    "managed_paths": [
                        relative.as_posix() for relative in managed_paths
                    ],
                    "originally_existing": [
                        relative.as_posix()
                        for relative in managed_paths
                        if relative in originally_existing
                    ],
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        try:
            for relative in managed_paths:
                staged = staging_root / relative
                target = output_root / relative
                backup = backup_root / relative
                if not staged.is_file():
                    raise RuntimeError(f"staged artifact is missing: {staged}")
                target.parent.mkdir(parents=True, exist_ok=True)
                if target.exists():
                    backup.parent.mkdir(parents=True, exist_ok=True)
                    os.replace(target, backup)
                os.replace(staged, target)
        except BaseException:
            _restore_promotion_backup(
                output_root,
                backup_root,
                managed_paths=managed_paths,
                originally_existing=originally_existing,
            )
            raise


def _build_staged(
    *,
    source_step: Path,
    staging_root: Path,
    source_sha256: str,
    product_name: str,
    y_min_mm: float,
    y_max_mm: float,
    max_nodes: int,
    visual_section_review: str,
) -> tuple[list[Path], dict[str, object]]:
    timings: dict[str, float] = {}
    command_started = time.perf_counter()
    source_stat = source_step.stat()
    source_before_hash = _sha256(source_step)
    if source_before_hash != source_sha256:
        raise RuntimeError(
            "source SHA-256 mismatch before CAD processing: "
            f"expected {source_sha256}, got {source_before_hash}"
        )

    tank_body = _timed(
        timings,
        "product_selection",
        lambda: load_named_product(
            source_step,
            product_name=product_name,
            expected_sha256=source_sha256,
        ),
    )
    fluid = _timed(
        timings,
        "fluid_construction",
        lambda: build_fluid_domain(
            tank_body,
            y_min_mm=y_min_mm,
            y_max_mm=y_max_mm,
            plane_tolerance_mm=PLANE_TOLERANCE_MM,
        ),
    )
    step_path = staging_root / STEP_RELATIVE
    base_audit_path = staging_root / AUDIT_RELATIVE
    imported, round_trip = _timed(
        timings,
        "step_round_trip",
        lambda: write_step_round_trip(
            fluid,
            step_path,
            base_audit_path,
        ),
    )
    _timed(
        timings,
        "step_timestamp_normalization",
        lambda: _normalize_step_header_timestamp(step_path),
    )
    normalized_step_hash = _sha256(step_path)
    round_trip = replace(
        round_trip,
        output_sha256=normalized_step_hash,
    )
    imported = _timed(
        timings,
        "normalized_step_reimport",
        lambda: _reimport_single_step(step_path),
    )
    exact_brep_metrics = _snapshot_brep_metrics(imported)
    fluid_step_hash = normalized_step_hash
    metadata = GeometryMetadata(
        schema_version=1,
        geometry_id=GEOMETRY_ID,
        source_step_sha256=source_sha256,
        fluid_step_sha256=fluid_step_hash,
        axis="+Y",
        gravity_direction="-Y",
        length_unit="ft",
        area_unit="ft^2",
        volume_unit="ft^3",
        y_min_mm=y_min_mm,
        y_max_mm=y_max_mm,
    )
    kernel = _timed(
        timings,
        "adaptive_measurement_and_coefficients",
        lambda: build_geometry_kernel(
            imported,
            metadata=metadata,
            max_nodes=max_nodes,
        ),
    )
    package_path = staging_root / PACKAGE_RELATIVE
    _timed(
        timings,
        "package_save",
        lambda: save_geometry_package(kernel, package_path),
    )
    independent_measurement = _timed(
        timings,
        "independent_cad_table_measurement",
        lambda: independent_cad_table_metrics(imported, kernel),
    )
    if (
        independent_measurement["max_midpoint_volume_relative"]
        > REFINEMENT_RELATIVE_TOLERANCE
        or independent_measurement["max_midpoint_sidewall_relative"]
        > REFINEMENT_RELATIVE_TOLERANCE
        or independent_measurement["cad_endpoint_volume_relative"]
        > REFINEMENT_RELATIVE_TOLERANCE
        or independent_measurement["integrated_dv_dh_relative"]
        > REFINEMENT_RELATIVE_TOLERANCE
        or independent_measurement["max_direct_area_relative"] > 2e-3
    ):
        raise RuntimeError(
            "independent CAD/table errors exceed tolerance: "
            f"{independent_measurement}"
        )
    section_images, section_render_modes, endpoint_mosaic_face_counts = _timed(
        timings,
        "section_rendering",
        lambda: _render_sections(
            imported,
            base_audit_path.parent / "sections",
            y_min_mm=y_min_mm,
            y_max_mm=y_max_mm,
        ),
    )

    source_after_hash = _sha256(source_step)
    source_after_stat = source_step.stat()
    source_unchanged = (
        source_after_hash == source_before_hash
        and source_after_stat.st_size == source_stat.st_size
        and source_after_stat.st_mtime_ns == source_stat.st_mtime_ns
    )
    if not source_unchanged:
        raise RuntimeError("immutable source STEP changed during the build")

    artifact_hashes = {
        "fluid_step": _sha256(step_path),
        "geometry_npz": _sha256(package_path),
        "geometry_json": _sha256(package_path.with_suffix(".json")),
        "geometry_csv": _sha256(package_path.with_suffix(".csv")),
    }
    for relative in section_images:
        artifact_hashes[relative] = _sha256(
            base_audit_path.parent / relative
        )

    audit_payload = {
        "schema": "liqlev.geometry.nasa_tank_audit",
        "version": 1,
        "passed": visual_section_review == "passed",
        "technical_validation_passed": True,
        "source_sha256": source_sha256,
        "source_size_bytes": source_stat.st_size,
        "source_path": str(source_step),
        "product_name": product_name,
        "closure_planes_mm": [y_min_mm, y_max_mm],
        "cap_plane_max_error_mm": round_trip.cap_plane_error_mm,
        "round_trip_relative_volume_error": (
            round_trip.relative_volume_difference
        ),
        "round_trip_absolute_volume_error_mm3": (
            round_trip.absolute_volume_difference_mm3
        ),
        **exact_brep_metrics,
        "node_count": len(kernel.height_ft),
        "max_nodes": max_nodes,
        "total_height_ft": kernel.total_height_ft,
        "total_volume_ft3": kernel.total_volume_ft3,
        "refinement_relative_tolerance": REFINEMENT_RELATIVE_TOLERANCE,
        "direct_area_relative_tolerance": 2e-3,
        "independent_measurement": independent_measurement,
        "axis": "+Y",
        "gravity_direction": "-Y",
        "artifact_sha256": artifact_hashes,
        "section_images": section_images,
        "section_render_modes": section_render_modes,
        "endpoint_mosaic_face_counts": endpoint_mosaic_face_counts,
        "visual_section_review": visual_section_review,
        "visual_review_evidence": {
            "review_scope": section_images,
            "reviewed_image_sha256": (
                {
                    relative: artifact_hashes[relative]
                    for relative in section_images
                }
                if visual_section_review == "passed"
                else {}
            ),
            "acceptance": (
                "No baffles, internal bodies, open ports, or outer-wall "
                "regions are present in the fluid solid."
                if visual_section_review == "passed"
                else "Awaiting human visual inspection."
            ),
        },
        "source_preservation": {
            "before_sha256": source_before_hash,
            "after_sha256": source_after_hash,
            "before_size_bytes": source_stat.st_size,
            "after_size_bytes": source_after_stat.st_size,
            "before_mtime_ns": source_stat.st_mtime_ns,
            "after_mtime_ns": source_after_stat.st_mtime_ns,
            "unchanged": source_unchanged,
        },
        "round_trip_audit": asdict(round_trip),
        "tool_versions": {
            "python": platform.python_version(),
            "cadquery": getattr(cq, "__version__", "unknown"),
            "ocp": getattr(OCP, "__version__", "unknown"),
            "numpy": np.__version__,
            "matplotlib": matplotlib.__version__,
            "platform": platform.platform(),
        },
        "build_command": {
            "executable": sys.executable,
            "argv": [sys.executable, *sys.argv],
            "working_directory": str(Path.cwd()),
        },
        "commands": [
            " ".join([sys.executable, *sys.argv]),
            "python scripts/check_nasa_tank_geometry.py",
            "python -m pytest tests/cad/test_nasa_tank_artifacts.py -q",
        ],
        "timings_seconds": timings,
    }
    base_audit_path.write_text(
        json.dumps(
            audit_payload,
            allow_nan=False,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    validation_started = time.perf_counter()
    checks = verify_geometry_root(
        source_step=source_step,
        geometry_root=staging_root,
        expected_source_sha256=source_sha256,
        require_visual_review=visual_section_review == "passed",
        precomputed_measurement=independent_measurement,
    )
    failed = [check for check in checks if not check.passed]
    if failed:
        details = "; ".join(f"{check.name}: {check.detail}" for check in failed)
        raise RuntimeError(f"staged artifact validation failed: {details}")
    timings["staged_validation"] = round(
        time.perf_counter() - validation_started,
        6,
    )
    timings["total_before_promotion"] = round(
        time.perf_counter() - command_started,
        6,
    )
    audit_payload["timings_seconds"] = timings
    base_audit_path.write_text(
        json.dumps(
            audit_payload,
            allow_nan=False,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    final_checks = verify_geometry_root(
        source_step=source_step,
        geometry_root=staging_root,
        expected_source_sha256=source_sha256,
        require_visual_review=visual_section_review == "passed",
        precomputed_measurement=independent_measurement,
    )
    final_failed = [check for check in final_checks if not check.passed]
    if final_failed:
        details = "; ".join(
            f"{check.name}: {check.detail}" for check in final_failed
        )
        raise RuntimeError(f"final staged artifact validation failed: {details}")
    return _managed_paths(section_images), audit_payload


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-step", type=Path, default=DEFAULT_SOURCE_STEP)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument(
        "--source-sha256",
        default=DEFAULT_SOURCE_SHA256,
    )
    parser.add_argument("--product-name", default=DEFAULT_PRODUCT_NAME)
    parser.add_argument("--y-min-mm", type=float, default=DEFAULT_Y_MIN_MM)
    parser.add_argument("--y-max-mm", type=float, default=DEFAULT_Y_MAX_MM)
    parser.add_argument("--max-nodes", type=int, default=DEFAULT_MAX_NODES)
    parser.add_argument(
        "--visual-section-review",
        choices=("pending", "passed"),
        default="pending",
        help=(
            "Use 'passed' only after visually inspecting an identical staged "
            "section set."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    source_sha256 = args.source_sha256.upper()
    if (
        not math.isclose(
            args.y_min_mm,
            DEFAULT_Y_MIN_MM,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        or not math.isclose(
            args.y_max_mm,
            DEFAULT_Y_MAX_MM,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
    ):
        print(
            "ERROR only the approved NASA closure planes are supported by "
            "the reviewed round-trip audit.",
            file=sys.stderr,
        )
        return 1
    output_root = args.output_root.resolve()
    output_root.parent.mkdir(parents=True, exist_ok=True)
    try:
        recovered = _recover_interrupted_promotions(output_root)
        if recovered:
            print(
                f"RECOVERY restored {recovered} interrupted promotion(s)",
                file=sys.stderr,
                flush=True,
            )
        with TemporaryDirectory(
            prefix=".nasa-tank-staging-",
            dir=output_root.parent,
        ) as staging_text:
            staging_root = Path(staging_text) / "geometry"
            managed_paths, audit = _build_staged(
                source_step=args.source_step,
                staging_root=staging_root,
                source_sha256=source_sha256,
                product_name=args.product_name,
                y_min_mm=args.y_min_mm,
                y_max_mm=args.y_max_mm,
                max_nodes=args.max_nodes,
                visual_section_review=args.visual_section_review,
            )
            _promote_staged_outputs(
                staging_root,
                output_root,
                managed_paths,
            )
    except Exception as exc:
        print(f"ERROR NASA tank geometry build failed: {exc}", file=sys.stderr)
        return 1
    print(
        "NASA tank geometry artifacts built and staged validation passed: "
        f"nodes={audit['node_count']}, output_root={output_root}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
