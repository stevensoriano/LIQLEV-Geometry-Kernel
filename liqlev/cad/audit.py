"""Deterministic AP214 STEP round-trip validation for fluid-domain solids."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from tempfile import TemporaryDirectory

import cadquery as cq
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepTools import BRepTools
from OCP.GeomAbs import GeomAbs_Plane
from OCP.IFSelect import IFSelect_RetDone
from OCP.Interface import Interface_Static
from OCP.STEPControl import STEPControl_AsIs, STEPControl_Writer

from .fluid_domain import FluidDomainError
from .xcaf import sha256_file


SCHEMA = "liqlev.cad.audit"
VERSION = 2
Y_MIN_MM = -275.406791
Y_MAX_MM = 275.296791
PLANE_TOLERANCE_MM = 1e-5
ABSOLUTE_VOLUME_TOLERANCE_MM3 = 1e-3
RELATIVE_VOLUME_TOLERANCE = 1e-8
SOURCE_PROVENANCE = (
    Path(__file__).resolve().parents[2] / "geometry" / "source" / "PROVENANCE.json"
)


@dataclass(frozen=True)
class CadAudit:
    schema: str
    version: int
    passed: bool
    step_schema: str
    source_sha256: str
    source_size_bytes: int
    input_sha256: str
    output_sha256: str
    pre_volume_mm3: float
    post_volume_mm3: float
    pre_surface_area_mm2: float
    post_surface_area_mm2: float
    pre_bounding_box_mm: tuple[float, float, float, float, float, float]
    post_bounding_box_mm: tuple[float, float, float, float, float, float]
    pre_solid_count: int
    post_solid_count: int
    pre_shell_count: int
    post_shell_count: int
    pre_face_count: int
    post_face_count: int
    pre_valid: bool
    post_valid: bool
    pre_closed_outer_shell: bool
    post_closed_outer_shell: bool
    closure_planes_mm: tuple[float, float]
    cap_plane_error_mm: float
    absolute_volume_difference_mm3: float
    relative_volume_difference: float
    surface_area_difference_mm2: float
    bounding_box_differences_mm: tuple[float, float, float, float, float, float]


@dataclass(frozen=True)
class _ShapeMetrics:
    volume_mm3: float
    surface_area_mm2: float
    bounds_mm: tuple[float, float, float, float, float, float]
    solid_count: int
    shell_count: int
    face_count: int
    valid: bool
    closed_outer_shell: bool


def _bounds(shape: cq.Shape) -> tuple[float, float, float, float, float, float]:
    box = shape.BoundingBox()
    return box.xmin, box.xmax, box.ymin, box.ymax, box.zmin, box.zmax


def _metrics(shape: object, *, label: str) -> _ShapeMetrics:
    if not isinstance(shape, cq.Solid):
        shape_type = (
            shape.ShapeType() if isinstance(shape, cq.Shape) else type(shape).__name__
        )
        solid_count = len(shape.Solids()) if isinstance(shape, cq.Shape) else 0
        raise FluidDomainError(
            f"{label} must be one CadQuery solid; observed "
            f"shape_type={shape_type}, solids={solid_count}"
        )
    solids = shape.Solids()
    shells = shape.Shells()
    volume = shape.Volume()
    area = shape.Area()
    valid = shape.isValid() and BRepCheck_Analyzer(shape.wrapped).IsValid()
    closed = len(shells) == 1 and shells[0].wrapped.Closed()
    if (
        len(solids) != 1
        or len(shells) != 1
        or not closed
        or not valid
        or not math.isfinite(volume)
        or volume <= 0.0
        or not math.isfinite(area)
        or area <= 0.0
    ):
        raise FluidDomainError(
            f"{label} failed solid acceptance: solids={len(solids)}, "
            f"shells={len(shells)}, faces={len(shape.Faces())}, "
            f"closed={closed}, valid={valid}, volume_mm3={volume}, "
            f"surface_area_mm2={area}"
        )
    return _ShapeMetrics(
        volume_mm3=volume,
        surface_area_mm2=area,
        bounds_mm=_bounds(shape),
        solid_count=len(solids),
        shell_count=len(shells),
        face_count=len(shape.Faces()),
        valid=valid,
        closed_outer_shell=closed,
    )


def _closure_plane_error(shape: cq.Solid, *, label: str) -> float:
    y_normal_planes: list[float] = []
    for face in shape.Faces():
        surface = BRepAdaptor_Surface(face.wrapped, True)
        if surface.GetType() != GeomAbs_Plane:
            continue
        plane = surface.Plane()
        direction = plane.Axis().Direction()
        if (
            math.hypot(direction.X(), direction.Z()) > 1e-8
            or abs(abs(direction.Y()) - 1.0) > 1e-8
        ):
            continue
        y_normal_planes.append(plane.Location().Y())

    deviations: list[float] = []
    candidate_counts: list[int] = []
    for closure_y in (Y_MIN_MM, Y_MAX_MM):
        candidates = [
            abs(plane_y - closure_y)
            for plane_y in y_normal_planes
            if abs(plane_y - closure_y) <= PLANE_TOLERANCE_MM
        ]
        candidate_counts.append(len(candidates))
        deviations.extend(candidates)
    if any(count == 0 for count in candidate_counts):
        raise FluidDomainError(
            f"{label} requires planar assembly-Y-normal closure faces at "
            f"y=({Y_MIN_MM}, {Y_MAX_MM}); observed candidate_counts="
            f"{tuple(candidate_counts)}, Y-normal plane ordinates="
            f"{tuple(sorted(y_normal_planes))}"
        )
    maximum_deviation = max(deviations)
    if maximum_deviation > PLANE_TOLERANCE_MM:
        raise FluidDomainError(
            f"{label} closure-face plane deviation exceeds tolerance: "
            f"maximum_deviation_mm={maximum_deviation}, "
            f"tolerance_mm={PLANE_TOLERANCE_MM}"
        )
    return maximum_deviation


def _source_authority() -> tuple[str, int, str]:
    try:
        provenance = json.loads(SOURCE_PROVENANCE.read_text(encoding="utf-8"))
        source_sha256 = str(provenance["source_sha256"]).upper()
        source_size = int(provenance["source_size_bytes"])
        source_name = Path(provenance["source_path"]).name
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise FluidDomainError(
            f"source provenance is missing or invalid: {SOURCE_PROVENANCE}"
        ) from exc
    if len(source_sha256) != 64 or source_size <= 0:
        raise FluidDomainError(
            "source provenance authority is invalid: "
            f"sha256={source_sha256!r}, size_bytes={source_size}"
        )
    return source_sha256, source_size, source_name


def _input_brep_sha256(fluid: cq.Solid) -> str:
    with TemporaryDirectory(prefix="liqlev-fluid-audit-") as temporary:
        path = Path(temporary) / "input.brep"
        if not BRepTools.Write_s(fluid.wrapped, str(path)):
            raise FluidDomainError("failed to serialize input B-Rep for audit hashing")
        return sha256_file(path)


def _write_ap214(fluid: cq.Solid, output_path: Path) -> None:
    writer = STEPControl_Writer()
    if not Interface_Static.SetCVal_s("write.step.schema", "AP214IS"):
        raise FluidDomainError("OpenCascade rejected AP214IS STEP schema selection")
    transfer_status = writer.Transfer(fluid.wrapped, STEPControl_AsIs)
    if transfer_status != IFSelect_RetDone:
        raise FluidDomainError(
            f"AP214 STEP transfer failed with status {transfer_status}"
        )
    write_status = writer.Write(str(output_path))
    if write_status != IFSelect_RetDone:
        raise FluidDomainError(
            f"AP214 STEP write failed with status {write_status}: {output_path}"
        )


def _independent_import(output_path: Path) -> cq.Solid:
    try:
        imported = cq.importers.importStep(str(output_path))
        solids = imported.solids().vals()
    except Exception as exc:
        raise FluidDomainError(
            f"independent STEP re-import failed: {output_path}"
        ) from exc
    if len(solids) != 1:
        raise FluidDomainError(
            f"independent STEP re-import must contain one solid; observed "
            f"{len(solids)} at {output_path}"
        )
    solid = solids[0]
    if not isinstance(solid, cq.Solid):
        raise FluidDomainError(
            "independent STEP re-import did not return a CadQuery Solid: "
            f"{type(solid).__name__}"
        )
    return solid


def write_step_round_trip(
    fluid: cq.Solid,
    output_path: str | Path,
    audit_path: str | Path,
) -> tuple[cq.Solid, CadAudit]:
    """Export AP214, independently re-import it, and enforce all CAD gates."""

    pre = _metrics(fluid, label="pre-export fluid")
    pre_closure_plane_error = _closure_plane_error(
        fluid,
        label="pre-export fluid",
    )
    output = Path(output_path)
    audit_file = Path(audit_path)
    source_sha256, source_size, source_name = _source_authority()
    try:
        if output.resolve() == audit_file.resolve():
            raise FluidDomainError("STEP output_path and audit_path must differ")
    except OSError as exc:
        raise FluidDomainError("STEP or audit path cannot be resolved") from exc
    if output.name.casefold() == source_name.casefold():
        raise FluidDomainError(
            f"output STEP filename must differ from source filename: {output.name}"
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    audit_file.parent.mkdir(parents=True, exist_ok=True)
    input_sha256 = _input_brep_sha256(fluid)
    _write_ap214(fluid, output)
    output_sha256 = sha256_file(output)
    imported = _independent_import(output)
    post = _metrics(imported, label="post-import fluid")
    post_closure_plane_error = _closure_plane_error(
        imported,
        label="post-import fluid",
    )

    absolute_volume_difference = abs(post.volume_mm3 - pre.volume_mm3)
    relative_volume_difference = absolute_volume_difference / pre.volume_mm3
    bounding_box_differences = tuple(
        abs(after - before)
        for before, after in zip(pre.bounds_mm, post.bounds_mm)
    )
    surface_area_difference = abs(
        post.surface_area_mm2 - pre.surface_area_mm2
    )
    closure_bounds_error = max(
        abs(pre.bounds_mm[2] - Y_MIN_MM),
        abs(pre.bounds_mm[3] - Y_MAX_MM),
        abs(post.bounds_mm[2] - Y_MIN_MM),
        abs(post.bounds_mm[3] - Y_MAX_MM),
    )
    cap_plane_error = max(
        pre_closure_plane_error,
        post_closure_plane_error,
    )
    volume_tolerance = max(
        ABSOLUTE_VOLUME_TOLERANCE_MM3,
        pre.volume_mm3 * RELATIVE_VOLUME_TOLERANCE,
    )
    differences = {
        "absolute_volume_difference_mm3": absolute_volume_difference,
        "relative_volume_difference": relative_volume_difference,
        "volume_tolerance_mm3": volume_tolerance,
        "bounding_box_differences_mm": bounding_box_differences,
        "cap_plane_error_mm": cap_plane_error,
        "closure_bounds_error_mm": closure_bounds_error,
        "pre_counts": (
            pre.solid_count,
            pre.shell_count,
            pre.face_count,
        ),
        "post_counts": (
            post.solid_count,
            post.shell_count,
            post.face_count,
        ),
    }
    passed = (
        absolute_volume_difference <= volume_tolerance
        and max(bounding_box_differences) <= PLANE_TOLERANCE_MM
        and cap_plane_error <= PLANE_TOLERANCE_MM
        and closure_bounds_error <= PLANE_TOLERANCE_MM
        and pre.solid_count == post.solid_count == 1
        and pre.shell_count == post.shell_count == 1
        and pre.face_count == post.face_count
        and pre.valid
        and post.valid
        and pre.closed_outer_shell
        and post.closed_outer_shell
    )
    audit = CadAudit(
        schema=SCHEMA,
        version=VERSION,
        passed=passed,
        step_schema="AP214IS",
        source_sha256=source_sha256,
        source_size_bytes=source_size,
        input_sha256=input_sha256,
        output_sha256=output_sha256,
        pre_volume_mm3=pre.volume_mm3,
        post_volume_mm3=post.volume_mm3,
        pre_surface_area_mm2=pre.surface_area_mm2,
        post_surface_area_mm2=post.surface_area_mm2,
        pre_bounding_box_mm=pre.bounds_mm,
        post_bounding_box_mm=post.bounds_mm,
        pre_solid_count=pre.solid_count,
        post_solid_count=post.solid_count,
        pre_shell_count=pre.shell_count,
        post_shell_count=post.shell_count,
        pre_face_count=pre.face_count,
        post_face_count=post.face_count,
        pre_valid=pre.valid,
        post_valid=post.valid,
        pre_closed_outer_shell=pre.closed_outer_shell,
        post_closed_outer_shell=post.closed_outer_shell,
        closure_planes_mm=(Y_MIN_MM, Y_MAX_MM),
        cap_plane_error_mm=cap_plane_error,
        absolute_volume_difference_mm3=absolute_volume_difference,
        relative_volume_difference=relative_volume_difference,
        surface_area_difference_mm2=surface_area_difference,
        bounding_box_differences_mm=bounding_box_differences,
    )
    audit_file.write_text(
        json.dumps(
            asdict(audit),
            allow_nan=False,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    if not passed:
        raise FluidDomainError(f"STEP round-trip acceptance failed: {differences}")
    return imported, audit
