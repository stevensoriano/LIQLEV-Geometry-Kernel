"""Exact OpenCascade construction of the closure-bounded tank fluid domain."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

import cadquery as cq
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Splitter
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepClass3d import BRepClass3d_SolidClassifier
from OCP.TopAbs import TopAbs_IN, TopAbs_ON, TopAbs_OUT
from OCP.TopTools import TopTools_ListOfShape
from OCP.gp import gp_Pnt


NORMAL_PARALLEL_TOLERANCE = 1e-8
PRIMARY_PADDING_MM = 25.0
SECONDARY_PADDING_MM = 100.0
ABSOLUTE_MEASURE_TOLERANCE = 1e-3
RELATIVE_MEASURE_TOLERANCE = 1e-8


class FluidDomainError(RuntimeError):
    """Raised when exact fluid-domain construction or validation fails."""


@dataclass(frozen=True)
class _LoopInventory:
    wire: cq.Wire
    area_mm2: float
    center_mm: tuple[float, float, float]


@dataclass(frozen=True)
class _RimInventory:
    closure_y_mm: float
    candidate_face_count: int
    wire_count: int
    outer_loop_area_mm2: float
    outer_loop_center_mm: tuple[float, float, float]
    wet_loop_area_mm2: float
    wet_loop_center_mm: tuple[float, float, float]
    face: cq.Face
    wet_loop: cq.Wire
    loops: tuple[_LoopInventory, ...]


@dataclass(frozen=True)
class _DirectCutDiagnostics:
    padding_mm: float
    result_solid_count: int
    result_shell_count: int
    result_face_count: int
    result_valid: bool
    seed_classifications: tuple[str, ...]
    selected_solid_count: int
    selected_shell_count: int
    selected_face_count: int
    selected_valid: bool
    volume_mm3: float
    surface_area_mm2: float
    bounds_mm: tuple[float, float, float, float, float, float]
    carrier_bounds_mm: tuple[float, float, float, float, float, float]
    carrier_xz_contact: bool
    exterior_probe_classifications: tuple[str, str, str, str]
    shared_face_count: int
    shared_edge_count: int
    shared_vertex_count: int
    volume_difference_mm3: float
    surface_area_difference_mm2: float
    bounding_box_differences_mm: tuple[float, float, float, float, float, float]
    selected: cq.Solid
    rejected: tuple[cq.Solid, ...]


@dataclass(frozen=True)
class _SplitterDiagnostics:
    result_solid_count: int
    seed_classifications: tuple[str, ...]
    selected_face_count: int
    selected_valid: bool
    volume_difference_mm3: float
    surface_area_difference_mm2: float
    bounding_box_differences_mm: tuple[float, float, float, float, float, float]


@dataclass(frozen=True)
class _ConstructionDiagnostics:
    fluid: cq.Solid
    seed_mm: tuple[float, float, float]
    rims: tuple[_RimInventory, _RimInventory]
    primary: _DirectCutDiagnostics
    secondary: _DirectCutDiagnostics
    splitter: _SplitterDiagnostics
    minimum_closure_face_count: int
    maximum_closure_face_count: int
    maximum_closure_plane_error_mm: float


def _bounds(shape: cq.Shape) -> tuple[float, float, float, float, float, float]:
    box = shape.BoundingBox()
    return box.xmin, box.xmax, box.ymin, box.ymax, box.zmin, box.zmax


def _measure_tolerance(value: float) -> float:
    return max(ABSOLUTE_MEASURE_TOLERANCE, abs(value) * RELATIVE_MEASURE_TOLERANCE)


def _state_name(state: object) -> str:
    if state == TopAbs_IN:
        return "IN"
    if state == TopAbs_OUT:
        return "OUT"
    if state == TopAbs_ON:
        return "ON"
    return f"UNKNOWN({state})"


def _classify(solid: cq.Solid, point_mm: tuple[float, float, float]) -> str:
    classifier = BRepClass3d_SolidClassifier(
        solid.wrapped,
        gp_Pnt(*point_mm),
        1e-9,
    )
    return _state_name(classifier.State())


def _validate_inputs(
    tank_body: object,
    *,
    y_min_mm: float,
    y_max_mm: float,
    plane_tolerance_mm: float,
) -> cq.Solid:
    if not math.isfinite(y_min_mm) or not math.isfinite(y_max_mm):
        raise FluidDomainError(
            f"y_min_mm and y_max_mm must be finite; observed "
            f"{y_min_mm!r}, {y_max_mm!r}"
        )
    if y_min_mm >= y_max_mm:
        raise FluidDomainError(
            f"y_min_mm must be less than y_max_mm; observed "
            f"{y_min_mm!r} >= {y_max_mm!r}"
        )
    if not math.isfinite(plane_tolerance_mm) or plane_tolerance_mm <= 0.0:
        raise FluidDomainError(
            "plane_tolerance_mm must be finite and positive; observed "
            f"{plane_tolerance_mm!r}"
        )
    if not isinstance(tank_body, cq.Shape):
        raise FluidDomainError(
            "tank_body must be one CadQuery Solid; observed "
            f"{type(tank_body).__name__}"
        )
    solids = tank_body.Solids()
    if not isinstance(tank_body, cq.Solid) or len(solids) != 1:
        raise FluidDomainError(
            "expected one input tank solid; observed "
            f"shape_type={tank_body.ShapeType()}, solids={len(solids)}"
        )
    if not tank_body.isValid() or not BRepCheck_Analyzer(tank_body.wrapped).IsValid():
        raise FluidDomainError(
            "input tank solid has an invalid source B-Rep; observed "
            f"cadquery_valid={tank_body.isValid()}"
        )
    volume = tank_body.Volume()
    if not math.isfinite(volume) or volume <= 0.0:
        raise FluidDomainError(
            f"input tank solid volume must be finite and positive; observed {volume}"
        )
    return tank_body


def _wire_inventory(wire: cq.Wire, *, closure_y_mm: float) -> _LoopInventory:
    if not wire.IsClosed():
        raise FluidDomainError(
            f"rim face at y={closure_y_mm} contains an open wire"
        )
    try:
        planar_face = cq.Face.makeFromWires(wire)
        area = planar_face.Area()
        center = planar_face.Center().toTuple()
    except Exception as exc:
        raise FluidDomainError(
            f"failed to evaluate exact planar wire at y={closure_y_mm}"
        ) from exc
    if not math.isfinite(area) or area <= 0.0:
        raise FluidDomainError(
            f"rim wire area at y={closure_y_mm} is not finite and positive: {area}"
        )
    return _LoopInventory(
        wire=wire,
        area_mm2=area,
        center_mm=tuple(float(value) for value in center),
    )


def _select_wet_loop(
    loops: Iterable[_LoopInventory],
    *,
    plane_tolerance_mm: float,
    closure_y_mm: float,
) -> tuple[_LoopInventory, _LoopInventory]:
    ordered = tuple(sorted(loops, key=lambda item: item.area_mm2, reverse=True))
    if len(ordered) < 2:
        raise FluidDomainError(
            f"wet-loop selection at y={closure_y_mm} requires at least two "
            f"closed wires; observed {len(ordered)}"
        )
    outer = ordered[0]
    concentric = tuple(
        loop
        for loop in ordered[1:]
        if math.hypot(
            loop.center_mm[0] - outer.center_mm[0],
            loop.center_mm[2] - outer.center_mm[2],
        )
        <= plane_tolerance_mm
    )
    if len(concentric) != 1:
        details = [
            {
                "area_mm2": loop.area_mm2,
                "center_mm": loop.center_mm,
            }
            for loop in ordered
        ]
        raise FluidDomainError(
            f"wet-loop candidates at y={closure_y_mm}: observed "
            f"{len(concentric)} concentric remaining loops; loops={details}"
        )
    return outer, concentric[0]


def _find_rim(
    tank_body: cq.Solid,
    *,
    closure_y_mm: float,
    plane_tolerance_mm: float,
) -> _RimInventory:
    candidates: list[cq.Face] = []
    candidate_details: list[dict[str, object]] = []
    for face in tank_body.Faces():
        if face.geomType() != "PLANE":
            continue
        center = face.Center()
        if abs(center.y - closure_y_mm) > plane_tolerance_mm:
            continue
        normal = face.normalAt()
        transverse = math.hypot(normal.x, normal.z)
        candidate_details.append(
            {
                "center_y_mm": center.y,
                "normal": normal.toTuple(),
                "wire_count": len(face.Wires()),
            }
        )
        if transverse <= NORMAL_PARALLEL_TOLERANCE and (
            abs(abs(normal.y) - 1.0) <= NORMAL_PARALLEL_TOLERANCE
        ):
            candidates.append(face)
    if len(candidates) != 1:
        raise FluidDomainError(
            f"rim candidate faces at y={closure_y_mm}: observed "
            f"{len(candidates)}; planar_near_plane={candidate_details}"
        )
    face = candidates[0]
    loops = tuple(
        _wire_inventory(wire, closure_y_mm=closure_y_mm)
        for wire in face.Wires()
    )
    outer, wet = _select_wet_loop(
        loops,
        plane_tolerance_mm=plane_tolerance_mm,
        closure_y_mm=closure_y_mm,
    )
    return _RimInventory(
        closure_y_mm=closure_y_mm,
        candidate_face_count=len(candidates),
        wire_count=len(loops),
        outer_loop_area_mm2=outer.area_mm2,
        outer_loop_center_mm=outer.center_mm,
        wet_loop_area_mm2=wet.area_mm2,
        wet_loop_center_mm=wet.center_mm,
        face=face,
        wet_loop=wet.wire,
        loops=loops,
    )


def _seed_from_rims(
    minimum: _RimInventory,
    maximum: _RimInventory,
    *,
    plane_tolerance_mm: float,
) -> tuple[float, float, float]:
    min_center = minimum.wet_loop_center_mm
    max_center = maximum.wet_loop_center_mm
    x_difference = abs(min_center[0] - max_center[0])
    z_difference = abs(min_center[2] - max_center[2])
    if x_difference > plane_tolerance_mm or z_difference > plane_tolerance_mm:
        raise FluidDomainError(
            "wet-loop centres do not define one assembly-Y axis: "
            f"x_difference_mm={x_difference}, z_difference_mm={z_difference}, "
            f"minimum={min_center}, maximum={max_center}"
        )
    return tuple(
        (minimum_value + maximum_value) / 2.0
        for minimum_value, maximum_value in zip(min_center, max_center)
    )


def _make_carrier(
    source_bounds: tuple[float, float, float, float, float, float],
    *,
    y_min_mm: float,
    y_max_mm: float,
    padding_mm: float,
) -> cq.Solid:
    xmin, xmax, _, _, zmin, zmax = source_bounds
    return cq.Solid.makeBox(
        xmax - xmin + 2.0 * padding_mm,
        y_max_mm - y_min_mm,
        zmax - zmin + 2.0 * padding_mm,
        cq.Vector(xmin - padding_mm, y_min_mm, zmin - padding_mm),
    )


def _shared_topology_counts(
    selected: cq.Solid,
    rejected: Iterable[cq.Solid],
) -> tuple[int, int, int]:
    def count_shared(
        selected_items: Iterable[cq.Shape],
        rejected_items: Iterable[cq.Shape],
    ) -> int:
        return sum(
            1
            for selected_item in selected_items
            if any(
                selected_item.wrapped.IsSame(rejected_item.wrapped)
                for rejected_item in rejected_items
            )
        )

    rejected_tuple = tuple(rejected)
    rejected_faces = tuple(
        face for solid in rejected_tuple for face in solid.Faces()
    )
    rejected_edges = tuple(
        edge for solid in rejected_tuple for edge in solid.Edges()
    )
    rejected_vertices = tuple(
        vertex for solid in rejected_tuple for vertex in solid.Vertices()
    )
    return (
        count_shared(selected.Faces(), rejected_faces),
        count_shared(selected.Edges(), rejected_edges),
        count_shared(selected.Vertices(), rejected_vertices),
    )


def _validate_selected(
    selected: cq.Solid,
    *,
    source_bounds: tuple[float, float, float, float, float, float],
    carrier_bounds: tuple[float, float, float, float, float, float],
    y_min_mm: float,
    y_max_mm: float,
    plane_tolerance_mm: float,
    seed_mm: tuple[float, float, float],
    rejected: tuple[cq.Solid, ...],
) -> tuple[
    tuple[float, float, float, float, float, float],
    bool,
    tuple[str, str, str, str],
    tuple[int, int, int],
]:
    bounds = _bounds(selected)
    solid_count = len(selected.Solids())
    shell_count = len(selected.Shells())
    valid = selected.isValid() and BRepCheck_Analyzer(selected.wrapped).IsValid()
    volume = selected.Volume()
    if (
        not isinstance(selected, cq.Solid)
        or solid_count != 1
        or shell_count != 1
        or not selected.Shells()[0].wrapped.Closed()
        or not valid
        or not math.isfinite(volume)
        or volume <= 0.0
    ):
        raise FluidDomainError(
            "selected fluid failed solid validation: "
            f"type={type(selected).__name__}, solids={solid_count}, "
            f"shells={shell_count}, closed={selected.Shells()[0].wrapped.Closed()}, "
            f"valid={valid}, volume_mm3={volume}, faces={len(selected.Faces())}"
        )
    if (
        abs(bounds[2] - y_min_mm) > plane_tolerance_mm
        or abs(bounds[3] - y_max_mm) > plane_tolerance_mm
    ):
        raise FluidDomainError(
            "selected fluid closure bounds failed: "
            f"observed_y=({bounds[2]}, {bounds[3]}), "
            f"required_y=({y_min_mm}, {y_max_mm}), "
            f"tolerance_mm={plane_tolerance_mm}"
        )
    containment_differences = tuple(
        fluid_value - source_value
        for fluid_value, source_value in zip(bounds, source_bounds)
    )
    if any(
        (
            bounds[index] < source_bounds[index] - plane_tolerance_mm
            if index % 2 == 0
            else bounds[index] > source_bounds[index] + plane_tolerance_mm
        )
        for index in range(6)
    ):
        raise FluidDomainError(
            "selected fluid bounds are not contained in the source tank: "
            f"fluid={bounds}, source={source_bounds}, "
            f"differences={containment_differences}"
        )
    carrier_contact = any(
        abs(bounds[index] - carrier_bounds[index]) <= plane_tolerance_mm
        for index in (0, 1, 4, 5)
    )
    if carrier_contact:
        raise FluidDomainError(
            "selected fluid contacts a carrier X/Z boundary: "
            f"fluid={bounds}, carrier={carrier_bounds}"
        )
    source_xmid = (source_bounds[0] + source_bounds[1]) / 2.0
    source_zmid = (source_bounds[4] + source_bounds[5]) / 2.0
    mid_y = (y_min_mm + y_max_mm) / 2.0
    probes = (
        ((source_bounds[0] + carrier_bounds[0]) / 2.0, mid_y, source_zmid),
        ((source_bounds[1] + carrier_bounds[1]) / 2.0, mid_y, source_zmid),
        (source_xmid, mid_y, (source_bounds[4] + carrier_bounds[4]) / 2.0),
        (source_xmid, mid_y, (source_bounds[5] + carrier_bounds[5]) / 2.0),
    )
    probe_classifications = tuple(_classify(selected, point) for point in probes)
    if probe_classifications != ("OUT", "OUT", "OUT", "OUT"):
        raise FluidDomainError(
            "selected fluid exterior probes are not all OUT: "
            f"probes={probes}, classifications={probe_classifications}"
        )
    shared_counts = _shared_topology_counts(selected, rejected)
    if shared_counts != (0, 0, 0):
        raise FluidDomainError(
            "selected fluid shares topology with rejected direct-cut "
            f"components: faces={shared_counts[0]}, edges={shared_counts[1]}, "
            f"vertices={shared_counts[2]}"
        )
    if _classify(selected, seed_mm) != "IN":
        raise FluidDomainError(
            f"selected fluid does not strictly contain seed {seed_mm}"
        )
    return bounds, carrier_contact, probe_classifications, shared_counts


def _direct_cut(
    tank_body: cq.Solid,
    *,
    source_bounds: tuple[float, float, float, float, float, float],
    y_min_mm: float,
    y_max_mm: float,
    padding_mm: float,
    plane_tolerance_mm: float,
    seed_mm: tuple[float, float, float],
    reference: cq.Solid | None = None,
) -> _DirectCutDiagnostics:
    carrier = _make_carrier(
        source_bounds,
        y_min_mm=y_min_mm,
        y_max_mm=y_max_mm,
        padding_mm=padding_mm,
    )
    operation = BRepAlgoAPI_Cut(carrier.wrapped, tank_body.wrapped)
    operation.Build()
    if not operation.IsDone():
        raise FluidDomainError(
            f"direct cut failed for padding_mm={padding_mm}: "
            "OpenCascade did not complete the operation"
        )
    result = cq.Shape.cast(operation.Shape())
    solids = tuple(result.Solids())
    classifications = tuple(_classify(solid, seed_mm) for solid in solids)
    matching_indices = [
        index for index, state in enumerate(classifications) if state == "IN"
    ]
    if "ON" in classifications or len(matching_indices) != 1:
        raise FluidDomainError(
            "direct-cut seed selection failed: "
            f"padding_mm={padding_mm}, components={len(solids)}, "
            f"seed={seed_mm}, classifications={classifications}"
        )
    selected = solids[matching_indices[0]]
    rejected = tuple(
        solid for index, solid in enumerate(solids) if index != matching_indices[0]
    )
    carrier_bounds = _bounds(carrier)
    selected_bounds, carrier_contact, probes, shared = _validate_selected(
        selected,
        source_bounds=source_bounds,
        carrier_bounds=carrier_bounds,
        y_min_mm=y_min_mm,
        y_max_mm=y_max_mm,
        plane_tolerance_mm=plane_tolerance_mm,
        seed_mm=seed_mm,
        rejected=rejected,
    )
    reference_shape = selected if reference is None else reference
    box_differences = tuple(
        abs(value - expected)
        for value, expected in zip(selected_bounds, _bounds(reference_shape))
    )
    return _DirectCutDiagnostics(
        padding_mm=padding_mm,
        result_solid_count=len(solids),
        result_shell_count=len(result.Shells()),
        result_face_count=len(result.Faces()),
        result_valid=result.isValid() and BRepCheck_Analyzer(result.wrapped).IsValid(),
        seed_classifications=classifications,
        selected_solid_count=len(selected.Solids()),
        selected_shell_count=len(selected.Shells()),
        selected_face_count=len(selected.Faces()),
        selected_valid=selected.isValid(),
        volume_mm3=selected.Volume(),
        surface_area_mm2=selected.Area(),
        bounds_mm=selected_bounds,
        carrier_bounds_mm=carrier_bounds,
        carrier_xz_contact=carrier_contact,
        exterior_probe_classifications=probes,
        shared_face_count=shared[0],
        shared_edge_count=shared[1],
        shared_vertex_count=shared[2],
        volume_difference_mm3=abs(selected.Volume() - reference_shape.Volume()),
        surface_area_difference_mm2=abs(selected.Area() - reference_shape.Area()),
        bounding_box_differences_mm=box_differences,
        selected=selected,
        rejected=rejected,
    )


def _require_invariant(
    *,
    label: str,
    reference: _DirectCutDiagnostics,
    candidate: _DirectCutDiagnostics,
) -> None:
    differences = {
        "solid_count": (
            reference.selected_solid_count,
            candidate.selected_solid_count,
        ),
        "shell_count": (
            reference.selected_shell_count,
            candidate.selected_shell_count,
        ),
        "face_count": (
            reference.selected_face_count,
            candidate.selected_face_count,
        ),
        "valid": (reference.selected_valid, candidate.selected_valid),
        "volume_difference_mm3": candidate.volume_difference_mm3,
        "surface_area_difference_mm2": candidate.surface_area_difference_mm2,
        "bounding_box_differences_mm": candidate.bounding_box_differences_mm,
    }
    if (
        candidate.selected_solid_count != reference.selected_solid_count
        or candidate.selected_shell_count != reference.selected_shell_count
        or candidate.selected_face_count != reference.selected_face_count
        or candidate.selected_valid != reference.selected_valid
        or candidate.volume_difference_mm3
        > _measure_tolerance(reference.volume_mm3)
        or candidate.surface_area_difference_mm2
        > _measure_tolerance(reference.surface_area_mm2)
        or max(candidate.bounding_box_differences_mm) > 1e-5
    ):
        raise FluidDomainError(f"{label} invariance failed: {differences}")


def _splitter_oracle(
    tank_body: cq.Solid,
    *,
    source_bounds: tuple[float, float, float, float, float, float],
    y_min_mm: float,
    y_max_mm: float,
    seed_mm: tuple[float, float, float],
    reference: _DirectCutDiagnostics,
) -> _SplitterDiagnostics:
    carrier = _make_carrier(
        source_bounds,
        y_min_mm=y_min_mm,
        y_max_mm=y_max_mm,
        padding_mm=PRIMARY_PADDING_MM,
    )
    splitter = BRepAlgoAPI_Splitter()
    arguments = TopTools_ListOfShape()
    arguments.Append(carrier.wrapped)
    tools = TopTools_ListOfShape()
    tools.Append(tank_body.wrapped)
    splitter.SetArguments(arguments)
    splitter.SetTools(tools)
    splitter.Build()
    if not splitter.IsDone():
        raise FluidDomainError(
            "exact splitter failed: OpenCascade did not complete the operation"
        )
    result = cq.Shape.cast(splitter.Shape())
    solids = tuple(result.Solids())
    classifications = tuple(_classify(solid, seed_mm) for solid in solids)
    matching_indices = [
        index for index, state in enumerate(classifications) if state == "IN"
    ]
    if "ON" in classifications or len(matching_indices) != 1:
        raise FluidDomainError(
            "splitter seed selection failed: "
            f"components={len(solids)}, seed={seed_mm}, "
            f"classifications={classifications}"
        )
    selected = solids[matching_indices[0]]
    bounds_differences = tuple(
        abs(value - expected)
        for value, expected in zip(_bounds(selected), reference.bounds_mm)
    )
    volume_difference = abs(selected.Volume() - reference.volume_mm3)
    area_difference = abs(selected.Area() - reference.surface_area_mm2)
    diagnostics = _SplitterDiagnostics(
        result_solid_count=len(solids),
        seed_classifications=classifications,
        selected_face_count=len(selected.Faces()),
        selected_valid=selected.isValid(),
        volume_difference_mm3=volume_difference,
        surface_area_difference_mm2=area_difference,
        bounding_box_differences_mm=bounds_differences,
    )
    if (
        not diagnostics.selected_valid
        or diagnostics.selected_face_count != reference.selected_face_count
        or volume_difference > _measure_tolerance(reference.volume_mm3)
        or area_difference > _measure_tolerance(reference.surface_area_mm2)
        or max(bounds_differences) > 1e-5
    ):
        raise FluidDomainError(
            "splitter oracle disagrees with direct cut: "
            f"result_solids={len(solids)}, "
            f"selected_faces={diagnostics.selected_face_count}, "
            f"volume_difference_mm3={volume_difference}, "
            f"surface_area_difference_mm2={area_difference}, "
            f"bounding_box_differences_mm={bounds_differences}"
        )
    return diagnostics


def _closure_mosaic(
    fluid: cq.Solid,
    *,
    y_min_mm: float,
    y_max_mm: float,
    plane_tolerance_mm: float,
) -> tuple[int, int, float]:
    counts = [0, 0]
    maximum_error = 0.0
    for face in fluid.Faces():
        if face.geomType() != "PLANE":
            continue
        face_bounds = _bounds(face)
        for index, closure_y in enumerate((y_min_mm, y_max_mm)):
            errors = (
                abs(face_bounds[2] - closure_y),
                abs(face_bounds[3] - closure_y),
            )
            if max(errors) <= plane_tolerance_mm:
                counts[index] += 1
                maximum_error = max(maximum_error, *errors)
    if counts[0] == 0 or counts[1] == 0 or maximum_error > plane_tolerance_mm:
        raise FluidDomainError(
            "closure-plane mosaic failed: "
            f"minimum_faces={counts[0]}, maximum_faces={counts[1]}, "
            f"maximum_error_mm={maximum_error}"
        )
    return counts[0], counts[1], maximum_error


def _diagnose_exact_construction(
    tank_body: cq.Shape,
    *,
    y_min_mm: float,
    y_max_mm: float,
    plane_tolerance_mm: float = 1e-5,
) -> _ConstructionDiagnostics:
    """Construct the direct-cut solid and retain exact validation diagnostics."""

    tank = _validate_inputs(
        tank_body,
        y_min_mm=y_min_mm,
        y_max_mm=y_max_mm,
        plane_tolerance_mm=plane_tolerance_mm,
    )
    minimum_rim = _find_rim(
        tank,
        closure_y_mm=y_min_mm,
        plane_tolerance_mm=plane_tolerance_mm,
    )
    maximum_rim = _find_rim(
        tank,
        closure_y_mm=y_max_mm,
        plane_tolerance_mm=plane_tolerance_mm,
    )
    seed = _seed_from_rims(
        minimum_rim,
        maximum_rim,
        plane_tolerance_mm=plane_tolerance_mm,
    )
    source_bounds = _bounds(tank)
    primary = _direct_cut(
        tank,
        source_bounds=source_bounds,
        y_min_mm=y_min_mm,
        y_max_mm=y_max_mm,
        padding_mm=PRIMARY_PADDING_MM,
        plane_tolerance_mm=plane_tolerance_mm,
        seed_mm=seed,
    )
    secondary = _direct_cut(
        tank,
        source_bounds=source_bounds,
        y_min_mm=y_min_mm,
        y_max_mm=y_max_mm,
        padding_mm=SECONDARY_PADDING_MM,
        plane_tolerance_mm=plane_tolerance_mm,
        seed_mm=seed,
        reference=primary.selected,
    )
    _require_invariant(
        label="25/100 mm padding",
        reference=primary,
        candidate=secondary,
    )
    splitter = _splitter_oracle(
        tank,
        source_bounds=source_bounds,
        y_min_mm=y_min_mm,
        y_max_mm=y_max_mm,
        seed_mm=seed,
        reference=primary,
    )
    minimum_faces, maximum_faces, closure_error = _closure_mosaic(
        primary.selected,
        y_min_mm=y_min_mm,
        y_max_mm=y_max_mm,
        plane_tolerance_mm=plane_tolerance_mm,
    )
    return _ConstructionDiagnostics(
        fluid=primary.selected,
        seed_mm=seed,
        rims=(minimum_rim, maximum_rim),
        primary=primary,
        secondary=secondary,
        splitter=splitter,
        minimum_closure_face_count=minimum_faces,
        maximum_closure_face_count=maximum_faces,
        maximum_closure_plane_error_mm=closure_error,
    )


def build_fluid_domain(
    tank_body: cq.Shape,
    *,
    y_min_mm: float,
    y_max_mm: float,
    plane_tolerance_mm: float = 1e-5,
) -> cq.Solid:
    """Build the exact closure-slab-minus-tank seed-selected fluid solid."""

    return _diagnose_exact_construction(
        tank_body,
        y_min_mm=y_min_mm,
        y_max_mm=y_max_mm,
        plane_tolerance_mm=plane_tolerance_mm,
    ).fluid
