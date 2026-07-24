from __future__ import annotations

from dataclasses import dataclass

import cadquery as cq
import numpy as np

from liqlev.geometry.coefficients import pchip_coefficients
from liqlev.geometry.jit import eval_ppoly, eval_ppoly_derivative
from liqlev.geometry.package import GeometryPackageError, validate_geometry_kernel
from liqlev.geometry.schema import GeometryKernel, GeometryMetadata


MM_TO_FT = 1.0 / 304.8
MM2_TO_FT2 = MM_TO_FT**2
MM3_TO_FT3 = MM_TO_FT**3
_REFINEMENT_RELATIVE_TOLERANCE = 5e-4
_VOLUME_VALIDATION_RELATIVE_TOLERANCE = 5e-4
_AREA_VALIDATION_RELATIVE_TOLERANCE = 2e-3
_INITIAL_UNIFORM_NODE_COUNT = 33
_MAXIMUM_NODE_LIMIT = 1025


class GeometryMeasurementError(ValueError):
    """Raised when exact CAD geometry cannot satisfy the table contract."""


@dataclass(frozen=True)
class GeometrySamplesMM:
    height_mm: np.ndarray
    volume_mm3: np.ndarray
    section_area_mm2: np.ndarray
    perimeter_mm: np.ndarray
    sidewall_area_mm2: np.ndarray
    total_wetted_area_mm2: np.ndarray


def _plane_faces_at_y(
    shape: cq.Shape,
    y_mm: float,
    normal_y: float,
    tolerance_mm: float,
) -> list[cq.Face]:
    result: list[cq.Face] = []
    for face in shape.Faces():
        if face.geomType() != "PLANE":
            continue
        center = face.Center()
        if abs(center.y - y_mm) > tolerance_mm:
            continue
        normal = face.normalAt(center)
        if normal_y * normal.y <= 1.0 - 1e-9:
            continue
        if abs(normal.x) > 1e-9 or abs(normal.z) > 1e-9:
            continue
        result.append(face)
    return result


def _section_face(
    shape: cq.Shape,
    y_mm: float,
    normal_y: float,
    tolerance_mm: float,
) -> cq.Face | None:
    candidates = _plane_faces_at_y(shape, y_mm, normal_y, tolerance_mm)
    if not candidates:
        return None
    if len(candidates) != 1:
        raise GeometryMeasurementError(
            f"section at Y={y_mm:.12g} mm has {len(candidates)} outer faces; "
            "exactly one is supported"
        )
    face = candidates[0]
    if len(face.Wires()) != 1:
        raise GeometryMeasurementError(
            f"section at Y={y_mm:.12g} mm contains an inner wire"
        )
    return face


def _require_single_solid(shape: cq.Shape, y_mm: float) -> cq.Solid:
    solids = shape.Solids()
    if len(solids) != 1:
        raise GeometryMeasurementError(
            f"clip at Y={y_mm:.12g} mm produced {len(solids)} solids; "
            "exactly one connected fluid region is supported"
        )
    solid = solids[0]
    if not solid.isValid():
        raise GeometryMeasurementError(
            f"clip at Y={y_mm:.12g} mm produced an invalid solid"
        )
    return solid


def _clip_below(
    fluid: cq.Solid,
    y_mm: float,
    *,
    y_min_mm: float,
) -> cq.Solid:
    bounds = fluid.BoundingBox()
    x_length = bounds.xlen + 2.0
    z_length = bounds.zlen + 2.0
    y_low = y_min_mm - 1.0
    y_length = y_mm - y_low
    clip_box = (
        cq.Workplane("XY")
        .box(x_length, y_length, z_length)
        .translate(
            (
                0.5 * (bounds.xmin + bounds.xmax),
                y_low + 0.5 * y_length,
                0.5 * (bounds.zmin + bounds.zmax),
            )
        )
        .val()
    )
    return _require_single_solid(fluid.intersect(clip_box), y_mm)


def _as_float64(values: list[float] | np.ndarray) -> np.ndarray:
    return np.ascontiguousarray(values, dtype=np.float64)


def measure_geometry(
    fluid: cq.Solid,
    absolute_y_mm: np.ndarray,
    *,
    y_min_mm: float,
    y_max_mm: float,
) -> GeometrySamplesMM:
    """Measure exact cumulative geometry at ordered absolute assembly-Y nodes."""

    if not isinstance(fluid, cq.Solid) or not fluid.isValid():
        raise GeometryMeasurementError("fluid must be a valid CadQuery Solid")
    nodes = np.asarray(absolute_y_mm, dtype=np.float64)
    if (
        nodes.ndim != 1
        or len(nodes) < 2
        or not np.all(np.isfinite(nodes))
        or np.any(np.diff(nodes) <= 0.0)
    ):
        raise GeometryMeasurementError(
            "absolute_y_mm must be a finite, strictly increasing 1D array"
        )
    if not np.isfinite(y_min_mm) or not np.isfinite(y_max_mm):
        raise GeometryMeasurementError("Y bounds must be finite")
    if y_min_mm >= y_max_mm:
        raise GeometryMeasurementError("y_min_mm must be less than y_max_mm")

    scale_mm = max(1.0, abs(y_min_mm), abs(y_max_mm))
    tolerance_mm = 1e-8 * scale_mm
    if (
        abs(nodes[0] - y_min_mm) > tolerance_mm
        or abs(nodes[-1] - y_max_mm) > tolerance_mm
    ):
        raise GeometryMeasurementError(
            "absolute_y_mm must start at y_min_mm and end at y_max_mm"
        )

    bounds = fluid.BoundingBox()
    if (
        abs(bounds.ymin - y_min_mm) > tolerance_mm
        or abs(bounds.ymax - y_max_mm) > tolerance_mm
    ):
        raise GeometryMeasurementError(
            "fluid Y bounds do not match y_min_mm and y_max_mm"
        )

    minimum_face = _section_face(
        fluid,
        y_min_mm,
        -1.0,
        tolerance_mm,
    )
    maximum_face = _section_face(
        fluid,
        y_max_mm,
        1.0,
        tolerance_mm,
    )
    minimum_cap_area = 0.0 if minimum_face is None else minimum_face.Area()
    maximum_cap_area = 0.0 if maximum_face is None else maximum_face.Area()

    volume: list[float] = []
    section_area: list[float] = []
    perimeter: list[float] = []
    sidewall_area: list[float] = []
    total_wetted_area: list[float] = []

    for index, y_mm in enumerate(nodes):
        if index == 0:
            volume.append(0.0)
            section_area.append(minimum_cap_area)
            perimeter.append(
                0.0 if minimum_face is None else minimum_face.outerWire().Length()
            )
            sidewall_area.append(0.0)
            total_wetted_area.append(minimum_cap_area)
            continue
        if index == len(nodes) - 1:
            volume.append(fluid.Volume())
            section_area.append(maximum_cap_area)
            perimeter.append(
                0.0 if maximum_face is None else maximum_face.outerWire().Length()
            )
            total_area = fluid.Area()
            total_wetted_area.append(total_area)
            sidewall_area.append(
                total_area - minimum_cap_area - maximum_cap_area
            )
            continue

        clipped = _clip_below(fluid, float(y_mm), y_min_mm=y_min_mm)
        cut_face = _section_face(
            clipped,
            float(y_mm),
            1.0,
            tolerance_mm,
        )
        if cut_face is None:
            raise GeometryMeasurementError(
                f"clip at Y={y_mm:.12g} mm has no planar +Y cut face"
            )
        cut_area = cut_face.Area()
        wetted_area = clipped.Area() - cut_area
        volume.append(clipped.Volume())
        section_area.append(cut_area)
        perimeter.append(cut_face.outerWire().Length())
        total_wetted_area.append(wetted_area)
        sidewall_area.append(wetted_area - minimum_cap_area)

    return GeometrySamplesMM(
        height_mm=_as_float64(nodes - y_min_mm),
        volume_mm3=_as_float64(volume),
        section_area_mm2=_as_float64(section_area),
        perimeter_mm=_as_float64(perimeter),
        sidewall_area_mm2=_as_float64(sidewall_area),
        total_wetted_area_mm2=_as_float64(total_wetted_area),
    )


def _metadata_interval(
    fluid: cq.Solid,
    metadata: GeometryMetadata,
) -> tuple[float, float, float]:
    if not isinstance(metadata, GeometryMetadata):
        raise GeometryMeasurementError("metadata must be GeometryMetadata")
    contract = (
        metadata.axis,
        metadata.gravity_direction,
        metadata.length_unit,
        metadata.area_unit,
        metadata.volume_unit,
    )
    if contract != ("+Y", "-Y", "ft", "ft^2", "ft^3"):
        raise GeometryMeasurementError(
            f"metadata axis/unit contract is invalid: {contract}"
        )
    y_min_mm = float(metadata.y_min_mm)
    y_max_mm = float(metadata.y_max_mm)
    if (
        not np.isfinite(y_min_mm)
        or not np.isfinite(y_max_mm)
        or y_min_mm >= y_max_mm
    ):
        raise GeometryMeasurementError("metadata Y bounds are invalid")
    bounds = fluid.BoundingBox()
    tolerance_mm = 1e-8 * max(1.0, abs(y_min_mm), abs(y_max_mm))
    if (
        abs(bounds.ymin - y_min_mm) > tolerance_mm
        or abs(bounds.ymax - y_max_mm) > tolerance_mm
    ):
        raise GeometryMeasurementError(
            "metadata Y bounds are inconsistent with the fluid bounds"
        )
    return y_min_mm, y_max_mm, tolerance_mm


def _initial_absolute_nodes(
    fluid: cq.Solid,
    y_min_mm: float,
    y_max_mm: float,
    tolerance_mm: float,
) -> tuple[np.ndarray, np.ndarray]:
    topology_y = sorted(
        {
            float(vertex.Center().y)
            for vertex in fluid.Vertices()
            if y_min_mm + tolerance_mm
            < float(vertex.Center().y)
            < y_max_mm - tolerance_mm
        }
    )
    nodes = list(
        np.linspace(
            y_min_mm,
            y_max_mm,
            _INITIAL_UNIFORM_NODE_COUNT,
        )
    )
    for topology_node in topology_y:
        distances = np.abs(np.asarray(nodes) - topology_node)
        nearest = int(np.argmin(distances))
        if distances[nearest] <= tolerance_mm:
            nodes[nearest] = topology_node
        else:
            nodes.append(topology_node)
    return (
        _as_float64(sorted(nodes)),
        _as_float64([y_min_mm, *topology_y, y_max_mm]),
    )


def _select_samples(
    samples: GeometrySamplesMM,
    indices: np.ndarray,
) -> GeometrySamplesMM:
    return GeometrySamplesMM(
        **{
            name: np.ascontiguousarray(getattr(samples, name)[indices])
            for name in (
                "height_mm",
                "volume_mm3",
                "section_area_mm2",
                "perimeter_mm",
                "sidewall_area_mm2",
                "total_wetted_area_mm2",
            )
        }
    )


def _midpoint_area_masks(
    midpoint_height: np.ndarray,
    direct_area: np.ndarray,
    breaks: np.ndarray,
    volume_coefficients: np.ndarray,
    absolute_y_mm: np.ndarray,
    topology_y_mm: np.ndarray,
    topology_tolerance_mm: float,
    degenerate_endpoint_y_mm: np.ndarray,
    degenerate_endpoint_clearance_mm: float,
) -> tuple[np.ndarray, np.ndarray]:
    derivative_area = np.asarray(
        [
            eval_ppoly_derivative(
                float(height),
                breaks,
                volume_coefficients,
            )
            for height in midpoint_height
        ]
    )
    topology_coincidence = np.asarray(
        [
            np.min(np.abs(topology_y_mm - absolute_y))
            <= topology_tolerance_mm
            for absolute_y in absolute_y_mm
        ]
    )
    if len(degenerate_endpoint_y_mm):
        degenerate_endpoint_neighborhood = np.asarray(
            [
                np.min(
                    np.abs(degenerate_endpoint_y_mm - absolute_y)
                )
                <= degenerate_endpoint_clearance_mm
                for absolute_y in absolute_y_mm
            ]
        )
    else:
        degenerate_endpoint_neighborhood = np.zeros(
            len(absolute_y_mm),
            dtype=bool,
        )
    excluded = topology_coincidence | degenerate_endpoint_neighborhood
    positive_area = direct_area > 1e-14 * max(
        1.0,
        float(np.max(direct_area)),
    )
    eligible = ~excluded & positive_area
    failing = (
        eligible
        & (
            np.abs(derivative_area - direct_area)
            > _AREA_VALIDATION_RELATIVE_TOLERANCE * direct_area
        )
    )
    return failing, eligible


def _degenerate_endpoint_y_mm(
    absolute_nodes_mm: np.ndarray,
    section_area_mm2: np.ndarray,
) -> np.ndarray:
    area_scale = max(1.0, float(np.max(section_area_mm2)))
    endpoint_indices = np.asarray([0, len(absolute_nodes_mm) - 1])
    degenerate = (
        section_area_mm2[endpoint_indices]
        <= 1e-14 * area_scale
    )
    return np.ascontiguousarray(
        absolute_nodes_mm[endpoint_indices][degenerate],
        dtype=np.float64,
    )


def _adaptive_measurements(
    fluid: cq.Solid,
    *,
    y_min_mm: float,
    y_max_mm: float,
    max_nodes: int,
    tolerance_mm: float,
) -> tuple[GeometrySamplesMM, np.ndarray, GeometrySamplesMM]:
    nodes, topology_y = _initial_absolute_nodes(
        fluid,
        y_min_mm,
        y_max_mm,
        tolerance_mm,
    )
    if len(nodes) > max_nodes:
        raise GeometryMeasurementError(
            f"initial topology grid requires {len(nodes)} nodes, "
            f"exceeding max_nodes={max_nodes}"
        )

    while True:
        midpoints = 0.5 * (nodes[:-1] + nodes[1:])
        combined = np.empty(len(nodes) + len(midpoints), dtype=np.float64)
        combined[0::2] = nodes
        combined[1::2] = midpoints
        combined_samples = measure_geometry(
            fluid,
            combined,
            y_min_mm=y_min_mm,
            y_max_mm=y_max_mm,
        )
        node_indices = np.arange(0, len(combined), 2)
        midpoint_indices = np.arange(1, len(combined), 2)
        node_samples = _select_samples(combined_samples, node_indices)
        midpoint_samples = _select_samples(
            combined_samples,
            midpoint_indices,
        )
        volume_coefficients = pchip_coefficients(
            node_samples.height_mm,
            node_samples.volume_mm3,
        )
        sidewall_coefficients = pchip_coefficients(
            node_samples.height_mm,
            node_samples.sidewall_area_mm2,
        )
        predicted_volume = np.asarray(
            [
                eval_ppoly(
                    float(height_mm),
                    node_samples.height_mm,
                    volume_coefficients,
                )
                for height_mm in combined_samples.height_mm[midpoint_indices]
            ]
        )
        predicted_sidewall = np.asarray(
            [
                eval_ppoly(
                    float(height_mm),
                    node_samples.height_mm,
                    sidewall_coefficients,
                )
                for height_mm in combined_samples.height_mm[midpoint_indices]
            ]
        )
        volume_scale = abs(float(node_samples.volume_mm3[-1]))
        sidewall_scale = abs(float(node_samples.sidewall_area_mm2[-1]))
        if volume_scale <= 0.0 or sidewall_scale <= 0.0:
            raise GeometryMeasurementError(
                "fluid volume and total sidewall area must be positive"
            )
        failing = (
            np.abs(
                combined_samples.volume_mm3[midpoint_indices]
                - predicted_volume
            )
            > _REFINEMENT_RELATIVE_TOLERANCE * volume_scale
        ) | (
            np.abs(
                combined_samples.sidewall_area_mm2[midpoint_indices]
                - predicted_sidewall
            )
            > _REFINEMENT_RELATIVE_TOLERANCE * sidewall_scale
        )
        degenerate_endpoint_y = _degenerate_endpoint_y_mm(
            nodes,
            node_samples.section_area_mm2,
        )
        endpoint_clearance_mm = max(
            tolerance_mm,
            (y_max_mm - y_min_mm)
            / (_INITIAL_UNIFORM_NODE_COUNT - 1),
        )
        area_failing, area_eligible = _midpoint_area_masks(
            midpoint_samples.height_mm,
            midpoint_samples.section_area_mm2,
            node_samples.height_mm,
            volume_coefficients,
            midpoints,
            topology_y,
            tolerance_mm,
            degenerate_endpoint_y,
            endpoint_clearance_mm,
        )
        if not np.any(area_eligible):
            raise GeometryMeasurementError(
                "direct area validation has no eligible interval midpoints"
            )
        failing |= area_failing
        if np.any(area_failing):
            topology_refinement_neighborhood = np.asarray(
                [
                    np.min(np.abs(topology_y - midpoint))
                    <= endpoint_clearance_mm
                    for midpoint in midpoints
                ]
            )
            failing |= topology_refinement_neighborhood
        if not np.any(failing):
            return node_samples, topology_y, midpoint_samples
        if len(nodes) + int(np.count_nonzero(failing)) > max_nodes:
            raise GeometryMeasurementError(
                f"adaptive measurement did not converge within "
                f"max_nodes={max_nodes}"
            )
        nodes = np.sort(np.concatenate((nodes, midpoints[failing])))


def _validate_measured_kernel(
    fluid: cq.Solid,
    kernel: GeometryKernel,
    topology_y_mm: np.ndarray,
    tolerance_mm: float,
    midpoint_samples: GeometrySamplesMM,
) -> None:
    expected_volume_ft3 = fluid.Volume() * MM3_TO_FT3
    volume_scale = max(abs(expected_volume_ft3), np.finfo(np.float64).tiny)
    if (
        abs(kernel.volume_ft3[-1] - expected_volume_ft3) / volume_scale
        > _VOLUME_VALIDATION_RELATIVE_TOLERANCE
    ):
        raise GeometryMeasurementError(
            "final table volume disagrees with the exact CAD solid by more "
            "than 0.05%"
        )

    evaluation_height_ft = np.linspace(
        kernel.height_ft[0],
        kernel.height_ft[-1],
        1025,
    )
    derivative_area_ft2 = np.asarray(
        [
            eval_ppoly_derivative(
                float(height_ft),
                kernel.height_ft,
                kernel.volume_coefficients,
            )
            for height_ft in evaluation_height_ft
        ]
    )
    integrated_volume_ft3 = float(
        np.trapezoid(derivative_area_ft2, evaluation_height_ft)
    )
    if (
        abs(integrated_volume_ft3 - kernel.volume_ft3[-1]) / volume_scale
        > _VOLUME_VALIDATION_RELATIVE_TOLERANCE
    ):
        raise GeometryMeasurementError(
            "integrated dV/dh disagrees with final volume by more than 0.05%"
        )

    absolute_nodes_mm = (
        kernel.height_ft / MM_TO_FT
        + kernel.metadata.y_min_mm
    )
    degenerate_endpoint_y = _degenerate_endpoint_y_mm(
        absolute_nodes_mm,
        kernel.section_area_ft2 / MM2_TO_FT2,
    )
    area_failing, area_eligible = _midpoint_area_masks(
        midpoint_samples.height_mm * MM_TO_FT,
        midpoint_samples.section_area_mm2 * MM2_TO_FT2,
        kernel.height_ft,
        kernel.volume_coefficients,
        midpoint_samples.height_mm + kernel.metadata.y_min_mm,
        topology_y_mm,
        tolerance_mm,
        degenerate_endpoint_y,
        max(
            tolerance_mm,
            (
                kernel.metadata.y_max_mm
                - kernel.metadata.y_min_mm
            )
            / (_INITIAL_UNIFORM_NODE_COUNT - 1),
        ),
    )
    if not np.any(area_eligible):
        raise GeometryMeasurementError(
            "direct area validation has no eligible interval midpoints"
        )
    if np.any(area_failing):
        raise GeometryMeasurementError(
            "direct CAD area and dV/dh disagree by more than 0.2% away "
            "from topology nodes"
        )


def build_geometry_kernel(
    fluid: cq.Solid,
    *,
    metadata: GeometryMetadata,
    max_nodes: int = _MAXIMUM_NODE_LIMIT,
) -> GeometryKernel:
    """Adaptively measure a BRep fluid solid and build a validated kernel."""

    if not isinstance(fluid, cq.Solid) or not fluid.isValid():
        raise GeometryMeasurementError("fluid must be a valid CadQuery Solid")
    if (
        isinstance(max_nodes, bool)
        or not isinstance(max_nodes, int)
        or not _INITIAL_UNIFORM_NODE_COUNT
        <= max_nodes
        <= _MAXIMUM_NODE_LIMIT
    ):
        raise GeometryMeasurementError(
            f"max_nodes must be an integer from "
            f"{_INITIAL_UNIFORM_NODE_COUNT} through "
            f"{_MAXIMUM_NODE_LIMIT}"
        )
    y_min_mm, y_max_mm, tolerance_mm = _metadata_interval(fluid, metadata)
    samples, topology_y_mm, midpoint_samples = _adaptive_measurements(
        fluid,
        y_min_mm=y_min_mm,
        y_max_mm=y_max_mm,
        max_nodes=max_nodes,
        tolerance_mm=tolerance_mm,
    )
    height_ft = _as_float64(samples.height_mm * MM_TO_FT)
    volume_ft3 = _as_float64(samples.volume_mm3 * MM3_TO_FT3)
    sidewall_area_ft2 = _as_float64(samples.sidewall_area_mm2 * MM2_TO_FT2)
    kernel = GeometryKernel(
        metadata=metadata,
        height_ft=height_ft,
        volume_ft3=volume_ft3,
        volume_coefficients=pchip_coefficients(height_ft, volume_ft3),
        section_area_ft2=_as_float64(
            samples.section_area_mm2 * MM2_TO_FT2
        ),
        perimeter_ft=_as_float64(samples.perimeter_mm * MM_TO_FT),
        sidewall_area_ft2=sidewall_area_ft2,
        sidewall_coefficients=pchip_coefficients(
            height_ft,
            sidewall_area_ft2,
        ),
        total_wetted_area_ft2=_as_float64(
            samples.total_wetted_area_mm2 * MM2_TO_FT2
        ),
    )
    try:
        validate_geometry_kernel(kernel)
    except GeometryPackageError as exc:
        raise GeometryMeasurementError(
            f"measured GeometryKernel is invalid: {exc}"
        ) from exc
    _validate_measured_kernel(
        fluid,
        kernel,
        topology_y_mm,
        tolerance_mm,
        midpoint_samples,
    )
    return kernel
