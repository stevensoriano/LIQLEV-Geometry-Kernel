from __future__ import annotations

from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from liqlev.cad import measure as measure_module
from liqlev.cad.measure import (
    GeometryMeasurementError,
    build_geometry_kernel,
    measure_geometry,
)
from liqlev.geometry.jit import eval_ppoly_derivative
from liqlev.geometry.package import validate_geometry_kernel
from liqlev.geometry.schema import GeometryMetadata


cq = pytest.importorskip("cadquery", reason="CadQuery is not installed")

MM_TO_FT = 1.0 / 304.8
MM2_TO_FT2 = MM_TO_FT**2
MM3_TO_FT3 = MM_TO_FT**3


def cylinder_expected(
    radius: float,
    height_mm: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    area = np.full_like(height_mm, np.pi * radius**2)
    volume = area * height_mm
    perimeter = np.full_like(height_mm, 2.0 * np.pi * radius)
    sidewall = 2.0 * np.pi * radius * height_mm
    return area, volume, perimeter, sidewall


def sphere_expected(
    radius: float,
    height_mm: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    area = np.pi * (2.0 * radius * height_mm - height_mm**2)
    volume = np.pi * height_mm**2 * (radius - height_mm / 3.0)
    perimeter = 2.0 * np.pi * np.sqrt(
        np.maximum(0.0, 2.0 * radius * height_mm - height_mm**2)
    )
    sidewall = 2.0 * np.pi * radius * height_mm
    return area, volume, perimeter, sidewall


def metadata(y_min_mm: float, y_max_mm: float) -> GeometryMetadata:
    return GeometryMetadata(
        schema_version=1,
        geometry_id="analytic-test",
        source_step_sha256="0" * 64,
        fluid_step_sha256="1" * 64,
        axis="+Y",
        gravity_direction="-Y",
        length_unit="ft",
        area_unit="ft^2",
        volume_unit="ft^3",
        y_min_mm=y_min_mm,
        y_max_mm=y_max_mm,
    )


def test_measures_exact_y_aligned_cylinder() -> None:
    radius = 17.0
    y_min_mm = -23.0
    y_max_mm = 59.0
    height_mm = np.linspace(0.0, y_max_mm - y_min_mm, 33)
    absolute_y_mm = y_min_mm + height_mm
    fluid = cq.Solid.makeCylinder(
        radius,
        y_max_mm - y_min_mm,
        cq.Vector(0.0, y_min_mm, 0.0),
        cq.Vector(0.0, 1.0, 0.0),
    )

    samples = measure_geometry(
        fluid,
        absolute_y_mm,
        y_min_mm=y_min_mm,
        y_max_mm=y_max_mm,
    )

    area, volume, perimeter, sidewall = cylinder_expected(radius, height_mm)
    assert samples.volume_mm3[0] == 0.0
    np.testing.assert_array_equal(samples.height_mm, height_mm)
    np.testing.assert_allclose(samples.volume_mm3[1:], volume[1:], rtol=1e-8)
    np.testing.assert_allclose(samples.section_area_mm2, area, rtol=1e-8)
    np.testing.assert_allclose(samples.perimeter_mm, perimeter, rtol=1e-8)
    np.testing.assert_allclose(samples.sidewall_area_mm2, sidewall, rtol=1e-8)

    total_wetted = area + sidewall
    total_wetted[-1] += area[-1]
    np.testing.assert_allclose(
        samples.total_wetted_area_mm2,
        total_wetted,
        rtol=1e-8,
    )
    for value in (
        samples.height_mm,
        samples.volume_mm3,
        samples.section_area_mm2,
        samples.perimeter_mm,
        samples.sidewall_area_mm2,
        samples.total_wetted_area_mm2,
    ):
        assert value.dtype == np.float64
        assert value.flags.c_contiguous
    with pytest.raises(FrozenInstanceError):
        samples.height_mm = samples.height_mm.copy()


def test_measures_exact_y_aligned_sphere() -> None:
    radius = 31.0
    y_min_mm = -11.0
    y_max_mm = y_min_mm + 2.0 * radius
    height_mm = np.linspace(0.0, 2.0 * radius, 33)
    absolute_y_mm = y_min_mm + height_mm
    fluid = cq.Solid.makeSphere(
        radius,
        cq.Vector(0.0, y_min_mm + radius, 0.0),
        cq.Vector(0.0, 1.0, 0.0),
        -90.0,
        90.0,
        360.0,
    )

    samples = measure_geometry(
        fluid,
        absolute_y_mm,
        y_min_mm=y_min_mm,
        y_max_mm=y_max_mm,
    )

    area, volume, perimeter, sidewall = sphere_expected(radius, height_mm)
    assert samples.volume_mm3[0] == 0.0
    np.testing.assert_array_equal(samples.height_mm, height_mm)
    np.testing.assert_allclose(samples.volume_mm3[1:], volume[1:], rtol=1e-8)
    np.testing.assert_allclose(
        samples.section_area_mm2[1:-1],
        area[1:-1],
        rtol=1e-8,
    )
    np.testing.assert_allclose(
        samples.perimeter_mm[1:-1],
        perimeter[1:-1],
        rtol=1e-8,
    )
    np.testing.assert_allclose(
        samples.sidewall_area_mm2[1:],
        sidewall[1:],
        rtol=1e-8,
    )
    np.testing.assert_allclose(
        samples.total_wetted_area_mm2[1:],
        sidewall[1:],
        rtol=1e-8,
    )


def test_builds_validated_sphere_geometry_kernel() -> None:
    radius = 31.0
    y_min_mm = -11.0
    y_max_mm = y_min_mm + 2.0 * radius
    fluid = cq.Solid.makeSphere(
        radius,
        cq.Vector(0.0, y_min_mm + radius, 0.0),
        cq.Vector(0.0, 1.0, 0.0),
        -90.0,
        90.0,
        360.0,
    )

    kernel = build_geometry_kernel(
        fluid,
        metadata=metadata(y_min_mm, y_max_mm),
        max_nodes=1025,
    )

    validate_geometry_kernel(kernel)
    assert 33 < len(kernel.height_ft) <= 1025
    height_mm = kernel.height_ft / MM_TO_FT
    area, volume, perimeter, sidewall = sphere_expected(radius, height_mm)
    np.testing.assert_allclose(kernel.volume_ft3[1:], volume[1:] * MM3_TO_FT3)
    np.testing.assert_allclose(
        kernel.section_area_ft2[1:-1],
        area[1:-1] * MM2_TO_FT2,
        rtol=1e-8,
    )
    np.testing.assert_allclose(
        kernel.perimeter_ft[1:-1],
        perimeter[1:-1] * MM_TO_FT,
        rtol=1e-8,
    )
    np.testing.assert_allclose(
        kernel.sidewall_area_ft2[1:],
        sidewall[1:] * MM2_TO_FT2,
        rtol=1e-8,
    )
    evaluation_height_ft = np.linspace(
        kernel.height_ft[0],
        kernel.height_ft[-1],
        1025,
    )
    evaluation_area_ft2 = np.asarray(
        [
            eval_ppoly_derivative(
                float(height_ft),
                kernel.height_ft,
                kernel.volume_coefficients,
            )
            for height_ft in evaluation_height_ft
        ]
    )
    integrated_volume_ft3 = np.trapezoid(
        evaluation_area_ft2,
        evaluation_height_ft,
    )
    assert (
        abs(integrated_volume_ft3 - kernel.volume_ft3[-1])
        / kernel.volume_ft3[-1]
        <= 5e-4
    )
    midpoint_height_ft = 0.5 * (
        kernel.height_ft[:-1] + kernel.height_ft[1:]
    )
    midpoint_height_mm = midpoint_height_ft / MM_TO_FT
    midpoint_area_mm2 = sphere_expected(radius, midpoint_height_mm)[0]
    derivative_area_mm2 = np.asarray(
        [
            eval_ppoly_derivative(
                float(height_ft),
                kernel.height_ft,
                kernel.volume_coefficients,
            )
            / MM2_TO_FT2
            for height_ft in midpoint_height_ft
        ]
    )
    topology_clearance_mm = 2.0 * radius / 32.0
    away_from_topology = (
        (midpoint_height_mm > topology_clearance_mm)
        & (
            midpoint_height_mm
            < 2.0 * radius - topology_clearance_mm
        )
    )
    np.testing.assert_array_less(
        np.abs(
            derivative_area_mm2[away_from_topology]
            - midpoint_area_mm2[away_from_topology]
        ),
        2e-3 * midpoint_area_mm2[away_from_topology],
    )


def test_dense_topology_cannot_make_area_validation_vacuous() -> None:
    topology_y_mm = np.linspace(0.0, 32.0, 17)
    half_width_mm = np.where(
        np.arange(len(topology_y_mm)) % 2 == 0,
        8.0,
        12.0,
    )
    depth_mm = 5.0
    right = list(zip(half_width_mm, topology_y_mm, strict=True))
    left = list(
        zip(
            -half_width_mm[::-1],
            topology_y_mm[::-1],
            strict=True,
        )
    )
    fluid = (
        cq.Workplane("XY")
        .polyline([*right, *left])
        .close()
        .extrude(depth_mm)
        .val()
    )

    kernel = build_geometry_kernel(
        fluid,
        metadata=metadata(0.0, 32.0),
        max_nodes=1025,
    )

    validate_geometry_kernel(kernel)
    assert len(kernel.height_ft) > 33
    midpoint_height_ft = 0.5 * (
        kernel.height_ft[:-1] + kernel.height_ft[1:]
    )
    midpoint_y_mm = midpoint_height_ft / MM_TO_FT
    coincidence_tolerance_mm = 1e-8 * 32.0
    eligible = np.asarray(
        [
            np.min(np.abs(topology_y_mm - y_mm))
            > coincidence_tolerance_mm
            for y_mm in midpoint_y_mm
        ]
    )
    assert np.any(eligible)
    direct_area_mm2 = (
        2.0
        * depth_mm
        * np.interp(
            midpoint_y_mm,
            topology_y_mm,
            half_width_mm,
        )
    )
    derivative_area_mm2 = np.asarray(
        [
            eval_ppoly_derivative(
                float(height_ft),
                kernel.height_ft,
                kernel.volume_coefficients,
            )
            / MM2_TO_FT2
            for height_ft in midpoint_height_ft
        ]
    )
    maximum_relative_error = np.max(
        np.abs(
            derivative_area_mm2[eligible]
            - direct_area_mm2[eligible]
        )
        / direct_area_mm2[eligible]
    )
    assert maximum_relative_error <= 2e-3


def test_builds_validated_cylinder_geometry_kernel() -> None:
    radius = 17.0
    y_min_mm = -23.0
    y_max_mm = 59.0
    fluid = cq.Solid.makeCylinder(
        radius,
        y_max_mm - y_min_mm,
        cq.Vector(0.0, y_min_mm, 0.0),
        cq.Vector(0.0, 1.0, 0.0),
    )
    expected_metadata = metadata(y_min_mm, y_max_mm)

    kernel = build_geometry_kernel(
        fluid,
        metadata=expected_metadata,
        max_nodes=33,
    )

    validate_geometry_kernel(kernel)
    assert kernel.metadata == expected_metadata
    assert len(kernel.height_ft) == 33
    height_mm = kernel.height_ft / MM_TO_FT
    area, volume, perimeter, sidewall = cylinder_expected(radius, height_mm)
    np.testing.assert_allclose(kernel.volume_ft3, volume * MM3_TO_FT3, rtol=1e-8)
    np.testing.assert_allclose(
        kernel.section_area_ft2,
        area * MM2_TO_FT2,
        rtol=1e-8,
    )
    np.testing.assert_allclose(
        kernel.perimeter_ft,
        perimeter * MM_TO_FT,
        rtol=1e-8,
    )
    np.testing.assert_allclose(
        kernel.sidewall_area_ft2,
        sidewall * MM2_TO_FT2,
        rtol=1e-8,
    )


def test_accepts_exact_endpoint_face_mosaics() -> None:
    width_mm = 8.0
    depth_mm = 6.0
    total_height_mm = 10.0
    left = cq.Solid.makeBox(
        width_mm / 2.0,
        total_height_mm,
        depth_mm,
        cq.Vector(-width_mm / 2.0, 0.0, -depth_mm / 2.0),
    )
    right = cq.Solid.makeBox(
        width_mm / 2.0,
        total_height_mm,
        depth_mm,
        cq.Vector(0.0, 0.0, -depth_mm / 2.0),
    )
    fluid = left.fuse(right, glue=False).Solids()[0]
    closure_counts: list[int] = []
    for y_mm, normal_y in ((0.0, -1.0), (total_height_mm, 1.0)):
        closure_counts.append(
            sum(
                face.geomType() == "PLANE"
                and abs(face.Center().y - y_mm) <= 1e-9
                and face.normalAt(face.Center()).y * normal_y > 1.0 - 1e-9
                for face in fluid.Faces()
            )
        )
    assert fluid.isValid()
    assert closure_counts == [2, 2]

    height_mm = np.linspace(0.0, total_height_mm, 33)
    samples = measure_geometry(
        fluid,
        height_mm,
        y_min_mm=0.0,
        y_max_mm=total_height_mm,
    )
    section_area_mm2 = np.full_like(
        height_mm,
        width_mm * depth_mm,
    )
    volume_mm3 = section_area_mm2 * height_mm
    perimeter_mm = np.full_like(
        height_mm,
        2.0 * (width_mm + depth_mm),
    )
    sidewall_area_mm2 = perimeter_mm * height_mm
    total_wetted_area_mm2 = section_area_mm2 + sidewall_area_mm2
    total_wetted_area_mm2[-1] += section_area_mm2[-1]
    np.testing.assert_allclose(samples.volume_mm3, volume_mm3, rtol=1e-8)
    np.testing.assert_allclose(
        samples.section_area_mm2,
        section_area_mm2,
        rtol=1e-8,
    )
    np.testing.assert_allclose(
        samples.perimeter_mm,
        perimeter_mm,
        rtol=1e-8,
    )
    np.testing.assert_allclose(
        samples.sidewall_area_mm2,
        sidewall_area_mm2,
        rtol=1e-8,
    )
    np.testing.assert_allclose(
        samples.total_wetted_area_mm2,
        total_wetted_area_mm2,
        rtol=1e-8,
    )

    kernel = build_geometry_kernel(
        fluid,
        metadata=metadata(0.0, total_height_mm),
        max_nodes=33,
    )
    validate_geometry_kernel(kernel)
    assert len(kernel.height_ft) == 33
    np.testing.assert_allclose(
        kernel.volume_ft3,
        volume_mm3 * MM3_TO_FT3,
        rtol=1e-8,
    )
    np.testing.assert_allclose(
        kernel.section_area_ft2,
        section_area_mm2 * MM2_TO_FT2,
        rtol=1e-8,
    )
    np.testing.assert_allclose(
        kernel.perimeter_ft,
        perimeter_mm * MM_TO_FT,
        rtol=1e-8,
    )
    np.testing.assert_allclose(
        kernel.sidewall_area_ft2,
        sidewall_area_mm2 * MM2_TO_FT2,
        rtol=1e-8,
    )


def test_rejects_endpoint_boundary_loops_touching_at_one_vertex() -> None:
    lower_left = cq.Solid.makeBox(
        1.0,
        1.0,
        1.0,
        cq.Vector(-1.0, 0.0, -1.0),
    )
    lower_right = cq.Solid.makeBox(
        1.0,
        1.0,
        1.0,
        cq.Vector(0.0, 0.0, 0.0),
    )
    upper_bridge = cq.Solid.makeBox(
        2.0,
        1.5,
        2.0,
        cq.Vector(-1.0, 0.5, -1.0),
    )
    fluid = lower_left.fuse(
        lower_right,
        upper_bridge,
        glue=False,
    ).Solids()[0]
    bottom_faces = [
        face
        for face in fluid.Faces()
        if face.geomType() == "PLANE"
        and abs(face.Center().y) <= 1e-9
        and face.normalAt(face.Center()).y < -1.0 + 1e-9
    ]

    assert fluid.isValid()
    assert len(bottom_faces) == 2
    assert [len(face.Edges()) for face in bottom_faces] == [4, 4]
    shared_vertices = [
        left_vertex
        for left_vertex in bottom_faces[0].Vertices()
        if any(
            left_vertex.isSame(right_vertex)
            for right_vertex in bottom_faces[1].Vertices()
        )
    ]
    assert len(shared_vertices) == 1

    with pytest.raises(
        GeometryMeasurementError,
        match="ambiguous endpoint boundary",
    ):
        measure_module._endpoint_section_metrics(
            fluid,
            0.0,
            -1.0,
            1e-8,
        )


@pytest.mark.parametrize("max_nodes", [True, 32, 1026, 33.0])
def test_rejects_invalid_adaptive_node_limits(max_nodes: object) -> None:
    fluid = cq.Solid.makeCylinder(
        2.0,
        4.0,
        cq.Vector(),
        cq.Vector(0.0, 1.0, 0.0),
    )

    with pytest.raises(GeometryMeasurementError, match="max_nodes"):
        build_geometry_kernel(
            fluid,
            metadata=metadata(0.0, 4.0),
            max_nodes=max_nodes,
        )


def test_max_node_failure_reports_each_refinement_driver() -> None:
    radius = 31.0
    y_min_mm = -11.0
    y_max_mm = y_min_mm + 2.0 * radius
    fluid = cq.Solid.makeSphere(
        radius,
        cq.Vector(0.0, y_min_mm + radius, 0.0),
        cq.Vector(0.0, 1.0, 0.0),
        -90.0,
        90.0,
        360.0,
    )

    with pytest.raises(
        GeometryMeasurementError,
        match=(
            r"current_nodes=33, proposed_nodes=\d+, "
            r"volume_failures=\d+, sidewall_failures=\d+, "
            r"direct_area_failures=\d+, "
            r"direct_area_neighbor_additions=\d+"
        ),
    ):
        build_geometry_kernel(
            fluid,
            metadata=metadata(y_min_mm, y_max_mm),
            max_nodes=33,
        )


@pytest.mark.parametrize(
    ("area_failing", "expected"),
    [
        ([True, False, False], [True, True, False]),
        ([False, False, True], [False, True, True]),
        (
            [False, False, True, False, False, True, False, False],
            [False, True, True, True, True, True, True, False],
        ),
        (
            [False, True, False, True, False],
            [True, True, True, True, True],
        ),
    ],
)
def test_direct_area_refinement_uses_one_ring_interval_dilation(
    area_failing: list[bool],
    expected: list[bool],
) -> None:
    area_failing_array = np.asarray(area_failing)

    neighborhood = measure_module._one_ring_refinement_mask(
        area_failing_array,
    )

    np.testing.assert_array_equal(
        neighborhood,
        np.asarray(expected),
    )


def test_adaptive_measurement_cache_reuses_prior_interior_ordinates(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fluid = cq.Solid.makeCylinder(
        7.0,
        10.0,
        cq.Vector(),
        cq.Vector(0.0, 1.0, 0.0),
    )
    first_nodes = np.linspace(0.0, 10.0, 5)
    second_nodes = np.sort(
        np.concatenate(
            (
                first_nodes,
                0.5 * (first_nodes[:-1] + first_nodes[1:]),
            )
        )
    )
    requests: list[np.ndarray] = []
    exact_measure_geometry = measure_module.measure_geometry

    def record_request(fluid, nodes, *, y_min_mm, y_max_mm):
        requests.append(np.asarray(nodes).copy())
        return exact_measure_geometry(
            fluid,
            nodes,
            y_min_mm=y_min_mm,
            y_max_mm=y_max_mm,
        )

    monkeypatch.setattr(
        measure_module,
        "measure_geometry",
        record_request,
    )
    cache: dict[float, tuple[float, ...]] = {}

    measure_module._measure_geometry_cached(
        fluid,
        first_nodes,
        y_min_mm=0.0,
        y_max_mm=10.0,
        cache=cache,
    )
    cached = measure_module._measure_geometry_cached(
        fluid,
        second_nodes,
        y_min_mm=0.0,
        y_max_mm=10.0,
        cache=cache,
    )
    direct = exact_measure_geometry(
        fluid,
        second_nodes,
        y_min_mm=0.0,
        y_max_mm=10.0,
    )

    assert len(requests) == 2
    assert set(requests[0][1:-1]).isdisjoint(requests[1][1:-1])
    for name in (
        "height_mm",
        "volume_mm3",
        "section_area_mm2",
        "perimeter_mm",
        "sidewall_area_mm2",
        "total_wetted_area_mm2",
    ):
        np.testing.assert_array_equal(
            getattr(cached, name),
            getattr(direct, name),
        )


def test_rejects_kernel_metadata_contract_mismatch() -> None:
    fluid = cq.Solid.makeCylinder(
        2.0,
        4.0,
        cq.Vector(),
        cq.Vector(0.0, 1.0, 0.0),
    )
    invalid_metadata = GeometryMetadata(
        **{
            **metadata(0.0, 4.0).as_dict(),
            "axis": "+Z",
        }
    )

    with pytest.raises(GeometryMeasurementError, match="metadata"):
        build_geometry_kernel(fluid, metadata=invalid_metadata)


def test_rejects_kernel_metadata_fluid_bound_mismatch() -> None:
    fluid = cq.Solid.makeCylinder(
        2.0,
        4.0,
        cq.Vector(),
        cq.Vector(0.0, 1.0, 0.0),
    )

    with pytest.raises(GeometryMeasurementError, match="bounds"):
        build_geometry_kernel(fluid, metadata=metadata(-1.0, 4.0))
