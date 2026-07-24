from __future__ import annotations

from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from liqlev.cad.measure import (
    GeometryMeasurementError,
    build_geometry_kernel,
    measure_geometry,
)
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
