from __future__ import annotations

import numpy as np

from .coefficients import pchip_coefficients
from .schema import GeometryKernel, GeometryMetadata


def _metadata(geometry_id: str, height_ft: float) -> GeometryMetadata:
    return GeometryMetadata(
        schema_version=1,
        geometry_id=geometry_id,
        source_step_sha256="0" * 64,
        fluid_step_sha256="0" * 64,
        axis="+Y",
        gravity_direction="-Y",
        length_unit="ft",
        area_unit="ft^2",
        volume_unit="ft^3",
        y_min_mm=0.0,
        y_max_mm=height_ft * 304.8,
    )


def cylinder_kernel(
    diameter_ft: float, height_ft: float, node_count: int = 1025
) -> GeometryKernel:
    h = np.linspace(0.0, height_ft, node_count)
    area = np.pi * diameter_ft**2 / 4.0
    perimeter = np.pi * diameter_ft
    volume = area * h
    wall = perimeter * h
    cap = area
    total_wetted = cap + wall
    total_wetted[-1] += cap
    return GeometryKernel(
        metadata=_metadata("analytic-cylinder", height_ft),
        height_ft=np.ascontiguousarray(h),
        volume_ft3=np.ascontiguousarray(volume),
        volume_coefficients=pchip_coefficients(h, volume),
        section_area_ft2=np.full(node_count, area),
        perimeter_ft=np.full(node_count, perimeter),
        sidewall_area_ft2=np.ascontiguousarray(wall),
        sidewall_coefficients=pchip_coefficients(h, wall),
        total_wetted_area_ft2=np.ascontiguousarray(total_wetted),
    )


def sphere_kernel(radius_ft: float, node_count: int = 1025) -> GeometryKernel:
    height_ft = 2.0 * radius_ft
    h = np.linspace(0.0, height_ft, node_count)
    radial_term = np.maximum(0.0, 2.0 * radius_ft * h - h**2)
    area = np.pi * radial_term
    volume = np.pi * h**2 * (radius_ft - h / 3.0)
    perimeter = 2.0 * np.pi * np.sqrt(radial_term)
    sidewall = 2.0 * np.pi * radius_ft * h
    return GeometryKernel(
        metadata=_metadata("analytic-sphere", height_ft),
        height_ft=np.ascontiguousarray(h),
        volume_ft3=np.ascontiguousarray(volume),
        volume_coefficients=pchip_coefficients(h, volume),
        section_area_ft2=np.ascontiguousarray(area),
        perimeter_ft=np.ascontiguousarray(perimeter),
        sidewall_area_ft2=np.ascontiguousarray(sidewall),
        sidewall_coefficients=pchip_coefficients(h, sidewall),
        total_wetted_area_ft2=np.ascontiguousarray(sidewall),
    )
