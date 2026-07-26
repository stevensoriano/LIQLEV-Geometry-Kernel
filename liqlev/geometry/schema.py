from __future__ import annotations

from dataclasses import asdict, dataclass, fields
from typing import Any

import numpy as np

# Schema v2 sentinel defaults applied when loading on-disk v1 packages in memory.
V1_MIGRATED_AUDIT_STATUS = "unaudited-v1"
SUPPORTED_GEOMETRY_SCHEMA_VERSIONS = frozenset({1, 2})


@dataclass(frozen=True)
class GeometryMetadata:
    """Geometry package metadata.

    Schema v1 fields are required. Schema v2 adds optional ``tolerances`` (table
    construction tolerances) and ``audit_status`` (CAD audit posture). On-disk v1
    packages are migrated in memory by the loader to v2-shaped metadata with
    ``tolerances=None`` and ``audit_status="unaudited-v1"``; committed package
    files are never rewritten.
    """

    schema_version: int
    geometry_id: str
    source_step_sha256: str
    fluid_step_sha256: str
    axis: str
    gravity_direction: str
    length_unit: str
    area_unit: str
    volume_unit: str
    y_min_mm: float
    y_max_mm: float
    # Schema v2 fields (defaults keep all existing v1 constructors valid).
    tolerances: dict[str, float] | None = None
    audit_status: str = V1_MIGRATED_AUDIT_STATUS

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class GeometryKernel:
    metadata: GeometryMetadata
    height_ft: np.ndarray
    volume_ft3: np.ndarray
    volume_coefficients: np.ndarray
    section_area_ft2: np.ndarray
    perimeter_ft: np.ndarray
    sidewall_area_ft2: np.ndarray
    sidewall_coefficients: np.ndarray
    total_wetted_area_ft2: np.ndarray

    @staticmethod
    def array_names() -> tuple[str, ...]:
        return tuple(
            field.name for field in fields(GeometryKernel) if field.name != "metadata"
        )

    def as_mapping(self) -> dict[str, Any]:
        return {
            "metadata": self.metadata,
            **{name: getattr(self, name) for name in self.array_names()},
        }

    @property
    def total_height_ft(self) -> float:
        return float(self.height_ft[-1])

    @property
    def total_volume_ft3(self) -> float:
        return float(self.volume_ft3[-1])
