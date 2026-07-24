from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys

import numpy as np
import pytest

from liqlev.cad.measure import build_geometry_kernel
from liqlev.geometry import load_geometry_package
from liqlev.geometry.schema import GeometryMetadata


REPOSITORY = Path(__file__).resolve().parents[2]
GEOMETRY_ROOT = Path(
    os.environ.get(
        "LIQLEV_NASA_GEOMETRY_ROOT",
        REPOSITORY / "geometry",
    )
)
STEP_PATH = (
    GEOMETRY_ROOT
    / "output"
    / "nhq01-m21a-0201_LIQLEV_FLUID_DOMAIN.step"
)
PACKAGE_PATH = (
    GEOMETRY_ROOT
    / "tables"
    / "nhq01-m21a-0201_LIQLEV_GEOMETRY.npz"
)
AUDIT_PATH = (
    GEOMETRY_ROOT
    / "audit"
    / "nhq01-m21a-0201_LIQLEV_AUDIT.json"
)
BUILD_COMMAND = REPOSITORY / "scripts" / "build_nasa_tank_geometry.py"
CHECK_COMMAND = REPOSITORY / "scripts" / "check_nasa_tank_geometry.py"
EXPECTED_SOURCE_HASH = (
    "0193EC296B5B754FFDC7B1652052F5B28BA99345A5BC89922D026DF011C837D5"
)
EXPECTED_SECTION_IMAGES = {
    "sections/orthogonal_xy.png",
    "sections/orthogonal_yz.png",
    "sections/orthogonal_xz.png",
    "sections/section_xz_h000.png",
    "sections/section_xz_h010.png",
    "sections/section_xz_h025.png",
    "sections/section_xz_h050.png",
    "sections/section_xz_h075.png",
    "sections/section_xz_h090.png",
    "sections/section_xz_h100.png",
}
EXPECTED_FLUID_BOUNDS_MM = (
    -279.67000017112,
    279.67000017112,
    -275.4067913371,
    275.29679133713,
    297.81275947409995,
    857.15275961635,
)
ALLOW_PENDING_VISUAL_REVIEW = (
    os.environ.get("LIQLEV_ALLOW_PENDING_VISUAL_REVIEW") == "1"
)

sys.path.insert(0, str(REPOSITORY / "scripts"))
import build_nasa_tank_geometry as build_command  # noqa: E402
import check_nasa_tank_geometry as check_command  # noqa: E402


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file_obj:
        for block in iter(lambda: file_obj.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def test_committed_nasa_geometry_artifact_set_is_complete_and_valid() -> None:
    required = (
        BUILD_COMMAND,
        CHECK_COMMAND,
        STEP_PATH,
        PACKAGE_PATH,
        PACKAGE_PATH.with_suffix(".json"),
        PACKAGE_PATH.with_suffix(".csv"),
        AUDIT_PATH,
        *(
            AUDIT_PATH.parent / relative_path
            for relative_path in sorted(EXPECTED_SECTION_IMAGES)
        ),
    )
    missing = [str(path.relative_to(REPOSITORY)) for path in required if not path.is_file()]
    assert not missing, f"missing NASA tank artifact files: {missing}"

    audit = json.loads(AUDIT_PATH.read_text(encoding="utf-8"))
    kernel = load_geometry_package(PACKAGE_PATH)
    cq = pytest.importorskip("cadquery", reason="CadQuery is not installed")
    imported = cq.importers.importStep(str(STEP_PATH))
    solids = imported.solids().vals()

    if ALLOW_PENDING_VISUAL_REVIEW:
        assert audit["passed"] is False
        assert audit["technical_validation_passed"] is True
    else:
        assert audit["passed"] is True
    assert audit["source_sha256"] == EXPECTED_SOURCE_HASH
    assert audit["solid_count"] == 1
    assert audit["shell_count"] == 1
    assert audit["valid"] is True
    assert audit["cap_plane_max_error_mm"] <= 1e-5
    assert audit["round_trip_relative_volume_error"] <= 1e-8
    assert audit["visual_section_review"] == (
        "pending" if ALLOW_PENDING_VISUAL_REVIEW else "passed"
    )
    assert audit["node_count"] == len(kernel.height_ft)
    assert audit["bounding_box_mm"] == pytest.approx(
        EXPECTED_FLUID_BOUNDS_MM,
        abs=1e-5,
    )
    assert audit["endpoint_mosaic_face_counts"] == [27, 3]
    assert audit["section_render_modes"]["sections/section_xz_h000.png"] == (
        "endpoint_closure_mosaic"
    )
    assert audit["section_render_modes"]["sections/section_xz_h100.png"] == (
        "endpoint_closure_mosaic"
    )
    assert audit["independent_measurement"]["cad_endpoint_volume_relative"] <= 5e-4
    assert audit["independent_measurement"]["integrated_dv_dh_relative"] <= 5e-4
    assert audit["independent_measurement"]["max_midpoint_volume_relative"] <= 5e-4
    assert audit["independent_measurement"]["max_midpoint_sidewall_relative"] <= 5e-4
    assert audit["independent_measurement"]["max_direct_area_relative"] <= 2e-3
    assert set(audit["tool_versions"]) >= {
        "python",
        "cadquery",
        "ocp",
        "numpy",
        "matplotlib",
    }
    assert audit["build_command"]["argv"]
    assert audit["build_command"]["executable"]
    assert audit["artifact_sha256"]["fluid_step"] == _sha256(STEP_PATH)
    assert audit["artifact_sha256"]["geometry_npz"] == _sha256(PACKAGE_PATH)
    assert audit["artifact_sha256"]["geometry_json"] == _sha256(
        PACKAGE_PATH.with_suffix(".json")
    )
    assert audit["artifact_sha256"]["geometry_csv"] == _sha256(
        PACKAGE_PATH.with_suffix(".csv")
    )
    assert len(solids) == audit["solid_count"]
    assert solids[0].isValid()
    assert kernel.metadata.axis == "+Y"
    assert kernel.metadata.gravity_direction == "-Y"
    assert kernel.metadata.fluid_step_sha256 == _sha256(STEP_PATH)
    assert kernel.total_volume_ft3 > 0.0
    assert kernel.total_height_ft == pytest.approx(
        (275.296791 - (-275.406791)) / 304.8,
        abs=1e-10,
    )

    listed_images = set(audit["section_images"])
    assert listed_images == EXPECTED_SECTION_IMAGES
    for relative_path in listed_images:
        image_path = AUDIT_PATH.parent / relative_path
        assert image_path.is_file()
        assert image_path.stat().st_size > 0
        assert audit["artifact_sha256"][relative_path] == _sha256(image_path)

    with PACKAGE_PATH.with_suffix(".csv").open(
        newline="",
        encoding="utf-8",
    ) as file_obj:
        rows = list(csv.reader(file_obj))
    assert len(rows) == audit["node_count"] + 1
    csv_values = np.asarray(rows[1:], dtype=np.float64)
    np.testing.assert_array_equal(csv_values[:, 0], kernel.height_ft)
    np.testing.assert_array_equal(csv_values[:, 1], kernel.volume_ft3)
    np.testing.assert_array_equal(csv_values[:, 2], kernel.section_area_ft2)
    np.testing.assert_array_equal(csv_values[:, 3], kernel.perimeter_ft)
    np.testing.assert_array_equal(csv_values[:, 4], kernel.sidewall_area_ft2)
    np.testing.assert_array_equal(csv_values[:, 5], kernel.total_wetted_area_ft2)


def test_independent_cad_table_metrics_remeasure_exact_sections() -> None:
    cq = pytest.importorskip("cadquery", reason="CadQuery is not installed")
    radius_mm = 7.0
    height_mm = 20.0
    fluid = cq.Solid.makeCylinder(
        radius_mm,
        height_mm,
        cq.Vector(),
        cq.Vector(0.0, 1.0, 0.0),
    )
    metadata = GeometryMetadata(
        schema_version=1,
        geometry_id="diagnostic-cylinder",
        source_step_sha256="0" * 64,
        fluid_step_sha256="1" * 64,
        axis="+Y",
        gravity_direction="-Y",
        length_unit="ft",
        area_unit="ft^2",
        volume_unit="ft^3",
        y_min_mm=0.0,
        y_max_mm=height_mm,
    )
    kernel = build_geometry_kernel(fluid, metadata=metadata, max_nodes=33)

    metrics = check_command.independent_cad_table_metrics(fluid, kernel)

    assert metrics["cad_endpoint_volume_relative"] <= 5e-4
    assert metrics["integrated_dv_dh_relative"] <= 5e-4
    assert metrics["max_midpoint_volume_relative"] <= 5e-4
    assert metrics["max_midpoint_sidewall_relative"] <= 5e-4
    assert metrics["max_direct_area_relative"] <= 2e-3
    assert metrics["eligible_direct_area_midpoints"] == 32


@pytest.mark.parametrize(
    ("field", "tampered_value"),
    [
        ("face_count", 39),
        ("surface_area_mm2", 1.0),
        ("cap_plane_max_error_mm", 0.0),
        ("total_volume_ft3", 1.0),
        ("total_height_ft", 1.0),
        ("axis", "-Y"),
        ("gravity_direction", "+Y"),
        ("refinement_relative_tolerance", 1.0),
        ("direct_area_relative_tolerance", 1.0),
    ],
)
def test_checker_rejects_tampered_authoritative_audit_fields(
    tmp_path: Path,
    field: str,
    tampered_value: object,
) -> None:
    source_step = check_command.DEFAULT_SOURCE_STEP
    if not source_step.is_file():
        source_step = tmp_path / "source.step"
        source_step.write_bytes(b"portable checker regression source")
        expected_source_hash = _sha256(source_step)
    else:
        expected_source_hash = EXPECTED_SOURCE_HASH
    geometry_root = tmp_path / "geometry"
    shutil.copytree(GEOMETRY_ROOT, geometry_root)
    audit_path = geometry_root / check_command.AUDIT_RELATIVE
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    audit[field] = tampered_value
    audit_path.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    checks = check_command.verify_geometry_root(
        source_step=source_step,
        geometry_root=geometry_root,
        expected_source_sha256=expected_source_hash,
        precomputed_measurement=audit["independent_measurement"],
    )

    authoritative = next(
        (
            check
            for check in checks
            if check.name == "authoritative audit fields"
        ),
        None,
    )
    assert authoritative is not None
    assert not authoritative.passed, field


def test_ordered_wire_render_points_follow_topology() -> None:
    cq = pytest.importorskip("cadquery", reason="CadQuery is not installed")
    wire = cq.Workplane("XZ").rect(8.0, 6.0).wire().val()

    points = build_command._ordered_wire_points(wire, samples_per_edge=8)

    assert points.shape[1] == 2
    np.testing.assert_allclose(points[0], points[-1], atol=1e-10)
    step_lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
    assert np.max(step_lengths) <= 2.0


def test_endpoint_renderer_uses_exact_closure_face_mosaic(
    tmp_path: Path,
) -> None:
    cq = pytest.importorskip("cadquery", reason="CadQuery is not installed")
    left = cq.Solid.makeBox(
        4.0,
        10.0,
        6.0,
        cq.Vector(-4.0, 0.0, -3.0),
    )
    right = cq.Solid.makeBox(
        4.0,
        10.0,
        6.0,
        cq.Vector(0.0, 0.0, -3.0),
    )
    fluid = left.fuse(right, glue=False).Solids()[0]

    faces = build_command._endpoint_closure_faces(
        fluid,
        y_mm=0.0,
        normal_y=-1.0,
    )
    output = tmp_path / "endpoint.png"
    mode = build_command._render_xz_section(
        fluid,
        output,
        normalized_height=0.0,
        y_min_mm=0.0,
        y_max_mm=10.0,
    )

    assert len(faces) == 2
    assert mode == "endpoint_closure_mosaic"
    assert output.is_file()
    assert output.stat().st_size > 0


def test_promotion_rolls_back_on_keyboard_interrupt(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    staging = tmp_path / "staging"
    output = tmp_path / "geometry"
    managed = [Path("output/fluid.step"), Path("tables/geometry.npz")]
    for relative in managed:
        (staging / relative).parent.mkdir(parents=True, exist_ok=True)
        (output / relative).parent.mkdir(parents=True, exist_ok=True)
        (staging / relative).write_text("new", encoding="utf-8")
        (output / relative).write_text("old", encoding="utf-8")
    real_replace = build_command.os.replace
    calls = 0

    def interrupt_third_replace(source, target) -> None:
        nonlocal calls
        calls += 1
        if calls == 3:
            raise KeyboardInterrupt
        real_replace(source, target)

    monkeypatch.setattr(build_command.os, "replace", interrupt_third_replace)

    with pytest.raises(KeyboardInterrupt):
        build_command._promote_staged_outputs(staging, output, managed)

    for relative in managed:
        assert (output / relative).read_text(encoding="utf-8") == "old"


def test_promotion_preflights_all_targets_before_first_move(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    staging = tmp_path / "staging"
    output = tmp_path / "geometry"
    regular = Path("output/fluid.step")
    non_file = Path("tables/geometry.npz")
    managed = [regular, non_file]
    for relative in managed:
        (staging / relative).parent.mkdir(parents=True, exist_ok=True)
        (staging / relative).write_text("new", encoding="utf-8")
    (output / regular).parent.mkdir(parents=True, exist_ok=True)
    (output / regular).write_text("old", encoding="utf-8")
    (output / non_file).mkdir(parents=True)
    directory_payload = output / non_file / "must-survive.txt"
    directory_payload.write_text("directory content", encoding="utf-8")
    real_replace = build_command.os.replace
    replace_calls: list[tuple[object, object]] = []

    def record_replace(source, target) -> None:
        replace_calls.append((source, target))
        real_replace(source, target)

    monkeypatch.setattr(build_command.os, "replace", record_replace)

    with pytest.raises(RuntimeError, match="regular file"):
        build_command._promote_staged_outputs(staging, output, managed)

    assert replace_calls == []
    assert not list(tmp_path.glob(".nasa-tank-backup-*"))
    assert (output / regular).read_text(encoding="utf-8") == "old"
    assert (staging / regular).read_text(encoding="utf-8") == "new"
    assert (output / non_file).is_dir()
    assert directory_payload.read_text(encoding="utf-8") == (
        "directory content"
    )
    assert (staging / non_file).read_text(encoding="utf-8") == "new"


@pytest.mark.parametrize("escape_kind", ["parent", "absolute"])
def test_promotion_rejects_unsafe_managed_paths_before_first_move(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    escape_kind: str,
) -> None:
    staging = tmp_path / "staging-area" / "staging"
    output = tmp_path / "output-area" / "geometry"
    if escape_kind == "parent":
        managed = Path("../escaped.step")
    else:
        managed = (tmp_path / "absolute-escaped.step").resolve()
    staged = staging / managed
    target = output / managed
    staged.parent.mkdir(parents=True, exist_ok=True)
    staged.write_text("new", encoding="utf-8")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("old", encoding="utf-8")
    real_replace = build_command.os.replace
    replace_calls: list[tuple[object, object]] = []

    def record_replace(source, destination) -> None:
        replace_calls.append((source, destination))
        real_replace(source, destination)

    monkeypatch.setattr(build_command.os, "replace", record_replace)

    with pytest.raises(RuntimeError, match="safe relative path"):
        build_command._promote_staged_outputs(
            staging,
            output,
            [managed],
        )

    assert replace_calls == []
    assert not list(output.parent.glob(".nasa-tank-backup-*"))
    assert staged.read_text(encoding="utf-8") == (
        "old" if staged == target else "new"
    )
    assert target.read_text(encoding="utf-8") == "old"


def test_promotion_rejects_linked_destination_component_before_first_move(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    staging = tmp_path / "staging"
    output = tmp_path / "geometry"
    external = tmp_path / "external"
    output.mkdir()
    external.mkdir()
    linked_parent = output / "linked"
    if sys.platform == "win32":
        import _winapi

        _winapi.CreateJunction(str(external), str(linked_parent))
    else:
        linked_parent.symlink_to(external, target_is_directory=True)
    managed = Path("linked/fluid.step")
    staged = staging / managed
    target = external / "fluid.step"
    staged.parent.mkdir(parents=True, exist_ok=True)
    staged.write_text("new", encoding="utf-8")
    target.write_text("old", encoding="utf-8")
    real_replace = build_command.os.replace
    replace_calls: list[tuple[object, object]] = []

    def record_replace(source, destination) -> None:
        replace_calls.append((source, destination))
        real_replace(source, destination)

    monkeypatch.setattr(build_command.os, "replace", record_replace)

    with pytest.raises(RuntimeError, match="linked destination component"):
        build_command._promote_staged_outputs(
            staging,
            output,
            [managed],
        )

    assert replace_calls == []
    assert not list(tmp_path.glob(".nasa-tank-backup-*"))
    assert staged.read_text(encoding="utf-8") == "new"
    assert target.read_text(encoding="utf-8") == "old"


def test_promotion_retains_backup_when_promotion_and_rollback_fail(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    staging = tmp_path / "staging"
    output = tmp_path / "geometry"
    first = Path("output/fluid.step")
    second = Path("tables/geometry.npz")
    managed = [first, second]
    for relative in managed:
        (staging / relative).parent.mkdir(parents=True, exist_ok=True)
        (output / relative).parent.mkdir(parents=True, exist_ok=True)
        (staging / relative).write_text("new", encoding="utf-8")
        (output / relative).write_text("old", encoding="utf-8")
    real_replace = build_command.os.replace
    calls = 0

    def fail_promotion_and_rollback(source, target) -> None:
        nonlocal calls
        calls += 1
        if calls == 3:
            raise OSError("promotion failed")
        if calls == 4:
            raise OSError("rollback failed")
        real_replace(source, target)

    monkeypatch.setattr(
        build_command.os,
        "replace",
        fail_promotion_and_rollback,
    )

    with pytest.raises(OSError, match="rollback failed"):
        build_command._promote_staged_outputs(staging, output, managed)

    backups = list(tmp_path.glob(".nasa-tank-backup-*"))
    assert len(backups) == 1
    backup = backups[0]
    assert (backup / "promotion-manifest.json").is_file()
    assert (backup / first).read_text(encoding="utf-8") == "old"
    assert (output / first).read_text(encoding="utf-8") == "new"
    assert (output / second).read_text(encoding="utf-8") == "old"

    monkeypatch.setattr(build_command.os, "replace", real_replace)
    assert build_command._recover_interrupted_promotions(output) == 1
    assert (output / first).read_text(encoding="utf-8") == "old"
    assert (output / second).read_text(encoding="utf-8") == "old"
    assert not backup.exists()


def test_next_build_recovers_hard_interrupted_promotion(
    tmp_path: Path,
) -> None:
    output = tmp_path / "geometry"
    backup = tmp_path / ".nasa-tank-backup-interrupted"
    existing = Path("output/fluid.step")
    new_only = Path("tables/geometry.npz")
    (backup / existing).parent.mkdir(parents=True, exist_ok=True)
    (backup / existing).write_text("old", encoding="utf-8")
    (output / existing).parent.mkdir(parents=True, exist_ok=True)
    (output / existing).write_text("new", encoding="utf-8")
    (output / new_only).parent.mkdir(parents=True, exist_ok=True)
    (output / new_only).write_text("new", encoding="utf-8")
    (backup / "promotion-manifest.json").write_text(
        json.dumps(
            {
                "output_root": str(output.resolve()),
                "managed_paths": [existing.as_posix(), new_only.as_posix()],
                "originally_existing": [existing.as_posix()],
            }
        ),
        encoding="utf-8",
    )

    recovered = build_command._recover_interrupted_promotions(output)

    assert recovered == 1
    assert (output / existing).read_text(encoding="utf-8") == "old"
    assert not (output / new_only).exists()
    assert not backup.exists()


def test_exact_brep_metric_snapshot_precedes_visual_tessellation() -> None:
    cq = pytest.importorskip("cadquery", reason="CadQuery is not installed")
    if not STEP_PATH.is_file():
        pytest.skip("provisional exact fluid STEP is not present")
    fluid = cq.importers.importStep(str(STEP_PATH)).solids().val()

    snapshot = build_command._snapshot_brep_metrics(fluid)
    fluid.tessellate(0.8, 0.08)
    tessellated_bounds = build_command._bounds(fluid)

    assert snapshot["bounding_box_mm"] == pytest.approx(
        EXPECTED_FLUID_BOUNDS_MM,
        abs=1e-5,
    )
    assert max(
        abs(before - after)
        for before, after in zip(
            snapshot["bounding_box_mm"],
            tessellated_bounds,
            strict=True,
        )
    ) > 0.5


def test_step_header_timestamp_is_normalized_deterministically(
    tmp_path: Path,
) -> None:
    step_path = tmp_path / "fluid.step"
    step_path.write_bytes(
        b"ISO-10303-21;\nHEADER;\n"
        b"FILE_NAME('Open CASCADE Shape Model','2026-07-24T00:39:52',"
        b"('Author'),('Open CASCADE'),'processor','system','Unknown');\n"
        b"ENDSEC;\nDATA;\nENDSEC;\nEND-ISO-10303-21;\n"
    )

    build_command._normalize_step_header_timestamp(step_path)
    first = step_path.read_bytes()
    build_command._normalize_step_header_timestamp(step_path)

    assert b"'1970-01-01T00:00:00'" in first
    assert step_path.read_bytes() == first


def test_generated_text_artifacts_preserve_exact_bytes_on_checkout() -> None:
    attributes = {
        line.strip()
        for line in (REPOSITORY / ".gitattributes").read_text(
            encoding="utf-8"
        ).splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }

    assert {
        "geometry/output/*.step -text",
        "geometry/tables/*.csv -text",
        "geometry/tables/*.json -text",
        "geometry/audit/*.json -text",
    } <= attributes
