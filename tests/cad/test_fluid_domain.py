from __future__ import annotations

from dataclasses import asdict
import json
import math
from pathlib import Path

import pytest


SOURCE = Path(
    r"C:\Users\sasorian\Documents\Eta_Space\geometry"
    r"\nhq01-m21a- 0201_TankAssy_NASA.STEP"
)
EXPECTED_HASH = "0193EC296B5B754FFDC7B1652052F5B28BA99345A5BC89922D026DF011C837D5"
EXPECTED_SIZE = 36_844_537
PRODUCT_NAME = "nhq01-m21a- 0202_short"
Y_MIN_MM = -275.406791
Y_MAX_MM = 275.296791
PLANE_TOLERANCE_MM = 1e-5
EXPECTED_VOLUME_MM3 = 98_109_377.7478651
EXPECTED_AREA_MM2 = 1_031_668.462706133
EXPECTED_BOUNDS_MM = (
    -279.6700000711255,
    279.6700000711255,
    -275.4067913371334,
    275.2967913371334,
    297.8127594740996,
    857.1527596163507,
)
REPOSITORY = Path(__file__).resolve().parents[2]
FINAL_STEP = (
    REPOSITORY
    / "geometry"
    / "output"
    / "nhq01-m21a-0201_LIQLEV_FLUID_DOMAIN.step"
)
FINAL_AUDIT = (
    REPOSITORY
    / "geometry"
    / "audit"
    / "nhq01-m21a-0201_LIQLEV_AUDIT.json"
)


@pytest.fixture(scope="module")
def cad_modules():
    cq = pytest.importorskip("cadquery", reason="CadQuery is not installed")
    pytest.importorskip("OCP", reason="OpenCascade Python bindings are not installed")
    import liqlev.cad.audit as audit
    import liqlev.cad.fluid_domain as fluid_domain
    import liqlev.cad.xcaf as xcaf

    return cq, audit, fluid_domain, xcaf


@pytest.fixture(scope="module")
def source_snapshot(cad_modules):
    if not SOURCE.is_file():
        pytest.skip("NASA source STEP is not present")
    _, _, _, xcaf = cad_modules
    stat = SOURCE.stat()
    return stat.st_size, stat.st_mtime_ns, xcaf.sha256_file(SOURCE)


@pytest.fixture(scope="module")
def tank_body(cad_modules, source_snapshot):
    _, _, _, xcaf = cad_modules
    return xcaf.load_named_product(
        SOURCE,
        product_name=PRODUCT_NAME,
        expected_sha256=EXPECTED_HASH,
    )


@pytest.fixture(scope="module")
def fluid(cad_modules, tank_body):
    _, _, fluid_domain, _ = cad_modules
    return fluid_domain.build_fluid_domain(
        tank_body,
        y_min_mm=Y_MIN_MM,
        y_max_mm=Y_MAX_MM,
        plane_tolerance_mm=PLANE_TOLERANCE_MM,
    )


@pytest.fixture(scope="module")
def construction_diagnostics(cad_modules, tank_body):
    _, _, fluid_domain, _ = cad_modules
    return fluid_domain._diagnose_exact_construction(
        tank_body,
        y_min_mm=Y_MIN_MM,
        y_max_mm=Y_MAX_MM,
        plane_tolerance_mm=PLANE_TOLERANCE_MM,
    )


@pytest.fixture(scope="module")
def round_trip(cad_modules, fluid, tmp_path_factory):
    _, audit_module, _, _ = cad_modules
    output_dir = tmp_path_factory.mktemp("fluid-domain-round-trip")
    output_path = output_dir / "fluid-domain.step"
    audit_path = output_dir / "fluid-domain.audit.json"
    imported, audit = audit_module.write_step_round_trip(
        fluid,
        output_path,
        audit_path,
    )
    return imported, audit, output_path, audit_path


@pytest.fixture(scope="module")
def final_round_trip(cad_modules, fluid):
    _, audit_module, _, _ = cad_modules
    imported, audit = audit_module.write_step_round_trip(
        fluid,
        FINAL_STEP,
        FINAL_AUDIT,
    )
    return imported, audit


def _bounds(shape) -> tuple[float, float, float, float, float, float]:
    box = shape.BoundingBox()
    return box.xmin, box.xmax, box.ymin, box.ymax, box.zmin, box.zmax


def test_builds_exact_reference_fluid_solid(cad_modules, tank_body, fluid) -> None:
    cq, _, _, _ = cad_modules

    assert isinstance(fluid, cq.Solid)
    assert fluid.isValid()
    assert len(fluid.Solids()) == 1
    assert len(fluid.Shells()) == 1
    assert fluid.Shells()[0].wrapped.Closed()
    assert len(fluid.Faces()) == 40
    assert fluid.Volume() == pytest.approx(EXPECTED_VOLUME_MM3, abs=1e-3)
    assert fluid.Area() == pytest.approx(EXPECTED_AREA_MM2, abs=1e-5)
    assert _bounds(fluid) == pytest.approx(EXPECTED_BOUNDS_MM, abs=1e-5)

    source_bounds = _bounds(tank_body)
    fluid_bounds = _bounds(fluid)
    for low_index, high_index in ((0, 1), (2, 3), (4, 5)):
        assert source_bounds[low_index] - PLANE_TOLERANCE_MM <= fluid_bounds[
            low_index
        ]
        assert fluid_bounds[high_index] <= source_bounds[
            high_index
        ] + PLANE_TOLERANCE_MM


def test_direct_cut_has_two_components_and_unique_strict_seed(
    construction_diagnostics,
) -> None:
    diagnostics = construction_diagnostics
    assert diagnostics.primary.result_solid_count == 2
    assert diagnostics.primary.result_shell_count == 2
    assert diagnostics.primary.result_face_count == 56
    assert diagnostics.primary.result_valid is True
    assert diagnostics.primary.seed_classifications.count("IN") == 1
    assert "ON" not in diagnostics.primary.seed_classifications
    assert diagnostics.primary.selected_solid_count == 1
    assert diagnostics.primary.selected_shell_count == 1
    assert diagnostics.primary.selected_face_count == 40
    assert diagnostics.primary.selected_valid is True
    assert diagnostics.seed_mm == pytest.approx(
        (0.0, -0.055, 577.4827595452252),
        abs=PLANE_TOLERANCE_MM,
    )


def test_exact_padding_and_splitter_oracles_agree(
    construction_diagnostics,
) -> None:
    diagnostics = construction_diagnostics
    primary = diagnostics.primary
    secondary = diagnostics.secondary
    splitter = diagnostics.splitter

    assert primary.padding_mm == 25.0
    assert secondary.padding_mm == 100.0
    assert secondary.selected_face_count == primary.selected_face_count == 40
    assert secondary.selected_solid_count == primary.selected_solid_count == 1
    assert secondary.selected_shell_count == primary.selected_shell_count == 1
    assert secondary.volume_difference_mm3 <= max(
        1e-3,
        primary.volume_mm3 * 1e-8,
    )
    assert max(secondary.bounding_box_differences_mm) <= 1e-5

    assert splitter.result_solid_count == 3
    assert splitter.seed_classifications.count("IN") == 1
    assert splitter.volume_difference_mm3 <= max(
        1e-3,
        primary.volume_mm3 * 1e-8,
    )
    assert max(splitter.bounding_box_differences_mm) <= 1e-5


def test_accepts_closure_mosaic_and_rejects_carrier_artifacts(
    construction_diagnostics,
) -> None:
    diagnostics = construction_diagnostics
    assert diagnostics.minimum_closure_face_count == 27
    assert diagnostics.maximum_closure_face_count == 3
    assert diagnostics.maximum_closure_plane_error_mm <= PLANE_TOLERANCE_MM
    assert diagnostics.primary.carrier_xz_contact is False
    assert diagnostics.primary.exterior_probe_classifications == (
        "OUT",
        "OUT",
        "OUT",
        "OUT",
    )
    assert diagnostics.primary.shared_face_count == 0
    assert diagnostics.primary.shared_edge_count == 0
    assert diagnostics.primary.shared_vertex_count == 0


def test_rims_select_exact_concentric_wet_loops(construction_diagnostics) -> None:
    minimum, maximum = construction_diagnostics.rims
    assert minimum.candidate_face_count == maximum.candidate_face_count == 1
    assert minimum.wire_count == 26
    assert maximum.wire_count == 2
    assert minimum.wet_loop_area_mm2 == pytest.approx(
        11432.587497137267,
        abs=1e-8,
    )
    assert maximum.wet_loop_area_mm2 == pytest.approx(
        11432.587497137267,
        abs=1e-8,
    )
    assert minimum.wet_loop_center_mm[0::2] == pytest.approx(
        maximum.wet_loop_center_mm[0::2],
        abs=PLANE_TOLERANCE_MM,
    )


def test_round_trip_preserves_exact_solid_geometry(
    cad_modules,
    fluid,
    round_trip,
) -> None:
    cq, _, _, _ = cad_modules
    imported, audit, _, _ = round_trip

    assert isinstance(imported, cq.Solid)
    assert imported.isValid()
    assert imported.Volume() > 0.0
    assert len(imported.Solids()) == 1
    assert len(imported.Shells()) == 1
    assert imported.Shells()[0].wrapped.Closed()
    assert len(imported.Faces()) == 40
    assert _bounds(imported) == pytest.approx(EXPECTED_BOUNDS_MM, abs=1e-5)
    assert audit.passed is True
    assert audit.absolute_volume_difference_mm3 <= max(
        1e-3,
        fluid.Volume() * 1e-8,
    )
    assert max(audit.bounding_box_differences_mm) <= 1e-5


def test_round_trip_writes_complete_deterministic_audit(round_trip) -> None:
    imported, audit, output_path, audit_path = round_trip

    assert output_path.is_file()
    assert audit_path.is_file()
    assert audit_path.read_bytes().endswith(b"\n")
    payload = json.loads(audit_path.read_text(encoding="utf-8"))
    assert payload == json.loads(
        json.dumps(asdict(audit), sort_keys=True, allow_nan=False)
    )
    assert payload["schema"] == "liqlev.cad.audit"
    assert payload["version"] == 2
    assert payload["passed"] is True
    assert payload["source_sha256"] == EXPECTED_HASH
    assert payload["source_size_bytes"] == EXPECTED_SIZE
    assert len(payload["input_sha256"]) == 64
    assert len(payload["output_sha256"]) == 64
    assert payload["input_sha256"] == payload["input_sha256"].upper()
    assert payload["output_sha256"] == payload["output_sha256"].upper()
    assert payload["pre_solid_count"] == payload["post_solid_count"] == 1
    assert payload["pre_shell_count"] == payload["post_shell_count"] == 1
    assert payload["pre_face_count"] == payload["post_face_count"] == 40
    assert payload["pre_valid"] is payload["post_valid"] is True
    assert payload["closure_planes_mm"] == [Y_MIN_MM, Y_MAX_MM]
    assert payload["cap_plane_error_mm"] <= PLANE_TOLERANCE_MM
    assert payload["post_face_count"] == len(imported.Faces())


def test_final_inspection_step_is_independently_reimported(
    cad_modules,
    final_round_trip,
) -> None:
    cq, _, _, _ = cad_modules
    imported, audit = final_round_trip
    independent = cq.importers.importStep(str(FINAL_STEP))
    independent_solids = independent.solids().vals()

    assert FINAL_STEP.is_file()
    assert FINAL_AUDIT.is_file()
    assert len(independent_solids) == 1
    inspected = independent_solids[0]
    assert inspected.isValid()
    assert len(inspected.Shells()) == 1
    assert len(inspected.Faces()) == 40
    assert inspected.Volume() == pytest.approx(imported.Volume(), abs=1e-3)
    assert _bounds(inspected) == pytest.approx(_bounds(imported), abs=1e-5)
    assert audit.passed is True


@pytest.mark.parametrize(
    ("y_min_mm", "y_max_mm", "tolerance", "message"),
    [
        (1.0, 1.0, PLANE_TOLERANCE_MM, "y_min_mm"),
        (2.0, 1.0, PLANE_TOLERANCE_MM, "y_min_mm"),
        (0.0, 1.0, 0.0, "plane_tolerance_mm"),
        (0.0, 1.0, -1e-5, "plane_tolerance_mm"),
        (0.0, 1.0, math.inf, "plane_tolerance_mm"),
    ],
)
def test_rejects_invalid_plane_arguments(
    cad_modules,
    y_min_mm: float,
    y_max_mm: float,
    tolerance: float,
    message: str,
) -> None:
    cq, _, fluid_domain, _ = cad_modules
    tank = cq.Workplane("XY").box(2.0, 2.0, 2.0).val()

    with pytest.raises(fluid_domain.FluidDomainError, match=message):
        fluid_domain.build_fluid_domain(
            tank,
            y_min_mm=y_min_mm,
            y_max_mm=y_max_mm,
            plane_tolerance_mm=tolerance,
        )


@pytest.mark.parametrize("invalid", [None, object()])
def test_rejects_invalid_non_shape_input(cad_modules, invalid) -> None:
    _, _, fluid_domain, _ = cad_modules
    with pytest.raises(fluid_domain.FluidDomainError, match="tank_body"):
        fluid_domain.build_fluid_domain(
            invalid,
            y_min_mm=Y_MIN_MM,
            y_max_mm=Y_MAX_MM,
        )


def test_rejects_non_solid_and_missing_rim_faces(cad_modules) -> None:
    cq, _, fluid_domain, _ = cad_modules
    face = cq.Face.makePlane(1.0, 1.0)
    with pytest.raises(fluid_domain.FluidDomainError, match="one input tank solid"):
        fluid_domain.build_fluid_domain(
            face,
            y_min_mm=Y_MIN_MM,
            y_max_mm=Y_MAX_MM,
        )

    box = cq.Workplane("XY").box(2.0, 2.0, 2.0).val()
    with pytest.raises(
        fluid_domain.FluidDomainError,
        match=r"rim candidate.*observed 0",
    ):
        fluid_domain.build_fluid_domain(
            box,
            y_min_mm=-3.0,
            y_max_mm=3.0,
            plane_tolerance_mm=PLANE_TOLERANCE_MM,
        )


def test_rejects_ambiguous_concentric_wet_loop_candidates(cad_modules) -> None:
    cq, _, fluid_domain, _ = cad_modules
    wires = (
        fluid_domain._LoopInventory(
            cq.Wire.makeCircle(3.0, cq.Vector(), cq.Vector(0.0, 1.0, 0.0)),
            30.0,
            (0.0, 0.0, 0.0),
        ),
        fluid_domain._LoopInventory(
            cq.Wire.makeCircle(2.0, cq.Vector(), cq.Vector(0.0, 1.0, 0.0)),
            20.0,
            (0.0, 0.0, 0.0),
        ),
        fluid_domain._LoopInventory(
            cq.Wire.makeCircle(1.0, cq.Vector(), cq.Vector(0.0, 1.0, 0.0)),
            10.0,
            (0.0, 0.0, 0.0),
        ),
    )

    with pytest.raises(fluid_domain.FluidDomainError, match=r"wet-loop.*observed 2"):
        fluid_domain._select_wet_loop(
            wires,
            plane_tolerance_mm=PLANE_TOLERANCE_MM,
            closure_y_mm=0.0,
        )


def test_round_trip_rejects_non_solid_input(cad_modules, tmp_path: Path) -> None:
    cq, audit_module, fluid_domain, _ = cad_modules
    face = cq.Face.makePlane(1.0, 1.0)

    with pytest.raises(fluid_domain.FluidDomainError, match="solid"):
        audit_module.write_step_round_trip(
            face,
            tmp_path / "invalid.step",
            tmp_path / "invalid.audit.json",
        )


def test_live_source_remains_byte_for_byte_unchanged(
    cad_modules,
    source_snapshot,
    tank_body,
    construction_diagnostics,
    round_trip,
    final_round_trip,
) -> None:
    _, _, _, xcaf = cad_modules
    del tank_body, construction_diagnostics, round_trip, final_round_trip

    size, mtime_ns, source_hash = source_snapshot
    assert (size, source_hash) == (EXPECTED_SIZE, EXPECTED_HASH)
    assert SOURCE.stat().st_size == size
    assert SOURCE.stat().st_mtime_ns == mtime_ns
    assert xcaf.sha256_file(SOURCE) == source_hash
