from __future__ import annotations

import json
from pathlib import Path

import pytest


SOURCE = Path(
    r"C:\Users\sasorian\Documents\Eta_Space\geometry"
    r"\nhq01-m21a- 0201_TankAssy_NASA.STEP"
)
EXPECTED_HASH = "0193EC296B5B754FFDC7B1652052F5B28BA99345A5BC89922D026DF011C837D5"
PROVENANCE = (
    Path(__file__).resolve().parents[2] / "geometry" / "source" / "PROVENANCE.json"
)
requires_source = pytest.mark.skipif(
    not SOURCE.is_file(),
    reason="NASA source STEP is not present",
)


@pytest.fixture(scope="module")
def cad_modules():
    cq = pytest.importorskip("cadquery", reason="CadQuery is not installed")
    pytest.importorskip("OCP", reason="OpenCascade Python bindings are not installed")
    import liqlev.cad.xcaf as xcaf

    return cq, xcaf


@requires_source
def test_selects_resolved_tank_body(cad_modules) -> None:
    _, xcaf = cad_modules
    source_stat_before = SOURCE.stat()
    sibling_names_before = {path.name for path in SOURCE.parent.iterdir()}
    assert xcaf.sha256_file(SOURCE) == EXPECTED_HASH

    shape = xcaf.load_named_product(
        SOURCE,
        product_name="nhq01-m21a- 0202_short",
        expected_sha256=EXPECTED_HASH,
    )

    box = shape.BoundingBox()
    assert box.xmin == pytest.approx(-280.94, abs=1e-6)
    assert box.xmax == pytest.approx(280.94, abs=1e-6)
    assert box.ymin == pytest.approx(-293.966791, abs=1e-6)
    assert box.ymax == pytest.approx(293.966791, abs=1e-6)
    assert box.zmin == pytest.approx(296.542759, abs=1e-6)
    assert box.zmax == pytest.approx(858.422760, abs=1e-6)
    assert len(shape.Faces()) == 342
    source_stat_after = SOURCE.stat()
    assert source_stat_after.st_size == source_stat_before.st_size
    assert source_stat_after.st_mtime_ns == source_stat_before.st_mtime_ns
    assert {path.name for path in SOURCE.parent.iterdir()} == sibling_names_before
    assert xcaf.sha256_file(SOURCE) == EXPECTED_HASH


def test_sha256_file_returns_uppercase_hexadecimal(
    tmp_path: Path,
    cad_modules,
) -> None:
    _, xcaf = cad_modules
    source = tmp_path / "source.step"
    source.write_bytes(b"abc")

    assert xcaf.sha256_file(source) == (
        "BA7816BF8F01CFEA414140DE5DAE2223"
        "B00361A396177A9CB410FF61F20015AD"
    )


def test_hash_mismatch_fails_before_step_reader_is_constructed(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    cad_modules,
) -> None:
    _, xcaf = cad_modules
    source = tmp_path / "source.step"
    source.write_bytes(b"not the approved source")

    def fail_if_constructed() -> None:
        raise AssertionError("STEP reader must not be constructed after a hash mismatch")

    monkeypatch.setattr(xcaf, "STEPCAFControl_Reader", fail_if_constructed)

    with pytest.raises(xcaf.StepProductError, match="SHA-256 mismatch"):
        xcaf.load_named_product(
            source,
            product_name="nhq01-m21a- 0202_short",
            expected_sha256="0" * 64,
        )


def test_reader_constructor_failure_is_wrapped_after_valid_hash(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    cad_modules,
) -> None:
    _, xcaf = cad_modules
    source = tmp_path / "source.step"
    source.write_bytes(b"approved source")
    expected_hash = xcaf.sha256_file(source)
    constructor_error = RuntimeError("reader construction failed")

    def fail_during_construction() -> None:
        raise constructor_error

    monkeypatch.setattr(xcaf, "STEPCAFControl_Reader", fail_during_construction)

    with pytest.raises(
        xcaf.StepProductError,
        match="STEP read, transfer, or XCAF traversal failed",
    ) as caught:
        xcaf.load_named_product(
            source,
            product_name="nhq01-m21a- 0202_short",
            expected_sha256=expected_hash,
        )

    assert caught.value.__cause__ is constructor_error


def test_missing_source_raises_step_product_error(
    tmp_path: Path,
    cad_modules,
) -> None:
    _, xcaf = cad_modules
    missing_source = tmp_path / "missing.step"

    with pytest.raises(xcaf.StepProductError, match="missing or unreadable"):
        xcaf.load_named_product(
            missing_source,
            product_name="nhq01-m21a- 0202_short",
            expected_sha256="0" * 64,
        )


@requires_source
def test_missing_product_name_reports_zero_exact_matches(cad_modules) -> None:
    _, xcaf = cad_modules
    with pytest.raises(xcaf.StepProductError, match="zero exact product matches"):
        xcaf.load_named_product(
            SOURCE,
            product_name="nhq01-m21a- 0202_short_DOES_NOT_EXIST",
            expected_sha256=EXPECTED_HASH,
        )


def test_multiple_product_matches_report_each_instance(
    tmp_path: Path,
    cad_modules,
) -> None:
    cq, xcaf = cad_modules
    shape = cq.Workplane("XY").box(1.0, 1.0, 1.0).val()
    matches = (
        xcaf._ProductMatch("0:1:1:90", "0:1:1:1:5", "duplicate", shape),
        xcaf._ProductMatch("0:1:1:90", "0:1:1:1:6", "duplicate", shape),
    )

    with pytest.raises(
        xcaf.StepProductError,
        match=r"multiple exact product matches \(2\).*0:1:1:1:5, 0:1:1:1:6",
    ):
        xcaf._select_unique_match(
            matches,
            product_name="duplicate",
            source=tmp_path / "duplicate.step",
        )


@requires_source
def test_provenance_identifies_the_read_only_source() -> None:
    provenance = json.loads(PROVENANCE.read_text(encoding="utf-8"))

    assert Path(provenance["source_path"]) == SOURCE
    assert provenance["source_sha256"] == EXPECTED_HASH
    assert provenance["source_size_bytes"] == SOURCE.stat().st_size
    assert provenance["tank_product"] == "nhq01-m21a- 0202_short"
    assert provenance["height_axis"] == "+Y"
    assert provenance["gravity_direction"] == "-Y"
    assert provenance["closure_planes_mm"] == [-275.406791, 275.296791]
