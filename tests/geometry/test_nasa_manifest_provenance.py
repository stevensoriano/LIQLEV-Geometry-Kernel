"""Four-criteria tests for NASA manifest WRITER F8 provenance hardening.

Call chain (real entry points):
  write_result_manifest
    -> _git_worktree_is_dirty  (refuse dirty tree)
    -> build_result_manifest
         -> _git_describe_dirty / solver_describe kwarg
         -> sha256_file(NASA_HARNESS_MODULE_PATH)
    -> Path.write_text (target = explicit path or MANIFEST_PATH)

All tests write only under tmp_path; the frozen D1 manifest is never opened for write.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from validation.custom_geometry_cases import (
    GEOMETRY_NPZ_PATH,
    NASA_HARNESS_MODULE_PATH,
    run_nasa_tank_validation,
    sha256_file,
    write_result_manifest,
)


@pytest.fixture(scope="module")
def nasa_validation():
    """Real NasaTankValidation so clean-write exercises the live build path."""

    return run_nasa_tank_validation()


def test_nasa_manifest_writer_refuses_dirty_tree_via_write_entry(
    monkeypatch, tmp_path: Path, nasa_validation
) -> None:
    """1 FIRES on bad (dirty) · 2 recovery-form raise · 3/4 real write_result_manifest."""

    monkeypatch.setattr(
        "validation.custom_geometry_cases._git_worktree_is_dirty",
        lambda: True,
    )
    target = tmp_path / "nasa_tank_geometry_manifest.json"

    with pytest.raises(RuntimeError, match="dirty git worktree"):
        write_result_manifest(nasa_validation, target)

    # Recovered state: no file written under the explicit path.
    assert not target.exists()


def test_nasa_manifest_writer_clean_write_records_describe_and_harness_hash(
    monkeypatch, tmp_path: Path, nasa_validation
) -> None:
    """1 clean case · 2 asserts recovered payload fields · 3/4 real write path to tmp."""

    monkeypatch.setattr(
        "validation.custom_geometry_cases._git_worktree_is_dirty",
        lambda: False,
    )
    monkeypatch.setattr(
        "validation.custom_geometry_cases._git_describe_dirty",
        lambda: "test-nasa-describe-dirty",
    )
    monkeypatch.setattr(
        "validation.custom_geometry_cases._git_head",
        lambda: "deadbeefcafefeed",
    )

    target = tmp_path / "nasa_tank_geometry_manifest.json"
    payload = write_result_manifest(nasa_validation, target)

    assert target.is_file()
    on_disk = json.loads(target.read_text(encoding="utf-8"))
    assert on_disk == payload

    assert payload["schema"] == "liqlev.validation.nasa_tank_geometry"
    assert payload["solver_describe"] == "test-nasa-describe-dirty"
    assert payload["solver_commit"] == "deadbeefcafefeed"
    assert payload["harness_module"] == "validation/custom_geometry_cases.py"
    assert payload["harness_module_sha256"] == sha256_file(NASA_HARNESS_MODULE_PATH)
    assert payload["geometry_npz_sha256"] == sha256_file(GEOMETRY_NPZ_PATH)
    assert "python" in payload["versions"]
    assert "numpy" in payload["versions"]
    assert "numba" in payload["versions"]
    assert payload["maximum_convergence_count"] == (
        nasa_validation.maximum_convergence_count
    )
    assert payload["passed"] is nasa_validation.passed
