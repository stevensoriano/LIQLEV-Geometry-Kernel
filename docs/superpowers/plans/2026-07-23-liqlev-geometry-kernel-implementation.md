# LIQLEV Geometry Kernel Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a reproducible CAD-to-Numba geometry pipeline, generate the cleaned NASA tank fluid-domain STEP, and add a volume-conserving custom-geometry path to LIQLEV while preserving the analytic cylinder path.

**Architecture:** OpenCascade runs only in offline preprocessing: it selects the tank body, constructs a single capped fluid solid, and exports exact section, volume, perimeter, and wall-area measurements. The solver loads a versioned numeric package and passes contiguous `float64` arrays to Numba functions for interpolation, volume inversion, and perimeter-aware boundary-layer integration. `core.py` keeps its current analytic cylinder branch and selects the new branch only when an explicit geometry package is supplied.

**Tech Stack:** Python 3.13, NumPy 2.3, Pandas, Numba 0.64, SciPy 1.16, CadQuery 2.7/OpenCascade 7.8, CoolProp, Pytest 8.4, STEP AP214.

## Global Constraints

- Keep every created or modified file under `C:\Users\sasorian\Documents\Eta_Space\LIQLEV-Geometry-Kernel`.
- Treat `C:\Users\sasorian\Documents\Eta_Space\geometry\nhq01-m21a- 0201_TankAssy_NASA.STEP` as read-only.
- Require source STEP SHA-256 `0193EC296B5B754FFDC7B1652052F5B28BA99345A5BC89922D026DF011C837D5`.
- Preserve the source assembly coordinate system; solver height increases along assembly `+Y`, and nominal gravity points along `-Y`.
- Use closure planes `Y=-275.406791 mm` and `Y=+275.296791 mm`.
- Keep CAD calculations in millimetres and export solver arrays in feet, square feet, and cubic feet.
- Produce exactly one valid, closed, positive-volume fluid solid without baffles, PMD, spray hardware, fasteners, or exterior viewport hardware.
- Keep OpenCascade and SciPy outside the JIT solver. Numba receives only scalars and contiguous `float64` arrays.
- Preserve the existing analytic cylinder path and require `python scripts/check_physics_baseline.py` to pass at every physics-critical commit.
- Use test-driven development: add a failing focused test, observe the expected failure, implement the smallest complete change, then run focused and regression tests.
- Do not include Codex or OpenAI as authors, co-authors, commit authors, or PR authors.

---

## Planned File Structure

### Preserved LIQLEV baseline

- `core.py` — legacy wrapper and JIT loop; later gains an explicit geometry-mode branch.
- `thermo_utils.py` — unchanged thermophysical-property authority.
- `liqlev/model/config.py` — adds an optional geometry-package path.
- `liqlev/model/builder.py` — builds cylinder or custom solver inputs.
- `liqlev/model/validation.py` — validates geometry configuration before execution.
- `liqlev/runner/single.py` — loads one geometry package per run.
- `liqlev/runner/sweep.py` — reports the actual custom tank height.
- `validation/` and `tests/` — preserved legacy baselines and tests.

### Numeric geometry kernel

- `liqlev/geometry/__init__.py` — supported public imports.
- `liqlev/geometry/schema.py` — immutable metadata and array container.
- `liqlev/geometry/package.py` — package validation, hashing, NPZ/JSON/CSV persistence.
- `liqlev/geometry/coefficients.py` — offline monotone-cubic coefficient construction.
- `liqlev/geometry/jit.py` — Numba interpolation, derivative, inverse-volume, and boundary-layer functions.
- `liqlev/geometry/fixtures.py` — analytic cylinder and sphere packages for tests.

### Offline CAD pipeline

- `liqlev/cad/__init__.py` — CAD public imports.
- `liqlev/cad/xcaf.py` — STEP/XCAF loading and named-product selection.
- `liqlev/cad/fluid_domain.py` — inner-wall traversal, port capping, sewing, and solid construction.
- `liqlev/cad/measure.py` — clipped-volume, section-area, perimeter, and wall-area measurements.
- `liqlev/cad/audit.py` — STEP round-trip and acceptance report.
- `scripts/build_nasa_tank_geometry.py` — deterministic end-to-end artifact command.
- `scripts/check_nasa_tank_geometry.py` — independent artifact verification command.

### Tests and artifacts

- `tests/geometry/test_package.py`
- `tests/geometry/test_jit_interpolation.py`
- `tests/geometry/test_boundary_layer.py`
- `tests/geometry/test_core_custom_geometry.py`
- `tests/cad/test_xcaf_selection.py`
- `tests/cad/test_fluid_domain.py`
- `tests/cad/test_analytic_measurements.py`
- `tests/cad/test_nasa_tank_artifacts.py`
- `geometry/source/PROVENANCE.json`
- `geometry/output/nhq01-m21a-0201_LIQLEV_FLUID_DOMAIN.step`
- `geometry/tables/nhq01-m21a-0201_LIQLEV_GEOMETRY.npz`
- `geometry/tables/nhq01-m21a-0201_LIQLEV_GEOMETRY.json`
- `geometry/tables/nhq01-m21a-0201_LIQLEV_GEOMETRY.csv`
- `geometry/audit/nhq01-m21a-0201_LIQLEV_AUDIT.json`
- `geometry/audit/sections/` — orthogonal validation images.
- `docs/geometry-kernel.md` — regeneration and solver-use instructions.

---

### Task 1: Import and Freeze the Current LIQLEV Baseline

**Files:**
- Create from source: `core.py`, `thermo_utils.py`, `gui.py`, `liqlev/`, `tests/`, `validation/`, `scripts/`, `packaging/`, `data/`
- Create from source: `AGENTS.md`, `.gitignore`, `LICENSE`, `README.md`, `requirements.txt`, `requirements-modernization.txt`, `LIQLEV.spec`, `LIQLEV-Modern.spec`, `MODERNIZATION_GOAL.md`
- Create: `validation/imported_solver_manifest.json`

**Interfaces:**
- Consumes: the read-only build at `C:\Users\sasorian\Documents\Cryo Vent LLR\LIQLEV_current_build_team_share_2026-05-13\github_repo_source`.
- Produces: a committed, testable LIQLEV baseline whose physics hashes and result baseline are fixed before modification.

- [ ] **Step 1: Copy the approved source baseline without generated state**

Run from the repository root:

```powershell
$source = 'C:\Users\sasorian\Documents\Cryo Vent LLR\LIQLEV_current_build_team_share_2026-05-13\github_repo_source'
$names = @(
  'core.py', 'thermo_utils.py', 'gui.py', 'liqlev', 'tests', 'validation',
  'scripts', 'packaging', 'data', 'AGENTS.md', 'LICENSE',
  'README.md', 'requirements.txt', 'requirements-modernization.txt',
  'LIQLEV.spec', 'LIQLEV-Modern.spec', 'MODERNIZATION_GOAL.md'
)
foreach ($name in $names) {
  Copy-Item -LiteralPath (Join-Path $source $name) -Destination . -Recurse
}
```

Expected: the imported files appear in `git status`; `.git`, caches, build
outputs, and environments are absent.

Merge the source `.gitignore` rules into the existing repository file without
removing its required `.worktrees/` entry.

- [ ] **Step 2: Record the immutable import manifest**

Create `validation/imported_solver_manifest.json` with exactly:

```json
{
  "source_directory": "C:\\Users\\sasorian\\Documents\\Cryo Vent LLR\\LIQLEV_current_build_team_share_2026-05-13\\github_repo_source",
  "import_date": "2026-07-23",
  "files": {
    "core.py": "B715978F408C7B4D66D499486F3A1FCD0F9177B12A54B981F9785A15D3124F29",
    "thermo_utils.py": "5D5817D58DBBC50F8E332A9CE34D9B84EEAAF8CEDC5F45FD5E37C82C6546C2EA",
    "validation/baselines/physics_baseline.json": "CD45206BECF128A28C55C7B00DAF494AC9FC8BD43FF356F724F95E2CF1A7E747"
  }
}
```

- [ ] **Step 3: Run the legacy physics gate**

Run:

```powershell
python scripts/check_physics_baseline.py
```

Expected final line:

```text
Physics baseline check passed.
```

- [ ] **Step 4: Run the imported headless tests**

Run:

```powershell
python -m pytest tests -q
```

Expected: all imported tests pass with zero failures.

- [ ] **Step 5: Commit the baseline**

```powershell
git add AGENTS.md .gitignore LICENSE README.md requirements.txt requirements-modernization.txt LIQLEV.spec LIQLEV-Modern.spec MODERNIZATION_GOAL.md core.py thermo_utils.py gui.py liqlev tests validation scripts packaging data
git commit -m "chore: import LIQLEV solver baseline"
```

---

### Task 2: Define and Validate the Geometry Package

**Files:**
- Create: `liqlev/geometry/__init__.py`
- Create: `liqlev/geometry/schema.py`
- Create: `liqlev/geometry/package.py`
- Create: `tests/geometry/test_package.py`

**Interfaces:**
- Consumes: sampled arrays in solver British units and a JSON-compatible metadata dictionary.
- Produces: `GeometryKernel`, `GeometryMetadata`, `load_geometry_package(path)`, `save_geometry_package(kernel, path)`, and `validate_geometry_kernel(kernel)`.

- [ ] **Step 1: Write failing package-validation tests**

Create `tests/geometry/test_package.py`:

```python
from __future__ import annotations

import json

import numpy as np
import pytest

from liqlev.geometry.package import (
    GeometryPackageError,
    load_geometry_package,
    save_geometry_package,
    validate_geometry_kernel,
)
from liqlev.geometry.schema import GeometryKernel, GeometryMetadata


def valid_kernel() -> GeometryKernel:
    h = np.array([0.0, 1.0, 2.0], dtype=np.float64)
    volume = np.array([0.0, 3.0, 6.0], dtype=np.float64)
    area = np.array([3.0, 3.0, 3.0], dtype=np.float64)
    perimeter = np.array([2.0, 2.0, 2.0], dtype=np.float64)
    wall = np.array([0.0, 2.0, 4.0], dtype=np.float64)
    coeff = np.array(
        [[0.0, 0.0], [0.0, 0.0], [3.0, 3.0], [0.0, 3.0]],
        dtype=np.float64,
    )
    wall_coeff = np.array(
        [[0.0, 0.0], [0.0, 0.0], [2.0, 2.0], [0.0, 2.0]],
        dtype=np.float64,
    )
    return GeometryKernel(
        metadata=GeometryMetadata(
            schema_version=1,
            geometry_id="unit-cylinder",
            source_step_sha256="0" * 64,
            fluid_step_sha256="1" * 64,
            axis="+Y",
            gravity_direction="-Y",
            length_unit="ft",
            area_unit="ft^2",
            volume_unit="ft^3",
            y_min_mm=0.0,
            y_max_mm=609.6,
        ),
        height_ft=h,
        volume_ft3=volume,
        volume_coefficients=coeff,
        section_area_ft2=area,
        perimeter_ft=perimeter,
        sidewall_area_ft2=wall,
        sidewall_coefficients=wall_coeff,
        total_wetted_area_ft2=np.array([3.0, 5.0, 7.0]),
    )


def test_package_round_trip_preserves_arrays_and_metadata(tmp_path) -> None:
    target = tmp_path / "geometry.npz"
    source = valid_kernel()
    save_geometry_package(source, target)
    loaded = load_geometry_package(target)
    assert loaded.metadata == source.metadata
    for name in source.array_names():
        np.testing.assert_array_equal(getattr(loaded, name), getattr(source, name))
    metadata = json.loads(target.with_suffix(".json").read_text(encoding="utf-8"))
    assert metadata["npz_sha256"]


def test_validation_rejects_nonmonotone_volume() -> None:
    kernel = valid_kernel()
    bad = GeometryKernel(
        **{
            **kernel.as_mapping(),
            "volume_ft3": np.array([0.0, 4.0, 3.0], dtype=np.float64),
        }
    )
    with pytest.raises(GeometryPackageError, match="volume_ft3"):
        validate_geometry_kernel(bad)


def test_load_rejects_changed_binary_hash(tmp_path) -> None:
    target = tmp_path / "geometry.npz"
    save_geometry_package(valid_kernel(), target)
    with target.open("ab") as file_obj:
        file_obj.write(b"changed")
    with pytest.raises(GeometryPackageError, match="SHA-256"):
        load_geometry_package(target)
```

- [ ] **Step 2: Run the test and observe the missing-module failure**

Run:

```powershell
python -m pytest tests/geometry/test_package.py -q
```

Expected: collection fails because `liqlev.geometry` does not exist.

- [ ] **Step 3: Implement the immutable schema**

Create `liqlev/geometry/schema.py` with:

```python
from __future__ import annotations

from dataclasses import asdict, dataclass, fields
from typing import Any

import numpy as np


@dataclass(frozen=True)
class GeometryMetadata:
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
```

Create `liqlev/geometry/__init__.py`:

```python
from .package import load_geometry_package, save_geometry_package
from .schema import GeometryKernel, GeometryMetadata

__all__ = [
    "GeometryKernel",
    "GeometryMetadata",
    "load_geometry_package",
    "save_geometry_package",
]
```

- [ ] **Step 4: Implement persistence and strict validation**

Create `liqlev/geometry/package.py` with functions that perform these exact
checks before save and after load:

```python
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from .schema import GeometryKernel, GeometryMetadata


class GeometryPackageError(ValueError):
    pass


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file_obj:
        for block in iter(lambda: file_obj.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def _float64_contiguous(name: str, value: np.ndarray) -> np.ndarray:
    array = np.asarray(value)
    if array.dtype != np.float64:
        raise GeometryPackageError(f"{name} must use float64")
    if array.ndim not in (1, 2) or not array.flags.c_contiguous:
        raise GeometryPackageError(f"{name} must be a contiguous 1D or 2D array")
    if not np.all(np.isfinite(array)):
        raise GeometryPackageError(f"{name} contains non-finite values")
    return array


def validate_geometry_kernel(kernel: GeometryKernel) -> None:
    if kernel.metadata.schema_version != 1:
        raise GeometryPackageError("schema_version must equal 1")
    expected_units = ("+Y", "-Y", "ft", "ft^2", "ft^3")
    actual_units = (
        kernel.metadata.axis,
        kernel.metadata.gravity_direction,
        kernel.metadata.length_unit,
        kernel.metadata.area_unit,
        kernel.metadata.volume_unit,
    )
    if actual_units != expected_units:
        raise GeometryPackageError(f"axis/unit contract mismatch: {actual_units}")
    for name in kernel.array_names():
        _float64_contiguous(name, getattr(kernel, name))
    count = len(kernel.height_ft)
    for name in (
        "volume_ft3",
        "section_area_ft2",
        "perimeter_ft",
        "sidewall_area_ft2",
        "total_wetted_area_ft2",
    ):
        if len(getattr(kernel, name)) != count:
            raise GeometryPackageError(f"{name} length must match height_ft")
    if kernel.volume_coefficients.shape != (4, count - 1):
        raise GeometryPackageError("volume_coefficients must have shape (4, N-1)")
    if kernel.sidewall_coefficients.shape != (4, count - 1):
        raise GeometryPackageError("sidewall_coefficients must have shape (4, N-1)")
    if count < 3 or np.any(np.diff(kernel.height_ft) <= 0.0):
        raise GeometryPackageError("height_ft must contain at least 3 increasing nodes")
    for name in ("volume_ft3", "sidewall_area_ft2", "total_wetted_area_ft2"):
        if np.any(np.diff(getattr(kernel, name)) < 0.0):
            raise GeometryPackageError(f"{name} must be non-decreasing")
    for name in ("section_area_ft2", "perimeter_ft"):
        if np.any(getattr(kernel, name) < 0.0):
            raise GeometryPackageError(f"{name} must be non-negative")
    if abs(kernel.height_ft[0]) > 1e-12 or abs(kernel.volume_ft3[0]) > 1e-12:
        raise GeometryPackageError("height_ft and volume_ft3 must begin at zero")
    if kernel.total_height_ft <= 0.0 or kernel.total_volume_ft3 <= 0.0:
        raise GeometryPackageError("geometry totals must be positive")


def save_geometry_package(kernel: GeometryKernel, path: str | Path) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    validate_geometry_kernel(kernel)
    np.savez_compressed(
        target,
        **{name: getattr(kernel, name) for name in kernel.array_names()},
    )
    payload = kernel.metadata.as_dict()
    payload["npz_sha256"] = _sha256(target)
    target.with_suffix(".json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    csv_path = target.with_suffix(".csv")
    with csv_path.open("w", newline="", encoding="utf-8") as file_obj:
        writer = csv.writer(file_obj)
        writer.writerow(
            [
                "height_ft",
                "volume_ft3",
                "section_area_ft2",
                "perimeter_ft",
                "sidewall_area_ft2",
                "total_wetted_area_ft2",
            ]
        )
        writer.writerows(
            zip(
                kernel.height_ft,
                kernel.volume_ft3,
                kernel.section_area_ft2,
                kernel.perimeter_ft,
                kernel.sidewall_area_ft2,
                kernel.total_wetted_area_ft2,
                strict=True,
            )
        )


def load_geometry_package(path: str | Path) -> GeometryKernel:
    target = Path(path)
    metadata_path = target.with_suffix(".json")
    payload = json.loads(metadata_path.read_text(encoding="utf-8"))
    expected_hash = payload.pop("npz_sha256").upper()
    if _sha256(target) != expected_hash:
        raise GeometryPackageError("NPZ SHA-256 does not match metadata")
    metadata = GeometryMetadata(**payload)
    with np.load(target, allow_pickle=False) as archive:
        arrays = {
            name: np.ascontiguousarray(archive[name], dtype=np.float64)
            for name in GeometryKernel.array_names()
        }
    kernel = GeometryKernel(metadata=metadata, **arrays)
    validate_geometry_kernel(kernel)
    return kernel
```

- [ ] **Step 5: Run focused and imported tests**

Run:

```powershell
python -m pytest tests/geometry/test_package.py tests/test_smoke.py -q
```

Expected: all selected tests pass.

- [ ] **Step 6: Commit**

```powershell
git add liqlev/geometry tests/geometry/test_package.py
git commit -m "feat: define versioned geometry package"
```

---

### Task 3: Add Monotone Coefficients and Numba Geometry Evaluation

**Files:**
- Create: `liqlev/geometry/coefficients.py`
- Create: `liqlev/geometry/jit.py`
- Create: `tests/geometry/test_jit_interpolation.py`

**Interfaces:**
- Consumes: strictly increasing nodes and non-decreasing cumulative values.
- Produces: `pchip_coefficients(x, y) -> ndarray[4,N-1]`, `eval_ppoly`, `eval_ppoly_derivative`, `interp_linear_nonnegative`, and `invert_monotone_volume`, all runtime evaluators compatible with Numba nopython mode.

- [ ] **Step 1: Write failing interpolation and inverse tests**

Create `tests/geometry/test_jit_interpolation.py`:

```python
from __future__ import annotations

import numpy as np
import pytest

from liqlev.geometry.coefficients import pchip_coefficients
from liqlev.geometry.jit import (
    eval_ppoly,
    eval_ppoly_derivative,
    interp_linear_nonnegative,
    invert_monotone_volume,
)


def test_linear_volume_is_exact_and_invertible() -> None:
    height = np.linspace(0.0, 4.0, 9)
    volume = 3.25 * height
    coefficients = pchip_coefficients(height, volume)
    for h in np.linspace(0.0, 4.0, 33):
        assert eval_ppoly(h, height, coefficients) == pytest.approx(3.25 * h)
        assert eval_ppoly_derivative(h, height, coefficients) == pytest.approx(3.25)
        assert invert_monotone_volume(
            3.25 * h, height, volume, coefficients
        ) == pytest.approx(h, abs=1e-11)


def test_inverse_rejects_out_of_domain_volume() -> None:
    height = np.array([0.0, 1.0, 2.0])
    volume = height**3
    coefficients = pchip_coefficients(height, volume)
    assert np.isnan(invert_monotone_volume(-1.0, height, volume, coefficients))
    assert np.isnan(invert_monotone_volume(9.0, height, volume, coefficients))


def test_perimeter_interpolation_is_nonnegative() -> None:
    nodes = np.array([0.0, 1.0, 2.0])
    values = np.array([0.0, 2.0, 0.0])
    assert interp_linear_nonnegative(0.5, nodes, values) == pytest.approx(1.0)
    assert interp_linear_nonnegative(3.0, nodes, values) == pytest.approx(0.0)
```

- [ ] **Step 2: Run the tests and observe missing imports**

Run:

```powershell
python -m pytest tests/geometry/test_jit_interpolation.py -q
```

Expected: collection fails because the coefficient and JIT modules do not
exist.

- [ ] **Step 3: Implement coefficient generation**

Create `liqlev/geometry/coefficients.py`. Use the Fritsch-Carlson/PCHIP
endpoint and interior slope rules, then export local-power coefficients in
Numba's expected order:

```python
from __future__ import annotations

import numpy as np


def pchip_coefficients(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y) or len(x) < 3:
        raise ValueError("x and y must be equal-length 1D arrays with N >= 3")
    h = np.diff(x)
    if np.any(h <= 0.0):
        raise ValueError("x must be strictly increasing")
    delta = np.diff(y) / h
    slopes = np.zeros_like(y)
    for index in range(1, len(y) - 1):
        if delta[index - 1] * delta[index] > 0.0:
            w1 = 2.0 * h[index] + h[index - 1]
            w2 = h[index] + 2.0 * h[index - 1]
            slopes[index] = (w1 + w2) / (
                w1 / delta[index - 1] + w2 / delta[index]
            )
    slopes[0] = ((2.0 * h[0] + h[1]) * delta[0] - h[0] * delta[1]) / (
        h[0] + h[1]
    )
    slopes[-1] = (
        (2.0 * h[-1] + h[-2]) * delta[-1] - h[-1] * delta[-2]
    ) / (h[-1] + h[-2])
    if slopes[0] * delta[0] <= 0.0:
        slopes[0] = 0.0
    elif abs(slopes[0]) > 3.0 * abs(delta[0]):
        slopes[0] = 3.0 * delta[0]
    if slopes[-1] * delta[-1] <= 0.0:
        slopes[-1] = 0.0
    elif abs(slopes[-1]) > 3.0 * abs(delta[-1]):
        slopes[-1] = 3.0 * delta[-1]
    coefficients = np.empty((4, len(x) - 1), dtype=np.float64)
    coefficients[0] = (slopes[:-1] + slopes[1:] - 2.0 * delta) / h**2
    coefficients[1] = (3.0 * delta - 2.0 * slopes[:-1] - slopes[1:]) / h
    coefficients[2] = slopes[:-1]
    coefficients[3] = y[:-1]
    return np.ascontiguousarray(coefficients)
```

- [ ] **Step 4: Implement the Numba evaluators and safeguarded inverse**

Create `liqlev/geometry/jit.py` with `@njit(cache=True)` functions. Interval
selection must clamp function evaluation at endpoints, while inverse-volume
must return `NaN` outside `[0, V_total]`. The inverse uses the `volume_ft3`
nodes to select an interval, then performs at most 64 Newton/bisection
iterations with a height tolerance of `1e-12 * max(1, total_height)`.
For a near-flat segment, the inverse contract is conditioning-aware: its
height acceptance bound is
`max(base_height_tolerance, 2 * volume_ulp / max(abs(local_dVdh), tiny))`,
where `volume_ulp = spacing(max(abs(target_volume), 1.0))`. The returned
height must remain inside the final monotone bracket (its midpoint), and its
forward interpolated volume must agree with the target within
`max(volume_ulp, 2 * abs(local_dVdh) * height_ulp)`, where
`height_ulp = spacing(max(abs(candidate_height), 1.0))`. One-volume-ULP
agreement remains required whenever representable height resolution permits it.

Use this complete implementation:

```python
from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True)
def _find_interval(x: float, nodes: np.ndarray) -> int:
    if x <= nodes[0]:
        return 0
    if x >= nodes[-1]:
        return len(nodes) - 2
    low = 0
    high = len(nodes) - 1
    while low < high - 1:
        middle = (low + high) >> 1
        if nodes[middle] <= x:
            low = middle
        else:
            high = middle
    return low


@njit(cache=True)
def eval_ppoly(x: float, breaks: np.ndarray, coefficients: np.ndarray) -> float:
    index = _find_interval(x, breaks)
    bounded = min(max(x, breaks[0]), breaks[-1])
    dx = bounded - breaks[index]
    return (
        (
            coefficients[0, index] * dx
            + coefficients[1, index]
        )
        * dx
        + coefficients[2, index]
    ) * dx + coefficients[3, index]


@njit(cache=True)
def eval_ppoly_derivative(
    x: float, breaks: np.ndarray, coefficients: np.ndarray
) -> float:
    index = _find_interval(x, breaks)
    bounded = min(max(x, breaks[0]), breaks[-1])
    dx = bounded - breaks[index]
    return (
        3.0 * coefficients[0, index] * dx
        + 2.0 * coefficients[1, index]
    ) * dx + coefficients[2, index]


@njit(cache=True)
def interp_linear_nonnegative(
    x: float, nodes: np.ndarray, values: np.ndarray
) -> float:
    index = _find_interval(x, nodes)
    bounded = min(max(x, nodes[0]), nodes[-1])
    width = nodes[index + 1] - nodes[index]
    fraction = (bounded - nodes[index]) / width
    value = values[index] + fraction * (values[index + 1] - values[index])
    return max(0.0, value)


@njit(cache=True)
def invert_monotone_volume(
    target_volume: float,
    height: np.ndarray,
    volume: np.ndarray,
    volume_coefficients: np.ndarray,
) -> float:
    if target_volume < volume[0] or target_volume > volume[-1]:
        return np.nan
    if target_volume == volume[0]:
        return height[0]
    if target_volume == volume[-1]:
        return height[-1]

    low_index = 0
    high_index = len(volume) - 1
    while low_index < high_index - 1:
        middle = (low_index + high_index) >> 1
        if volume[middle] <= target_volume:
            low_index = middle
        else:
            high_index = middle

    if target_volume == volume[low_index]:
        return height[low_index]

    lower = height[low_index]
    upper = height[low_index + 1]
    volume_width = volume[low_index + 1] - volume[low_index]
    if volume_width <= 0.0:
        return np.nan
    guess = lower + (
        (target_volume - volume[low_index])
        / volume_width
        * (upper - lower)
    )
    tolerance = 1e-12 * max(1.0, height[-1])
    volume_ulp = np.spacing(max(abs(target_volume), 1.0))
    result = guess

    for _ in range(64):
        residual = (
            eval_ppoly(guess, height, volume_coefficients)
            - target_volume
        )
        if residual <= 0.0:
            lower = guess
        else:
            upper = guess
        derivative = eval_ppoly_derivative(
            guess, height, volume_coefficients
        )
        candidate = guess - residual / derivative if derivative > 0.0 else np.nan
        if not np.isfinite(candidate) or candidate <= lower or candidate >= upper:
            candidate = 0.5 * (lower + upper)
        guess = candidate
        if upper - lower <= tolerance:
            result = 0.5 * (lower + upper)
            result_residual = abs(
                eval_ppoly(result, height, volume_coefficients)
                - target_volume
            )
            lower_residual = abs(
                eval_ppoly(lower, height, volume_coefficients)
                - target_volume
            )
            upper_residual = abs(
                eval_ppoly(upper, height, volume_coefficients)
                - target_volume
            )
            if lower_residual < result_residual:
                result = lower
                result_residual = lower_residual
            if upper_residual < result_residual:
                result = upper
                result_residual = upper_residual
            local_derivative = abs(
                eval_ppoly_derivative(result, height, volume_coefficients)
            )
            height_ulp = np.spacing(max(abs(result), 1.0))
            allowed_volume_error = max(
                volume_ulp, 2.0 * local_derivative * height_ulp
            )
            if result_residual <= allowed_volume_error:
                return result
    return result
```

The function bodies implement binary interval search, local cubic Horner
evaluation, and a bracket-preserving Newton step. A Newton candidate outside
the active bracket is replaced with its midpoint. Return the exact first or
last height for endpoint volumes. Duplicate-volume nodes use the explicit
rightmost-node inverse convention. Near a plateau, test height recovery with
the conditioning-aware bound above rather than the unconditional base height
tolerance. The midpoint is preferred; a final bracket endpoint is used only
when it is the more precise in-bracket point under the feasible residual
contract.

- [ ] **Step 5: Verify values and nopython compilation**

Add this assertion to `tests/geometry/test_jit_interpolation.py`:

```python
def test_runtime_functions_compile_without_object_mode() -> None:
    height = np.array([0.0, 1.0, 2.0])
    volume = height**2
    coefficients = pchip_coefficients(height, volume)
    invert_monotone_volume(1.0, height, volume, coefficients)
    assert invert_monotone_volume.nopython_signatures
    assert eval_ppoly.nopython_signatures
```

Run:

```powershell
python -m pytest tests/geometry/test_jit_interpolation.py -q
```

Expected: all tests pass.

- [ ] **Step 6: Commit**

```powershell
git add liqlev/geometry/coefficients.py liqlev/geometry/jit.py tests/geometry/test_jit_interpolation.py
git commit -m "feat: add JIT geometry interpolation and inversion"
```

---

### Task 4: Build Analytic Fixtures and the Perimeter-Aware Boundary Layer

**Files:**
- Create: `liqlev/geometry/fixtures.py`
- Modify: `liqlev/geometry/jit.py`
- Create: `tests/geometry/test_boundary_layer.py`

**Interfaces:**
- Consumes: `AK3`, top height, volume cubic coefficients, height nodes, and perimeter nodes.
- Produces: `integrate_boundary_layer(ak3, top_height, height, volume_coefficients, perimeter, substeps) -> (delta_top, vbl, normalized_exit_flow, status)`.

- [ ] **Step 1: Write the failing analytic-cylinder test**

Create `tests/geometry/test_boundary_layer.py` with:

```python
from __future__ import annotations

import numpy as np
import pytest

from liqlev.geometry.fixtures import cylinder_kernel
from liqlev.geometry.jit import integrate_boundary_layer


def legacy_cylinder_profile(
    diameter_ft: float, ak3: float, height_ft: float
) -> tuple[float, float, float]:
    powers = np.arange(1.0, 11.0)
    coeff_s1 = (
        4.0 ** (powers - 1.0)
        / diameter_ft**powers
        / (2.0**powers + 1.0)
    )
    delta = (0.375 * diameter_ft * ak3 * height_ft) ** (2.0 / 3.0)
    for _ in range(30):
        residual = (
            8.0
            * np.sum(coeff_s1 * delta ** (powers + 0.5))
            / ak3
            - height_ft
        )
        derivative = (
            4.0
            * np.sum(
                4.0 ** (powers - 1.0)
                / diameter_ft**powers
                * delta ** (powers - 0.5)
            )
            / ak3
        )
        delta -= residual / derivative
    coeff_vbl = (
        (2.0 * powers + 1.0)
        / (powers + 1.5)
        / diameter_ft ** (powers - 1.0)
    )
    vbl = np.pi * np.sum(coeff_vbl * delta ** (powers + 1.5)) / ak3
    normalized_exit = (2.0 / 3.0) * np.pi * diameter_ft * delta**1.5
    return float(delta), float(vbl), float(normalized_exit)


@pytest.mark.parametrize("fill", [0.1, 0.25, 0.5, 0.8, 0.95])
def test_numeric_cylinder_boundary_layer_matches_heritage(fill: float) -> None:
    diameter = 4.0
    tank_height = 8.0
    kernel = cylinder_kernel(diameter, tank_height, node_count=1025)
    top_height = fill * tank_height
    expected = legacy_cylinder_profile(diameter, 0.015, top_height)
    actual = integrate_boundary_layer(
        0.015,
        top_height,
        kernel.height_ft,
        kernel.volume_coefficients,
        kernel.perimeter_ft,
        4,
    )
    assert actual[3] == 0
    np.testing.assert_allclose(actual[:3], expected, rtol=1e-3, atol=1e-10)
```

- [ ] **Step 2: Run and observe the missing fixture/function failure**

Run:

```powershell
python -m pytest tests/geometry/test_boundary_layer.py -q
```

Expected: collection fails for `cylinder_kernel` or
`integrate_boundary_layer`.

- [ ] **Step 3: Implement analytic fixture builders**

Create `liqlev/geometry/fixtures.py` with:

```python
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
    return GeometryKernel(
        metadata=_metadata("analytic-cylinder", height_ft),
        height_ft=np.ascontiguousarray(h),
        volume_ft3=np.ascontiguousarray(volume),
        volume_coefficients=pchip_coefficients(h, volume),
        section_area_ft2=np.full(node_count, area),
        perimeter_ft=np.full(node_count, perimeter),
        sidewall_area_ft2=np.ascontiguousarray(wall),
        sidewall_coefficients=pchip_coefficients(h, wall),
        total_wetted_area_ft2=np.ascontiguousarray(cap + wall),
    )
```

Add `sphere_kernel(radius_ft, node_count=1025)` using exact spherical-segment
relations:

```python
area = np.pi * (2.0 * radius_ft * h - h**2)
volume = np.pi * h**2 * (radius_ft - h / 3.0)
perimeter = 2.0 * np.pi * np.sqrt(np.maximum(0.0, 2.0 * radius_ft * h - h**2))
sidewall = 2.0 * np.pi * radius_ft * h
```

- [ ] **Step 4: Implement the Numba boundary-layer integrator**

In `liqlev/geometry/jit.py`, implement the normalized flow state
`q=(2/3)P*delta^(3/2)` and boundary-layer volume state `vbl`. For each table
interval below `top_height`, use four fixed RK4 substeps. Evaluate:

```python
delta = (1.5 * max(q, 0.0) / perimeter) ** (2.0 / 3.0)
dq_dh = ak3 * (area - perimeter * delta)
dvbl_dh = perimeter * delta
```

where `area=eval_ppoly_derivative(h, height, volume_coefficients)` and
`perimeter=interp_linear_nonnegative(h, height, perimeter_values)`. Clamp only
the integrated state `q` and `vbl` to zero; do not clamp a negative derivative.
Return status `1` for an out-of-domain height, non-positive `AK3`, or non-finite
state, and status `0` otherwise. At the top, return
`normalized_exit_flow=q_top`.

- [ ] **Step 5: Run boundary-layer and interpolation tests**

Run:

```powershell
python -m pytest tests/geometry/test_boundary_layer.py tests/geometry/test_jit_interpolation.py -q
```

Expected: all tests pass and the cylinder comparison stays within `0.1%`.

- [ ] **Step 6: Commit**

```powershell
git add liqlev/geometry/fixtures.py liqlev/geometry/jit.py tests/geometry/test_boundary_layer.py
git commit -m "feat: generalize boundary layer by local perimeter"
```

---

### Task 5: Integrate Custom Geometry into the JIT Solver

**Files:**
- Modify: `core.py`
- Create: `tests/geometry/test_core_custom_geometry.py`

**Interfaces:**
- Consumes: solver input key `GeometryMode` and the seven numeric geometry arrays.
- Produces: the existing 29-column DataFrame with custom `Height`, `dh/dt`, `eps`, `VBL vol`, `BL thick`, and `BL Vap Out`; cylinder inputs remain unchanged.

- [ ] **Step 1: Write a failing cylinder-equivalence integration test**

Create `tests/geometry/test_core_custom_geometry.py`:

```python
from __future__ import annotations

import numpy as np
import pandas as pd

from core import liqlev_simulation
from liqlev.geometry.fixtures import cylinder_kernel
from validation.physics_cases import build_case_inputs, get_case


def attach_geometry(
    inputs: dict[str, object],
    diameter: float,
    height: float,
    fill_fraction: float,
) -> None:
    kernel = cylinder_kernel(diameter, height, node_count=1025)
    inputs.update(
        {
            "GeometryMode": 1,
            "GeomHeight": kernel.height_ft,
            "GeomVolume": kernel.volume_ft3,
            "GeomVolumeCoefficients": kernel.volume_coefficients,
            "GeomArea": kernel.section_area_ft2,
            "GeomPerimeter": kernel.perimeter_ft,
            "GeomSidewallArea": kernel.sidewall_area_ft2,
            "GeomSidewallCoefficients": kernel.sidewall_coefficients,
            "FillFraction": fill_fraction,
        }
    )


def test_custom_cylinder_matches_legacy_solver() -> None:
    case = get_case("hydrogen_height_dep_mid_fill")
    legacy_inputs = build_case_inputs(case)
    custom_inputs = build_case_inputs(case)
    attach_geometry(
        custom_inputs,
        case.dtank_ft,
        case.htank_ft,
        case.fill_fraction,
    )
    custom_inputs["Dtank"] = 0.75 * case.dtank_ft
    legacy = liqlev_simulation(legacy_inputs, verbose=False)
    custom = liqlev_simulation(custom_inputs, verbose=False)
    columns = ["Height", "eps", "VBL vol", "BL thick", "BL Vap Out"]
    pd.testing.assert_frame_equal(
        custom[columns],
        legacy[columns],
        rtol=1e-3,
        atol=1e-10,
        check_exact=False,
    )
    assert np.all(custom["Conv Failed"].to_numpy() == 0.0)
```

- [ ] **Step 2: Run the test and confirm custom inputs are not yet used**

Run:

```powershell
python -m pytest tests/geometry/test_core_custom_geometry.py -q
```

Expected: the equivalence assertion fails because `core.py` still uses the
deliberately inconsistent `Dtank` instead of the attached geometry arrays.

- [ ] **Step 3: Add geometry arguments to `_solver_loop`**

Append these arguments after the gravity arguments:

```python
geometry_mode,
geom_height,
geom_volume,
geom_volume_coefficients,
geom_area_samples,
geom_perimeter,
geom_sidewall_area,
geom_sidewall_coefficients,
```

Import the Numba functions from `liqlev.geometry.jit`. The wrapper supplies
two-element contiguous dummy arrays and `(4, 1)` zero coefficients when
`GeometryMode` is absent, preserving one compiled signature for legacy calls.

- [ ] **Step 4: Branch the geometry-dependent calculations**

Inside `_solver_loop`, make these exact substitutions only when
`geometry_mode == 1`:

```python
interface_area = max(
    0.0,
    eval_ppoly_derivative(h1, geom_height, geom_volume_coefficients),
)
wetted_side_area = max(
    0.0,
    eval_ppoly(h1, geom_height, geom_sidewall_coefficients),
)
eps_denom = wetted_side_area + interface_area
eps = wetted_side_area / eps_denom if eps_denom != 0.0 else 0.0
```

In each outer `AK3` trial, call `integrate_boundary_layer` at the predicted
height. Use:

```python
vbl2 = custom_vbl
custom_exit_rate = ak1 * normalized_exit_flow
fvbl = (
    vbl2
    - ak2 * xml1 * delta / rhol
    + custom_exit_rate * delta
    - vbl1
)
```

After the boundary-layer solve, update custom height from occupied volume:

```python
occupied_volume = xml2 / rhol + vbl2
h2 = invert_monotone_volume(
    occupied_volume,
    geom_height,
    geom_volume,
    geom_volume_coefficients,
)
dhdt = (h2 - h1) / delta
zht2 = h2
xmdtbl = custom_exit_rate * delta * rhov
```

If inverse-volume or the boundary-layer integrator fails, set
`Conv Failed=1.0`, retain the bounded prior height, and stop the timestep loop
after writing the diagnostic row. Leave every legacy expression in its
existing branch.

- [ ] **Step 5: Initialize custom fill by volume in `liqlev_simulation`**

When `geometry_mode == 1`, require `Volt` to match `GeomVolume[-1]` within
`1e-10` relative and set:

```python
initial_occupied_volume = inputs["fillfraction"] * volt
htzero = invert_monotone_volume(
    initial_occupied_volume,
    geom_height,
    geom_volume,
    geom_volume_coefficients,
)
```

The wrapper accepts explicit `FillFraction`; legacy callers continue deriving
fill from `Htzero * Ac / Volt`.

- [ ] **Step 6: Run focused tests and the strict physics baseline**

Run:

```powershell
python -m pytest tests/geometry/test_core_custom_geometry.py -q
python scripts/check_physics_baseline.py
```

Expected: custom-cylinder equivalence passes within `0.1%`, and the legacy
physics baseline passes at `rtol=1e-9`, `atol=1e-8`.

- [ ] **Step 7: Commit**

```powershell
git add core.py tests/geometry/test_core_custom_geometry.py
git commit -m "feat: add custom geometry branch to JIT solver"
```

---

### Task 6: Connect Geometry Packages to Configuration and Headless Runs

**Files:**
- Modify: `liqlev/model/config.py`
- Modify: `liqlev/model/builder.py`
- Modify: `liqlev/model/validation.py`
- Modify: `liqlev/io/config_json.py`
- Modify: `liqlev/runner/single.py`
- Modify: `liqlev/runner/sweep.py`
- Modify: `tests/test_model_helpers.py`
- Modify: `tests/test_runner.py`

**Interfaces:**
- Consumes: `TankConfig.geometry_path`.
- Produces: a validated `SimulationConfig` that loads one geometry package and builds the exact custom arrays expected by `core.liqlev_simulation`.

- [ ] **Step 1: Write failing configuration and runner tests**

Add tests that:

```python
def test_geometry_path_round_trips_in_schema_v2(tmp_path) -> None:
    config = SimulationConfig(
        tank=TankConfig(
            diameter_ft=4.0,
            height_ft=8.0,
            fill_fractions=(0.5,),
            geometry_path="tank.npz",
        )
    )
    path = tmp_path / "config.json"
    save_simulation_config(config, path)
    assert load_simulation_config(path).tank.geometry_path == "tank.npz"


def test_geometry_path_must_exist() -> None:
    config = SimulationConfig(
        tank=TankConfig(geometry_path="missing.npz")
    )
    with pytest.raises(InputValidationError) as exc_info:
        validate_simulation_config(config)
    assert "tank.geometry_path" in {
        issue.field for issue in exc_info.value.issues
    }
```

Add a runner test that saves `cylinder_kernel(4.0, 8.0)` to `tmp_path`, points
`TankConfig.geometry_path` to it, runs one short case, and asserts
`result.inputs["GeometryMode"] == 1` and `result.htank_ft == 8.0`.

- [ ] **Step 2: Run the focused tests and observe constructor failures**

Run:

```powershell
python -m pytest tests/test_model_helpers.py tests/test_runner.py -q
```

Expected: tests fail because `TankConfig` has no `geometry_path`.

- [ ] **Step 3: Add the versioned configuration field**

Add:

```python
@dataclass(frozen=True)
class TankConfig:
    diameter_ft: float = 21.670
    height_ft: float = 28.18
    fill_fractions: tuple[float, ...] = (0.5116,)
    geometry_path: str = ""
```

Set `CONFIG_SCHEMA_VERSION = 2`. `simulation_config_from_dict` accepts versions
1 and 2; version 1 receives `geometry_path=""`, while save always writes
version 2.

- [ ] **Step 4: Load and attach the package in the runner**

In `run_single_case`, call `load_geometry_package` once when
`config.tank.geometry_path` is non-empty and pass the returned
`GeometryKernel | None` to `build_inputs`. In `build_inputs`, custom mode sets:

```python
volt = geometry.total_volume_ft3
htzero = float(
    invert_monotone_volume(
        fill_fraction * volt,
        geometry.height_ft,
        geometry.volume_ft3,
        geometry.volume_coefficients,
    )
)
xmlzro = rhol * fill_fraction * volt
```

It then adds the exact eight `Geometry*`/`Geom*` keys from Task 5 plus
`FillFraction`. Legacy input construction remains byte-for-byte equivalent for
the existing builder test.

- [ ] **Step 5: Validate files and report actual tank height**

When `geometry_path` is non-empty, validation requires the `.npz` and adjacent
`.json` to exist and calls `load_geometry_package`; surface its
`GeometryPackageError` as `tank.geometry_path`. `SingleCaseResult.htank_ft` and
the sweep scenario's `htank` use `geometry.total_height_ft`.

- [ ] **Step 6: Run configuration, runner, and physics tests**

Run:

```powershell
python -m pytest tests/test_model_helpers.py tests/test_runner.py tests/geometry -q
python scripts/check_physics_baseline.py
```

Expected: all tests and the physics baseline pass.

- [ ] **Step 7: Commit**

```powershell
git add liqlev/model liqlev/io/config_json.py liqlev/runner tests/test_model_helpers.py tests/test_runner.py
git commit -m "feat: load geometry packages in headless runs"
```

---

### Task 7: Load the STEP Assembly and Select the Tank Product

**Files:**
- Create: `requirements-geometry.txt`
- Create: `liqlev/cad/__init__.py`
- Create: `liqlev/cad/xcaf.py`
- Create: `tests/cad/test_xcaf_selection.py`
- Create: `geometry/source/PROVENANCE.json`

**Interfaces:**
- Consumes: a STEP assembly path, expected SHA-256, and product name `nhq01-m21a- 0202_short`.
- Produces: `load_named_product(path, product_name, expected_sha256) -> cq.Shape` with the resolved assembly placement preserved.

- [ ] **Step 1: Add the optional CAD dependencies**

Create `requirements-geometry.txt`:

```text
-r requirements.txt
cadquery>=2.7,<2.8
scipy>=1.16,<1.17
pytest>=8.4,<9
```

- [ ] **Step 2: Write the failing source-product test**

Create `tests/cad/test_xcaf_selection.py`:

```python
from __future__ import annotations

from pathlib import Path

import pytest

from liqlev.cad.xcaf import load_named_product, sha256_file


SOURCE = Path(
    r"C:\Users\sasorian\Documents\Eta_Space\geometry"
    r"\nhq01-m21a- 0201_TankAssy_NASA.STEP"
)
EXPECTED_HASH = "0193EC296B5B754FFDC7B1652052F5B28BA99345A5BC89922D026DF011C837D5"


@pytest.mark.skipif(not SOURCE.exists(), reason="NASA source STEP is not present")
def test_selects_resolved_tank_body() -> None:
    assert sha256_file(SOURCE) == EXPECTED_HASH
    shape = load_named_product(
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
```

- [ ] **Step 3: Run and observe the missing CAD module**

Run:

```powershell
python -m pytest tests/cad/test_xcaf_selection.py -q
```

Expected: collection fails because `liqlev.cad.xcaf` does not exist.

- [ ] **Step 4: Implement XCAF assembly traversal**

Create `liqlev/cad/xcaf.py` using `STEPCAFControl_Reader`,
`XCAFDoc_DocumentTool.ShapeTool`, and `TDataStd_Name`. Traverse free shapes and
component labels recursively, accumulate each `TopLoc_Location`, and match the
exact product name. Resolve references before applying the accumulated
location. Require exactly one match and return a CadQuery `Shape` wrapper.

Expose `sha256_file(path: str | Path) -> str` and
`load_named_product(path: str | Path, *, product_name: str,
expected_sha256: str) -> cq.Shape`.

Raise `StepProductError` for a hash mismatch, STEP transfer failure, zero
matches, or multiple matches. Include the matched label entry and product name
in the successful diagnostic object returned by the internal traversal.

- [ ] **Step 5: Record source provenance**

Create `geometry/source/PROVENANCE.json`:

```json
{
  "source_path": "C:\\Users\\sasorian\\Documents\\Eta_Space\\geometry\\nhq01-m21a- 0201_TankAssy_NASA.STEP",
  "source_sha256": "0193EC296B5B754FFDC7B1652052F5B28BA99345A5BC89922D026DF011C837D5",
  "source_size_bytes": 36844537,
  "tank_product": "nhq01-m21a- 0202_short",
  "height_axis": "+Y",
  "gravity_direction": "-Y",
  "closure_planes_mm": [-275.406791, 275.296791]
}
```

- [ ] **Step 6: Run the source-selection test**

Run:

```powershell
python -m pytest tests/cad/test_xcaf_selection.py -q
```

Expected: one test passes and confirms the exact tank bounds and face count.

- [ ] **Step 7: Commit**

```powershell
git add requirements-geometry.txt liqlev/cad tests/cad/test_xcaf_selection.py geometry/source/PROVENANCE.json
git commit -m "feat: select tank body from STEP assembly"
```

---

### Task 8: Construct and Round-Trip the Watertight Fluid Domain

**Files:**
- Create: `liqlev/cad/fluid_domain.py`
- Create: `liqlev/cad/audit.py`
- Create: `tests/cad/test_fluid_domain.py`

**Interfaces:**
- Consumes: the placed tank-body shape and the two exact closure-plane ordinates.
- Produces: `build_fluid_domain(tank_body, y_min_mm, y_max_mm) -> cq.Solid`, `write_step_round_trip(fluid, output_path, audit_path) -> tuple[cq.Solid, CadAudit]`, and a validated one-solid STEP.

- [ ] **Step 1: Write failing topology and round-trip tests**

Create tests that load the named tank product, call:

```python
fluid = build_fluid_domain(
    tank,
    y_min_mm=-275.406791,
    y_max_mm=275.296791,
    plane_tolerance_mm=1e-5,
)
```

and assert:

```python
assert isinstance(fluid, cq.Solid)
assert fluid.isValid()
assert fluid.Volume() > 0.0
assert len(fluid.Solids()) == 1
assert len(fluid.Shells()) == 1
assert fluid.BoundingBox().ymin == pytest.approx(-275.406791, abs=1e-5)
assert fluid.BoundingBox().ymax == pytest.approx(275.296791, abs=1e-5)
```

Export to `tmp_path`, re-import, and assert the relative volume difference is
at most `1e-8` and every bounding-box coordinate differs by at most `1e-5 mm`.

- [ ] **Step 2: Run and observe the missing construction functions**

Run:

```powershell
python -m pytest tests/cad/test_fluid_domain.py -q
```

Expected: collection fails for `build_fluid_domain`.

- [ ] **Step 3: Implement deterministic inner-loop and inner-wall selection**

In `liqlev/cad/fluid_domain.py`:

1. Find planar faces whose centroids lie at either closure ordinate within
   `1e-5 mm` and whose normals are parallel to `Y` within `1e-8`.
2. From each matching annular rim, choose the closed wire with the smaller
   planar enclosed area; this is the wet-side opening loop.
3. Require exactly one wet-side loop per closure plane.
4. Build an edge-to-face adjacency map for all tank faces.
5. Start from the non-planar face adjacent to each wet-side loop edge and
   traverse connected faces without crossing either rim face.
6. Intersect the two traversed face sets and require one connected inner face
   network whose only free boundary wires are the two selected loops.
7. Build one planar face from each loop, sew the inner faces and caps with
   tolerance `1e-5 mm`, and convert the closed shell to a positive-volume
   solid.
8. Reject any repair that moves a vertex by more than `1e-5 mm`.

Expose `build_fluid_domain(tank_body: cq.Shape, *, y_min_mm: float,
y_max_mm: float, plane_tolerance_mm: float = 1e-5) -> cq.Solid`.

All topology-count failures raise `FluidDomainError` with observed face, wire,
shell, and solid counts.

- [ ] **Step 4: Implement the CAD audit and STEP round trip**

`CadAudit` stores source/output hash, volume, area, bounding box, solid count,
shell count, face count, validity, cap-plane error, and round-trip differences.
`write_step_round_trip` exports AP214, re-imports it, performs every acceptance
check from the approved specification, writes JSON, and returns the re-imported
solid plus audit.

- [ ] **Step 5: Run topology and round-trip tests**

Run:

```powershell
python -m pytest tests/cad/test_fluid_domain.py -q
```

Expected: all tests pass; temporary output has one valid closed solid.

- [ ] **Step 6: Commit**

```powershell
git add liqlev/cad/fluid_domain.py liqlev/cad/audit.py tests/cad/test_fluid_domain.py
git commit -m "feat: construct capped tank fluid domain"
```

---

### Task 9: Extract Exact Geometry Tables from BRep Measurements

**Files:**
- Create: `liqlev/cad/measure.py`
- Create: `tests/cad/test_analytic_measurements.py`

**Interfaces:**
- Consumes: a valid fluid solid and ordered absolute `Y` coordinates.
- Produces: exact cumulative volume, section area, perimeter, sidewall area, and total wetted area samples plus a validated `GeometryKernel`.

- [ ] **Step 1: Write cylinder and sphere measurement tests**

Create analytic CadQuery solids aligned with `Y`, cap them at their endpoints,
sample 33 nodes, and compare to:

```python
def cylinder_expected(radius: float, h: np.ndarray):
    area = np.full_like(h, np.pi * radius**2)
    volume = area * h
    perimeter = np.full_like(h, 2.0 * np.pi * radius)
    sidewall = 2.0 * np.pi * radius * h
    return area, volume, perimeter, sidewall


def sphere_expected(radius: float, h: np.ndarray):
    area = np.pi * (2.0 * radius * h - h**2)
    volume = np.pi * h**2 * (radius - h / 3.0)
    perimeter = 2.0 * np.pi * np.sqrt(
        np.maximum(0.0, 2.0 * radius * h - h**2)
    )
    sidewall = 2.0 * np.pi * radius * h
    return area, volume, perimeter, sidewall
```

Require volume/area relative error below `1e-8` away from degenerate endpoint
sections and exact zero cumulative volume at the bottom.

- [ ] **Step 2: Run and observe the missing measurement module**

Run:

```powershell
python -m pytest tests/cad/test_analytic_measurements.py -q
```

Expected: collection fails because `liqlev.cad.measure` does not exist.

- [ ] **Step 3: Implement exact clipped-solid measurements**

For each interior `Y`:

1. Intersect the fluid solid with a box spanning its full `X/Z` bounds plus
   `1 mm`, from `Y_min-1 mm` through the requested `Y`.
2. Use the clipped solid's volume as cumulative `V`.
3. Identify the new planar `+Y` cut face by centroid and normal; use its exact
   face area as `A`.
4. Use the total length of that face's outer wire as `P`; reject an inner wire
   or a second outer wire.
5. Compute total wetted area as clipped-solid surface area minus the cut-face
   area.
6. Compute sidewall area as total wetted area minus the minimum-`Y` cap area.
7. At the full endpoint, subtract both cap areas from the fluid-solid surface
   area to obtain sidewall area.

Create the immutable `GeometrySamplesMM` dataclass with `float64` arrays
`height_mm`, `volume_mm3`, `section_area_mm2`, `perimeter_mm`,
`sidewall_area_mm2`, and `total_wetted_area_mm2`. Expose
`measure_geometry(fluid: cq.Solid, absolute_y_mm: np.ndarray, *, y_min_mm:
float, y_max_mm: float) -> GeometrySamplesMM`.

- [ ] **Step 4: Implement adaptive refinement and package conversion**

Start with 33 uniform nodes plus all detected face-vertex `Y` ordinates inside
the domain. Measure interval midpoints and subdivide whenever the midpoint
volume or sidewall area differs from monotone-cubic prediction by more than
`5e-4` relative to the corresponding total. Stop after 1025 nodes; failure to
meet tolerance at 1025 nodes raises `GeometryMeasurementError`.

Convert with:

```python
MM_TO_FT = 1.0 / 304.8
MM2_TO_FT2 = MM_TO_FT**2
MM3_TO_FT3 = MM_TO_FT**3
```

Build volume and sidewall coefficients, validate the `GeometryKernel`, and
independently integrate `dV/dh` over 1025 evaluation points. Reject volume
disagreement above `0.05%` and direct CAD area/derivative disagreement above
`0.2%` away from topology nodes.

- [ ] **Step 5: Run analytic geometry tests**

Run:

```powershell
python -m pytest tests/cad/test_analytic_measurements.py tests/geometry -q
```

Expected: cylinder and sphere exact-function tests pass, along with all numeric
kernel tests.

- [ ] **Step 6: Commit**

```powershell
git add liqlev/cad/measure.py tests/cad/test_analytic_measurements.py
git commit -m "feat: extract exact BRep geometry tables"
```

---

### Task 10: Generate and Independently Verify the NASA Tank Artifacts

**Files:**
- Create: `scripts/build_nasa_tank_geometry.py`
- Create: `scripts/check_nasa_tank_geometry.py`
- Create: `tests/cad/test_nasa_tank_artifacts.py`
- Generate: `geometry/output/nhq01-m21a-0201_LIQLEV_FLUID_DOMAIN.step`
- Generate: `geometry/tables/nhq01-m21a-0201_LIQLEV_GEOMETRY.npz`
- Generate: `geometry/tables/nhq01-m21a-0201_LIQLEV_GEOMETRY.json`
- Generate: `geometry/tables/nhq01-m21a-0201_LIQLEV_GEOMETRY.csv`
- Generate: `geometry/audit/nhq01-m21a-0201_LIQLEV_AUDIT.json`
- Generate: `geometry/audit/sections/*.png`

**Interfaces:**
- Consumes: only the immutable source STEP and fixed design constants.
- Produces: reproducible CAD, solver tables, audit report, images, and a command that independently verifies them.

- [ ] **Step 1: Write the failing artifact test**

The test loads the committed output STEP and geometry package, then asserts:

```python
assert audit["passed"] is True
assert audit["source_sha256"] == EXPECTED_SOURCE_HASH
assert audit["solid_count"] == 1
assert audit["shell_count"] == 1
assert audit["valid"] is True
assert audit["cap_plane_max_error_mm"] <= 1e-5
assert audit["round_trip_relative_volume_error"] <= 1e-8
assert kernel.metadata.axis == "+Y"
assert kernel.metadata.gravity_direction == "-Y"
assert kernel.total_volume_ft3 > 0.0
assert kernel.total_height_ft == pytest.approx(
    (275.296791 - (-275.406791)) / 304.8,
    abs=1e-10,
)
```

It also checks that all section images listed in the audit exist and have
nonzero size.

- [ ] **Step 2: Run and observe missing artifacts**

Run:

```powershell
python -m pytest tests/cad/test_nasa_tank_artifacts.py -q
```

Expected: failure listing the absent STEP, NPZ, JSON, CSV, audit, and images.

- [ ] **Step 3: Implement the deterministic build command**

`scripts/build_nasa_tank_geometry.py` accepts:

```text
--source-step
--output-root
--source-sha256
--product-name
--y-min-mm
--y-max-mm
--max-nodes
```

Defaults are the approved NASA path, repository `geometry` directory, fixed
source hash, `nhq01-m21a- 0202_short`, both approved closure planes, and 1025
nodes. The command runs product selection, fluid construction, STEP round
trip, adaptive measurements, coefficient construction, and package save in
that order. It exits nonzero before replacing validated outputs if any stage
fails.

Generate orthogonal X-Y, Y-Z, and X-Z images plus X-Z sections at 0%, 10%,
25%, 50%, 75%, 90%, and 100% height. Each image labels assembly axes, absolute
`Y`, normalized `h`, and scale in millimetres.

- [ ] **Step 4: Implement the independent check command**

`scripts/check_nasa_tank_geometry.py` re-hashes the source and artifacts,
re-imports the STEP, reruns BRep checks, loads the numeric package, recomputes
volume from `dV/dh`, and verifies every specification tolerance without
calling the build command. It prints one line per check followed by:

```text
NASA tank geometry audit passed.
```

and exits `0`; any failed check exits `1`.

- [ ] **Step 5: Build the production artifacts**

Run:

```powershell
python scripts/build_nasa_tank_geometry.py
```

Expected: all six artifact groups are written below `geometry/`; the source
STEP remains unchanged with the approved hash.

- [ ] **Step 6: Verify independently and inspect sections**

Run:

```powershell
python scripts/check_nasa_tank_geometry.py
python -m pytest tests/cad/test_nasa_tank_artifacts.py -q
```

Expected: audit command and artifact test pass. Visually inspect every image
listed in the audit and record `"visual_section_review": "passed"` in the
audit only after confirming there are no baffles, internal bodies, open ports,
or outer-wall regions in the fluid solid.

- [ ] **Step 7: Commit**

```powershell
git add scripts/build_nasa_tank_geometry.py scripts/check_nasa_tank_geometry.py tests/cad/test_nasa_tank_artifacts.py geometry
git commit -m "feat: publish validated NASA tank geometry artifacts"
```

---

### Task 11: Run the NASA Tank Through LIQLEV and Check Refinement

**Files:**
- Create: `validation/custom_geometry_cases.py`
- Create: `tests/geometry/test_nasa_tank_solver.py`
- Create: `validation/results/nasa_tank_geometry_manifest.json`

**Interfaces:**
- Consumes: the production geometry package and a deterministic hydrogen validation case.
- Produces: solver convergence, physical-bound, refinement, and provenance evidence.

- [ ] **Step 1: Define the deterministic tank case**

Create a `SimulationConfig` using the production NPZ, fill fractions
`(0.10, 0.25, 0.50, 0.75, 0.90)`, the legacy hydrogen pressures and
temperature from `hydrogen_height_dep_mid_fill`, `height_dep` epsilon,
`0.001 g`, `10 s` timestep, and `300 s` duration.

- [ ] **Step 2: Write failing tank solver assertions**

For each fill, assert:

```python
assert not dataframe.empty
assert np.isfinite(dataframe.to_numpy(dtype=float)).all()
assert (dataframe["Height"] >= 0.0).all()
assert (dataframe["Height"] <= kernel.total_height_ft).all()
assert (dataframe["VBL vol"] >= 0.0).all()
assert (dataframe["BL thick"] >= 0.0).all()
assert dataframe["Conv Failed"].sum() == 0.0
```

Rebuild an in-memory 513-node package from the same exact sample authority and
require height and boundary-layer volume to differ by no more than `0.2%`
relative from the 1025-point evaluation.

- [ ] **Step 3: Run and classify any numerical failure**

Run:

```powershell
python -m pytest tests/geometry/test_nasa_tank_solver.py -q
```

Expected: all five fills and the refinement comparison pass. If a failure
occurs, retain the failing result manifest and adjust only integration
resolution or safeguards; do not change the governing equation or tolerance
without updating the approved specification.

- [ ] **Step 4: Write the result manifest**

Record geometry NPZ hash, fluid STEP hash, solver commit, Python/NumPy/Numba
versions, fill cases, maximum convergence count, maximum refinement
difference, and pass/fail in
`validation/results/nasa_tank_geometry_manifest.json`.

- [ ] **Step 5: Run all physics and geometry tests**

Run:

```powershell
python scripts/check_physics_baseline.py
python scripts/check_nasa_tank_geometry.py
python -m pytest tests -q
```

Expected: both independent checks pass and Pytest reports zero failures.

- [ ] **Step 6: Commit**

```powershell
git add validation/custom_geometry_cases.py validation/results/nasa_tank_geometry_manifest.json tests/geometry/test_nasa_tank_solver.py
git commit -m "test: validate LIQLEV with NASA tank geometry"
```

---

### Task 12: Document Operation and Perform Final Verification

**Files:**
- Create: `docs/geometry-kernel.md`
- Modify: `README.md`
- Modify: `.gitignore`

**Interfaces:**
- Consumes: completed commands, artifacts, configuration path, and validation results.
- Produces: reproducible operator instructions and a clean, auditable repository.

- [ ] **Step 1: Document the complete workflow**

`docs/geometry-kernel.md` must include:

- the supported phase-one geometry class and explicit unsupported cases;
- source STEP provenance and immutable hash;
- `+Y` height and `-Y` gravity convention;
- cap-plane locations and why caps are excluded from active sidewall area;
- commands to create the environment, build artifacts, verify artifacts, run
  cylinder equivalence, and run the NASA tank case;
- NPZ array names, shapes, units, and interpolation convention;
- configuration JSON example with `tank.geometry_path`;
- legacy-cylinder selection by leaving `geometry_path` empty;
- all acceptance tolerances and the distinction between numerical consistency
  and experimental validation.

- [ ] **Step 2: Add concise README entry points**

Add links to the approved specification, implementation plan, geometry guide,
fluid STEP, table metadata, CAD audit, and solver result manifest. Include:

```json
{
  "tank": {
    "diameter_ft": 1.0,
    "height_ft": 1.806769625984252,
    "fill_fractions": [0.5],
    "geometry_path": "geometry/tables/nhq01-m21a-0201_LIQLEV_GEOMETRY.npz"
  }
}
```

State that diameter and height fields remain for schema compatibility and are
not used by the custom-geometry solver.

- [ ] **Step 3: Ensure generated scratch state is ignored**

Keep validated STEP, NPZ, JSON, CSV, and audit images tracked. Ignore only:

```gitignore
__pycache__/
.pytest_cache/
.ruff_cache/
.numba_cache/
*.py[cod]
build/
dist/
geometry/work/
geometry/tmp/
```

- [ ] **Step 4: Run the final verification matrix**

Run:

```powershell
git diff --check
python scripts/check_physics_baseline.py
python scripts/check_nasa_tank_geometry.py
python -m pytest tests -q
python -m pytest tests/geometry/test_core_custom_geometry.py tests/geometry/test_nasa_tank_solver.py -q
git status --short
```

Expected:

- no whitespace errors;
- legacy physics baseline passes;
- NASA CAD/table audit passes;
- all tests pass;
- custom cylinder and NASA tank solver tests pass;
- `git status --short` lists only the intended documentation changes.

- [ ] **Step 5: Commit documentation**

```powershell
git add README.md .gitignore docs/geometry-kernel.md
git commit -m "docs: explain geometry kernel workflow"
```

- [ ] **Step 6: Verify authorship and final tree**

Run:

```powershell
git log --format="%h %an <%ae> %s"
git status --short --branch
```

Expected: every commit is authored by Steven Soriano using the configured
email, there are no assistant authors or co-authors, branch is `main`, and the
working tree is clean.
