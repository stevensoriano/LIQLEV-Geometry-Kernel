# Task 9 Report: Exact B-Rep Geometry Tables

## Status

Complete. Task 9 adds exact OpenCascade-backed measurements, adaptive table
construction, unit conversion, coefficient construction, and validated
`GeometryKernel` output. NASA production tables were not generated and Task 10
was not started.

## Controller-resolved interface

The brief explicitly named only:

```python
measure_geometry(
    fluid: cq.Solid,
    absolute_y_mm: np.ndarray,
    *,
    y_min_mm: float,
    y_max_mm: float,
) -> GeometrySamplesMM
```

It also required adaptive refinement and a validated `GeometryKernel`, but did
not name that public interface or explain how mandatory provenance metadata
would be supplied. Implementation paused before production coding and the
controller approved:

```python
build_geometry_kernel(
    fluid: cq.Solid,
    *,
    metadata: GeometryMetadata,
    max_nodes: int = 1025,
) -> GeometryKernel
```

The builder derives the absolute interval from metadata, requires the `+Y`,
`-Y`, `ft`, `ft^2`, and `ft^3` contract, and rejects invalid node limits or
fluid/metadata bound disagreement with `GeometryMeasurementError`.

## TDD evidence

The first test file was created before production code. The required initial
RED was observed:

```text
python -m pytest tests/cad/test_analytic_measurements.py -q

ModuleNotFoundError: No module named 'liqlev.cad.measure'
```

After implementing exact clipped-solid measurements, the first analytic GREEN
was:

```text
2 passed in 2.00s
```

Kernel-builder tests were then added and produced the expected second RED:

```text
ImportError: cannot import name 'build_geometry_kernel'
```

The completed focused suite passed:

```text
9 passed in 3.17s
```

## Exact measurement behavior

`GeometrySamplesMM` is a frozen dataclass containing contiguous `float64`
arrays for height, cumulative volume, section area, perimeter, cumulative
sidewall area, and cumulative total-wetted area.

For every interior node, the implementation:

1. intersects the fluid solid with an axis-aligned box padded by `1 mm` in
   `X/Z` and extending from `Y_min - 1 mm` through the requested plane;
2. requires one valid connected clipped solid;
3. identifies exactly one planar outward `+Y` cut face by centroid and normal;
4. rejects multiple cut faces or any section face with an inner wire;
5. uses exact B-Rep volume, cut-face area, and outer-wire length;
6. removes the artificial cut face from clipped surface area;
7. removes the minimum cap from sidewall area; and
8. handles both full-solid caps separately at the maximum endpoint.

Cylinder and sphere samples were checked at 33 nodes against analytic area,
volume, perimeter, and sidewall functions. Volume begins at exact zero.

## Adaptive kernel construction

The builder starts with 33 uniform nodes and every interior face-vertex
assembly-`Y` ordinate. At each pass it measures all interval midpoints and
subdivides intervals whose cumulative-volume or sidewall-area PCHIP prediction
differs by more than `5e-4` of the corresponding exact total. It enforces a
hard maximum of 1025 nodes and raises `GeometryMeasurementError` when the
requested node budget cannot meet the tolerance.

Measurements are converted using:

```text
MM_TO_FT  = 1 / 304.8
MM2_TO_FT2 = MM_TO_FT^2
MM3_TO_FT3 = MM_TO_FT^3
```

Volume and sidewall PCHIP coefficients are built through the existing
coefficient interface. The result is passed through
`validate_geometry_kernel`. Additional gates check:

- final table volume against the exact CAD solid at `0.05%`;
- the 1025-point trapezoidal integral of runtime `dV/dh` against final volume
  at `0.05%`; and
- direct CAD area against runtime `dV/dh` away from topology nodes at `0.2%`.

## Verification

Required and regression gates:

```text
python -m pytest tests/cad/test_analytic_measurements.py -q
9 passed

python -m pytest tests/cad/test_analytic_measurements.py tests/geometry -q
61 passed in 16.62s

python scripts/check_physics_baseline.py
Physics baseline check passed for all three cases.

python -m pytest -q
127 passed in 136.40s
```

The final Task 8 inspection pair was regenerated once through the normal
production API after the full suite, then independently re-imported:

```text
STEP SHA-256:
056D6E4B52270324E7F842708ADABCA274269E1DDD6EABEA12F1ADA12E44193F

audit JSON SHA-256:
D4182A02A8D65F61C4B378EBD8ECE49AF64B9A02D5A511289C171E4F23AC2E9F

audit passed:                     true
audit output hash matches STEP:  true
re-import:                        1 valid solid / 1 shell / 40 faces
```

The STEP and audit remain intentionally untracked and are excluded from the
Task 9 commit.

## Scope and concerns

- The legacy cylinder branch, `core.py`, `thermo_utils.py`, and physics
  baseline data were not changed.
- No NASA production geometry tables were created.
- Task 10 remains responsible for provenance-complete NASA metadata, artifact
  generation, images, and independent production verification.
