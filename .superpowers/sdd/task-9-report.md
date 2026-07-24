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

## Review correction: midpoint derivative refinement

Task 9 review identified that direct CAD area was compared with PCHIP
`dV/dh` only at table nodes after adaptive refinement had already stopped.
A valid Y-aligned sphere passed the cumulative-volume midpoint tolerance with
33 nodes, then failed terminally because endpoint-conditioned PCHIP nodal
slopes are not the analytic sphere derivatives.

The requested sphere kernel regression was written before the correction and
produced the expected RED:

```text
python -m pytest \
  tests/cad/test_analytic_measurements.py::\
test_builds_validated_sphere_geometry_kernel -q

GeometryMeasurementError:
direct CAD area and dV/dh disagree by more than 0.2% away from topology nodes

1 failed in 2.81s
```

The adaptive pass now compares exact CAD section area with PCHIP `dV/dh` at
the already measured interval midpoints. Midpoint disagreement above `0.2%`
drives subdivision together with the original cumulative-volume and
sidewall-area criteria. Degenerate endpoint and topology neighborhoods use
the fixed spacing of the initial 33-node grid; these neighborhoods are not
direct-area authorities, but are subdivided alongside adjacent area failures
so endpoint-conditioned slopes converge rather than contaminating the first
eligible midpoint. The final validator repeats the same midpoint criterion
from the accepted exact samples. Failure to converge within `max_nodes`
continues to raise `GeometryMeasurementError`.

The sphere regression independently checks:

- a validated `GeometryKernel`;
- more than 33 and no more than 1025 refined nodes;
- exact sphere volume, direct area, perimeter, and sidewall samples;
- 1025-point integrated `dV/dh` volume error no greater than `0.05%`; and
- midpoint direct-area/derivative disagreement below `0.2%` away from the
  fixed topology neighborhoods.

Post-correction evidence:

```text
python -m pytest \
  tests/cad/test_analytic_measurements.py::\
test_builds_validated_sphere_geometry_kernel -q

1 passed in 5.83s

python -m pytest tests/cad/test_analytic_measurements.py -q

10 passed in 7.23s

python -m pytest tests/cad/test_analytic_measurements.py tests/geometry -q

62 passed in 20.95s
```

Compared with the pre-correction gates, focused analytic runtime increased
from `3.26s` to `7.23s` and the combined analytic/geometry runtime increased
from `16.62s` to `20.95s`. The added time is the expected offline cost of exact
B-Rep midpoint sections for derivative-driven sphere refinement; runtime
solver behavior and the physics baseline are unchanged.

## Second review correction: non-vacuous dense-topology validation

Re-review identified that applying the initial-grid spacing as an exclusion
radius around every vertex ordinate could make the `0.2%` gate vacuous. A
deterministic valid `32 mm` zigzag prism was added with:

- 17 topology ordinates at `Y = 0, 2, ..., 32 mm`;
- alternating `8 mm` and `12 mm` half-widths; and
- a `5 mm` extrusion depth.

Before the correction, the overlapping `1 mm` exclusion neighborhoods covered
every interval midpoint. The builder certified the initial 33 nodes even
though the independently evaluated midpoint area/derivative disagreement was
`2.5%`.

The regression was written first and produced the expected RED:

```text
python -m pytest \
  tests/cad/test_analytic_measurements.py::\
test_dense_topology_cannot_make_area_validation_vacuous -q

assert len(kernel.height_ft) > 33
E assert 33 > 33

1 failed in 3.51s
```

The final area gate now:

- excludes an interval midpoint for topology only when it numerically
  coincides with an exact topology ordinate, using the scale-aware CAD
  tolerance;
- permits a fixed initial-grid neighborhood only for endpoints whose exact
  section area is numerically zero;
- requires at least one eligible positive-area midpoint in both adaptive and
  final validation;
- evaluates every other nearby vertex midpoint as authoritative; and
- uses topology neighborhoods only to drive symmetric subdivision when an
  authoritative midpoint fails, never to exclude those midpoints.

This preserves the sphere's necessary treatment at its two true zero-area
endpoints without allowing dense interior topology to suppress validation.
The `0.2%` area tolerance, `0.05%` integrated-volume tolerance, 1025-point
integration check, and 1025-node maximum remain unchanged.

Final exact metrics:

```text
sphere:
  nodes:                         79
  maximum eligible area error:  0.171749616878%
  standalone build time:        3.685330 s

dense zigzag prism:
  nodes:                         513
  eligible midpoints:            512
  maximum eligible area error:  0.192307692315%
  standalone build time:        37.478564 s
```

Post-correction gates:

```text
python -m pytest \
  tests/cad/test_analytic_measurements.py::\
test_dense_topology_cannot_make_area_validation_vacuous -q

1 passed in 40.76s

python -m pytest \
  tests/cad/test_analytic_measurements.py::\
test_builds_validated_sphere_geometry_kernel -q

1 passed in 5.53s

python -m pytest tests/cad/test_analytic_measurements.py -q

11 passed in 44.81s

python -m pytest tests/cad/test_analytic_measurements.py tests/geometry -q

63 passed in 57.33s
```

Focused analytic runtime increased from `7.23s` to `44.81s`, and the combined
analytic/geometry gate increased from `20.95s` to `57.33s`. The approximately
37-second increase is the deterministic offline cost of exact B-Rep clipping
through 513 nodes for the deliberately dense worst-case regression. Runtime
solver code and legacy physics remain untouched.

## Integration correction: exact endpoint closure mosaics

Task 10 integration exposed that the approved Task 8 fluid solid retains its
exact planar closure mosaics: 27 faces at minimum `Y` and 3 faces at maximum
`Y`. Task 9 incorrectly routed endpoint sections through the interior
single-cut-face rule and rejected the valid solid before any interior
measurement.

A deterministic analytic regression was created from two exactly fused
adjacent `4 × 10 × 6 mm` prisms with refinement disabled. The result is one
valid `8 × 10 × 6 mm` solid whose minimum and maximum closure planes each
retain two coplanar faces. Interior clipping still produces the expected
single rectangular cut face.

The regression was written first and produced the expected RED:

```text
python -m pytest \
  tests/cad/test_analytic_measurements.py::\
test_accepts_exact_endpoint_face_mosaics -q

GeometryMeasurementError:
section at Y=0 mm has 2 outer faces; exactly one is supported

1 failed in 1.99s
```

Endpoint metrics now preserve the exact B-Rep and:

1. collect every planar assembly-`Y`-normal face at the closure ordinate;
2. sum the exact face areas;
3. group all face-wire edges with OpenCascade `isSame` topological identity;
4. remove edges that occur exactly twice on different mosaic faces;
5. reject non-manifold or otherwise ambiguous shared-edge multiplicity;
6. partition the remaining edges by exact shared-vertex identity;
7. assemble each component into an exact wire; and
8. require exactly one closed boundary, thereby rejecting inner, open, or
   ambiguous boundaries.

The single-face endpoint path remains unchanged. No face unification, shape
healing, surface fitting, replacement solid, or approximate geometry is
performed. Both summed closure areas are subtracted from full-solid surface
area when computing total sidewall area.

Analytic GREEN:

```text
minimum/maximum mosaic faces:  2 / 2
section area:                  48 mm^2
full volume:                   480 mm^3
section perimeter:             28 mm
full sidewall area:            280 mm^2
kernel nodes:                  33

python -m pytest \
  tests/cad/test_analytic_measurements.py::\
test_accepts_exact_endpoint_face_mosaics -q

1 passed in 3.25s
```

A read-only focused probe then independently imported the current exact STEP
and measured minimum, midpoint, and maximum nodes:

```text
solid / valid / faces:         1 / true / 40
minimum/maximum mosaic faces:  27 / 3
nodes measured:                3
elapsed:                       0.443089 s

minimum cap:
  area:                        33883.637475439 mm^2
  perimeter:                   652.529038755 mm

midpoint:
  volume:                      49210880.430889934 mm^3
  section area:                245720.639963535 mm^2
  perimeter:                   1757.218435306 mm

maximum cap:
  area:                        33980.358091778 mm^2
  perimeter:                   653.459694976 mm

full volume:                   98109377.711638018 mm^3
full sidewall area:            963804.466720872 mm^2
```

The nonzero midpoint volume, area, and perimeter prove measurement progressed
through endpoint setup and completed an exact interior clip.

Final gates:

```text
python -m pytest tests/cad/test_analytic_measurements.py -q

12 passed in 47.21s

python -m pytest tests/cad/test_analytic_measurements.py tests/geometry -q

64 passed in 59.17s
```

Relative to the preceding Task 9 correction, focused runtime changed from
`44.81s` to `47.21s` and combined runtime changed from `57.33s` to `59.17s`.
The Task 10 scripts, Task 10 artifact test, and all generated artifacts remain
unmodified and uncommitted by this correction.

## Re-review correction: touching endpoint boundary loops

Closure-mosaic re-review identified an ambiguity in the endpoint boundary
assembly. Two otherwise separate closed loops that touch at one vertex form
one connected edge component. `Wire.assembleEdges` can return one closed wire
for that component while silently omitting the other loop, causing the
endpoint perimeter to be under-measured.

A deterministic valid-solid regression was written first. Two unit bottom
prisms have square closure faces that touch at exactly one OpenCascade
topological vertex, and a solid upper bridge joins the volumes above the
endpoint. The result is one valid solid with two bottom faces and eight
boundary edges.

The regression produced the expected RED because the endpoint helper returned
area `2.0 mm^2` and perimeter `4.0 mm` instead of rejecting the ambiguous
boundary:

```text
python -m pytest \
  tests/cad/test_analytic_measurements.py::\
test_rejects_endpoint_boundary_loops_touching_at_one_vertex -q

Failed: DID NOT RAISE <class 'liqlev.cad.measure.GeometryMeasurementError'>
1 failed in 1.99s
```

Endpoint mosaic validation now enforces two independent exact-topology
invariants:

1. Before assembly, every boundary vertex must have degree exactly two.
   Vertices are grouped with OpenCascade `isSame` identity, so the shared
   vertex in the regression has degree four and is rejected.
2. After assembly, the accepted wire must contain exactly the same number of
   edges as the input boundary and must match every input and assembled edge
   exactly once in both directions using `isSame`. This rejects dropped,
   duplicated, or added edges even if the returned wire is closed.

The pre-existing one-closed-boundary rule remains in place, retaining
rejection of open, multiple, and inner boundaries. No tolerances changed.

Focused GREEN:

```text
python -m pytest \
  tests/cad/test_analytic_measurements.py::\
test_rejects_endpoint_boundary_loops_touching_at_one_vertex -q

1 passed in 1.94s

python -m pytest \
  tests/cad/test_analytic_measurements.py::\
test_accepts_exact_endpoint_face_mosaics -q

1 passed in 3.51s
```

A read-only probe confirmed that the current exact STEP remains one valid
40-face solid with 27 minimum and 3 maximum closure faces. Its measurements
remain unchanged:

```text
minimum cap area/perimeter:
  33883.63747543928 mm^2 / 652.5290387545534 mm
midpoint volume/area/perimeter:
  49210880.43089363 mm^3
  245720.63996353454 mm^2
  1757.2184353059488 mm
maximum cap area/perimeter:
  33980.3580917779 mm^2 / 653.459694976168 mm
full volume/sidewall area:
  98109377.71163802 mm^3 / 963804.466720872 mm^2
```

Final gates:

```text
python -m pytest tests/cad/test_analytic_measurements.py -q

13 passed in 46.17s

python -m pytest tests/cad/test_analytic_measurements.py tests/geometry -q

65 passed in 60.70s
```

The Task 10 scripts, Task 10 artifact test, and all generated artifacts remain
unmodified and uncommitted by this correction.
