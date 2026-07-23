# Task 5 Report: Custom Geometry in the JIT Solver

## Status

Complete under the approved authority split:

- legacy cylinder formulas and saved baseline behavior remain unchanged;
- `GeometryMode == 1` uses only contiguous numeric geometry arrays in the
  Numba solver;
- the custom numeric cylinder is checked against NASA NTRS 19700017832
  Eq. 4-33, 4-37, and 4-38 rather than against the legacy transcription;
- custom inversion and boundary-layer failures produce one diagnostic row,
  retain the prior bounded height, and exit the timestep loop.

No Task 6 configuration work, CAD work, OpenCascade import, or Python geometry
object was added to the runtime path.

## TDD evidence

The Task 5 test file was authored before `core.py` was modified.

The first broad focused run was:

```text
python -m pytest tests/geometry/test_core_custom_geometry.py -q
```

The unmodified solver ignored `GeometryMode`; an invalid-custom diagnostic
case consequently entered the legacy iteration path and the run timed out.
The same unmodified source was rerun with the smallest concrete authority
assertion first:

```text
python -m pytest tests/geometry/test_core_custom_geometry.py -vv -x
```

Pytest collected 12 tests. The strict legacy baseline passed, then the custom
cylinder assertion failed as expected:

```text
1 passed, 1 failed
Obtained eps: 0.5136180928684225
Expected eps: 0.8
```

This RED showed that the inconsistent legacy `Dtank` still controlled epsilon
and the attached custom arrays were not used.

The first GREEN attempt remained red because the new custom AK3 bracket
searched in the wrong direction:

```text
1 passed, 1 failed
Conv Failed: 1.0
```

A numerical residual scan established the cause: the outer iteration varies
trial `AK3` while physical `AK2` remains fixed. Therefore, a positive custom
residual must search toward smaller `AK3`. The one-line bracket-direction
correction retained the existing test and produced:

```text
python -m pytest tests/geometry/test_core_custom_geometry.py -q -x
12 passed in 15.95s
```

## Implementation

`liqlev_simulation` now:

- accepts `GeometryMode`, explicit `FillFraction`, and the seven required
  geometry arrays;
- validates `float64`, C-contiguity, dimensionality, finiteness, shapes,
  monotonic arrays, non-negative samples, and geometry endpoints before JIT;
- requires `Volt` and `GeomVolume[-1]` to agree within `1e-10` relative;
- initializes custom height by monotone volume inversion;
- supplies two-element numeric dummies and `(4, 1)` coefficient dummies for
  legacy mode, keeping one compiled signature.

The compiled loop branches only for custom mode to use:

- `dV/dh` for interface area;
- cumulative sidewall interpolation for wetted area and epsilon;
- the perimeter-aware Task 4 boundary-layer integration in every AK3 trial;
- the exact normalized custom exit state multiplied by `AK1`;
- occupied-volume inversion for height and `dh/dt`;
- explicit diagnostic termination for integration or inversion failure.

The legacy cylinder series, rounded `2.1` exit coefficient, constant-area
kinematics, and legacy epsilon calculations remain in the legacy branch.
The 29-column DataFrame contract is unchanged.

## NASA report-cylinder results

For the required fills, absolute relative errors against the tested
report-correct analytic reference were:

| Fill | BL thickness | VBL volume | Normalized exit |
|---:|---:|---:|---:|
| 0.10 | `2.033857e-12` | `9.920125e-07` | `3.046889e-12` |
| 0.25 | `8.455544e-13` | `2.154182e-07` | `1.264228e-12` |
| 0.50 | `3.615640e-13` | `6.785231e-08` | `5.383455e-13` |
| 0.80 | `1.842305e-13` | `3.100018e-08` | `2.724252e-13` |
| 0.95 | `1.433915e-13` | `2.327955e-08` | `2.111632e-13` |

The maximum observed error was `9.920125e-07`, below the required
`rtol=1e-3`. These comparisons do not require or imply equality with the
preserved legacy transcription.

## Independent-review correction

Independent review found that the wrapper still unconditionally calculated
legacy `perim`, `ac`, `htank`, and diameter-series coefficients before calling
the branched JIT loop. Consequently, otherwise-valid custom geometry with
`Dtank=0.0` failed before entering custom mode.

A regression was added first and run against committed head
`6fdb920b651ab248750aaad28222cc22fed30697`:

```text
python -m pytest \
  tests/geometry/test_core_custom_geometry.py::test_custom_cylinder_does_not_require_a_nonzero_legacy_dtank \
  -q
1 failed
ZeroDivisionError: float division by zero
core.py: htank = volt / ac
```

The minimal correction moved the existing legacy scalar and coefficient
formulas under `geometry_mode == 0`. Custom mode supplies finite neutral
scalars plus contiguous ten-element `float64` dummy arrays, preserving the
single Numba signature without evaluating any diameter-based formula.

The unchanged regression then passed:

```text
1 passed in 8.93s
```

The zero-diameter custom outputs match the positive-placeholder custom outputs
for `Height`, `dh/dt`, `eps`, `VBL vol`, `BL thick`, and `BL Vap Out` within
`rtol=1e-12`, `atol=1e-14`. A separate SI custom probe with `Dtank=0.0` also
completed successfully and matched the equivalent British custom outputs at
the same tolerances. The diameter input key remains present only until the
planned Task 6 configuration work; its value no longer controls or crashes
custom mode.

## Verification

Task 5 focused:

```text
python -m pytest tests/geometry/test_core_custom_geometry.py -q
13 passed
```

Task 3 and Task 4 focused regressions:

```text
python -m pytest tests/geometry/test_boundary_layer.py \
  tests/geometry/test_jit_interpolation.py -q
27 passed
```

Strict legacy baseline, using the unchanged checker defaults
`rtol=1e-9`, `atol=1e-8`:

```text
python scripts/check_physics_baseline.py
as203_default_high_vent         PASS rows=19
hydrogen_height_dep_mid_fill    PASS rows=29
nitrogen_custom_epsilon         PASS rows=31
Physics baseline check passed.
```

Full suite:

```text
python -m pytest -q
54 passed
```

Whitespace:

```text
git diff --check
exit 0
```

## Legacy-preservation evidence

The strict saved baseline passed all three cases without changing its checker,
tolerances, or data. `thermo_utils.py`, `validation/physics_cases.py`, the
saved baseline directory, `liqlev/geometry/jit.py`, and
`liqlev/geometry/fixtures.py` have no Task 5 diff. The only physics-critical
file changed is `core.py`, where all custom behavior is visibly guarded by
`geometry_mode == 1`.

## Numba evidence

A fresh process ran one legacy case and a custom `Dtank=0.0` case, then
inspected `_solver_loop`:

```text
nopython_signature_count=1
compiled_signature_count=1
contains_pyobject=False
```

The signature contains only primitive scalars, contiguous `float64` one- and
two-dimensional arrays, and returns
`Tuple(array(float64, 2d, C), int64)`. It contains no `pyobject`.

## Files changed

- `core.py`
- `tests/geometry/test_core_custom_geometry.py`
- `.superpowers/sdd/task-5-report.md`

## Self-review

- Confirmed all seven arrays are validated before entering Numba.
- Confirmed legacy and custom wrapper calls share one nopython signature.
- Confirmed deliberately inconsistent legacy diameters do not change custom
  height, kinematics, epsilon, or boundary-layer outputs.
- Confirmed custom mode does not evaluate legacy diameter-derived scalars or
  series coefficients and runs with `Dtank=0.0` in British and SI inputs.
- Confirmed the legacy formulas are unchanged inside `geometry_mode == 0`.
- Confirmed explicit fill uses volume inversion rather than legacy
  `Htzero * Ac / Volt`.
- Confirmed interface area is the cumulative-volume derivative and sidewall
  area excludes modeled caps.
- Confirmed all required fills use the tested NASA report analytic reference.
- Confirmed normalized exit flow uses the exact `(2/3) * P * delta^(3/2)`
  state from Task 4 and applies `AK1` outside the geometry integrator.
- Confirmed boundary-layer and volume-inversion failures write one diagnostic
  row with `Conv Failed=1.0`, retain prior height, and stop cleanly.
- Confirmed no OpenCascade symbols, Python geometry objects, or equivalent
  diameter entered the JIT path.
- Confirmed the existing 29 output columns remain unchanged.
- Confirmed no Task 6 or CAD files were started.

## Commit

Implementation commit:

```text
878b44fe211af4d5dfa68227cf0140b5b14966e2
feat: add custom geometry branch to JIT solver
```

Independent-review correction:

```text
f5a142d62b17e620b1f747670070b19af0bef796
fix: decouple custom geometry from legacy diameter
```

Configured author:

```text
Steven Soriano <steven.a.soriano@nasa.gov>
```

No automated assistant attribution is included. The documentation follow-up
commit containing this report cannot include its own hash; its exact hash is
recorded in the controller handoff.

## Residual concerns

The custom mode inherits the approved report-model limitations: one
single-valued fill height, one usable total perimeter per section, monotone
volume, and adequate table/RK4 resolution. Passing the report cylinder and
regression gates establishes numerical consistency, not experimental
validation for an arbitrary tank. Custom-tank predictions still require later
experimental correlation or higher-fidelity analysis.
