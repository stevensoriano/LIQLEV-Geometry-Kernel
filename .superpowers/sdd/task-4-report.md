# Task 4 Report: Report-Based Perimeter-Aware Boundary Layer

## Status

Complete under the approved authority split:

- legacy cylinder mode remains governed by the unchanged `core.py` branch and
  saved physics baseline;
- custom geometry mode follows NASA NTRS 19700017832 Eq. 4-33, 4-37, and
  4-38;
- the custom numeric cylinder agrees with the report analytic cylinder within
  the required `rtol=1e-3` at every requested fill.

No equivalent diameter was introduced.

## Source authority and traceability

Official report:

`https://ntrs.nasa.gov/api/citations/19700017832/downloads/19700017832.pdf`

Local analysis copy:

`C:\Users\sasorian\Documents\Cryo Vent LLR\tmp\pdfs\19700017832.pdf`

SHA-256:

`F01404EE72A1EE8CE8CAAAEDD0812DB8A79ADF731452F7FDF7E35F12CA4313BB`

The relevant source is printed pages 4-22 through 4-24, PDF pages 72-74.
The Cryo Vent LLR PDF is only the analysis copy. All authored source, tests,
documentation, and project artifacts remain under
`C:\Users\sasorian\Documents\Eta_Space\LIQLEV-Geometry-Kernel`.

## TDD evidence

### Historical RED

`tests/geometry/test_boundary_layer.py` was created before the original
production files. The first run was:

```text
python -m pytest tests/geometry/test_boundary_layer.py -q
```

Collection failed as expected:

```text
ModuleNotFoundError: No module named 'liqlev.geometry.fixtures'
1 error during collection
```

This is the preserved RED evidence for the fixture and boundary-layer feature.

### Resumed fixture RED

The resumed work added exact total-wetted-area assertions before modifying the
fixtures. The focused run produced:

```text
2 failed, 14 passed
```

The failures showed that the cylinder omitted the top cap at the full endpoint
and the sphere incorrectly included its horizontal section area as wetted
wall. After the minimal fixture correction:

```text
16 passed in 6.48s
```

### Corrected-authority GREEN

The obsolete legacy-transcription acceptance oracle was replaced with the
published Eq. 4-33/Eq. 4-38 analytic cylinder. The required combined focused
run completed:

```text
26 passed in 6.63s
```

### Reviewer regression RED/GREEN

A reviewer found that Numba evaluates `max(0.0, np.nan)` as `0.0`. The
integrator previously clamped each raw RK4 candidate before checking the
completed state for finiteness, so nonfinite derivatives could be converted
to zero and returned with success status.

The regression injects `NaN` into the cumulative-volume derivative
coefficients and requires status `1` plus non-success (`NaN`) outputs. Before
the fix:

```text
python -m pytest \
  tests/geometry/test_boundary_layer.py::test_boundary_layer_rejects_nonfinite_rk4_candidate \
  -q
1 failed in 0.87s
E assert 0 == 1
```

The minimal correction computes raw `next_q` and `next_vbl`, checks both for
finiteness, and only then clamps finite completed states to zero. The
derivatives and RK4 formula are unchanged. After the fix:

```text
1 passed in 1.34s
```

## Governing equations

For a cylinder of diameter `D`, Eq. 4-33 is:

```text
delta^(1/2) d(delta) / (D/4 - delta) = K3 dZ.
```

Its integrated report height series is:

```text
Z = (8/K3) * sum(n=1..infinity) [
    4^(n-1) * delta^(n+1/2) / ((2*n+1) * D^n)
].
```

The denominators are `3, 5, 7, 9, ...`.

Eq. 4-37 is:

```text
VBL = pi*D * integral(0..h) delta(Z) dZ.
```

Eq. 4-38 is:

```text
VBL = sum(n=1..infinity) [
    (8*pi/K3)
    * 4^(n-1)
    / ((2*n+3) * D^(n-1))
    * delta^(n+3/2)
].
```

The report-correct normalized exit state is:

```text
q = (2/3) * pi*D * delta^(3/2).
```

Custom geometry uses the exactly corresponding local-perimeter state:

```text
q = (2/3) * P(h) * delta^(3/2)
delta = (1.5 * max(q, 0) / P(h))^(2/3)
dq/dh = K3 * (A(h) - P(h)*delta)
dVBL/dh = P(h)*delta.
```

`A(h)` is the derivative of cumulative-volume PCHIP and `P(h)` is the
nonnegative interpolated perimeter.

## Legacy/report discrepancy

The preserved legacy height series uses `(2**L + 1)`, giving denominators
`3, 5, 9, 17, ...`, rather than the report's `2*L+1` denominators
`3, 5, 7, 9, ...`. The first disagreement is at `L=3`.

The legacy VBL coefficient is:

```text
(pi/K3)
* (2*n+1)/(n+3/2)
* delta^(n+3/2)/D^(n-1).
```

The report coefficient is:

```text
(8*pi/K3)
* 4^(n-1)/(2*n+3)
* delta^(n+3/2)/D^(n-1).
```

At `n=1`, the report leading term is `8*pi/(5*K3)` and the legacy leading
term is `6*pi/(5*K3)`, a factor of `4/3`.

The legacy solver also keeps
`2.1*AK1*D*delta^(3/2)` unchanged; `2.1` is the rounded value of `2*pi/3`.
Custom mode uses the exact normalized `q` expression.

For `D=4 ft`, tank height `8 ft`, and `K3=0.015 ft^(-1/2)`, the analytic
results are:

| Fill | Report delta | Legacy delta | Report VBL | Legacy VBL | Report q | Legacy exit/AK1 | delta diff | VBL diff | exit diff |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.10 | 0.0668198927 | 0.0668396960 | 0.4062358108 | 0.2961674956 | 0.1447029102 | 0.1451546606 | -0.029628% | 37.164212% | -0.311220% |
| 0.25 | 0.1202618022 | 0.1203816565 | 1.8400908634 | 1.3104978031 | 0.3493897555 | 0.3508486089 | -0.099562% | 40.411595% | -0.415807% |
| 0.50 | 0.1852913972 | 0.1857513486 | 5.7193268083 | 3.9550009190 | 0.6681923347 | 0.6724767123 | -0.247617% | 44.609999% | -0.637104% |
| 0.80 | 0.2460556872 | 0.2471854662 | 12.2572577340 | 8.2368800420 | 1.0225127130 | 1.0323184371 | -0.457057% | 48.809472% | -0.949874% |
| 0.95 | 0.2722472500 | 0.2738109021 | 16.1679968417 | 10.7277844404 | 1.1900462974 | 1.2035257601 | -0.571070% | 50.711425% | -1.119998% |

The difference columns are `100*(report-legacy)/legacy`. These persistent
analytic differences explain why legacy output cannot be the custom-mode
acceptance reference.

## Report-cylinder tolerance evidence

Parameters:

```text
D = 4 ft
tank height = 8 ft
K3 = 0.015 ft^(-1/2)
node_count = 1025
RK4 substeps per table interval = 4
rtol = 1e-3
```

Absolute relative errors of the numeric custom integrator against the
report-correct analytic cylinder are:

| Fill | delta relative error | VBL relative error | q relative error |
|---:|---:|---:|---:|
| 0.10 | 2.650376e-08 | 9.440752e-07 | 3.975564e-08 |
| 0.25 | 1.034815e-08 | 1.964868e-07 | 1.552222e-08 |
| 0.50 | 5.010950e-09 | 5.854323e-08 | 7.516425e-09 |
| 0.80 | 3.030329e-09 | 2.527931e-08 | 4.545494e-09 |
| 0.95 | 2.513271e-09 | 1.849894e-08 | 3.769906e-09 |

The maximum observed relative error is `9.440752e-07`
(`0.00009440752%`), more than three orders of magnitude below the `0.1%`
limit. Error decreases rather than grows toward the high-fill end, and the
required 1025-node/four-substep discretization is sufficient.

## JIT boundary and safeguards

OpenCascade and Python objects are excluded from the runtime integration
boundary. The kernel receives only numeric scalars and arrays. The observed
nopython signature is:

```text
(float64, float64,
 Array(float64, 1, 'C'),
 Array(float64, 2, 'C'),
 Array(float64, 1, 'C'),
 int64) -> Tuple(float64, float64, float64, int64)
```

The integrator advances interval by interval with fixed RK4 substeps. Stage
derivatives are not clamped. Raw completed `q` and `VBL` candidates are
validated for finiteness before finite completed states are clamped to zero.
Invalid/nonfinite height, nonpositive/nonfinite `K3`, invalid substeps, or a
nonfinite state returns status `1`.

## Analytic fixtures

The tests verify exact arrays and metadata for:

- cylinder `A`, `V`, `P`, cumulative sidewall area, bottom cap, and the top
  cap added to total wetted area only at the full endpoint;
- sphere spherical-segment `A`, `V`, `P`, and
  `A_side=2*pi*R*h`, with no planar cap area added to total wetted area;
- `+Y` height, `-Y` gravity, and `ft`/`ft^2`/`ft^3` units.

## Legacy-baseline proof

The unchanged checker passed all saved cases:

```text
as203_default_high_vent         PASS rows=19
hydrogen_height_dep_mid_fill    PASS rows=29
nitrogen_custom_epsilon         PASS rows=31
Physics baseline check passed.
```

Neither `core.py`, the saved legacy baseline JSON, nor existing baseline
fixtures changed.

## Verification

Focused:

```text
python -m pytest tests/geometry/test_boundary_layer.py \
  tests/geometry/test_jit_interpolation.py -q
27 passed
```

Legacy baseline:

```text
python scripts/check_physics_baseline.py
3 cases passed
```

Full:

```text
python -m pytest -q
41 passed
```

Whitespace:

```text
git diff --check
exit 0
```

## Files changed

The reviewer fix changes `liqlev/geometry/jit.py`,
`tests/geometry/test_boundary_layer.py`, and this report. The complete Task 4
file set is:

- `liqlev/geometry/fixtures.py`
- `liqlev/geometry/jit.py`
- `tests/geometry/test_boundary_layer.py`
- `docs/physics/legacy-series-discrepancy.md`
- `docs/superpowers/specs/2026-07-23-liqlev-arbitrary-tank-geometry-design.md`
- `docs/superpowers/plans/2026-07-23-liqlev-geometry-kernel-implementation.md`
- `.superpowers/sdd/task-4-report.md`

## Self-review

- Confirmed the implementation uses local perimeter directly.
- Confirmed Eq. 4-33, 4-37, and 4-38 indexing and powers.
- Confirmed every requested fill is below the `0.1%` report-cylinder limit.
- Confirmed negative RK4 stage derivatives are not clamped.
- Confirmed raw RK4 candidates are rejected if nonfinite before clamping.
- Confirmed only finite completed `q` and `VBL` states are clamped.
- Confirmed invalid inputs and nonfinite states return status `1`.
- Confirmed normalized exit flow returns `q_top`.
- Confirmed analytic fixture arrays, axis, gravity, and unit metadata.
- Confirmed nopython compilation with no object fallback.
- Confirmed the unchanged legacy baseline passes.
- Confirmed `core.py`, baseline JSON, and baseline fixtures are untouched.
- Confirmed the design specification and implementation plan both encode the
  approved authority split and link the durable methodology.

## Commit

Original Task 4 commit:

```text
a2107769969eff1c8dd229d8857dd358b8ad0151
```

The report cannot contain the hash of a follow-up commit that contains the
updated report itself. The controller will record the reviewer-fix commit hash
in the handoff immediately after the commit.

Reviewer-fix subject:

```text
fix: reject nonfinite boundary layer states
```

Configured author:

```text
Steven Soriano <steven.a.soriano@nasa.gov>
```

No automated assistant attribution is included.

## Residual concerns

The custom-tank formulation assumes one usable total contact-line perimeter
per height, monotone cumulative volume, and adequate table/RK4 resolution. It
does not resolve baffles, disconnected liquid regions, multidimensional
circulation, or local thermal effects. Passing the analytic report cylinder
and table-refinement checks establishes numerical consistency, not
experimental validation. Custom-tank predictions require later experimental
correlation or higher-fidelity analysis.
