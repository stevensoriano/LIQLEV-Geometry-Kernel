# Legacy Cylinder Series and Report-Based Custom Geometry

## Purpose and authority

LIQLEV has two deliberately separate physics authorities:

1. **Legacy cylinder mode** is governed by the preserved implementation in
   `core.py` and the saved physics baseline. Its coefficients are unchanged so
   existing results remain reproducible.
2. **Custom geometry mode** is governed by NASA report NTRS 19700017832,
   printed pages 4-22 through 4-24 (PDF pages 72-74), especially Eq. 4-33,
   4-37, and 4-38.

The local report used for traceability is:

`C:\Users\sasorian\Documents\Cryo Vent LLR\tmp\pdfs\19700017832.pdf`

The official NTRS source is:

`https://ntrs.nasa.gov/api/citations/19700017832/downloads/19700017832.pdf`

Its SHA-256 is:

`F01404EE72A1EE8CE8CAAAEDD0812DB8A79ADF731452F7FDF7E35F12CA4313BB`

The Cryo Vent LLR file is only the local analysis copy. All authored source,
tests, documentation, and generated project artifacts remain under
`C:\Users\sasorian\Documents\Eta_Space\LIQLEV-Geometry-Kernel`.

The published equations take precedence for custom mode when they disagree
with the legacy transcription. This split prevents a correction to new
geometry physics from silently changing established legacy results.

## Published cylinder derivation

For a cylinder of diameter `D`, area `A=pi*D^2/4`, perimeter `P=pi*D`,
boundary-layer thickness `delta`, and report constant `K3`, Eq. 4-33 is

```text
delta^(1/2) d(delta) / (D/4 - delta) = K3 dZ.
```

Expanding `1/(D/4-delta)` as a geometric series and integrating gives

```text
Z = (8/K3) * sum(n=1..infinity) [
    4^(n-1) * delta^(n+1/2) / ((2*n+1) * D^n)
].
```

The height-series denominators are therefore `3, 5, 7, 9, ...`.

Eq. 4-37 defines boundary-layer volume:

```text
VBL = pi*D * integral(0..h) delta(Z) dZ.
```

Substituting Eq. 4-33 and integrating in `delta` gives Eq. 4-38:

```text
VBL = sum(n=1..infinity) [
    (8*pi/K3)
    * 4^(n-1)
    / ((2*n+3) * D^(n-1))
    * delta^(n+3/2)
].
```

The report-correct normalized exit state is

```text
q = (2/3) * pi*D * delta^(3/2).
```

In custom mode this cylinder result generalizes without an equivalent
diameter:

```text
q = (2/3) * P(h) * delta^(3/2)
delta = (1.5 * max(q, 0) / P(h))^(2/3)
dq/dh = K3 * (A(h) - P(h)*delta)
dVBL/dh = P(h)*delta.
```

Here `A(h)=dV/dh` is derived from the cumulative-volume PCHIP and `P(h)` is
the nonnegative interpolated contact-line perimeter.

## Documented legacy discrepancy

The legacy height coefficients use

```text
4^(n-1) / (D^n * (2^n + 1)),
```

so their denominators are `3, 5, 9, 17, ...`, rather than the report's
`3, 5, 7, 9, ...`. The discrepancy begins at `n=3` and is caused by the
preserved `(2**L + 1)` expression.

The legacy boundary-layer volume coefficient is

```text
(pi/K3)
* (2*n+1)/(n+3/2)
* delta^(n+3/2)/D^(n-1),
```

equivalently

```text
(2*pi/K3)
* (2*n+1)/(2*n+3)
* delta^(n+3/2)/D^(n-1).
```

That is not Eq. 4-38. At `n=1`, the report coefficient is
`8*pi/(5*K3)` while the legacy coefficient is `6*pi/(5*K3)`, so the report
leading term is `4/3` of the legacy leading term.

Legacy exit flow uses

```text
2.1 * AK1 * D * delta^(3/2).
```

The number `2.1` is the rounded form of `2*pi/3`. It remains unchanged in
legacy mode. Custom mode retains the exact normalized state
`q=(2/3)*P*delta^(3/2)` and applies `AK1` outside the geometry integration.

For `D=4 ft`, tank height `8 ft`, and `K3=0.015 ft^(-1/2)`, the two analytic
series give the following values. The final three columns are
`100*(report-legacy)/legacy`.

| Fill | Report delta | Legacy delta | Report VBL | Legacy VBL | Report q | Legacy exit/AK1 | delta diff | VBL diff | exit diff |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.10 | 0.0668198927 | 0.0668396960 | 0.4062358108 | 0.2961674956 | 0.1447029102 | 0.1451546606 | -0.029628% | 37.164212% | -0.311220% |
| 0.25 | 0.1202618022 | 0.1203816565 | 1.8400908634 | 1.3104978031 | 0.3493897555 | 0.3508486089 | -0.099562% | 40.411595% | -0.415807% |
| 0.50 | 0.1852913972 | 0.1857513486 | 5.7193268083 | 3.9550009190 | 0.6681923347 | 0.6724767123 | -0.247617% | 44.609999% | -0.637104% |
| 0.80 | 0.2460556872 | 0.2471854662 | 12.2572577340 | 8.2368800420 | 1.0225127130 | 1.0323184371 | -0.457057% | 48.809472% | -0.949874% |
| 0.95 | 0.2722472500 | 0.2738109021 | 16.1679968417 | 10.7277844404 | 1.1900462974 | 1.2035257601 | -0.571070% | 50.711425% | -1.119998% |

## Coordinate, units, and signs

Solver height `h=Y-Y_min` increases in the assembly `+Y` direction; supported
gravity points in `-Y`. Geometry packages explicitly carry this axis and
gravity metadata.

Runtime geometry arrays use feet, square feet, and cubic feet:

- `h`, `D`, and `delta`: ft;
- `A`: ft^2;
- `P`: ft;
- `VBL`: ft^3;
- normalized `q`: ft^(5/2);
- `K3`: ft^(-1/2), so `dq/dh` is dimensionally consistent.

The JIT kernel performs no hidden unit conversion.

## Numerical implementation

OpenCascade shapes and Python objects stop at the offline preprocessing
boundary. The compiled runtime receives only contiguous numeric arrays and
scalars.

The cumulative-volume table is represented by a monotone cubic PCHIP;
`A(h)` is its analytic piecewise-polynomial derivative. Perimeter is linearly
interpolated and bounded below by zero. The normalized `q` and `VBL` states
are integrated interval by interval with fixed RK4 substeps. RK4 stage
derivatives may be negative; only completed `q` and `VBL` states are clamped
to zero.

The integrator returns status `1` for an out-of-domain or nonfinite height,
nonpositive/nonfinite `K3`, invalid substep count, an unusable zero perimeter
with positive flow state, or a nonfinite integrated state. Otherwise it
returns status `0` together with `delta_top`, `VBL`, and `q_top`.

## Assumptions, validation, and limitations

The custom formulation assumes a single-valued fill height, one usable total
contact-line perimeter per section, a monotone cumulative-volume table, the
report's quasi-one-dimensional boundary-layer balance, and adequate table and
RK4 resolution. It does not substitute an equivalent diameter, represent
internal baffles, or resolve multidimensional circulation and local thermal
effects.

Acceptance separates reproducibility from model validation:

- the unchanged saved baseline proves legacy-cylinder regression stability;
- analytic cylinder and sphere arrays verify geometry definitions;
- the custom cylinder boundary layer must match the report-correct analytic
  solution for `delta`, `VBL`, and normalized exit flow within `0.1%`;
- custom-tank refinement checks establish numerical consistency only.

The report-cylinder comparison is not a claim of experimental validity for an
arbitrary tank. Custom-tank predictions remain report-based numerical
predictions and require later experimental correlation or higher-fidelity
analysis.
