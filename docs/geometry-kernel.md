# LIQLEV Geometry Kernel

## Purpose and phase-one scope

The geometry kernel extends LIQLEV beyond the analytic constant-diameter
cylinder while leaving the legacy cylinder branch unchanged. OpenCascade is
used offline to construct and measure an exact internal fluid volume. The
runtime solver receives only versioned, contiguous `float64` arrays and
primitive scalars, so no CAD object or Python geometry object enters the
Numba-compiled loop.

Phase one supports a fixed tank orientation with:

- one watertight, connected wetted fluid volume;
- one outer section loop and no inner holes at every interior height;
- a single-valued horizontal free surface;
- monotone cumulative volume and wetted-area functions; and
- assembly `+Y` as increasing liquid height with gravity in assembly `-Y`.

The preprocessor rejects disconnected sections, multiple cavities, internal
holes, baffles or PMD that partition the liquid, ambiguous port rims, open or
invalid solids, and geometries whose section topology cannot be represented by
one total contact perimeter. Phase one also does not support a different
gravity orientation, moving/deforming walls, slosh or multidimensional
circulation, local closure-hardware heat transfer, as-built deformation, or
cryogenic structural effects. Those cases require a new geometry/physics
model, not silent combination or an equivalent diameter.

The approved design and implementation records are:

- [geometry design specification](superpowers/specs/2026-07-23-liqlev-arbitrary-tank-geometry-design.md);
- [implementation plan](superpowers/plans/2026-07-23-liqlev-geometry-kernel-implementation.md); and
- [legacy/report equation discrepancy](physics/legacy-series-discrepancy.md).

## NASA tank authority and coordinate contract

The first production geometry is derived from this immutable source:

```text
C:\Users\sasorian\Documents\Eta_Space\geometry\nhq01-m21a- 0201_TankAssy_NASA.STEP
```

| Provenance field | Authoritative value |
| --- | --- |
| Source size | `36,844,537 bytes` |
| Source SHA-256 | `0193EC296B5B754FFDC7B1652052F5B28BA99345A5BC89922D026DF011C837D5` |
| STEP protocol | AP214 |
| Selected XCAF product | `nhq01-m21a- 0202_short` |
| Minimum closure | `Y = -275.406791 mm` |
| Maximum closure | `Y = +275.296791 mm` |
| Height coordinate | `h = Y - Y_min`, increasing in `+Y` |
| Supported gravity | `-Y` |

The source assembly remains outside the repository and is never written. Its
path, size, and hash are recorded in
[`geometry/source/PROVENANCE.json`](../geometry/source/PROVENANCE.json).

The exact OpenCascade construction selects the named vessel product, forms a
closure slab between the two wet-side rim planes, subtracts the tank material,
and selects the single component containing a strict interior seed. The
accepted result is checked with independent 25 mm and 100 mm carrier padding
and an OpenCascade splitter oracle. The two ports are capped at their wet-side
closure planes, deliberately excluding exterior port barrels and hardware.
The published fluid solid is
[`nhq01-m21a-0201_LIQLEV_FLUID_DOMAIN.step`](../geometry/output/nhq01-m21a-0201_LIQLEV_FLUID_DOMAIN.step).

The two caps close the computational domain; they are not active tank
sidewall. The runtime boundary-layer and wall heat-transfer calculations use
`A_w,side(h)`, which excludes both modeling caps because the omitted closure
hardware has no represented thermal model. `A_w,total(h)` includes the lower
cap and, only at the full endpoint, the upper cap. It is retained for audit and
future explicit closure-surface work.

## CAD-to-JIT data flow

1. OpenCascade selects the exact placed tank product and constructs the
   watertight capped fluid solid.
2. Exact BRep half-space intersections and planar sections measure cumulative
   volume, section area, contact perimeter, and wetted area. Tessellation is
   visualization only.
3. Adaptive refinement creates a versioned numeric package and independently
   checks every final interval midpoint.
4. Package loading verifies the adjacent metadata and NPZ hash before runtime
   arrays are passed to Numba.
5. The JIT branch interpolates/inverts the arrays and integrates the
   perimeter-aware boundary-layer equations without importing OpenCascade.

The production package has `N = 92` adaptively selected CAD-authority nodes:

| NPZ array | Shape | Units | Meaning and interpolation |
| --- | ---: | --- | --- |
| `height_ft` | `(N,)` | ft | Strictly increasing `h`, beginning at zero; interpolation breakpoints |
| `volume_ft3` | `(N,)` | ft³ | Monotone cumulative `V(h)`; cubic PCHIP |
| `volume_coefficients` | `(4, N-1)` | rows: `1`, ft, ft², ft³ | PCHIP `(a,b,c,d)`, evaluated as `((a*dx+b)*dx+c)*dx+d` |
| `section_area_ft2` | `(N,)` | ft² | Direct CAD samples used to validate, not override, runtime `dV/dh` |
| `perimeter_ft` | `(N,)` | ft | Total section contact perimeter; nonnegative linear interpolation |
| `sidewall_area_ft2` | `(N,)` | ft² | Monotone cumulative active sidewall area; cubic PCHIP |
| `sidewall_coefficients` | `(4, N-1)` | rows: ft⁻¹, `1`, ft, ft² | PCHIP `(a,b,c,d)` for `A_w,side(h)` |
| `total_wetted_area_ft2` | `(N,)` | ft² | Node-only audit area with the endpoint-cap convention; not a runtime solver input |

For the committed tank, `(N,)` is `(92,)` and `(4, N-1)` is `(4, 91)`.
Runtime interface area is the analytic derivative of the cumulative-volume
PCHIP, `A(h) = dV/dh`, so the height update and volume conservation use one
authority. The inverse level solve finds the PCHIP interval and performs a
safeguarded monotone solve of `V(h) = V_liquid`; values outside the tabulated
volume are errors, not extrapolations.

The adjacent
[`geometry metadata JSON`](../geometry/tables/nhq01-m21a-0201_LIQLEV_GEOMETRY.json)
declares hashes, units, axes, and closure locations. The
[`CSV table`](../geometry/tables/nhq01-m21a-0201_LIQLEV_GEOMETRY.csv) is a
human-readable copy of the node values. The NPZ remains the runtime authority.

## Physics authority: legacy versus custom mode

Leaving `tank.geometry_path` empty selects the preserved analytic cylinder
branch. Its constant area/perimeter, constant-area height update, coefficient
`2.1`, and historical boundary-layer series remain unchanged, and the saved
legacy baseline remains its regression authority.

Providing a geometry package selects custom mode. NASA report
[NTRS 19700017832](https://ntrs.nasa.gov/citations/19700017832), especially
Eqs. 4-33, 4-37, and 4-38, governs the new branch. The perimeter form is

```text
q = (2/3) P(h) delta^(3/2)
dq/dh = K3 [A(h) - P(h) delta]
dVBL/dh = P(h) delta
```

and uses the actual `A(h)` and `P(h)` rather than an equivalent diameter. The
legacy transcription has different higher-order series denominators and a
different boundary-layer-volume coefficient; this discrepancy is intentional
and documented in
[`docs/physics/legacy-series-discrepancy.md`](physics/legacy-series-discrepancy.md).
Custom-cylinder acceptance therefore compares with the published report
solution, while legacy results compare with the untouched saved baseline.

## Environment and reproducible artifact workflow

From the repository root in PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements-geometry.txt
```

The default build command expects the immutable source at the path and hash
listed above. Use explicit `--source-step` and `--source-sha256` arguments if
the same bytes are stored elsewhere.

Build a technically validated review set:

```powershell
python scripts/build_nasa_tank_geometry.py --visual-section-review pending
python scripts/check_nasa_tank_geometry.py --allow-pending-visual-review
```

Inspect all ten images under `geometry/audit/sections/` at original
resolution. Confirm one filled loop at every interior section, filled endpoint
caps, and no baffles, internal bodies, open ports, or outer-wall regions. Only
after the generated image hashes match the inspected set, publish the visual
attestation and run the independent checker:

```powershell
python scripts/build_nasa_tank_geometry.py --visual-section-review passed
python scripts/check_nasa_tank_geometry.py
```

The builder stages and verifies the entire set before promotion, retains a
recoverable backup if promotion/rollback is interrupted, and restores that
backup on the next build. It normalizes only the STEP `FILE_NAME` timestamp to
`1970-01-01T00:00:00`; STEP DATA entities are not changed. This makes the
generated STEP byte-reproducible.

Run the preserved legacy baseline and report-cylinder equivalence:

```powershell
python scripts/check_physics_baseline.py
python -m pytest tests/geometry/test_boundary_layer.py tests/geometry/test_core_custom_geometry.py -q
```

Run the committed NASA tank case and its 513/1025 evaluation-grid check:

```powershell
python -m pytest tests/geometry/test_nasa_tank_solver.py -q
```

To intentionally regenerate the tracked solver result manifest:

```powershell
python -m validation.custom_geometry_cases
```

That command records the current solver commit, dependency versions, artifact
hashes, five fill cases, convergence counts, and refinement results in
[`validation/results/nasa_tank_geometry_manifest.json`](../validation/results/nasa_tank_geometry_manifest.json).

## Solver configuration

Set `tank.geometry_path` to the NPZ path in a schema-version 2 configuration:

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

`diameter_ft` and `height_ft` are accepted legacy compatibility fields with
schema defaults; they are ignored by the custom-geometry solver. The loaded
package supplies total volume, actual height, interface area, perimeter, and
wetted sidewall area.

To select the legacy analytic cylinder instead:

```json
{
  "tank": {
    "diameter_ft": 4.0,
    "height_ft": 8.0,
    "fill_fractions": [0.5],
    "geometry_path": ""
  }
}
```

In legacy mode, `diameter_ft` and `height_ft` retain their historical meaning.

## Acceptance tolerances and current evidence

| Gate | Acceptance |
| --- | ---: |
| Source identity | Exact SHA-256 and size |
| Closure-plane and round-trip bound error | `1e-5 mm` maximum |
| STEP round-trip volume change | `max(1e-3 mm³, 1e-8 relative)` |
| Adaptive cumulative-volume and sidewall interpolation | `5e-4 relative` (`0.05%`) |
| Final table volume versus exact CAD | `5e-4 relative` |
| Integrated runtime `dV/dh` versus final volume | `5e-4 relative` |
| Independent midpoint cumulative volume and sidewall area | `5e-4 relative` |
| Direct CAD section area versus `dV/dh`, away from topology changes | `2e-3 relative` (`0.2%`) |
| Maximum adaptive CAD nodes | `1025`; failure is explicit |
| Custom cylinder versus report Eqs. 4-33, 4-37, 4-38 | `1e-3 relative` (`0.1%`) |
| NASA solver height and boundary-layer volume, 513 versus 1025 evaluation points | `2e-3 relative` (`0.2%`) |
| Legacy saved physics baseline | `rtol=1e-9`, `atol=1e-8` |

The committed continuous CAD authority contains 92 adaptive nodes. The
Task 11 `513` and `1025` grids are in-memory resamples of that same 92-node
PCHIP authority. They measure solver evaluation-grid sensitivity only; they
are not alternate CAD measurements or evidence that the CAD authority itself
contains 1025 nodes. The recorded maximum height/boundary-layer-volume
difference is `6.371214840947442e-05`, below the `0.2%` limit, with zero
reported convergence failures.

The exact CAD audit is
[`nhq01-m21a-0201_LIQLEV_AUDIT.json`](../geometry/audit/nhq01-m21a-0201_LIQLEV_AUDIT.json).
These CAD, table, analytic-cylinder, boundedness, and grid-sensitivity checks
establish reproducibility and numerical consistency. They are not experimental
validation of LIQLEV for this vessel. Predictions remain report-based
numerical results that require later experimental correlation or
higher-fidelity analysis before being treated as validated tank behavior.
