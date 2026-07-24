# Task 10 Report: Published NASA Tank Geometry Artifacts

## Status

Complete and ready for independent review.

The final repository artifact set is technically validated, independently
remeasured, visually reviewed, reproducible, and marked:

```text
technical_validation_passed = true
visual_section_review = "passed"
passed = true
```

Task 11 has not started.

## Authoritative input

- Source STEP:
  `C:\Users\sasorian\Documents\Eta_Space\geometry\nhq01-m21a- 0201_TankAssy_NASA.STEP`
- Source size: `36,844,537` bytes
- Source SHA-256:
  `0193EC296B5B754FFDC7B1652052F5B28BA99345A5BC89922D026DF011C837D5`
- Selected exact XCAF product: `nhq01-m21a- 0202_short`
- Approved closure planes:
  `Y = -275.406791 mm` and `Y = 275.296791 mm`
- Axis: `+Y`
- Gravity direction: `-Y`
- Node limit: `1025`

The build recorded identical source hashes, sizes, and modification times
before and after processing. `source_preservation.unchanged` is `true`.
The source STEP was never written.

## Implementation

### Deterministic staged build command

Created `scripts/build_nasa_tank_geometry.py`.

The command:

1. hashes the immutable source;
2. selects the exact named XCAF product;
3. constructs the reviewed exact closure-slab fluid solid;
4. exports and independently reimports AP214 STEP;
5. normalizes only the AP214 `FILE_NAME` timestamp to the fixed value
   `1970-01-01T00:00:00`;
6. reimports the normalized exact STEP;
7. adaptively measures the exact BRep and builds the geometry kernel;
8. saves NPZ, JSON, and CSV;
9. independently remeasures every final interval midpoint;
10. renders three orthogonal images and seven X-Z height sections;
11. verifies every staged file and recorded hash; and
12. promotes the complete set with rollback and interrupted-promotion
    recovery.

The exact STEP geometry is not simplified or substituted. Timestamp
normalization makes repeated STEP builds byte-reproducible without changing
DATA entities.

The promotion path:

- never touches the existing validated files before all staged checks pass;
- restores prior files on ordinary exceptions and `KeyboardInterrupt`;
- leaves a persistent manifest/backup during promotion; and
- recovers a hard-interrupted promotion before the next build.

### Independent verification command

Created `scripts/check_nasa_tank_geometry.py`.

It does not import or call the build command. It independently:

- re-hashes the source and every published artifact;
- reimports the STEP;
- reruns BRep validity, topology, closed-shell, fixed-bound, closure-plane,
  and 27/3 endpoint-mosaic checks;
- loads and validates the NPZ/JSON package;
- requires fixed metadata and closure constants;
- compares all CSV numeric values bit-for-bit with the NPZ arrays;
- remeasures exact CAD volume, sidewall area, and section area at all 91
  final interval midpoints;
- recomputes integrated volume from `dV/dh`;
- checks the 0.05% cumulative/refinement gates and 0.2% direct-area gate;
- verifies every image hash and render mode;
- verifies recorded tool/command provenance; and
- requires visual review to be passed by default.

The opt-in `--allow-pending-visual-review` flag exists only for inspecting a
technically validated review set before human visual acceptance.

### Adaptive refinement integration

The NASA run exposed that a fixed topology-neighborhood expansion could
consume the 1025-node budget even after cumulative volume and sidewall
errors had converged. The first captured failure was:

```text
adaptive measurement did not converge within max_nodes=1025
```

Per-driver diagnostics then showed:

```text
iteration=8
current_nodes=792
proposed_nodes=1161
volume_failures=0
sidewall_failures=0
direct_area_failures=5
topology_neighborhood_additions=364
```

The final rule follows the local PCHIP derivative stencil:

- every direct-area-failing interval is refined;
- only its immediate left and right intervals are added;
- indices are clipped, so the first interval cannot wrap to the last;
- topology ordinates remain part of eligibility/exclusion metadata;
- cumulative volume and sidewall masks are unchanged; and
- all numeric tolerances and the 1025-node limit are unchanged.

Exact per-ordinate caching reuses prior interior BRep measurements between
adaptive iterations. Tests prove exact array equality with uncached
measurement and disjoint prior/new interior request sets. This reduced the
NASA adaptive stage from a captured 824-second failure to a successful
41.8-second stage.

### Exact BRep metric ordering

OpenCascade visualization tessellation attaches an approximate triangulation
to the in-memory shape. A diagnostic probe proved:

```text
exact max bound before tessellation: authoritative
max bound change after tessellate(0.8, 0.08): 0.9019344254967336 mm
fresh STEP reimport: authoritative exact bounds restored
```

The build now snapshots volume, area, bounds, topology, validity, and shell
closure immediately after normalized STEP reimport and before any rendering.
A regression proves visualization cannot change the recorded exact snapshot.

## Final artifacts

### Core files

| Artifact | Size (bytes) | SHA-256 |
|---|---:|---|
| `geometry/output/nhq01-m21a-0201_LIQLEV_FLUID_DOMAIN.step` | 186,451 | `648C9ACEEEF9CD0283617608E7AECA62F92A814A5F67E130A152A330C30BA5AE` |
| `geometry/tables/nhq01-m21a-0201_LIQLEV_GEOMETRY.npz` | 11,750 | `859EE059A5B0E3A095A6F62697EA3B2F936B6CA057532C5734C81D39D5085C9E` |
| `geometry/tables/nhq01-m21a-0201_LIQLEV_GEOMETRY.json` | 528 | `C4AC34354118A8A207846EF92B37DAF1B0B890A4D2453FFEEE195D7FBE8E448B` |
| `geometry/tables/nhq01-m21a-0201_LIQLEV_GEOMETRY.csv` | 10,510 | `0DB1865758795BD2C4B00A540402541306DA4296D3544E92429AA43EE5E4C5A2` |
| `geometry/audit/nhq01-m21a-0201_LIQLEV_AUDIT.json` | 9,402 | `3142569BE1003C2A48D2702882772038C5F81E5867D40D922811E3725E33F355` |

The normalized review STEP and production STEP hashes match. The NPZ, CSV,
and all ten production image hashes match the visually inspected review set.

### Section images

| Image | SHA-256 |
|---|---|
| `sections/orthogonal_xy.png` | `1E086CAE3CC368A5441408CA3BB1BB3A890B5AEAAE46648AACD2938E92102EBE` |
| `sections/orthogonal_yz.png` | `4640E4B8C644D86126928B15DABA0685D4E9D6CF752755A352D23FFF667F1768` |
| `sections/orthogonal_xz.png` | `2498A52F08BE87A82B17A4ED011BC1DE4211E63768D126DC0B9E04E28F7C53EB` |
| `sections/section_xz_h000.png` | `4FF9F70F6425C5C759642726A7F25E7B52181C61F7907E0A19D69D4865F48472` |
| `sections/section_xz_h010.png` | `0F927443AEEF2B7470627B829530C5981BC934FF98B5C9AD25E4238589EC7ADE` |
| `sections/section_xz_h025.png` | `7FDD5B03B6C00171F8AD1FFCA30826D1560824DC3A700CE573F342320719E30F` |
| `sections/section_xz_h050.png` | `C0ED51F4942EAB75330B46B5CA217FDDA039EC93BD6AE116215C7C2CC297F8A4` |
| `sections/section_xz_h075.png` | `1AE41693D956EA83CD160D288D7985C0CF313E6B5D3AFB2C29B1A44D97224E6D` |
| `sections/section_xz_h090.png` | `E59EDED3B9E0A1849E497912127982F8076F792FECFA6DBF27761D995626D930` |
| `sections/section_xz_h100.png` | `D773F69FADFEFE822F81FDF10EE7A85F25B46E1198D5FE03573A1B7B3F4AEB6B` |

Every image labels assembly axes, absolute `Y`, normalized `h`, and a
100-millimetre scale bar. The 0% and 100% images render the actual endpoint
closure-face mosaics, not inset sections. Interior sections use topology-
ordered wires rather than polar point sorting.

## Final CAD and kernel metrics

- Solids: `1`
- Shells: `1`
- Faces: `40`
- Valid: `true`
- Closed outer shell: `true`
- Endpoint closure mosaic faces: `27` minimum / `3` maximum
- Volume: `98,109,377.71163802 mm^3`
- Volume: `3.4646999761391375 ft^3`
- Surface area: recorded in the audit
- Bounds `(xmin, xmax, ymin, ymax, zmin, zmax)` in millimetres:

  ```text
  (-279.67000017112,
    279.67000017112,
   -275.4067913371,
    275.29679133713,
    297.81275947409995,
    857.15275961635)
  ```

- Cap-plane maximum error: `2.37133406244538e-07 mm`
- Round-trip absolute volume difference: `0.03622707724571228 mm^3`
- Round-trip relative volume difference: `3.69251931643207e-10`
- Nodes: `92`
- Height: `1.806770282152231 ft`
- CAD endpoint volume relative error: `0.0`
- Integrated `dV/dh` volume relative error: `1.10606540386985e-07`
- Maximum midpoint cumulative-volume relative error:
  `3.16442237326085e-05`
- Maximum midpoint sidewall relative error:
  `4.47756186517317e-05`
- Maximum eligible direct-area relative error:
  `0.0017377546808950138`
- Eligible direct-area midpoint checks: `91`

All values are below their fixed acceptance tolerances.

## Visual review

All ten review PNGs were inspected at original resolution before the final
production build was marked passed.

- Orthogonal views show the expected vessel envelope.
- The 10%, 25%, 50%, 75%, and 90% sections each show one filled exterior
  loop.
- No baffles, internal bodies, open ports, or outer-wall regions are present.
- The 0% and 100% sections are filled endpoint caps.
- Circular/annular lines visible at the minimum endpoint are the approved
  coplanar closure-mosaic seams. The exact endpoint reconstruction has one
  assembled exterior boundary and no inner boundaries, so those lines are
  not holes or ports.

The final production image hashes match the inspected review hashes exactly.

## Tool and command provenance

Recorded by the final audit:

- Python `3.13.5`
- CadQuery `2.7.0`
- OCP `7.8.1.1`
- NumPy `2.3.4`
- Matplotlib `3.10.7`
- Platform `Windows-11-10.0.22631-SP0`
- Interpreter:
  `C:\Users\sasorian\Miniconda3\python.exe`
- Final build arguments:
  `scripts/build_nasa_tank_geometry.py --visual-section-review passed`

Final build stage timings in seconds:

| Stage | Seconds |
|---|---:|
| Product selection | 32.485855 |
| Fluid construction | 2.504082 |
| STEP round trip | 0.455953 |
| STEP timestamp normalization | 0.001144 |
| Normalized STEP reimport | 0.044066 |
| Adaptive measurement and coefficients | 41.795401 |
| Package save | 0.005094 |
| Independent CAD/table measurement | 20.328567 |
| Section rendering | 2.832745 |
| Staged validation | 0.283906 |
| Total before promotion | 101.081469 |

## TDD evidence

The Task 10 artifact contract was written first. Its initial RED listed the
missing builder, checker, NPZ, JSON, CSV, and ten PNGs.

Additional regressions were observed RED before their fixes:

- max-node failures lacked per-driver diagnostics;
- topology-neighborhood expansion selected nonlocal intervals;
- immediate-topology marching collapsed an ordinate into the sphere pole;
- one-ring refinement helper was absent;
- adaptive iterations remeasured prior interior ordinates;
- independent CAD/table measurement helper was absent;
- topology-ordered wire rendering helper was absent;
- endpoint renderer did not use actual closure mosaics;
- `KeyboardInterrupt` left a partial publication;
- hard-interrupted promotion recovery was absent;
- visualization tessellation could contaminate recorded exact bounds; and
- STEP header timestamps were not deterministic.

All focused regressions are green.

## Verification evidence

### Independent production check

```text
python scripts/check_nasa_tank_geometry.py
```

Result: every independent check passed, ending with:

```text
NASA tank geometry audit passed.
```

### Task 10 artifact tests

```text
python -m pytest tests/cad/test_nasa_tank_artifacts.py -q
```

Result:

```text
9 passed in 3.62s
```

### Analytic CAD measurement tests

```text
python -m pytest tests/cad/test_analytic_measurements.py -q
```

Result:

```text
19 passed in 15.94s
```

### CAD and geometry regression

```text
python -m pytest tests/cad tests/geometry -q
```

Result:

```text
107 passed in 136.82s
```

The production STEP and audit hashes were captured before and after this
gate and remained unchanged.

### Physics baseline

```text
python scripts/check_physics_baseline.py
```

Result:

```text
Physics baseline check passed.
```

All three baseline cases passed with their preserved row counts and numeric
summaries.

### Full repository suite

```text
python -m pytest -q
```

Result:

```text
145 passed in 138.91s
```

## Scope and residual concern

- `core.py`, `thermo_utils.py`, and solver physics were not changed.
- No tolerance was weakened.
- The 1025-node maximum was not raised.
- The exact fluid STEP topology was not simplified.
- The legacy Task 8 “final inspection” pytest fixture now writes to pytest
  temporary storage so regression tests cannot overwrite the first-class
  Task 10 publication.
- Narrow `.gitattributes` rules disable newline conversion for published
  STEP/CSV/JSON artifacts, preserving their audited bytes and hashes across
  checkouts; the artifact test covers this contract.
- `topology_coincidence` remains a narrow eligibility exclusion. Inserted
  topology ordinates normally become table nodes, so midpoint coincidence is
  uncommon; the direct-area one-ring regression and non-vacuous eligible
  midpoint requirement remain the active protections.
- Task 11 must wait for clean independent review of this Task 10 commit.
