# Task 8 Report: Deterministic Fluid-Domain Topology Blocker

## Status

Blocked before production implementation, as required by the Task 8 brief.

The approved deterministic traversal does not produce an inner-wall network
whose only free boundaries are the two selected wet-side loops. Exact
OpenCascade sewing and free-boundary analysis reports 18 closed free-boundary
wires before capping: the two intended wet-side loops and 16 additional
micro-hole loops. Adding exactly the two permitted caps leaves one open shell
with 16 closed free-boundary wires and zero solids.

No selection criterion was weakened. The controller authorized one in-memory,
evidence-only diagnostic that capped the 16 unexpected loops; it was not
adopted as production behavior and wrote no STEP. No face removal, Boolean
fill, mesh, voxel, healing, contour simplification, or approximate
reconstruction was applied. No production STEP was written and Task 9 was not
started.

The complete machine-readable evidence is in
`.superpowers/sdd/task-8-topology-inventory.json`.

## TDD evidence

`tests/cad/test_fluid_domain.py` was created before either required production
module. It includes a module-scoped live tank fixture, a module-scoped fluid
fixture, a module-scoped temporary round-trip fixture, exact solid/bounds
checks, audit JSON checks, source preservation, invalid plane/tolerance cases,
no-candidate rejection, non-solid audit rejection, and optional-dependency /
external-source skips.

The required missing-construction RED was:

```text
python -m pytest tests/cad/test_fluid_domain.py -q

ERROR at setup
ModuleNotFoundError: No module named 'liqlev.cad.fluid_domain'
11 errors in 1.65s
shell wall time: 4.2s
```

There is intentionally no GREEN. The Task 8 brief says to stop with a
machine-readable inventory when the actual topology violates the approved
deterministic rule. Creating production modules after confirming that
violation would conceal the blocker or require an unapproved selection rule.

## Source and selected body

The live product was loaded only through
`liqlev.cad.xcaf.load_named_product`, using exact XCAF product name
`nhq01-m21a- 0202_short` and the approved SHA-256.

```text
CadQuery shape type:       Solid
valid:                     true
solids:                    1
shells:                    9
faces:                     342
CadQuery edges:            700
CadQuery wires:            404
OCCT edge identity map:    764
xmin:                     -280.9400000711255 mm
xmax:                      280.9400000711255 mm
ymin:                     -293.96679123713346 mm
ymax:                      293.9667912371334 mm
zmin:                      296.54275947409957 mm
zmax:                      858.4227596163507 mm
```

The `764`-entry OCCT map includes topological/degenerated edge entries that
CadQuery's public `Edges()` list does not expose as 700 ordinary edges. All
adjacency and traversal used OCCT topological identity through
`TopExp.MapShapesAndAncestors_s`; coordinates were not used to match edges.

## Candidate rim faces and wires

The plane tolerance was `1e-5 mm`; normal parallelism tolerance was `1e-8`.
Exactly one planar, Y-normal candidate face was found at each closure plane.
The face-map indices below are diagnostics only, not selection criteria.

### Minimum-Y closure

```text
requested Y:               -275.406791 mm
candidate faces:           1
actual centroid Y:         -275.4067912371334 mm
cap-plane offset:           2.37133427560865e-7 mm
normal:                    (-6.475400626674115e-17,
                             1.0,
                             6.585737027495154e-19)
face area:                  13688.085013581547 mm^2
face wires:                 26
face edges:                 52
OCCT face-map diagnostic:   88
```

The minimum-Y candidate is not a simple two-wire annulus. Its exact closed
wire inventory is:

```text
outer rim wire:             25721.900513867593 mm^2, 2 edges
selected central wet loop:  11432.587497137267 mm^2, 2 edges
auxiliary loops:             8 x 75.10671608922654 mm^2, 2 edges each
micro-hole loops:           16 x approximately
                              0.02339215218554 mm^2, 2 edges each
```

The outer rim and selected wet loop have coincident planar enclosed-area
centres within floating-point evaluation noise at approximately
`(0, -275.4067912371334, 577.4827595452252) mm`. The selected loop is the
smaller member of that central annular pair (the largest non-outer closed
wire), not one of the auxiliary holes.

Selected wet-loop OCCT edge-map diagnostics were `305` and `306`. Their
adjacent non-rim traversal starts were face-map diagnostics `161` and `166`.

### Maximum-Y closure

```text
requested Y:                275.296791 mm
candidate faces:           1
actual centroid Y:          275.2967912371334 mm
cap-plane offset:           2.37133370717446e-7 mm
normal:                    ( 6.475400626674115e-17,
                            -1.0,
                            -6.585737027495154e-19)
face area:                  14289.3130167304 mm^2
face wires:                 2
face edges:                 4
OCCT face-map diagnostic:   52
outer rim wire:             25721.900513867666 mm^2, 2 edges
selected wet loop:          11432.587497137267 mm^2, 2 edges
```

Selected wet-loop OCCT edge-map diagnostics were `183` and `187`. Their
adjacent non-rim traversal starts were face-map diagnostics `50` and `86`.

## Traversal and intersection

Both traversals excluded the two candidate rim faces and used only
edge-to-face adjacency by OCCT topological identity.

```text
minimum-Y start faces:      2
maximum-Y start faces:      2
minimum-Y reachable faces:  218
maximum-Y reachable faces:  218
reachable sets equal:       true
intersection faces:         218
intersection components:    1
component face count:       218
```

The selected network's exact face geometry distribution was:

```text
CONE:                       128
CYLINDER:                    74
PLANE:                        4
TORUS:                       12
total:                      218
```

Thus the reachable-set and connectivity portions of the approved rule pass.
The free-boundary condition does not.

## Exact free-boundary blocker

The 218 selected faces were added unchanged to
`BRepBuilderAPI_Sewing(1e-5)`. `ShapeAnalysis_FreeBounds` was then evaluated
on the sewn result so periodic seam occurrences were not misclassified as
physical boundaries.

Before caps:

```text
input/sewn faces:            218 / 218
sewn shape type:             Shell
shells:                      1
closed shells:               0
solids:                      0
sewing free edges:           36
sewing degenerated shapes:   32
contiguous edges:            0
multiple edges:              0
deleted faces:               0
closed free-boundary wires:  18
open free-boundary wires:    0
expected wet boundaries:      2
unexpected boundaries:       16
```

The 18 closed wires are exactly:

- two 2-edge wet-side loops, each area
  `11432.587497137265 mm^2`, at the two closure planes; and
- sixteen 2-edge micro-hole loops at the minimum-Y plane, each area
  approximately `0.02339215218554 mm^2`.

All 16 unexpected loop centres and precise evaluated areas are recorded in
the machine-readable inventory. Their Y coordinates are
`-275.4067912371334 mm` or `-275.40679123713335 mm`.

After adding planar faces from exactly the two selected wet-side loops:

```text
input/sewn faces:            220 / 220
sewn shape type:             Shell
shells:                      1
closed shells:               0
solids:                      0
sewing free edges:           32
sewing degenerated shapes:   32
closed free-boundary wires:  16
open free-boundary wires:    0
```

This proves the earlier raw 68-edge diagnostic included periodic seam /
occurrence artifacts, but it did not create the failure. The exact sewn
free-boundary result still has 16 physical closed loops after the only two
approved caps.

The mandatory Task 8 conditions therefore fail:

1. The inner-wall network's only free boundary wires are not the two selected
   wet-side loops.
2. Sewing the network plus those two caps does not make one closed shell.
3. A positive-volume solid cannot be made from that open shell.

## Evidence-only cap-all diagnostic

After the blocker was confirmed, the controller requested one temporary
diagnostic to quantify a possible future rule revision. In memory only, all
18 exact free loops were capped: the two selected main loops plus all 16
minimum-plane micro-hole loops. This does not authorize that construction for
Task 8.

The micro loops have a defensible exact-topology classification:

```text
count:                              16
closure plane:                      all at minimum-Y within 1e-5 mm
edges per loop:                      2
geometry of both edges:             CIRCLE
enclosed area range:                0.02339215218554371 to
                                     0.023392152185574534 mm^2
edge radius range:                  0.08628993741899293 to
                                     0.08628993741904845 mm
perimeter range:                    0.5421756669484624 to
                                     0.5421756669488111 mm
distance from central port axis:    78.31449058183607 to
                                     79.92303334527219 mm
```

They are therefore 16 distinct, two-semicircle, circular micro-hole loops on
the minimum closure flange, radially distributed around the central port.
Their exact centres are retained in the JSON inventory.

The cap-all sewing result was:

```text
selected network faces:             218
main cap faces:                       2
micro-hole cap faces:                16
total cap faces:                     18
maximum cap-plane error:              2.371334062445385e-7 mm
sewn faces:                          236
sewn shells:                           1
sewn solids:                           0
sewing free edges:                     0
remaining closed free wires:           0
remaining open free wires:             0
sewing degenerated shapes:             32
TopoDS shell Closed flag:             false
```

Although the shell's stored `TopoDS.Closed()` flag was not set by sewing,
`BRepBuilderAPI_MakeSolid` completed and both CadQuery validity and
`BRepCheck_Analyzer` validity passed:

```text
solids:                                1
shells:                                1
faces:                               236
valid (CadQuery / BRepCheck):          true / true
volume:                                99912238.01532683 mm^3
surface area:                          1087207.7669801908 mm^2
x bounds:                             -280.9400000711255 to
                                       280.9400000711255 mm
y bounds:                             -293.96679123713346 to
                                       293.9667912371334 mm
z bounds:                              296.54275947409957 to
                                       858.4227596163507 mm
```

This result supplies stronger blocker evidence. Capping all 18 free loops
removes the free boundaries, but the resulting valid B-Rep retains the entire
source tank-body bounds. Its minimum Y misses the approved closure by
`18.560000237133486 mm`, and its maximum Y misses by
`18.670000237133405 mm`. The 218-face traversal is thus the broader source
shell network, not the closure-bounded internal fluid wall.

The diagnostic took `43.345113700000184 s` total, including
`42.35172950000015 s` for the single XCAF load. It wrote no STEP or audit
artifact.

## Source preservation

The source was re-statted and re-hashed after the live topology gates:

```text
path:      C:\Users\sasorian\Documents\Eta_Space\geometry\
           nhq01-m21a- 0201_TankAssy_NASA.STEP
size:      36844537 bytes
modified:  2026-07-21T19:35:04.1500332Z
SHA-256:   0193EC296B5B754FFDC7B1652052F5B28BA99345A5BC89922D026DF011C837D5
```

The hash, size, and modification timestamp match the approved Task 7
authority. No source copy or modification was made.

## Verification and timing

Completed:

```text
missing-function RED:       11 errors in 1.65s (4.2s shell wall time)
candidate/wire inventory:   44.8s shell wall time
identity traversal probe:   46.9s shell wall time
sewn free-boundary probe:   46.6s shell wall time
source post-hash/stat:       pass, 0.6s shell wall time
```

Construction, STEP round-trip, audit, focused GREEN, combined CAD, geometry,
physics-baseline, and full-suite gates are not applicable after the mandated
stop. No production solid exists to audit or round-trip. Running a GREEN or
claiming acceptance would require changing the approved topology rule.

## Files and commit

Created as blocker evidence:

- `tests/cad/test_fluid_domain.py`
- `.superpowers/sdd/task-8-topology-inventory.json`
- `.superpowers/sdd/task-8-report.md`

Not created:

- `liqlev/cad/fluid_domain.py`
- `liqlev/cad/audit.py`
- any final or temporary STEP artifact in the repository

The production commit `feat: construct capped tank fluid domain` was not made.
The controller authorized preserving only this blocked report and its
machine-readable inventory in an evidence-only commit:

```text
docs: record fluid-domain topology blocker
```

The evidence commit excludes the intentionally RED test and all production
modules. Its sole author and committer is
`Steven Soriano <steven.a.soriano@nasa.gov>`. The commit cannot contain its own
hash; the hash is returned in the controller handoff.

## Self-review

- Confirmed the required RED occurred before any production module existed.
- Confirmed the assembly was loaded by exact XCAF product name and approved
  source hash.
- Confirmed face adjacency used OCCT topological identity, not coordinate
  matching.
- Confirmed one candidate rim face exists at each approved closure plane.
- Confirmed each plausible central wet loop has exact enclosed area
  `11432.587497137267 mm^2`.
- Confirmed both reachable sets contain the same connected 218-face network.
- Confirmed free boundaries with OCCT sewing and free-boundary analysis,
  eliminating periodic seam artifacts from the blocker decision.
- Confirmed exactly 16 unexpected closed loops remain after the two permitted
  caps.
- Confirmed the controller-authorized cap-all diagnostic remained in memory,
  was not adopted as Task 8 behavior, and independently failed the closure
  bounds by more than 18 mm at both ends.
- Confirmed no topology-changing repair was attempted.
- Confirmed the external source STEP remains byte-for-byte unchanged.
- Confirmed no production STEP was published and no Task 9 work began.

## Required controller decision

The deterministic selection/construction rule must be revised before Task 8
can continue. The exact evidence indicates that traversing the full connected
face network also reaches the walls of 16 minimum-flange micro-holes.

Any continuation needs explicit authority for a new exact-topology rule, such
as a deterministic way to exclude those micro-hole wall faces or to treat
their loops as additional closure boundaries. This report does not recommend
or implement either change because both alter the approved requirement that
the selected network have only the two wet-side free boundaries.

---

# Task 8 Revision Continuation: Exact Closure-Slab Fluid Domain

## Revised status

GREEN after explicit user approval of the exact slab-minus-tank construction.
This continuation preserves the blocker history above and supersedes only its
obsolete conclusion that Task 8 lacked an approved deterministic construction.

The normative result is the unique direct-cut component that strictly contains
the loop-derived seed. It is one valid 40-face solid with one closed outer
shell, exact closure-plane bounds, no carrier contact, no rejected-component
topology, and the approved reference volume and surface area. The independent
`100 mm` padding construction and exact splitter oracle agree.

Task 9 was not started. No solver physics, numeric geometry-table extraction,
mesh operation, healing, source-face ranking, equivalent primitive, or
approximate fill was introduced.

## Revised TDD evidence

The intentionally untracked original test remained the starting point. Before
either production module existed, its expectations were revised to cover the
approved construction, exact reference metrics, Boolean inventories,
padding/splitter oracles, closure mosaic, rejection gates, AP214 round trip,
final inspection artifact, and source preservation.

The fresh missing-production-module RED was:

```text
python -m pytest tests/cad/test_fluid_domain.py -q

19 errors in 1.83s
ModuleNotFoundError: No module named 'liqlev.cad.audit'
```

The first live implementation runs then exposed only installed OCP/CadQuery
binding details: OCP 7.8 provides `IsDone()` but not `HasErrors()` on these
Boolean wrappers, `BRepAlgoAPI_Splitter` accepts `SetArguments`/`SetTools`, the
closed flag belongs to the sole outer shell rather than the solid container,
and `cadquery.importers.importStep` returns a `Workplane`. These corrections
did not change the approved algorithm or numeric result.

The focused GREEN was:

```text
python -m pytest tests/cad/test_fluid_domain.py -q

19 passed in 57.56s
```

## Exact construction

The placed source solid is loaded only through
`liqlev.cad.xcaf.load_named_product` using the exact product name
`nhq01-m21a- 0202_short` and approved source SHA-256.

At `1e-5 mm` plane tolerance and `1e-8` normal-parallelism tolerance, rim
inventory produced:

```text
minimum-Y rim candidates:             1
maximum-Y rim candidates:             1
minimum-Y closed wires:              26
maximum-Y closed wires:               2
minimum-Y wet-loop area:  11432.587497137267 mm^2
maximum-Y wet-loop area:  11432.587497137267 mm^2
minimum wet-loop centre:
  ( 2.1283774576526403e-14,
   -275.40679123713335,
    577.4827595452252) mm
maximum wet-loop centre:
  (-1.4376488654129075e-14,
    275.2967912371334,
    577.4827595452252) mm
```

The exact loop areas and centres select the central wet loops rather than the
eight auxiliary loops or sixteen minimum-flange micro loops documented above.
The derived seed, with no hard-coded location, is:

```text
(3.453642961198664e-15,
 -0.0549999999999784,
 577.4827595452252) mm
```

The primary carrier uses exact closure `Y` faces and `25.0 mm` padding beyond
the tank's `X/Z` bounds. Direct `BRepAlgoAPI_Cut(carrier, tank_body)` produced:

```text
whole result solids / shells / faces:   2 / 2 / 56
whole result valid:                     true
seed classifications:                  OUT, IN
selected solids / shells / faces:       1 / 1 / 40
selected valid / closed outer shell:    true / true
volume:                                 98109377.7478651 mm^3
surface area:                           1031668.462706133 mm^2
bounds:
  (-279.6700000711255, 279.6700000711255,
   -275.4067913371334, 275.2967913371334,
    297.8127594740996, 857.1527596163507) mm
carrier X/Z contact:                    false
exterior probes:                        OUT, OUT, OUT, OUT
shared faces / edges / vertices:        0 / 0 / 0
minimum/maximum closure faces:          27 / 3
maximum closure-plane error:            3.371334287294303e-7 mm
```

No volume ranking is used. Every direct-cut solid is classified, exactly one
must contain the seed strictly `IN`, and `ON` is rejected.

The independently repeated `100.0 mm` construction selected the same valid
40-face topology:

```text
solid / shell / face counts:            1 / 1 / 40
volume difference:                      4.470348358154297e-8 mm^3
surface-area difference:                3.4924596548080444e-10 mm^2
six bounding-box differences:           all 0.0 mm
```

The independent exact splitter produced three solids with seed
classifications `OUT, OUT, IN`. Its selected valid 40-face fluid had zero
volume, area, and six-coordinate bounds difference from the direct cut.

The accepted 27/3 closure-face mosaic is exact planar Boolean topology. No
face-merging requirement or common-then-cut substitute is applied.

## AP214 round trip and deterministic audit

The final handoff STEP was exported explicitly as AP214IS, independently
re-imported, and checked again. The source authority is recorded separately
from the input B-Rep serialization and generated output hashes.

```text
pre volume:                       98109377.7478651 mm^3
post volume:                      98109377.71163802 mm^3
absolute volume difference:       0.03622707724571228 mm^3
relative volume difference:       3.6925193164320727e-10
pre area:                         1031668.462706133 mm^2
post area:                        1031668.4622880891 mm^2
surface-area difference:          0.0004180439282208681 mm^2
maximum bounds difference:        9.999450867326232e-8 mm
cap-plane error:                  3.371334287294303e-7 mm
pre/post solids:                  1 / 1
pre/post closed outer shells:     1 / 1
pre/post faces:                   40 / 40
pre/post validity:                true / true
audit passed:                     true
```

An additional independent CadQuery import of the final file reported one
solid, one closed outer shell, 40 faces, validity `true`, volume
`98109377.71163802 mm^3`, and area `1031668.4622880891 mm^2`.

Final artifact authorities:

```text
source SHA-256:
0193EC296B5B754FFDC7B1652052F5B28BA99345A5BC89922D026DF011C837D5

input B-Rep SHA-256:
0D1A0B60BDC2056DCF7714E1CBE0C79A4A049B96DAB29AC44CDD12C440AB637F

output STEP SHA-256:
F139D3EBB745DF7CEAD64998AF1B6CB953EEE303C89914A6EE2FA56759677BA9

audit JSON SHA-256:
93DB82CAFB6CF9AA574868CD96B554B91ABF71C6767C9D5FA0DFD4F0B2A1D246
```

Final artifact paths:

```text
geometry/output/nhq01-m21a-0201_LIQLEV_FLUID_DOMAIN.step
geometry/audit/nhq01-m21a-0201_LIQLEV_AUDIT.json
```

The final timed live gate took:

```text
XCAF load:                  34.91323830000056 s
exact construction:         2.6782555999998294 s
AP214 write/re-import/audit: 0.5387886000007711 s
```

## Source preservation

After focused, combined, CAD/geometry, full-suite, and final timed live gates:

```text
size:             36844537 bytes
modified time ns: 1784662504150033200
SHA-256:          0193EC296B5B754FFDC7B1652052F5B28BA99345A5BC89922D026DF011C837D5
```

All three values exactly match the pre-operation snapshot. The source file was
never copied into the repository, modified, overwritten, or used as an export
target.

## Verification

```text
python -m pytest tests/cad/test_fluid_domain.py -q
19 passed in 57.56s

python -m pytest tests/cad/test_xcaf_selection.py tests/cad/test_fluid_domain.py -q
27 passed in 158.79s

python -m pytest tests/cad tests/geometry -q
79 passed in 166.41s

python scripts/check_physics_baseline.py
Physics baseline check passed for all three cases.

python -m pytest -q
117 passed in 128.50s
```

The final `git diff --check`, status, artifact inspection, and authorship
checks are performed after this report is staged and are recorded in the
controller handoff.

## Files and artifact workflow

Committed Task 8 intent:

- `liqlev/cad/fluid_domain.py`
- `liqlev/cad/audit.py`
- `liqlev/cad/__init__.py`
- `tests/cad/test_fluid_domain.py`
- `docs/superpowers/specs/2026-07-23-liqlev-arbitrary-tank-geometry-design.md`
- `docs/superpowers/plans/2026-07-23-liqlev-geometry-kernel-implementation.md`
- `.superpowers/sdd/task-8-report.md`

The final generated STEP and audit JSON remain uncommitted inspection
artifacts for the later repository artifact workflow, as explicitly permitted
by the revised brief.

The intended commit message is:

```text
feat: construct exact tank fluid domain
```

Its sole author and committer are
`Steven Soriano <steven.a.soriano@nasa.gov>`. The commit cannot contain its own
hash; the controller handoff records that hash. No automated-assistant
authorship, co-authorship, or attribution is included.

## Assumptions and remaining validation limitations

- The exact selected XCAF product remains the authoritative tank-material
  solid. Separate PMD, baffle, spray-bar, fastener, and viewport-hardware
  assembly products are neither fused nor subtracted.
- The loop-derived central seed is safely interior for this source and is
  rejected unless exactly one Boolean component classifies it strictly `IN`.
- The `25/100 mm` padding comparison and splitter agreement demonstrate
  construction invariance for the two approved exact oracles; they are not
  proof against every possible future OpenCascade-version change.
- STEP AP214 translation preserves acceptance metrics within the approved
  tolerances but, as shown by the nonzero differences, is not byte-identical
  B-Rep serialization.
- This work establishes exact B-Rep numerical/topological consistency with
  the approved source. Experimental validation, as-built tolerances,
  deformation, cryogenic operating effects, and any higher-fidelity model
  comparison remain later engineering validation activities.
