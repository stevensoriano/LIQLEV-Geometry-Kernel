# Task 7 Report: Exact XCAF Tank Product Selection

## Status

Complete under the approved Task 7 boundary:

- the source STEP is verified before OpenCascade constructs a reader;
- the assembly is loaded through `STEPCAFControl_Reader` into an XCAF
  document;
- free shapes and component labels are traversed recursively;
- references are resolved, names are read from `TDataStd_Name`, and every
  assembly-level `TopLoc_Location` is accumulated;
- exactly one product named `nhq01-m21a- 0202_short` is selected with its
  complete assembly placement;
- no global solid index, approximate bounds, file-order selection, fluid
  construction, face removal, capping, or source-file copy is used.

## Installed CAD environment

```text
Python             3.13.5
cadquery           2.7.0
cadquery-ocp       7.8.1.1.post1
OCP                7.8.1.1
pytest             8.4.2
scipy              1.16.2
```

The installed OCP binding exposes the relevant static methods with `_s`
suffixes:

```text
XCAFDoc_DocumentTool.ShapeTool_s
XCAFDoc_ShapeTool.GetComponents_s
XCAFDoc_ShapeTool.GetReferredShape_s
XCAFDoc_ShapeTool.GetLocation_s
XCAFDoc_ShapeTool.GetShape_s
TDF_Tool.Entry_s
```

Only these Python binding spellings differ from the plan. The required XCAF
semantics are unchanged.

## TDD evidence

The required live-source test was written before either `liqlev/cad` module
existed. Its first run was:

```text
python -m pytest tests/cad/test_xcaf_selection.py -q

ERROR tests/cad/test_xcaf_selection.py
ModuleNotFoundError: No module named 'liqlev.cad'
1 error in 0.25s
```

The adversarial tests were also present before the first implementation. They
cover uppercase hashing, hash rejection before reader construction, missing
source attribution, zero exact product matches, provenance agreement, and
source immutability.

The first implementation run selected the correct shape and passed all error
and provenance assertions. The only failure was the test comparing the whole
Windows stat record:

```text
1 failed, 5 passed in 89.48s
```

Reading the STEP changed its last-access timestamp, while size, modification
time, directory contents, and hash remained unchanged. The test was corrected
to assert those write-sensitive fields and the post-load hash.

A separate duplicate-match RED/GREEN cycle established the ambiguity error:

```text
RED:
ImportError: cannot import name '_select_unique_match'

GREEN:
1 passed, 6 deselected in 1.58s
```

The final focused GREEN was:

```text
python -m pytest tests/cad/test_xcaf_selection.py -q
7 passed in 89.41s
```

## XCAF traversal and placement

`load_named_product` performs these operations in order:

1. stream-hash the source STEP and compare uppercase SHA-256 values;
2. create a `TDocStd_Document` and enable reader name transfer;
3. call `STEPCAFControl_Reader.ReadFile` and `Transfer`;
4. obtain `XCAFDoc_DocumentTool.ShapeTool_s(document.Main())`;
5. start at every label returned by `GetFreeShapes`;
6. recursively visit component labels from `GetComponents_s`;
7. accumulate placement as
   `parent_location.Multiplied(local_location)` at every level;
8. resolve each reference through `GetReferredShape_s`;
9. extract the exact product name from the resolved label's
   `TDataStd_Name`;
10. retrieve the resolved shape and apply only the complete accumulated
    location with `Located`;
11. require one and only one diagnostic match.

The internal traversal returns `_ProductMatch` diagnostics containing the
resolved matched label entry, component-instance label entry, exact product
name, and placed CadQuery shape.

The live assembly traversal reported:

```text
free shapes visited:          1
component/free labels visited: 552
exact matches:                1
instance label entry:         0:1:1:1:5
instance label name:          NAUO5
resolved product label entry: 0:1:1:90
resolved product name:        nhq01-m21a- 0202_short
```

## Selected shape identity

The placed CadQuery solid reported:

```text
xmin = -280.9400000711255 mm
xmax =  280.9400000711255 mm
ymin = -293.96679123713346 mm
ymax =  293.9667912371334 mm
zmin =  296.54275947409957 mm
zmax =  858.4227596163507 mm
faces = 342
shape type = Solid
```

All six bounds satisfy the approved rounded values within absolute tolerance
`1e-6 mm`.

## Source preservation and provenance

Before implementation:

```text
path:
C:\Users\sasorian\Documents\Eta_Space\geometry\nhq01-m21a- 0201_TankAssy_NASA.STEP
size:     36844537 bytes
SHA-256:  0193EC296B5B754FFDC7B1652052F5B28BA99345A5BC89922D026DF011C837D5
modified: 2026-07-21 19:35:04 UTC
```

After all focused, combined, baseline, and full-suite testing:

```text
size:     36844537 bytes
SHA-256:  0193EC296B5B754FFDC7B1652052F5B28BA99345A5BC89922D026DF011C837D5
modified: 2026-07-21 19:35:04 UTC
```

The live test also records the source directory's entry names before loading
and verifies that no sibling file appears, disappears, or changes the source
hash, byte count, or modification time. The 36.8 MB source assembly was not
copied into the repository.

`geometry/source/PROVENANCE.json` records the approved absolute path, hash,
byte count, exact product, `+Y` height axis, `-Y` gravity direction, and closure
planes. The closure planes are provenance only in Task 7; no cap or fluid
geometry is constructed.

## Verification

Pre-change repository baseline:

```text
python -m pytest -q
90 passed in 17.75s
```

Focused Task 7:

```text
python -m pytest tests/cad/test_xcaf_selection.py -q
7 passed in 89.41s
```

Task 7 plus all geometry tests:

```text
python -m pytest tests/cad/test_xcaf_selection.py tests/geometry -q
59 passed in 102.41s
```

Strict physics baseline:

```text
as203_default_high_vent         PASS rows=19 final_pressure=13.626865
hydrogen_height_dep_mid_fill    PASS rows=29 final_pressure=19.496498
nitrogen_custom_epsilon         PASS rows=31 final_pressure=29.823545
Physics baseline check passed.
```

Full suite:

```text
python -m pytest -q
97 passed in 120.62s
```

Whitespace and final status checks are recorded in the controller handoff
after this report is staged.

## Files

- `requirements-geometry.txt`
- `liqlev/cad/__init__.py`
- `liqlev/cad/xcaf.py`
- `tests/cad/test_xcaf_selection.py`
- `geometry/source/PROVENANCE.json`
- `.superpowers/sdd/task-7-report.md`

## Commit

The implementation and this report are committed together as:

```text
feat: select tank body from STEP assembly
```

Configured sole author and committer:

```text
Steven Soriano <steven.a.soriano@nasa.gov>
```

The commit cannot contain its own hash; the final hash is recorded in the
controller handoff. No automated-assistant authorship, co-authorship, or
attribution is included.

## Self-review

- Confirmed SHA-256 verification completes before the STEP reader is
  constructed or transfer is attempted.
- Confirmed read, transfer, zero-match, multiple-match, and null/unusable-shape
  failures cross the public boundary as `StepProductError`.
- Confirmed exact string equality is used for the resolved XCAF product name.
- Confirmed the matched product is resolved through its reference label before
  retrieving and placing the shape.
- Confirmed all parent/component locations are accumulated in assembly order.
- Confirmed successful internal diagnostics retain both resolved and instance
  label entries plus the exact product name.
- Confirmed the returned object is a CadQuery `Solid`/`Shape` wrapper with the
  required placed bounds and 342 faces.
- Confirmed no generic flattened STEP import, global shape index, approximate
  geometric search, or file-order fallback exists.
- Confirmed the source file remains external, read-only in practice, and
  byte-for-byte unchanged after repeated live tests.
- Confirmed `core.py`, `thermo_utils.py`, geometry kernels, saved baselines,
  and solver physics were not modified.
- Confirmed no Task 8 face-network selection, fluid construction, sewing,
  capping, healing, or output STEP work was started.

## Remaining concerns

No unresolved Task 7 contract concern remains. The selected body is the tank
material authority for Task 8 and is deliberately not an internal fluid-volume
solid.
