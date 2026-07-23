# Task 6 Report: Configuration and Headless Geometry Packages

## Status

Complete under the approved Task 6 boundary:

- saved configurations now use schema version 2 and include
  `TankConfig.geometry_path`;
- explicit schema-version-1 files still load with an empty geometry path;
- headless single and sweep runs load validated numeric geometry packages;
- custom inputs use the exact Task 5 primitive-array contract;
- legacy builder inputs and strict physics baselines remain unchanged;
- no solver-physics, CAD, OpenCascade, GUI, or baseline-data file was changed.

## TDD evidence

The Task 6 tests were authored before any production file was edited.

The initial focused RED was:

```text
python -m pytest tests/test_model_helpers.py tests/test_runner.py -q
7 failed, 7 passed in 2.78s
```

Six failures reported:

```text
TypeError: TankConfig.__init__() got an unexpected keyword argument
'geometry_path'
```

The explicit version-1 fixture failed with:

```text
AttributeError: 'TankConfig' object has no attribute 'geometry_path'
```

These failures established that configuration, validation, single-run, and
sweep geometry support was absent.

A second RED covered a corrupt archive whose metadata hash matched the corrupt
bytes. NumPy's loader error initially escaped structured validation:

```text
python -m pytest \
  tests/test_model_helpers.py::test_corrupt_geometry_archive_is_attributed_to_geometry_path \
  -q
1 failed in 2.41s
ValueError: This file contains pickled (object) data.
```

The minimal error-routing correction retained the test and produced:

```text
1 passed in 2.10s
```

The final focused GREEN was:

```text
python -m pytest tests/test_model_helpers.py tests/test_runner.py -q
15 passed in 2.57s
```

## Schema compatibility evidence

`CONFIG_SCHEMA_VERSION` and the default in-memory schema are version 2.
Serialization always replaces the payload version with 2. Loading accepts only
versions 1 and 2; a version-1 payload receives `geometry_path=""`, even if that
field is absent.

An explicit serialization/load probe reported:

```text
saved_schema_version=2
saved_geometry_path=tank.npz
legacy_schema_version=1
legacy_geometry_path=''
relative_geometry_height_ft=8.0
```

The relative-path probe changed the process working directory to the package
directory and validated `geometry_path="relative.npz"`. This preserves the
existing current-working-directory resolution convention used by other
configuration paths.

## Package validation and error routing

When `geometry_path` is non-empty, validation:

- requires the path to identify an existing `.npz` file;
- requires the adjacent `.json` metadata file;
- calls `load_geometry_package`, enforcing metadata, hash, schema, units,
  axis/gravity, array shape/dtype/contiguity, finiteness, and monotonicity;
- returns the loaded `GeometryKernel` so one direct single-case run loads the
  package once;
- surfaces missing files, `GeometryPackageError`, malformed metadata, and
  corrupt NumPy/ZIP data as `InputValidationError` issues on
  `tank.geometry_path`.

`run_sweep` validates and loads the package once, then passes the already
validated kernel into each `run_single_case` call. No repeated package load is
performed for each sweep scenario.

Tests cover a missing NPZ, missing adjacent JSON, SHA-256 mismatch, and a
hash-consistent but structurally corrupt archive.

## Builder and runtime integration

Custom builder mode uses:

```text
Volt = geometry.total_volume_ft3
Htzero = invert_monotone_volume(FillFraction * Volt, ...)
Xmlzro = saturated_liquid_density * FillFraction * Volt
```

It supplies `GeometryMode=1`, `FillFraction`, and exactly:

```text
GeomHeight
GeomVolume
GeomVolumeCoefficients
GeomAreaSamples
GeomPerimeter
GeomSidewallArea
GeomSidewallCoefficients
```

The arrays are passed directly from the validated frozen `GeometryKernel`;
builder and runner code do not copy, coerce, or mutate them. Only numeric
arrays and primitive scalars cross into `core.liqlev_simulation`.

The cylinder-package test deliberately uses legacy placeholders
`diameter_ft=25.0` and `height_ft=99.0` for a 4-foot-diameter,
8-foot-high numeric cylinder. It verifies:

```text
GeometryMode=1
FillFraction=0.5
Htzero=4.0 ft
SingleCaseResult.htank_ft=8.0 ft
```

It also verifies exact equality of all seven custom arrays. The sweep test
verifies scenario `htank=8.0`.

## Legacy and physics preservation

The legacy builder test now asserts equality of the complete key set as well
as all prior values. Therefore legacy inputs contain no new `GeometryMode`,
`FillFraction`, or `Geom*` entries.

The strict baseline checker passed without changing its data or tolerances:

```text
as203_default_high_vent         PASS rows=19
hydrogen_height_dep_mid_fill    PASS rows=29
nitrogen_custom_epsilon         PASS rows=31
Physics baseline check passed.
```

`core.py`, all saved baseline/tolerance files, and all geometry-package and
JIT physics modules have no Task 6 implementation diff.

## Verification

Task 6 focused:

```text
python -m pytest tests/test_model_helpers.py tests/test_runner.py -q
15 passed
```

Task 6 plus all geometry tests:

```text
python -m pytest tests/test_model_helpers.py tests/test_runner.py tests/geometry -q
58 passed
```

Strict physics baseline:

```text
python scripts/check_physics_baseline.py
3 cases passed; 79 total baseline rows
```

Full suite:

```text
python -m pytest -q
62 passed
```

Whitespace and worktree checks:

```text
git diff --check
exit 0
```

Final clean-status evidence is recorded in the controller handoff after this
report follow-up commit.

## Files changed

- `liqlev/model/config.py`
- `liqlev/model/builder.py`
- `liqlev/model/validation.py`
- `liqlev/io/config_json.py`
- `liqlev/runner/single.py`
- `liqlev/runner/sweep.py`
- `tests/test_model_helpers.py`
- `tests/test_runner.py`
- `.superpowers/sdd/task-6-report.md`

## Commits

Implementation:

```text
9f9dcae7efa0292c01e3f0ea77ab67a3a20dcc73
feat: load geometry packages in headless runs
```

Configured author and committer:

```text
Steven Soriano <steven.a.soriano@nasa.gov>
```

No automated-assistant authorship, co-authorship, or attribution is included.
The documentation follow-up commit containing this report cannot include its
own hash; its hash is recorded in the controller handoff.

## Self-review

- Confirmed saved files always use schema version 2 and version 1 remains
  loadable with an empty geometry path.
- Confirmed unsupported schema versions remain rejected.
- Confirmed direct and sweep runs reuse a validated package instead of loading
  it again during the same execution context.
- Confirmed custom `Volt`, `Htzero`, and `Xmlzro` use package volume and
  monotone inversion rather than legacy diameter/height placeholders.
- Confirmed the exact seven Task 5 arrays are passed without mutation or
  coercion.
- Confirmed custom reported height is package-derived in both single and sweep
  results.
- Confirmed legacy builder keys/values and strict physics baselines remain
  unchanged.
- Confirmed no OpenCascade, CAD, GUI, or Python geometry object enters the
  solver.
- Confirmed frozen dataclass semantics and existing positional call contracts
  remain intact; new arguments were appended or keyword-only.
- Confirmed all authored changes are confined to the approved worktree.

## Remaining concerns

No unresolved Task 6 contract concern remains. Runtime custom geometry retains
the approved Task 5 and design limitations: packages must satisfy the
single-region, monotone-height numeric schema, and successful numerical
consistency is not experimental validation of arbitrary-tank physics.
