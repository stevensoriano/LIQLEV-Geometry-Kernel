# Task 6 Report: Configuration and Headless Geometry Packages

## Status

Complete under the approved Task 6 boundary:

- saved configurations now use schema version 2 and include
  `TankConfig.geometry_path`;
- explicit schema-version-1 files still load with an empty geometry path;
- headless single, sweep, and Monte Carlo runs load validated numeric geometry
  packages;
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

## Independent review corrections

Independent review found four Task 6 integration and input-boundary issues.
Each regression group was written and run against committed head
`fb87766910a35b95f0512b18d252ea4d3e485220` before corrective production
edits.

### Monte Carlo geometry propagation

`run_monte_carlo` validated and loaded a package but discarded the returned
kernel. Its real `build_inputs` calls therefore stayed in legacy mode.

The RED solver-boundary capture was:

```text
python -m pytest tests/test_monte_carlo.py -q
1 failed in 2.43s
KeyError: 'FillFraction'
```

The correction captures the one validated kernel and passes it to every
Monte Carlo `build_inputs` call. The GREEN result was:

```text
1 passed in 3.26s
```

For three deterministic samples, the test verifies:

- one package load for the complete Monte Carlo execution;
- `GeometryMode=1` and package `Volt` for every sample;
- each sampled fill equals the input `FillFraction`;
- all seven arrays match the package exactly;
- all samples reuse the same seven array objects;
- the exposed geometry height endpoint is `8.0 ft`.

### Public single-run source of truth

The initial Task 6 reuse parameter made `run_single_case(..., geometry=...)`
publicly injectable. An empty-path config could activate custom mode, and a
config pointing at package A could execute a caller-supplied package B.

The RED source-of-truth group was:

```text
2 failed, 1 passed in 2.84s
```

The public signature still exposed `geometry`, and the package-B call
completed instead of raising `TypeError`. The pre-existing two-fill sweep
already demonstrated one package load.

The correction restores the entire pre-Task-6 public signature. Public
`run_single_case` now loads only from `config.tank.geometry_path` and delegates
to private `_run_single_case_prevalidated`. `run_sweep` validates once and
uses that private helper for package reuse.

GREEN:

```text
3 passed in 3.45s
```

The two-fill sweep loads package A once, executes two fills, and reports
`htank=8.0` for both scenarios. Public package-B injection is rejected by the
function signature.

### Package and path error boundary

The package loader leaked raw parsing, metadata, and archive exceptions.
Validation also constructed `Path` before checking the configured value type.

Direct package RED:

```text
python -m pytest tests/geometry/test_package.py -q
9 failed, 3 passed in 0.47s
```

Raw `TypeError`, `AttributeError`, `JSONDecodeError`, `KeyError`, and NumPy
`ValueError` exceptions escaped.

Structured-validation RED:

```text
4 failed, 4 passed in 2.58s
```

Numeric paths and scalar/null/numeric metadata cases escaped
`InputValidationError`.

`load_geometry_package` now raises `GeometryPackageError` at its boundary for
invalid path values, unreadable or corrupt JSON, non-object JSON roots,
missing/null/numeric hashes, malformed metadata fields, unreadable NPZ files,
and corrupt or incomplete archives. Wrapped parser/archive exceptions retain
their causes with `raise ... from exc`; no `BaseException` catch is used.
Validation checks `geometry_path` type before `Path` and catches only
`GeometryPackageError`.

GREEN:

```text
package boundary:        12 passed in 0.16s
structured validation:    8 passed in 2.33s
```

### Exact configuration schema versions

Configuration loading used `int(...)`, silently normalizing unsupported
values.

RED:

```text
4 failed, 3 passed in 2.43s
```

`2.9`, `True`, and `"2"` were accepted, while `None` leaked a `TypeError`.

The loader now requires `type(version) is int` and membership in `(1, 2)`.
Versions `0`, `3`, `2.9`, `True`, `"2"`, and `None` are all rejected with
`ValueError`. An absent version retains the established version-1 default and
empty geometry path.

GREEN:

```text
7 passed in 2.24s
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
validated kernel into the private prevalidated single-case helper. No repeated
package load is performed for each sweep scenario. Monte Carlo follows the
same one-load execution-context contract.

Tests cover a missing NPZ, missing adjacent JSON, SHA-256 mismatch, and a
hash-consistent but structurally corrupt archive, plus malformed path, JSON
root, hash, metadata-field, and archive types.

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

`core.py`, all saved baseline/tolerance files, and all JIT physics modules have
no corrective-review diff. The geometry-package change is confined to loader
error normalization and does not alter validated numeric arrays or physics.

## Verification

Corrective Task 6 focused:

```text
python -m pytest tests/test_model_helpers.py tests/test_runner.py \
  tests/test_monte_carlo.py tests/geometry/test_package.py -q
46 passed
```

Corrective Task 6 plus all geometry tests:

```text
python -m pytest tests/test_model_helpers.py tests/test_runner.py \
  tests/test_monte_carlo.py tests/geometry -q
86 passed
```

Schema/path adversarial probes:

```text
27 passed
```

Strict physics baseline:

```text
python scripts/check_physics_baseline.py
3 cases passed; 79 total baseline rows
```

Full suite:

```text
python -m pytest -q
90 passed
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
- `liqlev/geometry/package.py`
- `liqlev/runner/monte_carlo.py`
- `liqlev/runner/single.py`
- `liqlev/runner/sweep.py`
- `tests/geometry/test_package.py`
- `tests/test_model_helpers.py`
- `tests/test_monte_carlo.py`
- `tests/test_runner.py`
- `.superpowers/sdd/task-6-report.md`

## Commits

Implementation:

```text
9f9dcae7efa0292c01e3f0ea77ab67a3a20dcc73
feat: load geometry packages in headless runs
```

Initial Task 6 report:

```text
fb87766910a35b95f0512b18d252ea4d3e485220
docs: record task 6 verification
```

Independent-review correction:

```text
73563bc6e34e6dd67eb56443f455ca657456caeb
fix: harden headless geometry integration
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
- Confirmed Monte Carlo loads once per execution and reuses identical array
  objects for every sample while preserving each sampled fill fraction.
- Confirmed `config.tank.geometry_path` is the sole public custom-geometry
  source; the reuse helper is private and caller-supplied package replacement
  is unavailable.
- Confirmed schema versions are not coerced and only exact integer versions 1
  and 2 are accepted, while the absent-version compatibility default remains
  version 1.
- Confirmed malformed path, metadata root/hash/field, JSON, and NPZ cases
  become `GeometryPackageError` and then field-level `InputValidationError`.
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
