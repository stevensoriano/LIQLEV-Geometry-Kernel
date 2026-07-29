# Vent study re-run interface (`liqlev.analysis.vent_study`)

**Module:** `liqlev/analysis/vent_study.py`
**Preset:** `configs/lox_43L_40to35_ssm3.json`
**Tests:** `tests/analysis/test_vent_study.py`
**Case of record:** [lox-vent-test-definition.md](lox-vent-test-definition.md)

A config-driven front end for re-running the LOX vent acceptance case (or a
variation of it) without editing Python. It is a **driver, not a model**: every
physical quantity is produced by machinery that already carries the campaign's
acceptance evidence.

| Concern | Owner (reused unchanged) |
|---|---|
| Case structure (epsilon mode, ramp contract, dual termination, tank constants) | `validation.lox_vent_cases.build_lox_vent_config` |
| Execution | `liqlev.runner.single.run_single_case` |
| Summary, ullage-closure metric, physicality + classifications | `validation.lox_vent_cases.summarize_lox_vent_result` |
| 5% ullage-closure gate | `validation.lox_vent_cases.ULLAGE_CLOSURE_RELATIVE_TOLERANCE` |
| Dirty-tree refusal and `git describe` stamp (F8) | `validation.lox_vent_cases` guards, delegated |

What is new here is the **input surface**, the automatic dt-plateau check, the
study manifest, and a report that always quotes §5 of the test definition.

---

## 1. CLI

```
python -m liqlev.analysis.vent_study --config configs/lox_43L_40to35_ssm3.json
python -m liqlev.analysis.vent_study --config <path> --output-dir <dir>
```

| Flag | Effect |
|---|---|
| `--config` | study config JSON (required) |
| `--output-dir` | manifest/report directory (default `validation/results/vent_study/`) |
| `--gravity ROW` | run only these rows; repeatable or comma-separated (`--gravity G0,G3`) |
| `--duration-s` / `--timestep-s` | override the config for a bounded re-run |
| `--skip-dt-plateau` | skip the dt/2 and dt/4 re-runs |
| `--no-manifest` | print the report only; write nothing |

Exit codes: `0` success, `2` config error (message on stderr, no traceback),
`3` manifest refused because the git worktree is dirty.

Outputs written to the output directory: `<name>_manifest.json` and
`<name>_report.txt`. The report also goes to stdout; progress lines go to
stderr.

Importable equivalent (e.g. from Spyder):

```python
from liqlev.analysis.vent_study import load_study_config, run_study, format_report

result = run_study(load_study_config("configs/lox_43L_40to35_ssm3.json"))
print(format_report(result))
```

---

## 2. Config schema

```json
{
  "schema": "liqlev.analysis.vent_study",
  "version": 1,
  "name": "lox_43L_40to35_ssm3",
  "description": "free text, echoed into the report and manifest",
  "fluid": "Oxygen",
  "geometry_package": "geometry/tables/nhq01-m21a-0201_LIQLEV_GEOMETRY.npz",
  "fill":      {"fraction": 0.438286},
  "pressure":  {"initial_psia": 40.0, "final_psia": 35.0},
  "vent_rate": {"lbm_per_s": 0.026212963},
  "duration_s": 60.0,
  "timestep_s": 0.02,
  "gravity": {"levels": {"G0": 0.0, "G3": 5.826950e-04}},
  "dt_plateau": {"enabled": true, "row": "G3"}
}
```

Rules — violations raise `VentStudyConfigError` before any solver runs:

- `fluid` — one of `Nitrogen`, `Oxygen`, `Hydrogen`, `Methane`.
- `geometry_package` — optional; relative paths resolve against the **repo
  root**, and the default is the committed NPZ.
- `fill` — **exactly one** of `liters` or `fraction`. Litres are divided by the
  geometry package's own total volume (`total_volume_ft3 × 28.316846592`).
- `pressure` — `initial_psia` must exceed `final_psia`.
- `vent_rate` — **exactly one** of `g_per_s` or `lbm_per_s`
  (`g/s ÷ 453.59237`).
- `timestep_s` must be smaller than `duration_s`.
- `gravity` — **exactly one** of `levels` or `spaceflight` (below).
- `dt_plateau.row` must name a configured level; default is `G3` when present,
  otherwise the highest level. `enabled` defaults to `true`.
- Unknown keys anywhere are rejected rather than ignored, so a typo cannot
  silently change what is run.

### Fill and rate exactness

The preset states `fraction` and `lbm_per_s` directly because those are the
values pinned by the case definition. The convenience forms are close but not
bit-identical (43 L resolves to 0.4382863392 against 0.438286; 11.89 g/s to
0.0262129629738 against 0.026212963), so a config meant to reproduce a
committed manifest should state the pinned form.

### Spaceflight block

```json
"gravity": {
  "spaceflight": {
    "thrust_N": 1.60,
    "vehicle_mass_kg": 280.0,
    "on_ms": 120.0,
    "period_ms": 3400.0
  }
}
```

`duty_cycle` may be given directly instead of the `on_ms`/`period_ms` pair —
exactly one of the two forms. The matrix is derived as in §3 of the test
definition:

| Row | Value | Meaning |
|---|---|---|
| `G0` | `0.0` | zero-g bound (bound, not a prediction — §5.3) |
| `G1` | `peak × duty` | duty-averaged minimum |
| `G3` | `thrust_N / (vehicle_mass_kg × 9.80665)` | peak, nominal maximum |

SSM-3 (1.60 N, 280 kg, 120/3400 ms) reproduces the pinned levels exactly at the
precision the definition carries: `G3 = 5.826950e-04`, `G1 = 2.056571e-05`.
The even rows `G2`/`G4` of that table are the sensitivity-thrust variants and
are only produced by an explicit `levels` block — which is why the shipped
preset uses `levels` and reproduces all five production rows.

The block also accepts the optional `pulsed_row` / `phase_ms` pair, which adds
the separate `GP` row described in [§6](#6-pulsed-profile-row-gp). It never
changes the matrix above.

---

## 3. Results table, dt-plateau, assumptions

Every report prints the per-row table (rise in mm, `dh/h0`, max AK3, ullage
closure with its `under`/`EXCEEDS` verdict against the 5% gate, convergence
failures, physicality), the dt-plateau block, and then **§5 of
`docs/lox-vent-test-definition.md` quoted verbatim** — read from the committed
document at run time, never restated in code, so the report cannot drift from
the definition of record.

The dt-plateau check re-runs the nominal-maximum row at `dt/2` and `dt/4` and
reports `100 × (rise(dt_i) − rise(dt)) / rise(dt)` plus the maximum absolute
deviation.

---

## 4. Manifest provenance

`<name>_manifest.json` carries the same F8 discipline as the production LOX
manifest and refuses to be written from a dirty worktree:

- sha256 of the geometry NPZ, the fluid STEP, `validation/lox_vent_cases.py`,
  `liqlev/analysis/vent_study.py`, `docs/lox-vent-test-definition.md` and the
  config file itself;
- `solver_describe` = `git describe --dirty --always`;
- env pins (python / numpy / numba);
- the resolved case definition (including the derived gravity matrix and the
  spaceflight block when used);
- per-row summaries serialized by the same function as the production manifest,
  plus `rise_mm`, and the dt-plateau block.

---

## 5. Worked example — the production G3 row

Run from a clean worktree:

```
python -m liqlev.analysis.vent_study --config configs/lox_43L_40to35_ssm3.json \
    --gravity G3 --duration-s 60 --output-dir <evidence dir>
```

```
Row          g (std)  t_end (s)   rise (mm)    dh/h0   max AK3    ullage    vs 5%  ConvFail  physical
G3      5.826950e-04      58.58      31.304   0.1252    0.0992    1.237%    under         0  yes

dt-plateau check — row G3 (g = 5.826950e-04), base dt = 0.02 s
  dt = 0.02      rise =     31.304 mm   (base)
  dt = 0.01      rise =     31.338 mm   +0.110%
  dt = 0.005     rise =     31.364 mm   +0.193%
  max |delta| vs base = 0.193%
```

`rise_mm = 31.303825173274227` against `validation/results/lox_vent_manifest.json`
`results.G3.dh_ft × 304.8 = 31.3038251732747` — agreement to 5e-13 mm, and the
plateau matches the campaign evidence (+0.109% / +0.192%) to rounding.

Any report quoting these rises must restate §5 items 1–3 and the per-row
ullage-closure column (§9 of the test definition) — which is why the printed
report carries them by construction.

---

## 6. Pulsed-profile row (`GP`)

**Constant-g bracketing remains the agreed engineering method.** Gravity is
treated as steady state and bracketed with a maximum and a minimum level per
direction from the engineer (§3 of the test definition); `G1`–`G3` is what gets
reported. The `GP` row exists to answer a different, narrower question — *what
does this model do when the acceleration history itself is fed in* — and is a
**model-consistency experiment only**. It is not a validated prediction and
does not supersede anything.

### Enabling it

```json
"gravity": {
  "spaceflight": {
    "thrust_N": 1.60,
    "vehicle_mass_kg": 280.0,
    "on_ms": 120.0,
    "period_ms": 3400.0,
    "pulsed_row": true,
    "phase_ms": 0.0
  }
}
```

| Key | Type | Default | Meaning |
|---|---|---|---|
| `pulsed_row` | bool | `false` | run the extra `GP` row |
| `phase_ms` | float | `0.0` | offset of the pulse train, `0 ≤ phase_ms < period_ms` |

Rules, each rejected before any solver runs:

- `pulsed_row` requires the `on_ms`/`period_ms` form. `duty_cycle` alone does
  not define pulse timing, so it cannot synthesize a train.
- `phase_ms` is only meaningful with `pulsed_row: true`; supplying it otherwise
  is an error rather than a silent no-op.
- `phase_ms` must satisfy `0 ≤ phase_ms < period_ms`.
- `timestep_s` must put at least **two steps inside the ON interval**
  (`on_s ≥ 2·dt`); anything coarser does not represent the square wave at all
  and is refused. This is re-checked after `--timestep-s`.
- If `timestep_s` does not divide `on_s` or `period_s` to 1e-9 s, the run
  proceeds but the report carries an **aliasing warning** line. The production
  timesteps 0.02 / 0.01 / 0.005 s all divide 120 ms and 3400 ms exactly.

`GP` is selectable like any other row: `--gravity G3,GP` runs both, and
`--gravity G3` alone turns the pulsed row **off**. `dt_plateau.row` may name
`GP`; the synthesized table depends only on the pulse timing, so halving the
timestep merely resamples it.

### Waveform and phase convention

`pulsed_gravity_table` synthesizes the infinite pulse train

```
g(t) = peak_g   when  ((t - phase_s) mod period_s) < on_s
g(t) = 0.0      otherwise
```

sampled on `[0, duration_s]`, with `xggo` in ft/s² (`peak_g × 32.174`). Square
edges are emitted as **double breakpoints** — at a transition `t_e` the table
carries `(t_e, previous)` then `(t_e + 1e-9, new)` — the same riser convention
`gui.py` uses when extending a CSV gravity profile, so the solver's linear
interpolation between table points reproduces vertical edges. The first node is
`(0, g(0))`, the last is exactly `(duration_s, g(duration_s))`, times are
strictly increasing and both arrays are contiguous float64.

Phase shifts the train later in time. A `phase_ms` in
`(period_ms − on_ms, period_ms)` therefore leaves a pulse **already ON at
t = 0** — that is the intended reading of the definition, not an artifact. At
`on_ms == period_ms` (100% duty) the table collapses to the two-point constant
table bit for bit, which is the tripwire the test suite uses to prove the
profile path adds no arithmetic drift.

### What the solver actually integrates

The `GP` case is built by the same `simulation_config_for` path as every
constant row, at `gravity_g = peak_g × duty_cycle` — the true time mean of the
train. **That value is bookkeeping and reporting only.** The synthesized table
is handed to `run_single_case(..., prepared_gravity=PreparedGravity(...))`,
which overrides gravity entirely (`liqlev/runner/single.py`:
`gravity = prepared_gravity or prepare_gravity(config)`). The builder echoes
`Tggo`/`Xggo` into `result.inputs`, and that echo is the provenance truth for
the profile the solver saw; an integration test asserts it matches the
synthesized table exactly. `core.py` is untouched — its per-step gravity
interpolation already supports `g(t)`.

### Reporting and provenance

When a `GP` row is present the report gains a `mode` column (`constant` /
`pulsed`), a pulse block showing the ON/period, phase, peak and mean levels,
the table point count and its sha256, plus any aliasing warnings — and then
prints the **mandatory pulsed-profile caveat verbatim**. The caveat is a module
constant asserted by test, so it cannot be edited away from a report: it states
that `GP` feeds the square wave directly, that this departs from the §3 method
of record, that §5 item 2 applies with full force (the quasi-steady boundary
layer has no relaxation timescale, so the model jumps between a thin film and
the saturated maximum inside each 3.4 s cycle — representative of the
acceleration history but not more physical than the `G1`–`G3` bracket), and
that `GP` is not a validated prediction.

A constant-g study is unaffected: with no pulsed row the report has no `mode`
column and no caveat, exactly as before.

The manifest records `gravity_mode` per row and, for `GP`, a `pulse`
descriptor: `thrust_N`, `vehicle_mass_kg`, `on_ms`, `period_ms`, `phase_ms`,
`peak_g`, `mean_g`, `n_table_points` and `table_sha256` (sha256 over
`tggo.tobytes() + xggo.tobytes()`, uppercase hex, matching the digest style of
the other manifest hashes).
