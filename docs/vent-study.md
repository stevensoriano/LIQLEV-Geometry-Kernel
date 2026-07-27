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
