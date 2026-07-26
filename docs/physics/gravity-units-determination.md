# Gravity-Units Determination for the AK1 Correlation (Finding F1)

## Status and scope

- **Determination: FINAL.** For the calibrated LIQLEV lineage, gravity enters the AK1
  correlation as a dimensionless standard-gravity level, not as ft/s^2.
- **Correction: APPLIED (plan Phase 1 / path-forward Step 2).** Solver fix + heritage
  regression test + physics-baseline regeneration landed on `wt/step2` (see
  [Applied correction](#applied-correction-plan-phase-1--measured) for measured impact and
  commit SHAs). Pre-fix adjudication tables below are retained as evidence and labeled
  **pre-fix**.
- Determination date: 2026-07-24 (review verification round); recorded and independently
  re-verified after cluster transfer: 2026-07-25; Step-2 measured update: 2026-07-25.
  Author: Steven Soriano.
- Finding of record: F1 in `docs/2026-07-24-authoritative-review-findings.md`
  (py_files branch `feature/liqlev-capability-review`).

## The question

**Pre-fix** `core.py` line 197 (`_solver_loop`, boundary-layer coefficients) evaluated
the AK1 buoyancy correlation with gravity in ft/s^2:

```python
# PRE-FIX (defective) form at core.py:197
ak1_term = 10.8 * (1 + spacl) * (1 + spacv) * ggo_ft_s2 * (rhol - rhov) / rhol
```

with `ak1 = 1.089 * ak1_term ** 0.5`. The heritage VBA->Python reference
(`LIQLEV-Heritage-Python`, `src/liqlev_heritage/solver.py:40`) converts to standard-g
units at the correlation boundary first:

```python
gravity_g = gravity_ft_s2 / STD_GRAVITY_FT_S2      # STD_GRAVITY_FT_S2 = 32.174
```

The two conventions differ by a factor of 32.174 inside the square root, i.e. by
`sqrt(32.174) = 5.672213` on AK1. Only one can be what the calibrated FORTRAN/VBA lineage
constants `(10.8, 1.089)` were tuned for. Which one is a physics determination with
safety consequences: AK1 drives bubble removal from the wall boundary layer, so an AK1
that is too large thins the computed boundary layer and **under-predicts liquid rise** —
the non-conservative direction for a liquid-level-rise analysis.

## Determination

**The dimensionless standard-gravity convention (VBA `Ggo`) is correct for this lineage.
The kernel's ft/s^2 usage at `core.py:197` is a units defect (finding F1): AK1 is too
large by exactly `sqrt(32.174) = 5.672213`.**

The determination rests on empirical reproduction of the published AS-203 flight data,
because the printed report equation cannot adjudicate (retuned constants, below) and the
original FORTRAN listing, VBA source, and workbook are unavailable.

## Empirical adjudication — reproduction of the published AS-203 figure

### Pre-fix adjudication (retained evidence)

Probe `probes/probe_as203_adjudication.py` (py_files review branch) runs the documented
AS-203 reference case — orbital gravity level 3.0e-4 g, vent rates 3.3069 / 2.2046 /
1.1023 lbm/s, per `liqlev_heritage.reference_cases` — through three configurations and
compares peak dimensionless liquid rise `dh/h0` against the published curves.
**Table below is the pre-fix record** (defective kernel ft/s^2 path vs heritage and
g0-input control).

| Vent, lbm/s | Published peak dh/h0 | Heritage (g0 convention) | Kernel pre-fix (ft/s^2, as-is) | Kernel, g0-input control (pre-fix preview) |
|---:|---:|---:|---:|---:|
| 3.3069 | 0.1350 | 0.137863 (1.021x) | 0.035833 (0.265x) | 0.138265 (1.024x) |
| 2.2046 | 0.1120 | 0.114431 (1.022x) | 0.026122 (0.233x) | 0.114564 (1.023x) |
| 1.1023 | 0.0650 | 0.066821 (1.028x) | 0.016289 (0.251x) | 0.066371 (1.021x) |

Three independent facts, one conclusion (pre-fix):

1. The **heritage g0 convention reproduces the published figure within ~2.5%** at all
   three vent rates — inside the heritage acceptance bands (`docs/VALIDATION.md`:
   0.115–0.155 at 3.3069 lbm/s, 0.090–0.135 at 2.2046 lbm/s), which allow for
   digitization and reference-plot uncertainty.
2. The **kernel's ft/s^2 convention lands at ~25% of published** — a ~4x
   under-prediction of liquid rise, far outside any digitization uncertainty.
3. The **control** — the same kernel binary with `Xggo` entered per the VBA g0
   convention (`Xggo / 32.174`) — also matches published within 2.4%. The units
   convention is therefore the **sole** defect; the modernized solver structure is
   otherwise sound.

Post-fix native kernel measurements (no g0-input hack) are in
[Applied correction §B](#b-as-203-adjudication-probesprobe_as203_adjudicationpy).

### Published-data provenance

Figure 4-14, "Dimensionless Liquid Level Increase With Time for Various Vent Rates,"
report page 4-28 (PDF page 78) of R. D. Bradshaw, *Evaluation and Application of Data
from Low-Gravity Orbital Experiment, Phase I — Final Report*, NASA CR-109847 /
GDC-DDB-70-003, 1 April 1970 (NTRS 19700017832,
`https://ntrs.nasa.gov/api/citations/19700017832/downloads/19700017832.pdf`, SHA-256
`F01404EE72A1EE8CE8CAAAEDD0812DB8A79ADF731452F7FDF7E35F12CA4313BB`). Peaks digitized
from slide 16 of Matthew E. Moran, *Liquid Rise During Venting: SpaceX Tanker LOx Header
Tank*, briefing dated 28 January 2022, which reproduces that figure; the digitized points
are stored at heritage `data/moran_2022_fig4_14.csv`.

### Mechanism confirmation (pre-fix)

Probe `probes/probe_ak1_units.py`, identical legacy-cylinder inputs at g = 0.001 through
both implementations. **Table below is the pre-fix record.**

| Row | heritage AK1 | kernel AK1 (pre-fix) | ratio |
|---:|---:|---:|---:|
| 0 | 0.22342309 | 1.26730335 | 5.67221 |
| 1 | 0.22342313 | 1.26730355 | 5.67221 |

The ratio is exactly `sqrt(32.174) = 5.672213`. Downstream on that case (pre-fix):
`BL Vap Out` +461% (row 0); final `VBL vol` and `BL thick` both −41.3%. AK3 is nearly
unchanged (ratio 0.95–0.99) because the solver iterates AK3 to close the vapour balance
and partially compensates — which is why the small low-vent case shows a −41%
boundary-layer error while the AS-203 published case shows the full ~4x rise error.

Post-fix pure-AK1 measurements are in
[Applied correction §A](#a-pure-ak1-unit-probe-probesprobe_ak1_unitspy-g--0001).

## Why the printed report equation cannot adjudicate

**Do not re-litigate this determination from Eq. 4-30.** The report's Eq. 4-30 defines

```text
K1 = (3/2) * 0.73 * [ 6*(1+S_v)*(1+S_l)*(rho_l - rho_v)*g / (pi*rho_l) ]^(1/2)
```

with dimensional `g` inside a Davies-Taylor-type rise velocity (Eq. 4-27,
`u ~ sqrt(g*delta)`) — superficially favouring the kernel's dimensional form. But the
lineage constants `(10.8, 1.089)` match the report's `(6/pi, 1.095)` under **neither**
unit convention: the effective inner constant is 5.65x the report's dimensional form and
0.176x its g0-folded form. The FORTRAN/VBA lineage was evidently **retuned** relative to
the printed equation, so the printed equation cannot determine which gravity
representation the retuned constants absorb. Only reproduction of the published AS-203
curves can — and that test is unambiguous (table above).

A related dimensional argument also fails to settle it: `xmdtbl =
2.1*AK1*dtank*delta^1.5*delta*rhov` in lbm requires `AK1 = [ft^0.5/s]`, which
`sqrt(g[ft/s^2])` supplies and `sqrt(g/g0)` does not. That argument assumes the empirical
constant `10.8` is dimensionless — unfounded for a 1960s empirical correlation, where the
constant carries units and absorbs g0. The documented VBA mapping (below) is
authoritative over dimensional inference.

## Provenance chain (durable anchor: commit `114e295`)

- The defect entered this repository at commit `114e295`
  (`114e29589edfda5b89b5a4ab39400f61645c87be`, "chore: import LIQLEV solver baseline",
  2026-07-23), which imported `core.py` from the **2026-05-13 team-share build**. The
  imported `core.py` is **byte-identical** to that build's (hash-verified during the
  2026-07-24 review; the same AK1 line sits at the import source's `core.py:160`).
- The 2026-05-13 import source **predates the 2026-07-22 heritage correction** — that is
  how this kernel imported an uncorrected `core.py` while the heritage reference already
  carried the fix.
- Both team-share copies were **deleted 2026-07-24 with user confirmation** after
  preservation was hash-verified (`core.py` byte-identical at `114e295`; remaining files
  identical modulo CRLF except a trivially adapted `.gitignore`). **Commit `114e295` is
  therefore the durable provenance anchor**; the design-spec §2.2 citation of the
  team-share path remains valid as history.
- Bug-copy census (2026-07-24): **fixed** in `LIQLEV-Heritage-Python` and
  `LIQLEV-Optimized-Python`; **still carrying the defect** — this kernel (F1), archived
  `LIQLEV-Python-mp`, `liqlev_from_cluster`, and the vendored copy inside
  `cryovent_scaling` including its 2026-07-14 cluster upload (finding F19, remediated in
  that project after plan Phase 1 lands here).

## The heritage correction record this determination mirrors

`LIQLEV-Heritage-Python` (`main` @ `f84de57`, 2026-07-22; private GitHub
`stevensoriano/LIQLEV-Heritage-Python`) documents the same fault and its repair in
`docs/GRAVITY_UNITS_CORRECTION.md`:

> "An affected Python calculation passed the ft/s^2 value directly into the correlation
> position occupied by VBA `Ggo`. That made the gravity term inside the `AK1` square
> root too large by a factor of `32.174`. Because `AK1` is proportional to the square
> root of that term, `AK1` was too large by approximately sqrt(32.174) = 5.67. The
> stronger bubble-removal term changed the boundary-layer balance and underpredicted
> liquid rise in the affected Python results."

and fixes it at the correlation boundary (`src/liqlev_heritage/solver.py:13` defines
`STD_GRAVITY_FT_S2 = 32.174`; `:40` applies `gravity_g = gravity_ft_s2 /
STD_GRAVITY_FT_S2`), keeping the external input contract in physical ft/s^2.
`docs/VBA_MAPPING.md` states the convention plainly:

> "The VBA `Ggo` value is a number in standard-gravity units at the `Ak1` expression.
> The Python input schedule instead carries `ggo_ft_s2` in ft/s^2; `boundary_layer_ak1`
> obtains the corresponding conceptual `ggo_g` as `gravity_g = gravity_ft_s2 / 32.174`
> before evaluating the same correlation."

As the heritage record states, and this determination adopts: **the correction is a unit
repair, not a change to the heritage AK1 correlation.**

## Consequence for the committed baseline

**Pre-fix:** `validation/baselines/physics_baseline.json` was generated from the
uncorrected `core.py`, so it locked in the defective behaviour as the legacy authority;
`scripts/check_physics_baseline.py` could not detect F1.

**Post-fix (Step 2):** the baseline was **deliberately regenerated** via
`scripts/write_physics_baseline.py` in its own commit (see
[Applied correction §C](#c-physics-baseline-regeneration)), recorded as a documented
physics correction — AK1 was 5.672x too large; prior results under-predicted liquid rise
— not as a regression. `scripts/check_physics_baseline.py` PASSes against the new
baseline.

## Applied correction (plan Phase 1 — measured)

**Code:** `core.py` correlation boundary now uses
`ggo_g = ggo_ft_s2 / STD_GRAVITY_FT_S2` with
`STD_GRAVITY_FT_S2 = 32.174` (module constant). Fix+test commit:
`3f86b920839414a5b4d1640c2fc205dbb7dc9a30`. Baseline regen commit:
`757dfb4fd7a0cd554450ad5fafc4bbec103552d1`.
`res[step, 24]` (`Gravity_g`) was already correct and was not changed.

Post-fix form at the correlation boundary:

```python
# POST-FIX form
if rhol != 0:
    ggo_g = ggo_ft_s2 / STD_GRAVITY_FT_S2       # 32.174
    ak1_term = 10.8 * (1 + spacl) * (1 + spacv) * ggo_g * (rhol - rhov) / rhol
else:
    ak1_term = 0.0
ak1 = 1.089 * (ak1_term ** 0.5) if ak1_term > 0 else 0.0
```

Permanent regression: `tests/geometry/test_ak1_heritage_units.py` (vendored heritage
`boundary_layer_ak1`; RED pre-fix at ratio 5.672213, GREEN post-fix at rel ≤ 1e-9).

### Pre-fix vs post-fix (same pin: py 3.13.5 / numpy 2.3.4 / numba 0.64.0)

#### A. Pure AK1 unit probe (`probes/probe_ak1_units.py`, g = 0.001)

| Row | pre-fix kernel AK1 | post-fix kernel AK1 | heritage AK1 | ratio pre/post |
|---:|---:|---:|---:|---:|
| 0 | 1.26730335 | **0.22342309** | 0.22342309 | F = 5.672213 |
| 1 | 1.26730355 | **0.22342313** | 0.22342313 | F |

Kernel/heritage AK1 ratio post-fix = **1.000000** (to printed digits).

| Quantity | pre-fix kernel | post-fix kernel | heritage |
|---|---:|---:|---:|
| BL Vap Out (row 0) | 2.2517e-05 | **4.01717542e-06** | 4.01347387e-06 |
| VBL vol (final) | 4.3274e-01 | **7.26554521e-01** | 7.37097184e-01 |
| BL thick (final) | 1.0025e-03 | **1.68313304e-03** | 1.70755151e-03 |
| AK3 (row 0 ratio k/h) | ~0.95–0.99× heritage | **1.00092×** heritage | — |

#### B. AS-203 adjudication (`probes/probe_as203_adjudication.py`)

Post-fix native kernel column (exit gate: within ~2.5% of published, matching former
control column):

| Vent, lbm/s | Published | pre-fix kernel | post-fix kernel | heritage |
|---:|---:|---:|---:|---:|
| 3.3069 | 0.1350 | 0.035833 (0.265×) | **0.138265 (1.024×)** | 0.137863 |
| 2.2046 | 0.1120 | 0.026122 (0.233×) | **0.114564 (1.023×)** | 0.114431 |
| 1.1023 | 0.0650 | 0.016289 (0.251×) | **0.066371 (1.021×)** | 0.066821 |

All three post-fix kern/pub ratios fall in **1.021–1.024**, inside the ~2.5% band and
matching the pre-fix g0-input control column exactly (to printed digits).

#### C. Physics baseline regeneration

- Command: `python scripts/write_physics_baseline.py` (pinned env); verified by
  `scripts/check_physics_baseline.py` (PASS).
- Pre-fix summary anchors (defective):
  - as203_default_high_vent `max_dh_h0` = 0.13461627071131826
  - hydrogen_height_dep_mid_fill `max_dh_h0` = 0.00018668091840282378
  - nitrogen_custom_epsilon `max_dh_h0` = 0.0024743887496442418
- Post-fix summary (regenerated):
  - as203_default_high_vent `max_dh_h0` = **0.2619626166693403** (`final_height_ft` 18.193573704708815)
  - hydrogen_height_dep_mid_fill `max_dh_h0` = **0.0003463253782541249** (`final_height_ft` 14.0948797245796)
  - nitrogen_custom_epsilon `max_dh_h0` = **0.007927313661292758** (`final_height_ft` 4.535672911475817)
- Direction of record: prior AK1 was 5.672× too large; prior results under-predicted
  liquid rise. Baseline regen is a deliberate physics correction, not a regression.

#### D. LOX case AK3 band (plan projection)

- Pre-fix projection language: AK3 × 5.672 → ~0.07–0.79 on LOX.
- Measured post-fix LOX AK3 min/max: **not measured in Step 2** (deferred to LOX flight
  case / plan Phase 3 matrix when run).

## Independent re-verification after cluster transfer (2026-07-25)

The work package was transferred to a Linux cluster (bundles SHA-256-verified) and the
adjudication was re-run under the exact pinned environment — conda-forge Python 3.13.5,
numpy 2.3.4, numba 0.64.0 — against clones at the expected tips (kernel `67f082d`,
py_files review branch `580b43b`, heritage `f84de57`):

- `probe_as203_adjudication.py` reproduced **every tabulated dh/h0 value above to all
  printed digits** (heritage 1.021–1.028x, kernel 0.233–0.265x, control 1.021–1.024x).
  The adjudication is platform-independent, not a Windows artifact.
- Heritage suite 15/15; kernel suite (without `tests/cad` — OpenCascade unavailable)
  97 passed of 98 collected; `check_physics_baseline.py` passes unchanged.
- The single kernel failure is the result-manifest floating-point equality check
  (`test_nasa_tank_result_manifest_matches_current_evidence`): recomputed
  `maximum_refinement_difference` 6.371214840947442e-05 vs committed
  6.37121484250719e-05 — a 1.56e-14 Linux-vs-Windows floating-point delta against the
  manifest's 1e-15 absolute tolerance, with interpreter and library versions identical.
  That is the manifest consistency guard detecting a platform change, not a solver or
  physics discrepancy, and it does not touch the F1 evidence. The manifest is left
  unmodified.

---

## Post-fix docs map (Phase 5 cross-reference only)

Measured F1 content above landed in Step 2; this block does **not** rewrite that
history. Readers of later campaign docs should start here for units, then follow:

| Topic | Document |
|-------|----------|
| Low-g BL physics (δ ≤ A/P, custom g=0 saturation, rate-scaled residual gate, G0 bound-only caveat, production max AK3) | [`low-gravity-boundary-layer.md`](low-gravity-boundary-layer.md) |
| LOX 43 L / SSM-3 test definition (all eight §5 assumptions; production G0–G4 matrix; hydrogen stand-in F4 deprecation) | [`../lox-vent-test-definition.md`](../lox-vent-test-definition.md) |
| On-cluster gating convention (suite green EXCEPT frozen D1; F4 xfail; F10 Solver Status; F18 headless custom geometry) | Repository [`README.md`](../../README.md) (gating / F10 / custom-geometry sections) and [`../geometry-kernel.md`](../geometry-kernel.md) |

**D1 after F1 (pointer only):** the frozen manifest value `6.371214840947442e-05`
is unchanged by this determination. Live post-F1 recomputation on the cluster is
`4.977864341653149e-05`; the red is expected by lead ruling. See the README
two-layer D1 note — not a reopening of the pre-fix platform-FP appendix above.
