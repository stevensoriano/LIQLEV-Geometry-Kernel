# Low-Gravity Boundary Layer (Findings F2/F3 and the BL Residual Gate)

## Status and scope

- **Status: record of applied physics.** F2 (δ cap), F3 (custom-mode zero-g
  saturation path), and the custom-mode BL residual-gate fix are landed on the
  Step-6 documentation branch (`wt/step6-docs`) and earlier campaign branches.
- **Scope:** custom-geometry mode only for zero-g saturation and the rate-scaled
  residual gate. Legacy cylinder mode stays baseline-locked at g = 0 (lead
  decision); its absolute residual gate is byte-retained so
  `validation/baselines/physics_baseline.json` remains reproducible.
- **Companion docs:** [gravity-units determination (F1)](gravity-units-determination.md),
  [legacy series discrepancy](legacy-series-discrepancy.md),
  [LOX/SSM-3 test definition](../lox-vent-test-definition.md).
- **Numbers authority:** supervisor campaign report
  `/gen/home/sasorian/cfdmate/data/liqlev-f1-verify-a1/report.md` and the durable
  evidence tree under `…/evidence/`. Committed LOX matrix:
  `validation/results/lox_vent_manifest.json` at commit `7d6a744`
  (`solver_describe` = `f8fedae`).

## Physical bound: δ ≤ A(h)/P(h) (F2, commit `21521cc`)

The report boundary-layer ODE (NTRS 19700017832, Eq. 4-33 and custom-mode
generalization) advances thickness `δ` along the wall. At a free surface of
cross-section area `A(h)` and contact-line perimeter `P(h)`, the section cannot
hold a film thicker than the local hydraulic half-width:

```text
δ ≤ A(h) / P(h)
```

when `A` and `P` are positive. That bound is geometric, not a numerical
convenience: a larger film would leave the fluid domain of the section.

Commit `21521cc` (`fix(F2): cap boundary-layer thickness at the physical section
limit…`) applies the cap in `_boundary_layer_derivatives`, clamps the final
`delta_top`, and guards post-substep `q` overshoot to the local saturated exit
state `(2/3)·P·(A/P)^1.5`. It also treats `ak3 == 0` as the exact empty-film
solution and rejects only `ak3 < 0`. Post-integration self-consistency returns
status 1 with NaNs if `V_BL > V(top)·(1+1e-6)` or `δ_top > (A/P)·(1+1e-6)`.

### Analytic AK3 → ∞ saturation limit

As the vapour-balance coefficient AK3 → ∞, the film fills the section and the
integrator must land on the saturation state (δ → A/P, physical `V_BL ≤ V`), not
on a non-physical overshoot. Measured RED → GREEN on the LOX NASA tank at the
43 L liquid height (`h = 0.820444 ft`), campaign evidence
`5x_s3a_red_garbage_lox_geom.txt` / `5x_s3a_green_garbage_lox_geom.txt`:

| AK3 | Pre-F2 (status, δ/(A/P), V_BL/V) | Post-F2 (status, δ/(A/P), V_BL/V) |
|----:|----------------------------------|----------------------------------|
| 1e2 | 0, 0.999245, 0.992151 | 0, **0.999248**, 0.992151 |
| 1e3 | 0, 0.0, **3.23** | 0, **1.000000**, **0.996627** |
| 1e4 | 0, 0.0, **15.0** | 0, **1.000000**, **0.996644** |
| 1e6 | 0, 0.0, **323.2** | 0, **1.000000**, **0.996644** |

Extreme AK3 therefore resolves to **physical saturation (status 0)**, not forced
rejection of a still-garbage film. The self-consistency net remains for residual
overshoot classes; the δ cap is the primary fix.

## Custom-mode g = 0 takes the saturation path; legacy stays locked (F3)

Commit `7dbd5b8` (`fix(F3): custom-mode zero-gravity boundary layer takes the
saturation limit`):

| Path | At `ak1 == 0` | Intent |
|------|---------------|--------|
| **Custom** (`geometry_mode == 1`) | Seed AK3 positive (`hldak3` if > 0, else bootstrap 1.0); bracket floor at expansion | Close vapour balance with zero exit flow → film grows toward A/P |
| **Legacy** | `ak3 = 0.0` exactly (byte-locked) | Preserve `physics_baseline.json` and heritage cylinder results |

Lead decision: **legacy g = 0 remains baseline-locked.** Zero-g physics that is
intended as a conservative bound is exercised only on the custom path.

### Measured observation: legacy g = 0 is practically non-terminating

Campaign Step-2/3 diagnostics recorded that the **legacy** path at exactly
`g = 0` is **practically non-terminating** (>90 s wall time for a short 4-step
schedule), not merely “meaningless but finite.” That observation sharpens F3: the
legacy lock is not a validated zero-g prediction path; it is a reproducibility
rail. Custom-mode multi-step g = 0 on a cylinder completes in sub-second class
with physical bounds (evidence `5x_s3a_f3_zerog_key_numbers.txt`: 14 rows,
0.64 s wall, Conv Failed = 0, δ growing toward A/P).

## Measured AK3 ranges (LOX 43 L vent)

### Pre-F1-fix (review record)

From the review-branch LOX test definition
(`docs/2026-07-24-lox-vent-test-definition.md` §7) and probe
`probes/probe_lox_ak3_regime.py` (Bond-span gravity levels, dt = 0.02 s, all
`Conv Failed = 0`):

| gravity (probe) | AK3 max |
|---:|---:|
| 9.80e-06 | 0.140 |
| 1.96e-05 | 0.097 |
| 2.78e-04 | 0.025 |
| 5.55e-04 | 0.018 |
| 1.11e-03 | 0.012 |

**Pre-fix band: AK3 ≈ 0.012–0.14**, three to four orders of magnitude below the
~1e2 integrator-stiffness threshold (F2). Plan language projected that correcting
F1 would raise AK3 by `√32.174 ≈ 5.672` → roughly **0.07–0.79**, still inside the
stable regime.

### Post-fix production matrix (manifest `7d6a744`, rate-scaled gate)

Committed file `validation/results/lox_vent_manifest.json`
(`solver_describe` = `f8fedae`). Per-row **max AK3** (full precision in file;
four decimals below):

| Row | g (std) | max AK3 | Source field |
|-----|--------:|--------:|--------------|
| G0 | 0 | **0.4421** | `results.G0.max_ak3` = 0.44212902895860823 |
| G1 | 2.056571e-05 | **0.3375** | `results.G1.max_ak3` = 0.3375343785667821 |
| G2 | 4.113141e-05 | **0.2915** | `results.G2.max_ak3` = 0.29149092873142984 |
| G3 | 5.826950e-04 | **0.0992** | `results.G3.max_ak3` = 0.09916301940436156 |
| G4 | 1.165390e-03 | **0.0698** | `results.G4.max_ak3` = 0.06978242076524582 |

Cross-check: durable evidence `15x_s5b_primary_matrix.json` /
`13x_supervisor_lowg_rate_scaled_gate.txt` match the committed manifest digits.

**Note on stock-gate AK3:** under the superseded absolute residual gate
(manifest `74a0d7d` / evidence `8x_s4b_G*.txt`) G1 max AK3 was **0.2997** and G4
**0.0705**. Those rows are historical only; production numbers are the
rate-scaled table above.

## BL residual-gate correction (commit `f8fedae`)

### Discrete closure

At each step the AK3 bracket drives the residual of the discrete vapour balance

```text
fvbl = vbl2 - vbl1 - S·dt + E·dt
```

toward zero. `S·dt` and `E·dt` are the source and exit volume increments over the
step; both scale as `O(dt)`.

### Why the absolute gate produced false quasi-steady states

The stock custom-mode (and still-current legacy) gate accepted

```text
|fvbl| ≤ 0.001 · vbl2
```

Once `|S·dt|` fit inside `0.001 · vbl2`, the iteration could declare “converged”
with **E ≠ S** — a **false quasi-steady** film inventory frozen at the previous
step’s value. Campaign measurements under that gate (pre-`f8fedae`, LOX matrix
at dt = 0.02 unless noted):

| Failure class | Evidence | Measured effect |
|---------------|----------|-----------------|
| **dt-collapse of rise** | `10x_i1_*`, `12x_i3_intervention_4dt_g3.txt` | G3 rise ~31.8 mm at dt 0.02 → ~9.0 mm at dt 0.005 (~linear with dt) |
| **G1 > G0 inversion** | stock matrix `74a0d7d` / report | G1 93.32 mm > G0 79.17 mm — zero-g not the upper bound |
| **Ullage-closure blow-up** | same | G0 closure **61.4%**, G1 14.0%, G2 22.8% vs 5% report band |

Causation was proven by intervention (rate-scaled gate on stock code) before the
fix was approved; see evidence `12x` / `13x` and supervisor report gate-fix
section.

### Rate-scaled custom-mode gate

Commit `f8fedae` changes **custom mode only** (`geometry_mode == 1`, `rhol ≠ 0`):

```text
fvbl_tol = 0.001 · (|S·dt| + |E·dt|)
accept if |fvbl| ≤ fvbl_tol
```

Legacy keeps `0.001 · vbl2` byte-for-byte. If `rhol == 0`, custom falls back to
the absolute gate (no division path in the njit kernel).

### Production consequences after the fix

| Metric | Stock gate (`74a0d7d`) | Rate-scaled (`7d6a744`) |
|--------|------------------------:|-------------------------:|
| G0 rise | 79.17 mm | **184.43 mm** |
| G1 rise | 93.32 mm | **108.59 mm** |
| G3 rise | 31.82 mm | **31.30 mm** |
| Ordering | G1 > G0 inverted | **G0 ≥ G1 ≥ G2 ≥ G3 ≥ G4** |
| G0 ullage closure | 61.4% | **5.247%** (still reported EXCEEDS 5%) |

### dt-plateau proof (0.2%-class)

Evidence `15x_s5b_dt_plateau.txt` / `.json` — primary reference dt = 0.02 s:

| Row | rise @ 0.02 mm | vs 0.01 | vs 0.005 |
|-----|---------------:|--------:|---------:|
| G0 | 184.433 | +0.015% | +0.022% |
| G1 | 108.589 | +0.047% | +0.045% |
| G3 | 31.304 | +0.109% | +0.192% |

All checked rows lie in the **0.2%-class** across 0.02 / 0.01 / 0.005 s. Rise is
therefore a converged physical result of the quasi-steady model at the agreed
primary dt, not an artifact of residual freeze.

## g = 0 is a conservative upper bound on swelling, not a validated prediction

**Caveat class (must appear with any G0 number):** AK1 is a buoyancy-driven
rise-velocity coefficient. At zero gravity the real physics becomes diffusion-,
Marangoni-, and nucleation-dominated — none of which LIQLEV models. **G0 is a
conservative upper bound on swelling only**, outside the model’s validated
buoyancy-driven regime. Exact wording class is restated verbatim in
[LOX/SSM-3 test definition §5 item 3](../lox-vent-test-definition.md).

### G0 measured anchor (production)

From `validation/results/lox_vent_manifest.json` row G0 (and
`15x_s5b_primary_matrix.txt`):

| Quantity | Value | Notes |
|----------|------:|-------|
| Rise | **184.43 mm** (`dh_ft` = 0.605095 → 184.433 mm) | Absolute bound under the rate-scaled gate |
| Ullage closure max relative | **5.247%** | **EXCEEDS** the 5% reporting band — **reported, not relaxed** |
| Final V_BL | 1.5156 ft³ | ~full 43 L liquid volume scale — model saturation limit |
| max AK3 | 0.4421 | Well below F2 stiffness threshold |
| Conv Failed | 0 | Physical completion |

Use G0 only with the bound-only caveat and the 5.25%-class closure caveat. Prefer
G1 (duty-averaged minimum) and G3 (nominal peak) for operational brackets; see
the LOX test definition for the full gravity matrix.

## Summary for consumers of results

1. **δ ≤ A/P** is physics (`21521cc`); AK3 → ∞ saturates physically.
2. **Custom g = 0** takes the saturation path (`7dbd5b8`); **legacy g = 0** is
   baseline-locked and is practically non-terminating if forced.
3. **Absolute residual gate** falsified low-g quasi-steady states; **rate-scaled
   gate** (`f8fedae`) restores monotone rise, dt-plateau, and honest closures.
4. **Production max AK3** on LOX G0–G4 is **0.07–0.44**, not the pre-fix
   0.012–0.14 band and not the pre-F1 projection ceiling alone.
5. **G0 = 184.4 mm with 5.25% ullage-closure EXCEEDS** is a **bound**, not a
   flight prediction.

## Commit map (this topic)

| SHA | Role |
|-----|------|
| `21521cc` | F2 δ ≤ A/P + AK3=0 exact + self-consistency |
| `7dbd5b8` | F3 custom zero-g saturation; legacy locked |
| `f8fedae` | Custom-mode rate-scaled BL residual gate |
| `7d6a744` | LOX G0–G4 production manifest under rate-scaled gate |
| `d7f2a4f` | F1 gravity-units fix (raises AK3 × √32.174) — see gravity-units doc |
