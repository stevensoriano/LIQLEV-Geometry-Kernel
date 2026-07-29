# LOX Vent Test Definition — 43 L, 40→35 psia, SSM-3 Settling

**Date:** 2026-07-24 (definition); measured production anchors updated 2026-07-25–26  
**Status:** production acceptance case on the geometry kernel  
**Kernel path:** `validation/lox_vent_cases.py`  
**Production matrix manifest:** `validation/results/lox_vent_manifest.json`  
  (commit `7d6a744`, `solver_describe` = `f8fedae`, rate-scaled BL residual gate)  
**Supersedes:** the hydrogen `hydrogen_height_dep_mid_fill` stand-in used by
`validation/custom_geometry_cases.py` and the NASA tank refinement path as the
**ullage-acceptance story** (finding F4 — see §8). Geometry refinement evidence
for the NASA package remains a separate frozen-manifest topic (D1).

**Ported from:** review-branch
`docs/2026-07-24-lox-vent-test-definition.md` (py_files
`feature/liqlev-capability-review`), with pre-fix numerical anchors **replaced**
by post-F1-fix / post-BL-gate production measurements and explicit supersession
notes.

**Companion physics notes:**
[low-gravity boundary layer](physics/low-gravity-boundary-layer.md),
[gravity-units determination (F1)](physics/gravity-units-determination.md).

All geometry and CoolProp property anchors below were verified against the
committed geometry package. Definition-era probes (review branch):
`probes/probe_lox_and_zerog.py`, `probes/probe_lox_ak3_regime.py`. Production
matrix: campaign evidence `15x_s5b_primary_matrix.*`,
`15x_s5b_dt_plateau.*`, committed manifest `7d6a744`.

---

## 1. System

| Quantity | Value | Provenance |
|---|---|---|
| Fluid | **LOX** (CoolProp `Oxygen`) | engineer |
| Tank | committed NASA fluid domain `nhq01-m21a-0201` | geometry package |
| Tank total volume | **3.464700 ft³ = 98.109 L** | geometry package |
| Tank total height | **1.806770 ft** | geometry package |
| Liquid load | **43 L** at saturation | engineer |
| **Fill fraction** | **0.438286** | 43 L / 98.109 L |
| **Initial liquid height** | **0.820444 ft** (45.41% of tank height) | inverted from volume PCHIP |
| Initial LOX mass | ~102.9 lbm (46.7 kg) | ρ_l × 43 L |
| Ullage gas | GOX at saturation temperature (autogenous) | engineer |
| Thermal state | saturated — conservative max-rise bound | engineer |

### LOX saturation properties (CoolProp)

| Pressure | T_sat | ρ_l | ρ_v | σ |
|---:|---:|---:|---:|---:|
| 40 psia | 100.990 K = **181.783 R** | 1085.61 kg/m³ | 11.253 kg/m³ | 10.4976 mN/m |
| 35 psia | 99.394 K = 178.908 R | 1094.10 kg/m³ | 9.941 kg/m³ | 10.8820 mN/m |

`build_property_table("Oxygen", 35.0, 40.0)` verified working: 400 points,
T range 175.308–185.383 R.

---

## 2. Vent

| Quantity | Value | Provenance |
|---|---|---|
| Initial pressure | **40 psia** | engineer |
| Final pressure | **35 psia** | engineer |
| Vent mass flow | **11.89 g/s = 0.026212963 lbm/s** (constant, average) | engineer |
| Maximum vent-valve open time | **60 s** | engineer |
| Termination | pressure ≤ 35 psia **OR** 60 s elapsed | engineer |
| Vent valve Cv | 1.4, all tubing 1/2" | engineer — not yet used; see §6 |

`while p1 > pfinal` plus the `theta2 >= tvmdot[-1] - delta` break gives the dual
termination natively.

### Blowdown duration — supersession of the pre-fix cross-check

| Era | t_end (dual dt) | Evidence |
|-----|----------------:|----------|
| **Pre-F1-fix** (definition-era) | **59.10 s** at both dt = 0.10 and 0.02 | Review definition; engineer ~59 s analysis |
| **Post-F1-fix, stock BL gate** (manifest `74a0d7d`) | **58.70 s** at both dt = 0.10 and 0.02 (p_end 34.992 / 34.9997 psia) | Evidence `9x_s4c_comparison.txt` |
| **Post rate-scaled gate** (manifest `7d6a744`, production) | G3 **t_end = 58.58 s** at primary dt = 0.02 (p_end 34.9988 psia) | `validation/results/lox_vent_manifest.json` `results.G3` |

**Standing criterion:** blowdown completes in the **~59 s** class and remains
timestep-independent on the pressure path. The −0.40 s post-F1 shift vs 59.10 s
is denom-driven (F1 gravity-units correction) and **holds** the ~59 s criterion;
the further −0.12 s under the rate-scaled gate vs stock 58.70 s is recorded, not
a failure of the dual-termination contract.

---

## 3. Settling and gravity

### SSM-3 profile

| Quantity | Value |
|---|---|
| Profile | SSM-3, "Pulsed Thrust SSM, Mid-Range" |
| Thruster on/off | **120 ms ON / 3280 ms OFF** |
| Period | 3400 ms |
| **Duty cycle** | **3.5294%** |
| Pulses during a 60 s vent | 17.65 |
| Settling time before vent | 78 s |
| Total profile duration | 78 + 60 + 10 = **148 s** |

### Vehicle and thrust — resolved

From the Eta Space chart *"Model Outputs: Pulsed Thrust Settling"* (02/20/2024
Rev F, slide 19):

- **Vehicle mass: 280 kg**
- **Thrust: 1.60 N — both RCS engines combined** (not per engine)
- Chart acceleration plot peaks at ~0.0057 m/s²

Confirmed: `a = 1.60 N / 280 kg = 5.714286e-3 m/s²`, matching the plotted 0.0057.

**This resolves the earlier ambiguity.** The table header "@ 1.6 N Effective
Thrust" and the chart both mean 1.60 N *total*. The verbal description "two
thrusters, each with 1.6 N NET axially" (which would imply 3.2 N) is **not**
supported by the chart and is carried only as a sensitivity case.

### Bond-number cross-check

Using the committed fluid-domain half-width R = 0.28094 m and LOX properties at
40 psia:

```
Bo = (rho_l - rho_v) * a * R^2 / sigma
   = 1074.357 * 5.714286e-3 * 0.28094^2 / 0.0104976
   = 46.16
```

against the engineer's stated **Bo = 44** — agreement to **4.9%**. This confirms
the radius (not diameter) convention and that Bo = 44 describes the
**during-pulse peak** acceleration, not a duty-averaged value. Residual
difference is attributable to their choice of σ, T_sat or characteristic radius.

### Gravity matrix — steady-state constants

Per direction from the engineer: **treat gravity as steady state, and bracket
with a maximum and a minimum level.** The pulsed time history is deliberately
*not* fed to the solver (see §5).

| ID | g (standard g) | a (m/s²) | Interpretation |
|---|---:|---:|---|
| **G0** | **0** | 0 | absolute maximum-swelling bound; requires F2/F3 and the rate-scaled gate |
| **G1** | **2.056571e-05** | 2.016807e-04 | **minimum** — 1.60 N duty-averaged over SSM-3 |
| G2 | 4.113141e-05 | 4.033613e-04 | 3.2 N duty-averaged (sensitivity only) |
| **G3** | **5.826950e-04** | 5.714286e-03 | **nominal maximum** — 1.60 N peak, 280 kg |
| G4 | 1.165390e-03 | 1.1428571e-02 | 3.2 N peak (sensitivity only) |

g-column conversion uses g₀ = 9.80665 m/s² exactly. (A prior revision of this
table used 9.8084 — a uniform 0.017% slip, corrected in the verification
round; the m/s² column was always exact.)

G3 is the nominal prediction. G1 and G3 are the requested min/max bracket. G0 is
the absolute bound. G2 and G4 exist only to cover the unresolved verbal "3.2 N"
statement and can be dropped once that is settled.

**F1 units convention (applied):** gravity enters AK1 as a **dimensionless
standard-g level** at the correlation boundary. The solver still stores `Xggo`
in ft/s² for I/O; `prepare_gravity` multiplies the value in g by
`G_TO_FT_S2 = 32.174`, and the F1 fix divides by standard g inside AK1. See
[gravity-units-determination.md](physics/gravity-units-determination.md).

---

## 4. Timestep

Both dt = 0.02 s and dt = 0.01 s divide the SSM-3 pulse timing exactly
(120 ms and 3400 ms), which matters only if the pulsed profile is ever used.
For steady-state gravity the blowdown is already timestep-independent between
0.1 s and 0.02 s on the **pressure** path.

| dt | steps/pulse | steps for 60 s | note |
|---:|---:|---:|---|
| 0.10 | 1.20 | 600 | blowdown already converged |
| 0.02 | 6.00 | 3000 | **recommended primary** |
| 0.01 | 12.00 | 6000 | confirmation |
| 0.005 | 24.00 | 12000 | convergence check |

**Rise dt-convergence** is a separate question from blowdown. Under the stock
absolute BL residual gate, rise collapsed ~linearly with dt (false quasi-steady
film). Under the rate-scaled gate (`f8fedae`), G0/G1/G3 rises are
**0.2%-class plateaus** across 0.02 / 0.01 / 0.005
(evidence `15x_s5b_dt_plateau.txt`). Primary reporting dt remains **0.02 s**.

---

## 5. Assumptions and limitations — must be stated in any report

These eight items are **verbatim-faithful** to the review-branch definition §5
(wording preserved; only cross-references modernized where the original named
“finding F4” still applies).

1. **Steady-state gravity.** SSM-3 is a *pulsed* profile at 3.53% duty. The model
   is run at constant gravity levels bracketing the peak and duty-averaged
   values. **The model therefore does not capture interface dynamics between
   pulses**, including slosh, interface re-formation, or the settling transient
   itself. Bo = 44 and We = 41.2 describe the *settled* interface and are the
   justification for treating the interface as flat and settled during the vent.

2. **No boundary-layer time constant.** LIQLEV's boundary layer is quasi-steady —
   an instantaneous function of the current state. It has no relaxation
   timescale, so it cannot respond correctly to a 3.4 s pulse cycle. This is the
   technical reason the pulsed history is not fed directly: doing so would make
   the model jump instantaneously between a thin boundary layer (120 ms) and the
   saturated maximum (3280 ms), which is representative of the *acceleration* but
   **not** more physical.

3. **g = 0 is a bound, not a prediction.** AK1 is a buoyancy-driven rise-velocity
   coefficient. At zero gravity the real physics becomes diffusion-, Marangoni-
   and nucleation-dominated, none of which LIQLEV models. G0 is a conservative
   upper bound on swelling only.

4. **Saturated initial state** is assumed as a conservative maximum-rise bound;
   any subcool margin would reduce predicted rise.

5. **Constant average vent rate.** 11.89 g/s is used as a constant. The real
   valve (Cv 1.4, 1/2" tubing) will produce a pressure-dependent flow that
   decreases as the tank blows down.

6. **Heritage ullage formulation.** The accumulated vapour mass is never
   reconciled against the ullage equation of state (finding F4). This is
   authentic heritage behaviour and is preserved. Its magnitude must be reported
   per run, not assumed small.

7. **Numerical, not experimental, validation.** All geometry and solver
   acceptance criteria are internal-consistency checks. Custom-tank predictions
   remain report-based numerical predictions requiring later experimental or
   higher-fidelity correlation.

8. **No internal hardware.** The fluid domain is baffle-free; PMD, spray bar,
   fasteners and viewport hardware are excluded, and the two axial ports are
   capped at the wet-side closure planes. Thermal properties of the omitted
   closure hardware are not represented, which is why phase-one heat transfer
   uses sidewall area only (`A_w,side`), not total wetted area.

---

## 6. Deferred / not yet in scope

- Cv-derived pressure-dependent vent flow (Cv 1.4, 1/2" tubing) instead of a
  constant 11.89 g/s.
- The other settling profiles SSM-1, SSM-2, SSM-4. Note the source table carries
  a handwritten correction: **SSM-4 is 250/7040 ms, not 260/7040**. SSM-3 is
  unaffected.
- Concept-of-operations comparison (settle-then-vent vs concurrent vs vent
  cycling). The agreed baseline is concurrent: settle 78 s, vent 60 s while still
  thrusting, then 10 s post-vent thrusting.
- Helium pressurization. Current assumption is autogenous GOX only.

**2026-07-29 (`feature/pulsed-gravity`), additive note.** A pulsed-profile
capability now exists in the vent-study harness: an optional `GP` row feeds the
SSM-3 square-wave acceleration history directly to the solver. It is a
model-consistency experiment only — the agreed steady-state bracketing method
of §3 is unchanged, and §5 item 2 applies to anything it produces.

---

## 7. Production G0–G4 matrix (rate-scaled gate) — supersedes pre-fix anchors

### 7.1 Pre-fix anchors (historical only)

Review definition §7 and pre-F1 probe rows (not the definitive G1–G4 g levels):

- AK3 band **0.012–0.14**; F1 projected ~0.07–0.79.
- G3 pre-F1-fix level rise 0.82060 → 0.84013 ft → **5.95 mm (+2.4%)**.
- Blowdown **59.10 s** at dt 0.10 and 0.02.

### 7.2 Stock-gate matrix (superseded — do not adopt low-g rises)

Manifest `74a0d7d` under absolute residual gate `abs(fvbl) ≤ 0.001·vbl2` produced
false quasi-steady films at low g (dt-collapse, G1>G0 inversion, ullage blow-up).
**Historical record only** (supervisor report Steps 2–4 table). Do **not** adopt
G0–G2 rises from that matrix; G3/G4 were nearer true QS at dt = 0.02.

### 7.3 Production matrix — manifest `7d6a744` (authoritative)

Primary dt = 0.02 s, dual termination, rate-scaled custom BL gate
(`f8fedae`). File: `validation/results/lox_vent_manifest.json`.

| Row | g (std) | t_end (s) | rise (mm) | dh/h0 | max AK3 | ullage closure | vs 5% | Conv Failed |
|-----|--------:|----------:|----------:|------:|--------:|---------------:|-------|------------:|
| **G0** | 0 | 57.18 | **184.433** | 0.737 | 0.4421 | **5.247%** | **EXCEEDS (reported)** | 0 |
| **G1** | 2.056571e-05 | 57.72 | **108.589** | 0.434 | 0.3375 | 2.216% | under | 0 |
| G2 | 4.113141e-05 | 57.88 | 90.265 | 0.361 | 0.2915 | 1.895% | under | 0 |
| **G3** | 5.826950e-04 | 58.58 | **31.304** | 0.125 | 0.0992 | 1.237% | under | 0 |
| G4 | 1.165390e-03 | 58.72 | 22.860 | 0.092 | 0.0698 | 1.146% | under | 0 |

Source fields: `results.G*.dh_ft` (× 304.8 → mm), `ullage_closure_max_relative`,
`max_ak3`, `t_end_s`, `conv_failed_total`. Cross-check evidence
`15x_s5b_primary_matrix.txt`, `15x_s5b_manifest_write.txt`, supervisor report
production table.

**Ordering:** monotone **G0 ≥ G1 ≥ G2 ≥ G3 ≥ G4**.  
**G0 caveat:** rise is a **conservative upper bound**, not a validated zero-g
prediction (§5.3), and carries the **5.247% ullage-closure EXCEEDS** report.

**dt-plateau (evidence `15x_s5b_dt_plateau.txt`):**

| Row | vs dt=0.01 | vs dt=0.005 |
|-----|-----------:|------------:|
| G0 | +0.015% | +0.022% |
| G1 | +0.047% | +0.045% |
| G3 | +0.109% | +0.192% |

All **0.2%-class**.

**Operational quote (campaign):** minimum-gravity rise **G1 = 108.6 mm**,
nominal-maximum **G3 = 31.3 mm**, absolute zero-g bound **G0 = 184.4 mm** (with
the bound-only and 5.25% closure caveats).

### 7.4 Post-fix AK3 band

Production max AK3 across G0–G4 is **0.0698–0.4421** (manifest `7d6a744`), still
well below the ~1e2 F2 stiffness threshold. See
[low-gravity-boundary-layer.md](physics/low-gravity-boundary-layer.md).

---

## 8. Hydrogen stand-in: DEPRECATED (F4-quarantined)

The NASA hydrogen acceptance stand-in (`hydrogen_height_dep_mid_fill` / fill-0.90
ullage path) is **DEPRECATED as the ullage-acceptance case** for custom geometry.

| Item | Record |
|------|--------|
| Finding | **F4** — heritage ullage mass not reconciled to EOS; stand-in is ~3000× volume-mismatched vs the NASA fluid domain |
| Quarantine | `tests/geometry/test_ullage_guards.py::test_nasa_hydrogen_standin_ullage_passes_strict_guard` is **`@pytest.mark.xfail(strict=True)`** naming F4 (commit `5a6df96`) |
| Tripwire intent | Asserts the *recovered* state (strict ullage guard passes, no negative ullage). **Fails today → xfail satisfied.** XPASS later forces marker removal when bookkeeping is reconciled or the stand-in is retired |
| Attribution (supervisor-confirmed) | Pathology is the stand-in’s own geometry/bookkeeping mismatch, not campaign Step-3/4/5 code. Pre-campaign main `@67f082d` already showed large closure; F1 deepened fill-0.90 (~121% closure, negative ullage) |
| **Replacement** | **This LOX 43 L / SSM-3 case** is the production custom-geometry ullage and gravity-matrix acceptance path. Strict ullage guards remain **STRICT on LOX/custom** |

Pre-existing NASA **refinement** tests are not the F4 quarantine target and remain
separate from this deprecation note (D1 frozen manifest is a different ruling —
see README / geometry-kernel gating convention).

---

## 9. How to re-run / cite

- Module: `validation/lox_vent_cases.py`
- Manifest writer refuses dirty trees (F8); production hash-bound fields include
  geometry NPZ, fluid STEP, harness module sha256, `solver_describe`, pin versions.
- Any report that quotes rise numbers **must** restate §5 items 1–3 and the
  per-row ullage-closure column (especially G0 EXCEEDS).
