# Solver column contract

The authoritative reading of the solver result table.

`liqlev_simulation` returns the historical **29-column** DataFrame whose names
are `core._COL_NAMES` (`core.py:581-589`). Those names and every number under
them are **frozen**: they are pinned byte-identically by
`tests/geometry/test_f10_solver_status.py::test_f10_default_off_29_column_contract`
and consumed by name string across the GUI, the Qt app, the viz layer, the
validation harness and external study scripts. This document does not change
them — it is the precise statement of *what they mean*, and three of them do not
mean what they are called (see [Misleading names](#misleading-names--precise-semantics)).

Every row is the **END of a step**: row `k` holds the state after step `k`
advanced by `Delta`. The start of step `k` is therefore row `k-1`, and for step 0
it is the initial condition (`core.py:110` `vbl1 = 0.0`, `:113` `t1 = tinit`,
`:116` `xml1 = xmlzro`). Several columns mix the two ends of a step; those are
called out individually below.

Units are the kernel's heritage **British** set (`"Units": "British"`) —
ft, lbm, psia, Rankine, seconds. Columns 25/26/27 (`Superheat`,
`Vap Gen Rate (kg/s)`, `Total Vap Gen (kg)`) are the SI outliers, in K and kg.

## The 29 public columns

Line numbers in the last column are the **write site** in `core.py`; where the
value is *computed* elsewhere in the step that line is given too.

| # | Name | Units | Meaning | `core.py` |
|--:|------|-------|---------|-----------|
| 0 | `Time` | s | `theta2` — end-of-step time, `theta1 + delta` (`:147`) | `:495` |
| 1 | `Press` | psia | `p2` — end-of-step ullage pressure, `p1 + delp` (`:158`) | `:496` |
| 2 | `Temp` | °R | `t2` — end-of-step saturation temperature, `t1 + dtdps*delp` (`:159`) | `:497` |
| 3 | `Liq Mass` | lbm | `xml2` — end-of-step bulk liquid mass, `xml1 + delme` (`:163`) | `:498` |
| 4 | `Ullage Mass` | lbm | `xmvap2` — the **bookkeeping** ullage accumulator, `xmvap1 - (1-eps)*delme - delmv + xmdtbl` (`:478`). Not an EOS quantity; compare col 22 | `:499` |
| 5 | `Height` | ft | `h2` — end-of-step liquid level. Custom geometry inverts the occupied volume (`:414-420`); legacy integrates `dh/dt` (`:431-435`) | `:500` |
| 6 | `dh/dt` | ft/s | Level rise rate over the step (`:428` custom, `:431` legacy) | `:501` |
| 7 | `delh` | ft | `h2 - h1` — level change over the step (`:440`) | `:502` |
| 8 | `Hratio` | — | `(h2 - htzero)/htzero` — fractional rise against the initial level (`:441`) | `:503` |
| 9 | `dP/dtha` | psia/s | `dpdtha` — depressurization rate, `-vmdot/denom` (`:156`) | `:504` |
| 10 | `delp` | psia | `delp = dpdtha*delta` — pressure change over the step (`:157`) | `:505` |
| 11 | `eps` | — | Wall-contact fraction. Does **not** add a vapour source: it *routes* a share of the same bulk flash into the wall film (`:167-194`) | `:506` |
| 12 | `beta` | — | Heritage BL similarity parameter, `vbl2*rhov/xmlzro/beta_denom`, with `rhov` at **start**-of-step (`:455`) | `:507` |
| 13 | `VBL vol` | ft³ | `vbl2` — end-of-step wall-film (boundary-layer) volume (`:269` custom, `:301` legacy). Carried across steps as a *volume* (`:567`) | `:508` |
| 14 | `BL thick` | ft | `delblz` — converged boundary-layer thickness (`:268` custom, `:274-295` legacy Newton–Raphson) | `:509` |
| 15 | `AK3` | heritage† | The **converged** BL root `ak3` returned by the secant/bracket solver (`:245-405`) | `:510` |
| 16 | `Vent Rate` | lbm/s | `vmdot` interpolated at the step **midpoint** `thetav`, not at either endpoint (`:148-149`) | `:511` |
| 17 | `AK1` | heritage† | Buoyancy/gravity BL coefficient (`:207-210`) | `:512` |
| 18 | `AK2` | heritage† | Flash-source BL coefficient, `-eps*cs*rhol*dtdps*dpdtha/rhov/hfg` (`:211`) | `:513` |
| 19 | `AK2/AK1` | heritage† | The step's **initial** `ak3` estimate (`:216`), recomputed at the write site. This is the seed, *not* the converged root — that is col 15 | `:514` |
| 20 | `Vapor in BL` | lbm | `xmvbl2 = vbl2 * rhov` — end-of-step film volume × **start**-of-step density. **Misleading name**, see below (`:457`) | `:515` |
| 21 | `BL Vap Out` | lbm | `xmdtbl` — vapour mass released from the film to the ullage this step (`:459` custom, `:461-468` legacy) | `:516` |
| 22 | `Ullage from Calc` | lbm | `xmvap3 = volgas * rhov` — the EOS ullage mass, evaluated from **start**-of-step state (`:139-141`) | `:517` |
| 23 | `Conv Iterations` | count | `nconv` — BL solver iterations consumed on this step | `:518` |
| 24 | `Gravity_g` | g | Gravity at the step midpoint in standard-g, `ggo_ft_s2/32.174` (`:199-202`) | `:519` |
| 25 | `Superheat` | K | `(t2 - t_sat)/1.8` with `t_sat` back-projected along the pressure slope, `t2 - delp*dtdps`; reduces to `delp*dtdps/1.8` (`:481-486`) | `:520` |
| 26 | `Vap Gen Rate (kg/s)` | kg/s | Rate of vapour **delivered to the ullage**. **Misleading name**, see below (`:471-472`) | `:521` |
| 27 | `Total Vap Gen (kg)` | kg | Cumulative vapour **delivered to the ullage**, **omitting step 0**. **Misleading name**, see below (`:474-475`) | `:522` |
| 28 | `Conv Failed` | flag | `1.0` when the BL solver exhausted its iterations or a custom-mode failure fired, else `0.0`. The OR of the F10 failure codes | `:523-527` |

† `AK1`/`AK2`/`AK3`/`AK2/AK1` are heritage groups. The numeric constants that
surround them (`10.8`, `1.089`, `0.375`, `2.1`, `8.0`, `4.0`) absorb the
dimensional bookkeeping, so these four are **not** dimensionally consistent
across their uses and no unit is claimed for them here. `AK1`'s gravity
normalization is finding F1 — see `tests/geometry/test_ak1_heritage_units.py`.

### The opt-in 30th column

The solver always writes an internal 30th column, `Solver Status`, with
per-step status codes; it is sliced off unless explicitly enabled, so the public
contract stays 29 columns. Status code table, the three ways to opt in, and the
default-off contract are in the **F10 section of `README.md`**. Write site:
`core.py:554` (codes assigned `:532-553`).

## Misleading names — precise semantics

Found by the F4 vapour mass-balance audit. Run of record:
`LIQLEV_vent_study_spyder/results_vapor_balance/2026-07-29_204715/report.txt`
(the companion Spyder study project, outside this repository), whose ledger
closes to 6e-13 lbm. Locked by
`tests/geometry/test_column_semantics.py`.

### col 27 `Total Vap Gen (kg)` — it is not total vapour generation

It is the **cumulative vapour delivered to the ullage**, and it **omits the
first timestep**.

Per step the kernel forms

```
mass_gen_lbm = xmdtbl - (1.0 - eps) * delme        # core.py:471
```

— the film's release to the ullage (`xmdtbl`, col 21) plus the share of the bulk
flash that went straight to the ullage rather than into the film
(`-(1-eps)*delme`). That is a *delivery*, not a generation: the `eps` share of
the flash is real vapour that has been created and is sitting in the wall film,
and it is absent from this total.

The accumulator is then guarded:

```
if not is_first:
    cumulative_vap_kg += mass_gen_lbm * LBM_TO_KG   # core.py:474-475
```

so `col27[0] == 0.0` exactly and step 0's delivery never enters the running sum
at any later row either. The offset is permanent, not a startup transient.

### col 26 `Vap Gen Rate (kg/s)` — same quantity as a rate, on every step

```
vap_gen_rate = (mass_gen_lbm / delta) * LBM_TO_KG   # core.py:472
```

It is the col-27 quantity as a rate, so the "delivered, not generated"
correction applies identically.

**Code-truth correction.** The original audit finding text stated that the rate
"likewise" omits the first step. It does not. The write at `core.py:521` is
**unconditional** — it sits outside the `if not is_first` guard, which encloses
only the cumulative accumulator. `col26[0]` carries step 0's delivery rate
normally. Only col 27 (`core.py:474-475`) skips a step.

### col 20 `Vapor in BL` — a one-step-stale mixed state

```
xmvbl2 = vbl2 * rhov                                # core.py:457
```

`vbl2` is the **end**-of-step film volume, but `rhov` was interpolated at `t1`
(`core.py:130`), the **start-of-step** temperature. The product is therefore a
mixed state one step stale in density — not the film's current EOS inventory,
which would be `vbl2 * rhov(t2)`. On the 1 s G3 short run the two differ by
~4e-5 relative on every step; over a full blowdown the density drop is far
larger.

This is the same staleness family as the reported wall-film mass drift: the
solver carries the film as a *volume* across steps (`core.py:567`) while
saturated vapour density falls. The harness-side surfacing of that effect is the
`film_drift_lbm` metric documented in
[`docs/vent-study.md`](vent-study.md) under "Wall-film mass drift — reported,
not gated"; col 20 is where the same mixed-state convention shows up inside the
table itself.

### True total vapour generation

Total vapour generated is the **liquid-mass depletion**, and **no column carries
it**:

```
true_gen = -sum(delme) = xmlzro - LiqMass[-1]
```

`delme = xml1*cs*(t2-t1)/hfg` (`core.py:162`) and `xml2 = xml1 + delme`
(`core.py:163`), so the per-step diffs telescope exactly. Reconstruct `delme`
from the `Liq Mass` column seeded with `inputs["Xmlzro"]`; the timestep is
`inputs["Delta"]`.

The exact gap against col 27 is

```
true_gen - col27[-1]/LBM_TO_KG
    = delivered_0 + sum_k( eps_k*(-delme_k) - xmdtbl_k )
    = delivered_0 + film_retention + bl_residual
```

where `delivered_0` is the step-0 delivery col 27 omits, `film_retention =
sum((vbl_k - vbl_{k-1}) * rhov(t1_k))` is the flash routed into the film and not
released (the `film_retention_lbm` metric), and `bl_residual = sum(fvbl_k *
rhov(t1_k))` is the BL solver's unconverged remainder. The second form follows
from the film volume balance `core.py:306-311`, in which
`ak2*xml1*delta/rhol == -eps*delme/rhov` exactly (`core.py:211` with `:157-159`
and `:162`) and `xmdtbl == custom_exit_rate*delta*rhov` (`core.py:459`). On the
1 s G3 short run col 27 reports **35.8%** of true generation, and the identity
closes to 1.5e-14 lbm.

### `LBM_TO_KG` — heritage 6-digit rounding

Columns 26 and 27 are converted to kg with

```
LBM_TO_KG = 0.453592                                # core.py:122
```

which is the kernel's **6-digit heritage rounding**, not the exact NIST factor
`0.45359237`. The two differ by 8.2e-7 relative — negligible against the
semantic errors above, but it means cols 26/27 are not interchangeable with an
exactly-converted value. `LBM_TO_KG` is a function-local of the njit solver body
and cannot be imported.

## Why the names were not changed

Reviewer ruling on the F4 follow-up: **document precisely, lock the semantics
with tests, and apply no rename.** Three reasons.

1. **The names are a deliberate heritage contract.**
   `tests/geometry/test_f10_solver_status.py` pins all 29 name strings
   byte-identically (`test_f10_default_off_29_column_contract`). Renaming means
   breaking a contract that was written to be unbreakable.
2. **Every consumer keys by name string.** `gui.py:2572,2574`,
   `liqlev/ui_qt/app.py:1030`, `liqlev/viz/datasets.py:98`, `tests/test_viz.py:34`,
   `tests/geometry/test_conservation.py:11,189`, plus the external Spyder study
   scripts and the deferred-F19 downstream vendor copies. A rename breaks
   external consumers silently at runtime, for zero numerical benefit.
3. **Heritage discipline.** Outputs stay byte-identical; the truth is documented
   alongside them rather than encoded into a new interface. This document is
   that record, and `tests/geometry/test_column_semantics.py` is its executable
   half.
