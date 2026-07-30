"""Executable contract for the misleadingly-named solver columns 20/26/27.

Companion to ``docs/solver-columns.md``. The names in ``core._COL_NAMES`` are a
frozen heritage contract (pinned byte-identically by
``tests/geometry/test_f10_solver_status.py::test_f10_default_off_29_column_contract``)
and the F4 vapour mass-balance audit found three of them misleading. The
reviewer ruling was *document precisely, lock with tests, do not rename* — so
this file is the lock: every assertion characterises what the live solver
actually writes, and the readings the NAMES invite are kept as explicit negative
assertions so a future change that made a misreading true would fail here.

Reconstruction provenance: the closed mass-ledger audit ``run_vapor_balance.py``
in the companion Spyder study project (run of record
``results_vapor_balance/2026-07-29_204715/report.txt``), which balances the full
per-step ledger to 6e-13 lbm using exactly the seeding below —
``inputs["Xmlzro"]`` for the initial liquid mass and ``inputs["Delta"]`` for the
timestep.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Any

import numpy as np
import pytest

from liqlev.runner.single import SingleCaseResult, run_single_case
from validation.lox_vent_cases import (
    G3,
    ROOT,
    _saturated_vapor_density_table,
    build_lox_vent_config,
)

# core.py:122 defines this INSIDE the njit solver body, so it is a function-local
# and cannot be imported. It is mirrored here and the mirror is verified two ways
# by ``test_lbm_to_kg_is_the_heritage_six_digit_rounding``: against the core.py
# source line and against the factor recovered from the live run's own output.
LBM_TO_KG = 0.453592
NIST_LBM_TO_KG = 0.45359237

# Reconstructing ``delme`` as a difference of ~102.9 lbm liquid masses costs
# ~5e-11 relative precision on a ~4.7e-4 lbm per-step flash (catastrophic
# cancellation at the double-precision floor, not a physics disagreement). Every
# tolerance below that touches ``delme`` therefore sits at 1e-9 relative, which
# is still ~7 decades tighter than any of the misreadings this file refutes.
RECONSTRUCTION_REL = 1e-9


@lru_cache(maxsize=1)
def _short_g3_result() -> SingleCaseResult:
    """One shared 1 s G3 run — the fast pattern used by the smoke tests.

    Shared so the whole file costs one solver run and one property-table build.
    ``_primitives`` rebuilds its arrays from this result on every call, so no
    test can poison another through a shared buffer.
    """

    return run_single_case(
        build_lox_vent_config(G3, timestep_s=0.02, duration_s=1.0)
    )


def _primitives() -> dict[str, Any]:
    """Per-step ledger primitives reconstructed from the table and the inputs.

    Row ``k`` of the dataframe is the END of step ``k``, so the start of step
    ``k`` is row ``k-1`` and for step 0 it is the initial condition
    (``core.py:110`` ``vbl1 = 0.0``, ``:113`` ``t1 = tinit``, ``:116``
    ``xml1 = xmlzro``).
    """

    result = _short_g3_result()
    dataframe = result.dataframe
    inputs = result.inputs

    xmlzro = float(inputs["Xmlzro"])          # core.py:116  xml1 = xmlzro
    delta = float(inputs["Delta"])            # core.py:147  theta2 = theta1 + delta
    tinit = float(inputs["Tinit"])            # core.py:113  t1 = tinit

    xml = dataframe["Liq Mass"].to_numpy(float)          # col 3   xml2
    eps = dataframe["eps"].to_numpy(float)               # col 11  wall-contact
    vbl = dataframe["VBL vol"].to_numpy(float)           # col 13  vbl2, ft^3
    xmdtbl = dataframe["BL Vap Out"].to_numpy(float)     # col 21  film -> ullage
    temp = dataframe["Temp"].to_numpy(float)             # col 2   t2, Rankine

    xml1 = np.concatenate([[xmlzro], xml[:-1]])
    temp1 = np.concatenate([[tinit], temp[:-1]])
    delme = xml - xml1                        # core.py:162-163, NEGATIVE
    delivered = xmdtbl - (1.0 - eps) * delme  # core.py:471  mass_gen_lbm

    pt_temps, pt_rhov = _saturated_vapor_density_table(
        str(inputs["Liquid"]), float(inputs["Pfinal"]), float(inputs["Pinit"])
    )

    return {
        "xmlzro": xmlzro,
        "delta": delta,
        "xml": xml,
        "xml1": xml1,
        "eps": eps,
        "vbl": vbl,
        "xmdtbl": xmdtbl,
        "temp": temp,
        "temp1": temp1,
        "delme": delme,
        "delivered": delivered,
        # core.py:130 interpolates rho_v at t1; :129-133 at the same table.
        "rhov1": np.interp(temp1, pt_temps, pt_rhov),
        "rhov2": np.interp(temp, pt_temps, pt_rhov),
        "col20": dataframe["Vapor in BL"].to_numpy(float),
        "col26": dataframe["Vap Gen Rate (kg/s)"].to_numpy(float),
        "col27": dataframe["Total Vap Gen (kg)"].to_numpy(float),
    }


# ═══════════════════════════════════════════════════════════════════════════
#  (a) The reconstruction itself — seeding and the telescoping identity
# ═══════════════════════════════════════════════════════════════════════════
def test_liquid_mass_diffs_seed_from_xmlzro_and_telescope() -> None:
    """Contract: ``delme`` is recoverable from ``'Liq Mass'`` diffs seeded with
    ``inputs['Xmlzro']``, and the diffs telescope to the total depletion.

    This is the audit's reconstruction (run_vapor_balance.py:254, 259) and every
    assertion below depends on it, so it is locked first.
    """

    p = _primitives()

    # The seed is the initial condition, not row 0: liquid has already flashed
    # by the end of step 0, so xmlzro strictly exceeds Liq Mass[0].
    assert p["xmlzro"] > p["xml"][0]
    # Bulk flashing only removes liquid over this run (core.py:162, t2 < t1).
    assert (p["delme"] < 0.0).all()
    # Telescoping: sum of the diffs is the endpoint depletion, seed included.
    assert float(p["delme"].sum()) == pytest.approx(
        p["xml"][-1] - p["xmlzro"], rel=1e-12
    )


# ═══════════════════════════════════════════════════════════════════════════
#  (b) col 26 'Vap Gen Rate (kg/s)' — delivered-to-ullage rate, EVERY step
# ═══════════════════════════════════════════════════════════════════════════
def test_col26_rate_is_delivered_mass_on_every_step_including_the_first() -> None:
    """Contract: col 26 is ``mass_gen_lbm/delta * LBM_TO_KG`` (core.py:472) and
    is written UNCONDITIONALLY (core.py:521) — including step 0.

    This is the code-truth correction to the audit finding text, which claimed
    the rate "likewise" omits the first step. It does not; only col 27 does.
    """

    p = _primitives()
    expected_kg = p["delivered"] * LBM_TO_KG

    assert p["col26"] * p["delta"] == pytest.approx(
        expected_kg, rel=RECONSTRUCTION_REL
    )
    # Step 0 specifically: the rate is present and nonzero there.
    assert p["col26"][0] != 0.0
    assert p["col26"][0] * p["delta"] == pytest.approx(
        expected_kg[0], rel=RECONSTRUCTION_REL
    )
    # NEGATIVE: refutes the finding text's "the rate likewise omits step 1".
    assert p["delivered"][0] != 0.0


# ═══════════════════════════════════════════════════════════════════════════
#  (c) col 27 'Total Vap Gen (kg)' — cumulative delivered, SKIPPING step 0
# ═══════════════════════════════════════════════════════════════════════════
def test_col27_cumulative_omits_the_first_step() -> None:
    """Contract: col 27 accumulates only under ``if not is_first``
    (core.py:474-475), so step 0's delivered mass never enters it.
    """

    p = _primitives()
    delivered_kg = p["delivered"] * LBM_TO_KG
    cumulative_skipping_step0 = np.concatenate(
        [[0.0], np.cumsum(delivered_kg[1:])]
    )

    # The accumulator starts at its initial 0.0 (core.py:120) and is still 0.0
    # when row 0 is written — exactly, not approximately.
    assert p["col27"][0] == 0.0
    assert p["col27"] == pytest.approx(
        cumulative_skipping_step0, rel=RECONSTRUCTION_REL
    )

    # NEGATIVE, sharpened: the naive reading (cumulative over ALL steps) is
    # wrong by exactly step 0's delivered mass, at every row.
    cumulative_all_steps = np.cumsum(delivered_kg)
    assert cumulative_all_steps - p["col27"] == pytest.approx(
        np.full(p["col27"].shape, delivered_kg[0]), rel=RECONSTRUCTION_REL
    )
    # And that offset is a real physical shortfall, not float noise: >1% of the
    # reported endpoint value.
    assert abs(delivered_kg[0]) > 0.01 * abs(p["col27"][-1])


# ═══════════════════════════════════════════════════════════════════════════
#  (d) True vapour generation is liquid depletion — and is NOT col 27
# ═══════════════════════════════════════════════════════════════════════════
def test_true_vapor_generation_is_liquid_depletion_not_col27() -> None:
    """Contract: total vapour generated is ``-sum(delme)`` == initial liquid
    mass minus final ``'Liq Mass'``, and no column carries it.

    Derivation of the exact gap. With ``delivered_k = xmdtbl_k -
    (1-eps_k)*delme_k`` (core.py:471) and col 27 summing ``k >= 1`` only
    (core.py:474-475):

        true_gen - col27[-1]/LBM_TO_KG
            = delivered_0 + sum_k( eps_k*(-delme_k) - xmdtbl_k )

    The bracketed sum is the eps-share of the flash that was routed into the
    wall film and never released. Reading it through the film volume balance
    (core.py:306-311), where ``ak2*xml1*delta/rhol == -eps*delme/rhov`` exactly
    (core.py:211 with core.py:157-159, :162) and ``xmdtbl ==
    custom_exit_rate*delta*rhov`` (core.py:459), it splits as

        sum_k( eps_k*(-delme_k) - xmdtbl_k )
            = film_retention + bl_residual

    with ``film_retention = sum((vbl_k - vbl_{k-1}) * rhov(t1_k))`` (the metric
    reported by ``validation.lox_vent_cases.film_drift_metric``) and
    ``bl_residual = sum(fvbl_k * rhov(t1_k))``, the BL solver's unconverged
    remainder. Both forms are asserted; the residual term is measured rather
    than assumed small (it is ~-0.02% of true generation on this run).
    """

    p = _primitives()
    true_gen = p["xmlzro"] - p["xml"][-1]
    col27_lbm = p["col27"][-1] / LBM_TO_KG

    # Identity: depletion == -sum(delme).
    assert true_gen == pytest.approx(-float(p["delme"].sum()), rel=1e-12)
    assert true_gen > 0.0

    # NEGATIVE: col 27 is NOT total vapour generation. The shortfall is 64% of
    # true generation here — three orders above any reconstruction noise.
    assert abs(true_gen - col27_lbm) / true_gen > 0.01

    # The gap is exactly what the algebra above requires.
    film_route = float(np.sum(p["eps"] * (-p["delme"])))
    released = float(np.sum(p["xmdtbl"]))
    assert true_gen - col27_lbm == pytest.approx(
        p["delivered"][0] + film_route - released, rel=1e-10
    )

    # Same gap read as film retention plus the BL-solver residual.
    vbl1 = np.concatenate([[0.0], p["vbl"][:-1]])   # core.py:110  vbl1 = 0.0
    film_step = (p["vbl"] - vbl1) * p["rhov1"]
    film_retention = float(film_step.sum())
    bl_residual = float(
        np.sum(p["eps"] * (-p["delme"]) - p["xmdtbl"] - film_step)
    )
    assert true_gen - col27_lbm == pytest.approx(
        p["delivered"][0] + film_retention + bl_residual, rel=1e-10
    )
    # Reported, not assumed: the residual is a small but nonzero share of the gap.
    assert abs(bl_residual) < 0.01 * true_gen


# ═══════════════════════════════════════════════════════════════════════════
#  (e) col 20 'Vapor in BL' — end-of-step volume x START-of-step density
# ═══════════════════════════════════════════════════════════════════════════
def test_col20_is_end_of_step_volume_times_start_of_step_density() -> None:
    """Contract: col 20 is ``xmvbl2 = vbl2 * rhov`` (core.py:457) where ``rhov``
    was interpolated at ``t1`` (core.py:130) — a one-step-stale mixed state, not
    the film's current EOS inventory.
    """

    p = _primitives()
    stale_mixed_state = p["vbl"] * p["rhov1"]
    current_eos_inventory = p["vbl"] * p["rhov2"]

    assert p["col20"] == pytest.approx(stale_mixed_state, rel=1e-12)

    # NEGATIVE: the name's reading (current EOS inventory) is wrong, and not
    # marginally so. Temperature falls on every step of a blowdown, so the two
    # products separate by a measurable relative margin everywhere.
    assert (np.abs(p["temp"] - p["temp1"]) > 0.0).all()
    relative_gap = np.abs(p["col20"] - current_eos_inventory) / np.abs(
        current_eos_inventory
    )
    assert (relative_gap > 1e-5).all()


# ═══════════════════════════════════════════════════════════════════════════
#  (f) The kg conversion factor is the heritage 6-digit rounding
# ═══════════════════════════════════════════════════════════════════════════
def test_lbm_to_kg_is_the_heritage_six_digit_rounding() -> None:
    """Contract: cols 26/27 are in kg via ``LBM_TO_KG = 0.453592``
    (core.py:122), the kernel's 6-digit heritage rounding — NOT the exact NIST
    0.45359237.

    ``LBM_TO_KG`` is a function-local of the njit solver body and cannot be
    imported, so it is pinned twice: at the core.py source line, and by
    recovering the factor the live run actually applied.
    """

    source = (ROOT / "core.py").read_text(encoding="utf-8").splitlines()
    assert source[121].strip() == "LBM_TO_KG = 0.453592"
    assert "0.45359237" not in "\n".join(source)

    p = _primitives()
    # Summed rather than per-step: both terms of ``delivered`` are positive, so
    # the sum has no cancellation of its own and recovers the factor to ~1e-12.
    recovered = float(
        np.sum(p["col26"]) * p["delta"] / np.sum(p["delivered"])
    )
    assert recovered == pytest.approx(LBM_TO_KG, rel=1e-10)
    # NEGATIVE: the run is precise enough to tell the two constants apart, and
    # it is not using the exact one (they differ by 8.2e-7 relative).
    assert abs(recovered - NIST_LBM_TO_KG) / NIST_LBM_TO_KG > 1e-7


# ═══════════════════════════════════════════════════════════════════════════
#  Doc record
# ═══════════════════════════════════════════════════════════════════════════
def test_docs_record_the_column_contract() -> None:
    """``docs/solver-columns.md`` carries the contract, the cites and the ruling."""

    doc = (ROOT / "docs" / "solver-columns.md").read_text(encoding="utf-8")

    # All three misleading names, by their exact frozen strings.
    assert "Total Vap Gen (kg)" in doc
    assert "Vap Gen Rate (kg/s)" in doc
    assert "Vapor in BL" in doc
    # The precise semantics, in the words that carry them.
    assert "delivered" in doc.lower()
    assert "omits" in doc.lower()
    assert "start-of-step" in doc.lower()
    # Write-site and definition cites for the three corrections.
    assert "core.py:471" in doc      # mass_gen_lbm definition
    assert "core.py:474" in doc      # the `if not is_first` skip
    assert "core.py:521" in doc      # the unconditional rate write
    assert "core.py:457" in doc      # xmvbl2 = vbl2 * rhov(t1)
    # Heritage rounding note.
    assert "0.453592" in doc
    assert "0.45359237" in doc
    # The ruling and its three reasons.
    assert "no rename" in doc.lower()
    assert "test_f10_solver_status.py" in doc
    # Cross-references: the F10 opt-in column, the film-drift metric, provenance.
    assert "Solver Status" in doc
    assert "docs/vent-study.md" in doc
    assert "2026-07-29_204715" in doc
    # Names and numbers are frozen.
    assert "frozen" in doc.lower()
    # Every public column is tabulated, by name.
    from core import _COL_NAMES

    for name in _COL_NAMES:
        assert f"`{name}`" in doc
