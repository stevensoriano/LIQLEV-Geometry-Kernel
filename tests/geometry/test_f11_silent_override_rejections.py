"""F11: silent-override rejections (plan Phase 4.7).

(a) ``xmlzro_override`` + custom geometry previously ignored the override
    silently in ``build_inputs``; now raises ``ValueError``.
(b) ``Neps > 0`` + ``GeometryMode == 1`` previously discarded the epsilon
    schedule silently in the custom-geometry branch order; now
    ``liqlev_simulation`` raises ``ValueError``.

Four-criteria map
-----------------
| Guard                                      | C1 bad-case fires | C2 recovered state | C3 real path | C4 reachability |
|--------------------------------------------|-------------------|--------------------|--------------|-----------------|
| xmlzro_override + geometry raises          | ValueError msg    | no silent ignore   | build_inputs | ``test_f11a_xmlzro_override_with_geometry_raises`` |
| geometry without override still builds     | n/a (recovery)    | mass from fill·V   | build_inputs | ``test_f11a_geometry_without_override_builds`` |
| override without geometry still applies    | n/a (recovery)    | Xmlzro==override   | build_inputs | ``test_f11a_override_without_geometry_applies`` |
| Neps>0 + GeometryMode==1 raises            | ValueError msg    | schedule not used  | liqlev_simulation | ``test_f11b_neps_with_geometry_mode_raises`` |
| Neps>0 legacy mode still runs              | n/a (recovery)    | GeometryMode 0 ok  | liqlev_simulation | ``test_f11b_neps_legacy_mode_runs`` |
| GeometryMode==1 with Neps==0 still runs    | n/a (recovery)    | custom eps path    | liqlev_simulation / runner | ``test_f11b_geometry_mode_neps_zero_runs`` |

Pre-flight: no existing authorized test feeds the newly-rejected combinations
(physics_cases overrides are legacy-only; custom-geometry cases use
height_dep → Neps=0). Suite edit not required for green except D1.
"""

from __future__ import annotations

import numpy as np
import pytest

from core import liqlev_simulation
from liqlev.geometry.fixtures import cylinder_kernel
from liqlev.model.builder import build_inputs
from liqlev.runner.single import run_single_case
from validation.lox_vent_cases import G3, build_lox_vent_config


def _legacy_build_kwargs(**overrides):
    base = dict(
        fluid="Hydrogen",
        pinit_psia=19.5,
        pfinal_psia=13.8,
        dtank=4.0,
        htank=8.0,
        fill_fraction=0.5,
        duration=2.0,
        delta_t=0.5,
        vent_rate=0.1,
        neps=0,
        teps=None,
        xeps=None,
        ramp_duration=2.0,
        ramp_target_factor=1.0,
        nggo=2,
        tggo=np.array([0.0, 2.0]),
        # Non-zero: legacy g=0 can stall the Newton BL path; recovery only needs
        # GeometryMode absent + Neps>0, not a zero-g stress case.
        xggo=np.array([0.01, 0.01]),
    )
    base.update(overrides)
    return base


# ---------------------------------------------------------------------------
# F11(a) — xmlzro_override + geometry
# ---------------------------------------------------------------------------


def test_f11a_xmlzro_override_with_geometry_raises() -> None:
    """C1/C3/C4: conflicting override+geometry raises at the real builder entry."""

    geom = cylinder_kernel(4.0, 8.0, node_count=65)
    with pytest.raises(ValueError, match="xmlzro_override cannot be combined"):
        build_inputs(
            **_legacy_build_kwargs(
                geometry=geom,
                xmlzro_override=100.0,
            )
        )


def test_f11a_geometry_without_override_builds() -> None:
    """C2 recovery: geometry alone still builds (mass from fill·volume·ρ)."""

    geom = cylinder_kernel(4.0, 8.0, node_count=65)
    inputs = build_inputs(**_legacy_build_kwargs(geometry=geom))
    assert inputs["GeometryMode"] == 1
    assert inputs["Xmlzro"] > 0.0
    assert "xmlzro_override" not in {k.lower() for k in inputs}


def test_f11a_override_without_geometry_applies() -> None:
    """C2 recovery: override on legacy (no geometry) still sets Xmlzro."""

    inputs = build_inputs(
        **_legacy_build_kwargs(xmlzro_override=16300.0, geometry=None)
    )
    assert "GeometryMode" not in inputs
    assert inputs["Xmlzro"] == pytest.approx(16300.0)


# ---------------------------------------------------------------------------
# F11(b) — Neps > 0 + GeometryMode == 1
# ---------------------------------------------------------------------------


def test_f11b_neps_with_geometry_mode_raises() -> None:
    """C1/C3/C4: Neps>0 + GeometryMode==1 raises at liqlev_simulation entry."""

    config = build_lox_vent_config(G3, timestep_s=0.02, duration_s=0.04)
    # Real path: build a legal GeometryMode==1 input dict, then inject Neps>0.
    result = run_single_case(config)
    inputs = dict(result.inputs)
    assert inputs["GeometryMode"] == 1
    inputs["Neps"] = 2
    inputs["Teps"] = np.array([0.0, 1.0])
    inputs["Xeps"] = np.array([0.4, 0.4])

    with pytest.raises(ValueError, match="Neps > 0 is incompatible with GeometryMode"):
        liqlev_simulation(inputs, verbose=False)


def test_f11b_neps_legacy_mode_runs() -> None:
    """C2 recovery: Neps>0 on legacy GeometryMode (absent/0) still runs."""

    inputs = build_inputs(
        **_legacy_build_kwargs(
            neps=2,
            teps=np.array([0.0, 2.0]),
            xeps=np.array([0.4, 0.4]),
            tinit_override=38.3,
        )
    )
    assert int(inputs["Neps"]) > 0
    assert "GeometryMode" not in inputs
    df = liqlev_simulation(inputs, verbose=False)
    assert len(df) >= 1
    assert np.isfinite(df.to_numpy(dtype=float)).all()


def test_f11b_geometry_mode_neps_zero_runs() -> None:
    """C2 recovery: GeometryMode==1 with Neps==0 (height_dep) still runs."""

    config = build_lox_vent_config(G3, timestep_s=0.02, duration_s=0.04)
    result = run_single_case(config)
    assert result.inputs["GeometryMode"] == 1
    assert int(result.inputs["Neps"]) == 0
    assert len(result.dataframe) >= 1
    assert np.isfinite(result.dataframe.to_numpy(dtype=float)).all()
