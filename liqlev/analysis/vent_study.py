"""Config-driven re-run interface for the LOX vent study (``vent_study``).

This module is a *driver*, not a model. Every physical quantity is produced by
machinery that already carries the campaign's acceptance evidence:

* case structure — :func:`validation.lox_vent_cases.build_lox_vent_config`
  (epsilon mode, ramp contract, dual termination, tank constants);
* execution — :func:`liqlev.runner.single.run_single_case`;
* summary, ullage-closure metric and physicality assessment —
  :func:`validation.lox_vent_cases.summarize_lox_vent_result`;
* provenance discipline — the F8 hash-bound manifest fields, ``git describe``
  stamp and dirty-tree refusal used by
  :func:`validation.lox_vent_cases.write_lox_vent_manifest`.

What is new here is the *input surface*: a JSON config (fluid, geometry
package, fill in litres or fraction, pressures, vent rate in g/s or lbm/s,
duration, timestep, and either explicit gravity levels or an SSM-style
spaceflight block from which the levels are derived), an automatic dt-plateau
check on the nominal-maximum row, and a report that always quotes the eight
assumptions from ``docs/lox-vent-test-definition.md`` §5.

CLI::

    python -m liqlev.analysis.vent_study --config configs/lox_43L_40to35_ssm3.json

Importable::

    from liqlev.analysis.vent_study import load_study_config, run_study
    result = run_study(load_study_config("configs/lox_43L_40to35_ssm3.json"))
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, TextIO

import numba
import numpy as np

from liqlev.geometry.package import GeometryPackageError, load_geometry_package
from liqlev.model.builder import G_TO_FT_S2
from liqlev.model.config import SimulationConfig
from liqlev.runner.single import PreparedGravity, run_single_case
from validation import lox_vent_cases as lox


ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_MODULE_PATH = Path(__file__).resolve()
ASSUMPTIONS_DOC_PATH = (ROOT / "docs" / "lox-vent-test-definition.md").resolve()
ASSUMPTIONS_SECTION = "5"
DEFAULT_GEOMETRY_NPZ_PATH = lox.GEOMETRY_NPZ_PATH
DEFAULT_OUTPUT_DIR = (ROOT / "validation" / "results" / "vent_study").resolve()

SCHEMA = "liqlev.analysis.vent_study"
SCHEMA_VERSION = 1

# Unit conversions (exact definitions).
FT_TO_MM = 304.8
LITERS_PER_FT3 = 28.316846592
GRAMS_PER_LBM = 453.59237
STANDARD_GRAVITY_M_S2 = 9.80665

# Row names mirror docs/lox-vent-test-definition.md §3: G0 is the zero-g bound,
# G1 the duty-averaged minimum and G3 the peak (nominal-maximum) level. The
# even rows G2/G4 in that table are the sensitivity-thrust variants and are only
# produced by an explicit ``gravity.levels`` block.
ZERO_G_ROW = "G0"
DUTY_AVERAGE_ROW = "G1"
PEAK_ROW = "G3"

# Optional pulsed-profile row. It is *not* part of the constant-g matrix: it
# feeds the synthesized square wave to the solver and is a model-consistency
# experiment only (see PULSED_ROW_CAVEAT and §3 of the test definition).
PULSED_ROW = "GP"
CONSTANT_GRAVITY_MODE = "constant"
PULSED_GRAVITY_MODE = "pulsed"

# Printed verbatim by ``format_report`` whenever a GP row is present. ASCII
# only and plain hyphens so it survives any console encoding untouched.
PULSED_ROW_CAVEAT = """\
Pulsed-profile caveat (row GP) - mandatory wherever GP is quoted
------------------------------------------------------------------------------
Row GP feeds the SSM-3 square wave directly to the solver: g(t) is the pulse
train itself, not a constant level. This departs from the engineer-directed
method of record, which treats gravity as steady state and brackets it with a
maximum and a minimum constant level (docs/lox-vent-test-definition.md
section 3). GP is a model-consistency experiment only. The agreed engineering
method is unchanged and the G1-G3 constant-g bracket remains what is reported.

Assumption 2 (section 5, item 2) applies with full force. LIQLEV's boundary
layer is quasi-steady - an instantaneous function of the current state, with
no relaxation timescale - so it cannot respond correctly to the pulse cycle.
Driving it with the pulse train makes the model jump between a thin
boundary-layer state during the ON interval and the saturated maximum during
the OFF interval, within each 3.4 s cycle. That is representative of the
acceleration history, but it is not more physical than the G1-G3 bracket.
GP is not a validated prediction and must not be quoted as one.
"""

# Square edges are emitted as double breakpoints separated by this riser — the
# same 1e-9 s convention the GUI uses when it extends a CSV gravity profile.
PULSE_RISER_S = 1e-9

# A timestep must resolve the ON interval with at least this many samples, and
# should divide both the ON interval and the period to this tolerance.
PULSE_MIN_STEPS_PER_ON = 2.0
PULSE_DIVISION_TOLERANCE_S = 1e-9

# Plan §"Re-run interface spec": dt-plateau at dt, dt/2, dt/4 (0.02/0.01/0.005).
DT_PLATEAU_DIVISORS = (2, 4)

FLUIDS = ("Nitrogen", "Oxygen", "Hydrogen", "Methane")

_TOP_LEVEL_KEYS = frozenset(
    {
        "schema",
        "version",
        "name",
        "description",
        "fluid",
        "geometry_package",
        "tank_diameter_ft",
        "fill",
        "pressure",
        "vent_rate",
        "duration_s",
        "timestep_s",
        "gravity",
        "dt_plateau",
    }
)


class VentStudyConfigError(ValueError):
    """Raised for malformed, ambiguous or incomplete study configs."""


# --------------------------------------------------------------------------
# Config model
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class SpaceflightSpec:
    """Pulsed-settling block from which a gravity matrix is derived.

    ``peak_g = thrust_N / (vehicle_mass_kg * 9.80665)`` and
    ``duty_average_g = peak_g * duty_cycle`` — the arithmetic that produced the
    committed G1/G3 levels (1.60 N, 280 kg, SSM-3 120/3400 ms).

    ``pulsed_row`` opts the study into the additional ``GP`` row, which feeds
    the square wave itself to the solver; it never changes the constant-g
    matrix. ``phase_ms`` offsets the pulse train and is only meaningful there.
    """

    thrust_n: float
    vehicle_mass_kg: float
    duty_cycle: float
    on_ms: float | None = None
    period_ms: float | None = None
    pulsed_row: bool = False
    phase_ms: float = 0.0

    @property
    def peak_g(self) -> float:
        return self.thrust_n / (self.vehicle_mass_kg * STANDARD_GRAVITY_M_S2)

    @property
    def duty_average_g(self) -> float:
        return self.peak_g * self.duty_cycle

    @property
    def on_s(self) -> float | None:
        return None if self.on_ms is None else self.on_ms / 1000.0

    @property
    def period_s(self) -> float | None:
        return None if self.period_ms is None else self.period_ms / 1000.0

    @property
    def phase_s(self) -> float:
        return self.phase_ms / 1000.0

    def as_dict(self) -> dict[str, float | bool | None]:
        return {
            "thrust_N": self.thrust_n,
            "vehicle_mass_kg": self.vehicle_mass_kg,
            "duty_cycle": self.duty_cycle,
            "on_ms": self.on_ms,
            "period_ms": self.period_ms,
            "pulsed_row": self.pulsed_row,
            "phase_ms": self.phase_ms,
            "peak_g": self.peak_g,
            "duty_average_g": self.duty_average_g,
        }


@dataclass(frozen=True)
class VentStudyConfig:
    """Validated, fully resolved study inputs (no ambiguous alternatives left)."""

    name: str
    fluid: str
    geometry_path: Path
    fill_fraction: float
    initial_pressure_psia: float
    final_pressure_psia: float
    vent_rate_lbm_s: float
    duration_s: float
    timestep_s: float
    gravity_levels: tuple[tuple[str, float], ...]
    gravity_source: str
    tank_diameter_ft: float
    tank_height_ft: float
    dt_plateau_enabled: bool
    dt_plateau_row: str
    description: str = ""
    fill_liters: float | None = None
    spaceflight: SpaceflightSpec | None = None
    source_path: Path | None = None

    @property
    def vent_rate_g_s(self) -> float:
        return self.vent_rate_lbm_s * GRAMS_PER_LBM

    @property
    def gravity_matrix(self) -> dict[str, float]:
        return dict(self.gravity_levels)

    @property
    def has_pulsed_row(self) -> bool:
        """True when this study also runs the pulsed-profile ``GP`` row."""

        return self.spaceflight is not None and self.spaceflight.pulsed_row

    def row_names(self) -> tuple[str, ...]:
        return tuple(name for name, _ in self.gravity_levels)


def derive_gravity_levels(spec: SpaceflightSpec) -> dict[str, float]:
    """Derive the G-matrix from a spaceflight block (mirrors §3 of the doc).

    Returns the zero-g bound, the duty-averaged minimum and the peak
    nominal-maximum keyed ``G0`` / ``G1`` / ``G3``.
    """

    return {
        ZERO_G_ROW: 0.0,
        DUTY_AVERAGE_ROW: spec.duty_average_g,
        PEAK_ROW: spec.peak_g,
    }


# --------------------------------------------------------------------------
# Pulsed-profile synthesis (row GP)
# --------------------------------------------------------------------------


def pulsed_gravity_table(
    peak_g: float,
    on_s: float,
    period_s: float,
    phase_s: float,
    duration_s: float,
    riser_s: float = PULSE_RISER_S,
) -> tuple[np.ndarray, np.ndarray]:
    """Synthesize the SSM-style square wave as a solver gravity table.

    The train is the infinite pulse train sampled on ``[0, duration_s]``::

        g(t) = peak_g   when  ((t - phase_s) mod period_s) < on_s
        g(t) = 0.0      otherwise

    so a phase in ``(period_s - on_s, period_s)`` leaves a pulse already ON at
    ``t = 0`` — that is intended, not an artifact. Square edges are emitted as
    double breakpoints: at a transition ``t_e`` the table carries
    ``(t_e, previous)`` followed by ``(t_e + riser_s, new)``, the same 1e-9 s
    riser the GUI uses when extending a CSV gravity profile. The first node is
    ``(0.0, g(0))`` and the last is exactly ``(duration_s, g(duration_s))``.

    Returns ``(tggo, xggo)`` with ``xggo`` in **ft/s²** (``peak_g`` is a
    standard-g level, multiplied by ``G_TO_FT_S2``), both float64 and
    C-contiguous, with strictly increasing times.
    """

    peak_g = float(peak_g)
    on_s = float(on_s)
    period_s = float(period_s)
    phase_s = float(phase_s)
    duration_s = float(duration_s)
    riser_s = float(riser_s)

    if not np.isfinite([peak_g, on_s, period_s, phase_s, duration_s, riser_s]).all():
        raise ValueError("pulsed_gravity_table arguments must be finite")
    if peak_g < 0.0:
        raise ValueError(f"peak_g must be non-negative (got {peak_g})")
    if period_s <= 0.0:
        raise ValueError(f"period_s must be positive (got {period_s})")
    if on_s <= 0.0:
        raise ValueError(f"on_s must be positive (got {on_s})")
    if on_s > period_s:
        raise ValueError(
            f"on_s must not exceed period_s (got {on_s} > {period_s})"
        )
    if not 0.0 <= phase_s < period_s:
        raise ValueError(
            f"phase_s must satisfy 0 <= phase_s < period_s "
            f"(got {phase_s} with period_s {period_s})"
        )
    if duration_s <= 0.0:
        raise ValueError(f"duration_s must be positive (got {duration_s})")
    if riser_s <= 0.0:
        raise ValueError(f"riser_s must be positive (got {riser_s})")
    if riser_s >= on_s or riser_s >= duration_s:
        raise ValueError(
            f"riser_s must be smaller than on_s and duration_s (got {riser_s})"
        )

    peak_ft_s2 = peak_g * G_TO_FT_S2

    # Edge train built analytically from the pulse index rather than by a
    # modulo of the sample time: at an exact edge the modulo can land on either
    # side of the comparison, which would flip a breakpoint's value.
    first_index = int(np.floor((0.0 - phase_s) / period_s)) - 1
    last_index = int(np.ceil((duration_s - phase_s) / period_s)) + 1
    edges: list[tuple[float, float]] = []
    for index in range(first_index, last_index + 1):
        rise_s = phase_s + index * period_s
        edges.append((rise_s, peak_ft_s2))
        if on_s < period_s:
            edges.append((rise_s + on_s, 0.0))
    edges.sort(key=lambda edge: edge[0])

    def _value_at(time_s: float) -> float:
        value = 0.0
        for edge_time_s, after in edges:
            if edge_time_s > time_s:
                break
            value = after
        return value

    times: list[float] = [0.0]
    values: list[float] = [_value_at(0.0)]

    def _push(time_s: float, value: float) -> None:
        if time_s > times[-1]:
            times.append(time_s)
            values.append(value)

    for edge_time_s, after in edges:
        if edge_time_s <= 0.0 or edge_time_s > duration_s:
            continue
        before = values[-1]
        if after == before:
            continue
        lead_s = edge_time_s
        trail_s = edge_time_s + riser_s
        if trail_s >= duration_s:
            # Keep the last node exactly on the duration boundary.
            lead_s = duration_s - riser_s
            trail_s = duration_s
        _push(lead_s, before)
        _push(trail_s, after)

    if times[-1] < duration_s:
        _push(duration_s, _value_at(duration_s))
    else:
        values[-1] = _value_at(duration_s)

    tggo = np.ascontiguousarray(times, dtype=np.float64)
    xggo = np.ascontiguousarray(values, dtype=np.float64)
    return tggo, xggo


def pulse_table_sha256(tggo: np.ndarray, xggo: np.ndarray) -> str:
    """SHA-256 over the raw bytes of both table arrays (provenance binding)."""

    digest = hashlib.sha256()
    digest.update(np.ascontiguousarray(tggo, dtype=np.float64).tobytes())
    digest.update(np.ascontiguousarray(xggo, dtype=np.float64).tobytes())
    return digest.hexdigest().upper()


def _divides(value_s: float, timestep_s: float) -> bool:
    remainder = value_s % timestep_s
    return min(remainder, timestep_s - remainder) <= PULSE_DIVISION_TOLERANCE_S


def check_pulse_resolution(spec: SpaceflightSpec, timestep_s: float) -> None:
    """Hard guard: the timestep must resolve the ON interval.

    Fewer than two steps inside the ON interval means the square wave is not
    represented at all, so this refuses rather than reporting a warning.
    """

    on_s = spec.on_s
    if on_s is None:
        return
    timestep_s = float(timestep_s)
    if on_s < PULSE_MIN_STEPS_PER_ON * timestep_s:
        raise VentStudyConfigError(
            f"pulsed row is under-resolved: on_ms {spec.on_ms:g} "
            f"({on_s:g} s) needs at least "
            f"{PULSE_MIN_STEPS_PER_ON:g} timesteps but timestep_s is "
            f"{timestep_s:g} s; shorten timestep_s or drop pulsed_row"
        )


def pulse_timing_notes(
    spec: SpaceflightSpec, timestep_s: float
) -> tuple[str, ...]:
    """Aliasing notes when the timestep does not divide the pulse timing.

    The production timesteps 0.02 / 0.01 / 0.005 s all divide the SSM-3 120 ms
    and 3400 ms exactly, so this is normally empty.
    """

    on_s = spec.on_s
    period_s = spec.period_s
    if on_s is None or period_s is None:
        return ()
    timestep_s = float(timestep_s)
    notes: list[str] = []
    if not _divides(on_s, timestep_s):
        notes.append(
            f"WARNING (aliasing): timestep {timestep_s:g} s does not divide "
            f"on_s = {on_s:g} s; the pulse train is sampled unevenly."
        )
    if not _divides(period_s, timestep_s):
        notes.append(
            f"WARNING (aliasing): timestep {timestep_s:g} s does not divide "
            f"period_s = {period_s:g} s; the pulse train is sampled unevenly."
        )
    return tuple(notes)


def default_dt_plateau_row(levels: Sequence[tuple[str, float]]) -> str:
    """Nominal-maximum row for the dt-plateau check.

    ``G3`` when present (the §3 nominal-maximum name, which the spaceflight
    derivation also emits); otherwise the highest gravity level supplied.
    """

    names = [name for name, _ in levels]
    if PEAK_ROW in names:
        return PEAK_ROW
    return max(levels, key=lambda item: item[1])[0]


def _require_mapping(value: Any, field: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise VentStudyConfigError(f"{field} must be a JSON object")
    return value


def _positive_float(value: Any, field: str, *, allow_zero: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise VentStudyConfigError(f"{field} must be a number")
    numeric = float(value)
    if not np.isfinite(numeric):
        raise VentStudyConfigError(f"{field} must be finite")
    if numeric < 0.0 or (numeric == 0.0 and not allow_zero):
        limit = "non-negative" if allow_zero else "positive"
        raise VentStudyConfigError(f"{field} must be {limit} (got {numeric})")
    return numeric


def _exactly_one(
    payload: Mapping[str, Any], options: Sequence[str], field: str
) -> str:
    present = [name for name in options if name in payload]
    if len(present) != 1:
        choices = " or ".join(repr(name) for name in options)
        got = ", ".join(repr(name) for name in present) if present else "none"
        raise VentStudyConfigError(
            f"{field} requires exactly one of {choices} (got {got})"
        )
    return present[0]


def _resolve_geometry_path(value: Any) -> Path:
    """Resolve a geometry package path; relative paths are repo-root relative."""

    if value is None:
        return DEFAULT_GEOMETRY_NPZ_PATH
    if not isinstance(value, str) or not value.strip():
        raise VentStudyConfigError("geometry_package must be a non-empty string")
    candidate = Path(value)
    if not candidate.is_absolute():
        candidate = ROOT / candidate
    resolved = candidate.resolve()
    if not resolved.is_file():
        raise VentStudyConfigError(f"geometry package not found: {resolved}")
    return resolved


def _geometry_tank_constants(geometry_path: Path) -> tuple[float, float]:
    """Return ``(diameter_ft, total_height_ft)`` for a geometry package.

    The committed package keeps the pinned case constants byte-for-byte so the
    production case reproduces exactly; any other package gets its total height
    from the package and an equivalent diameter from the maximum section area
    (the same convention used to pin 1.84 ft for the NASA fluid domain).
    """

    if geometry_path == DEFAULT_GEOMETRY_NPZ_PATH:
        return lox.TANK_DIAMETER_FT, lox.TANK_TOTAL_HEIGHT_FT
    kernel = _load_geometry(geometry_path)
    max_area = float(np.max(kernel.section_area_ft2))
    return 2.0 * float(np.sqrt(max_area / np.pi)), float(kernel.total_height_ft)


def _load_geometry(geometry_path: Path):
    try:
        return load_geometry_package(geometry_path)
    except GeometryPackageError as exc:
        raise VentStudyConfigError(f"invalid geometry package: {exc}") from exc


def _fill_fraction(payload: Mapping[str, Any], geometry_path: Path) -> tuple[
    float, float | None
]:
    fill = _require_mapping(payload, "fill")
    unknown = set(fill) - {"liters", "fraction"}
    if unknown:
        raise VentStudyConfigError(
            f"fill has unknown key(s): {', '.join(sorted(unknown))}"
        )
    chosen = _exactly_one(fill, ("liters", "fraction"), "fill")
    if chosen == "fraction":
        fraction = _positive_float(fill["fraction"], "fill.fraction")
        liters: float | None = None
    else:
        liters = _positive_float(fill["liters"], "fill.liters")
        total_liters = _load_geometry(geometry_path).total_volume_ft3 * LITERS_PER_FT3
        fraction = liters / total_liters
    if not 0.0 < fraction <= 1.0:
        raise VentStudyConfigError(
            f"fill resolves to fraction {fraction!r}, outside (0, 1]"
        )
    return fraction, liters


def _vent_rate_lbm_s(payload: Mapping[str, Any]) -> float:
    vent = _require_mapping(payload, "vent_rate")
    unknown = set(vent) - {"g_per_s", "lbm_per_s"}
    if unknown:
        raise VentStudyConfigError(
            f"vent_rate has unknown key(s): {', '.join(sorted(unknown))}"
        )
    chosen = _exactly_one(vent, ("g_per_s", "lbm_per_s"), "vent_rate")
    if chosen == "lbm_per_s":
        return _positive_float(vent["lbm_per_s"], "vent_rate.lbm_per_s")
    return _positive_float(vent["g_per_s"], "vent_rate.g_per_s") / GRAMS_PER_LBM


def _spaceflight_spec(payload: Mapping[str, Any]) -> SpaceflightSpec:
    block = _require_mapping(payload, "gravity.spaceflight")
    allowed = {
        "thrust_N",
        "vehicle_mass_kg",
        "duty_cycle",
        "on_ms",
        "period_ms",
        "pulsed_row",
        "phase_ms",
    }
    unknown = set(block) - allowed
    if unknown:
        raise VentStudyConfigError(
            "gravity.spaceflight has unknown key(s): "
            f"{', '.join(sorted(unknown))}"
        )
    for required in ("thrust_N", "vehicle_mass_kg"):
        if required not in block:
            raise VentStudyConfigError(
                f"gravity.spaceflight requires {required!r}"
            )
    thrust = _positive_float(block["thrust_N"], "gravity.spaceflight.thrust_N")
    mass = _positive_float(
        block["vehicle_mass_kg"], "gravity.spaceflight.vehicle_mass_kg"
    )

    pulsed_row = block.get("pulsed_row", False)
    if not isinstance(pulsed_row, bool):
        raise VentStudyConfigError(
            "gravity.spaceflight.pulsed_row must be true or false"
        )

    has_duty = "duty_cycle" in block
    has_pulse = "on_ms" in block or "period_ms" in block
    if pulsed_row and not has_pulse:
        raise VentStudyConfigError(
            "gravity.spaceflight.pulsed_row requires the 'on_ms'/'period_ms' "
            "form: pulse timing is required to synthesize the square wave, and "
            "'duty_cycle' alone does not define it"
        )
    if has_duty == has_pulse:
        raise VentStudyConfigError(
            "gravity.spaceflight requires exactly one of 'duty_cycle' or the "
            "'on_ms'/'period_ms' pair"
        )
    if has_duty:
        duty = _positive_float(block["duty_cycle"], "gravity.spaceflight.duty_cycle")
        on_ms: float | None = None
        period_ms: float | None = None
    else:
        if "on_ms" not in block or "period_ms" not in block:
            raise VentStudyConfigError(
                "gravity.spaceflight needs both 'on_ms' and 'period_ms'"
            )
        on_ms = _positive_float(block["on_ms"], "gravity.spaceflight.on_ms")
        period_ms = _positive_float(block["period_ms"], "gravity.spaceflight.period_ms")
        if on_ms > period_ms:
            raise VentStudyConfigError(
                "gravity.spaceflight.on_ms must not exceed period_ms"
            )
        duty = on_ms / period_ms
    if not 0.0 < duty <= 1.0:
        raise VentStudyConfigError(
            f"gravity.spaceflight duty cycle {duty!r} is outside (0, 1]"
        )

    if "phase_ms" in block and not pulsed_row:
        raise VentStudyConfigError(
            "gravity.spaceflight.phase_ms requires pulsed_row true "
            "(the phase only shifts the synthesized pulse train)"
        )
    phase_ms = _positive_float(
        block.get("phase_ms", 0.0), "gravity.spaceflight.phase_ms", allow_zero=True
    )
    if period_ms is not None and not 0.0 <= phase_ms < period_ms:
        raise VentStudyConfigError(
            "gravity.spaceflight.phase_ms must satisfy 0 <= phase_ms < "
            f"period_ms (got {phase_ms:g} with period_ms {period_ms:g})"
        )

    return SpaceflightSpec(
        thrust_n=thrust,
        vehicle_mass_kg=mass,
        duty_cycle=duty,
        on_ms=on_ms,
        period_ms=period_ms,
        pulsed_row=pulsed_row,
        phase_ms=phase_ms,
    )


def _explicit_levels(payload: Any) -> dict[str, float]:
    levels = _require_mapping(payload, "gravity.levels")
    if not levels:
        raise VentStudyConfigError("gravity.levels must contain at least one row")
    resolved: dict[str, float] = {}
    for name, value in levels.items():
        if not isinstance(name, str) or not name.strip():
            raise VentStudyConfigError("gravity.levels keys must be non-empty strings")
        resolved[name] = _positive_float(
            value, f"gravity.levels[{name!r}]", allow_zero=True
        )
    return resolved


def _gravity_block(
    block: Any,
) -> tuple[tuple[tuple[str, float], ...], str, SpaceflightSpec | None]:
    gravity = _require_mapping(block, "gravity")
    unknown = set(gravity) - {"levels", "spaceflight"}
    if unknown:
        raise VentStudyConfigError(
            f"gravity has unknown key(s): {', '.join(sorted(unknown))}"
        )
    chosen = _exactly_one(gravity, ("levels", "spaceflight"), "gravity")
    if chosen == "levels":
        return tuple(_explicit_levels(gravity["levels"]).items()), "levels", None
    spec = _spaceflight_spec(gravity["spaceflight"])
    return tuple(derive_gravity_levels(spec).items()), "spaceflight", spec


def study_config_from_mapping(
    payload: Mapping[str, Any], *, source_path: str | Path | None = None
) -> VentStudyConfig:
    """Validate a config mapping into a :class:`VentStudyConfig`."""

    payload = _require_mapping(payload, "config")
    unknown = set(payload) - _TOP_LEVEL_KEYS
    if unknown:
        raise VentStudyConfigError(
            f"config has unknown key(s): {', '.join(sorted(unknown))}"
        )
    schema = payload.get("schema", SCHEMA)
    if schema != SCHEMA:
        raise VentStudyConfigError(
            f"config schema must be {SCHEMA!r} (got {schema!r})"
        )
    version = payload.get("version", SCHEMA_VERSION)
    if version != SCHEMA_VERSION:
        raise VentStudyConfigError(
            f"unsupported config version {version!r} (expected {SCHEMA_VERSION})"
        )

    for required in ("fluid", "fill", "pressure", "vent_rate", "duration_s",
                     "timestep_s", "gravity"):
        if required not in payload:
            raise VentStudyConfigError(f"config is missing required key {required!r}")

    fluid = payload["fluid"]
    if fluid not in FLUIDS:
        raise VentStudyConfigError(
            f"fluid must be one of {', '.join(FLUIDS)} (got {fluid!r})"
        )

    geometry_path = _resolve_geometry_path(payload.get("geometry_package"))
    fill_fraction, fill_liters = _fill_fraction(payload["fill"], geometry_path)

    pressure = _require_mapping(payload["pressure"], "pressure")
    unknown = set(pressure) - {"initial_psia", "final_psia"}
    if unknown:
        raise VentStudyConfigError(
            f"pressure has unknown key(s): {', '.join(sorted(unknown))}"
        )
    for required in ("initial_psia", "final_psia"):
        if required not in pressure:
            raise VentStudyConfigError(f"pressure requires {required!r}")
    initial_psia = _positive_float(pressure["initial_psia"], "pressure.initial_psia")
    final_psia = _positive_float(pressure["final_psia"], "pressure.final_psia")
    if initial_psia <= final_psia:
        raise VentStudyConfigError(
            "pressure.initial_psia must exceed pressure.final_psia "
            f"({initial_psia} <= {final_psia})"
        )

    vent_rate = _vent_rate_lbm_s(payload["vent_rate"])
    duration_s = _positive_float(payload["duration_s"], "duration_s")
    timestep_s = _positive_float(payload["timestep_s"], "timestep_s")
    if timestep_s >= duration_s:
        raise VentStudyConfigError(
            f"timestep_s ({timestep_s}) must be smaller than duration_s "
            f"({duration_s})"
        )

    levels, gravity_source, spaceflight = _gravity_block(payload["gravity"])
    if spaceflight is not None and spaceflight.pulsed_row:
        check_pulse_resolution(spaceflight, timestep_s)

    diameter_default, height_default = _geometry_tank_constants(geometry_path)
    if "tank_diameter_ft" in payload:
        diameter = _positive_float(payload["tank_diameter_ft"], "tank_diameter_ft")
    else:
        diameter = diameter_default

    plateau = _require_mapping(payload.get("dt_plateau", {}), "dt_plateau")
    unknown = set(plateau) - {"enabled", "row"}
    if unknown:
        raise VentStudyConfigError(
            f"dt_plateau has unknown key(s): {', '.join(sorted(unknown))}"
        )
    enabled = plateau.get("enabled", True)
    if not isinstance(enabled, bool):
        raise VentStudyConfigError("dt_plateau.enabled must be true or false")
    plateau_row = plateau.get("row", default_dt_plateau_row(levels))
    selectable = set(dict(levels))
    if spaceflight is not None and spaceflight.pulsed_row:
        selectable.add(PULSED_ROW)
    if plateau_row not in selectable:
        raise VentStudyConfigError(
            f"dt_plateau.row {plateau_row!r} is not a configured gravity level"
        )

    name = payload.get("name")
    if name is None:
        name = Path(source_path).stem if source_path is not None else "vent_study"
    if not isinstance(name, str) or not name.strip():
        raise VentStudyConfigError("name must be a non-empty string")
    description = payload.get("description", "")
    if not isinstance(description, str):
        raise VentStudyConfigError("description must be a string")

    return VentStudyConfig(
        name=name,
        description=description,
        fluid=fluid,
        geometry_path=geometry_path,
        fill_fraction=fill_fraction,
        fill_liters=fill_liters,
        initial_pressure_psia=initial_psia,
        final_pressure_psia=final_psia,
        vent_rate_lbm_s=vent_rate,
        duration_s=duration_s,
        timestep_s=timestep_s,
        gravity_levels=levels,
        gravity_source=gravity_source,
        spaceflight=spaceflight,
        tank_diameter_ft=diameter,
        tank_height_ft=height_default,
        dt_plateau_enabled=enabled,
        dt_plateau_row=plateau_row,
        source_path=Path(source_path).resolve() if source_path is not None else None,
    )


def load_study_config(path: str | Path) -> VentStudyConfig:
    """Load and validate a study config JSON file."""

    target = Path(path)
    try:
        text = target.read_text(encoding="utf-8")
    except OSError as exc:
        raise VentStudyConfigError(f"could not read config: {target}") from exc
    try:
        payload = json.loads(text)
    except json.JSONDecodeError as exc:
        raise VentStudyConfigError(f"config is not valid JSON: {exc}") from exc
    return study_config_from_mapping(payload, source_path=target)


def select_rows(config: VentStudyConfig, names: Iterable[str]) -> VentStudyConfig:
    """Return a config restricted to ``names`` (config order preserved).

    ``GP`` is selectable exactly when the config enables the pulsed row; not
    naming it turns that row off, so a narrowed run never smuggles it back in.
    """

    wanted = list(dict.fromkeys(names))
    available = dict(config.gravity_levels)
    defined = list(available)
    if config.has_pulsed_row:
        defined.append(PULSED_ROW)
    missing = [name for name in wanted if name not in defined]
    if missing:
        raise VentStudyConfigError(
            f"unknown gravity row(s): {', '.join(missing)}; "
            f"config defines {', '.join(defined)}"
        )
    levels = tuple(
        (name, value) for name, value in config.gravity_levels if name in set(wanted)
    )
    keeps_pulsed = config.has_pulsed_row and PULSED_ROW in wanted
    updates: dict[str, Any] = {"gravity_levels": levels}
    if config.has_pulsed_row and not keeps_pulsed:
        updates["spaceflight"] = replace(config.spaceflight, pulsed_row=False)
    plateau_row = config.dt_plateau_row
    if plateau_row not in dict(levels) and not (
        keeps_pulsed and plateau_row == PULSED_ROW
    ):
        plateau_row = (
            PULSED_ROW if keeps_pulsed and not levels
            else default_dt_plateau_row(levels)
        )
    updates["dt_plateau_row"] = plateau_row
    return replace(config, **updates)


# --------------------------------------------------------------------------
# Execution
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class VentStudyRow:
    """One gravity row: the reused summary plus the reported rise in mm.

    ``gravity_mode`` is ``"constant"`` for the bracketing matrix rows and
    ``"pulsed"`` for ``GP``; ``pulse`` carries the GP profile descriptor (and
    is ``None`` otherwise), ``notes`` any run-time timing warnings.
    """

    name: str
    gravity_g: float
    rise_mm: float
    summary: lox.LoxVentSummary
    gravity_mode: str = CONSTANT_GRAVITY_MODE
    pulse: dict[str, Any] | None = None
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class DtPlateauPoint:
    timestep_s: float
    rise_mm: float
    delta_percent: float


@dataclass(frozen=True)
class DtPlateauReport:
    row: str
    gravity_g: float
    base_timestep_s: float
    base_rise_mm: float
    points: tuple[DtPlateauPoint, ...]
    max_abs_delta_percent: float


@dataclass(frozen=True)
class VentStudyResult:
    config: VentStudyConfig
    rows: tuple[VentStudyRow, ...]
    dt_plateau: DtPlateauReport | None
    assumptions: str


def simulation_config_for(
    config: VentStudyConfig,
    gravity_g: float,
    *,
    timestep_s: float | None = None,
    duration_s: float | None = None,
) -> SimulationConfig:
    """Build the solver config for one row from the committed case builder.

    ``build_lox_vent_config`` stays the single definition of case structure
    (epsilon mode, ramp contract, gravity mode, run controls); only the
    config-driven fields are substituted.
    """

    timestep = config.timestep_s if timestep_s is None else float(timestep_s)
    duration = config.duration_s if duration_s is None else float(duration_s)
    base = lox.build_lox_vent_config(
        gravity_g, timestep_s=timestep, duration_s=duration
    )
    return replace(
        base,
        fluid=replace(
            base.fluid,
            name=config.fluid,
            initial_pressure_psia=config.initial_pressure_psia,
            final_pressure_psia=config.final_pressure_psia,
        ),
        tank=replace(
            base.tank,
            diameter_ft=config.tank_diameter_ft,
            height_ft=config.tank_height_ft,
            fill_fractions=(config.fill_fraction,),
            geometry_path=str(config.geometry_path),
        ),
        vent=replace(base.vent, rates_lbm_s=(config.vent_rate_lbm_s,)),
    )


def _summarize_row(
    config: VentStudyConfig,
    name: str,
    gravity_g: float,
    result,
    timestep_s: float,
    *,
    gravity_mode: str = CONSTANT_GRAVITY_MODE,
    pulse: dict[str, Any] | None = None,
    notes: Sequence[str] = (),
) -> VentStudyRow:
    """Summarise one executed case with the validation machinery."""

    summary = lox.summarize_lox_vent_result(
        result, gravity_g=gravity_g, timestep_s=timestep_s
    )
    if abs(config.tank_height_ft - lox.TANK_TOTAL_HEIGHT_FT) > 1e-9:
        # Non-default geometry: re-assess bounds against this tank's height
        # using the same assessment function (no metric is recomputed here).
        physical, classifications = lox.assess_lox_dataframe_physicality(
            result.dataframe, total_height_ft=config.tank_height_ft
        )
        summary = replace(
            summary, physical=physical, failure_classifications=classifications
        )
    return VentStudyRow(
        name=name,
        gravity_g=float(gravity_g),
        rise_mm=summary.dh_ft * FT_TO_MM,
        summary=summary,
        gravity_mode=gravity_mode,
        pulse=pulse,
        notes=tuple(notes),
    )


def run_row(
    config: VentStudyConfig,
    name: str,
    gravity_g: float,
    *,
    timestep_s: float | None = None,
    duration_s: float | None = None,
) -> VentStudyRow:
    """Run one constant-gravity row and summarise it."""

    timestep = config.timestep_s if timestep_s is None else float(timestep_s)
    simulation = simulation_config_for(
        config, gravity_g, timestep_s=timestep, duration_s=duration_s
    )
    result = run_single_case(simulation)
    return _summarize_row(config, name, gravity_g, result, timestep)


def prepare_pulsed_gravity(
    config: VentStudyConfig, *, duration_s: float | None = None
) -> tuple[PreparedGravity, dict[str, Any]]:
    """Synthesize the GP gravity profile and its provenance descriptor.

    The returned :class:`~liqlev.runner.single.PreparedGravity` is handed to
    ``run_single_case`` as ``prepared_gravity``, which overrides gravity
    entirely (``single.py``: ``gravity = prepared_gravity or prepare_gravity``),
    so the solver integrates this table and nothing else.
    """

    spec = config.spaceflight
    if spec is None or not spec.pulsed_row:
        raise VentStudyConfigError(
            "the pulsed row requires gravity.spaceflight.pulsed_row true"
        )
    duration = config.duration_s if duration_s is None else float(duration_s)
    tggo, xggo = pulsed_gravity_table(
        spec.peak_g, spec.on_s, spec.period_s, spec.phase_s, duration
    )
    descriptor: dict[str, Any] = {
        "thrust_N": spec.thrust_n,
        "vehicle_mass_kg": spec.vehicle_mass_kg,
        "on_ms": spec.on_ms,
        "period_ms": spec.period_ms,
        "phase_ms": spec.phase_ms,
        "peak_g": spec.peak_g,
        "mean_g": spec.duty_average_g,
        "n_table_points": int(len(tggo)),
        "table_sha256": pulse_table_sha256(tggo, xggo),
    }
    message = (
        f"Gravity: pulsed square wave {spec.on_ms:g}/{spec.period_ms:g} ms, "
        f"phase {spec.phase_ms:g} ms, peak {spec.peak_g:.6e} g, "
        f"mean {spec.duty_average_g:.6e} g, {len(tggo)} table points"
    )
    prepared = PreparedGravity(
        nggo=int(len(tggo)), tggo=tggo, xggo=xggo, messages=(message,)
    )
    return prepared, descriptor


def run_pulsed_row(
    config: VentStudyConfig,
    *,
    timestep_s: float | None = None,
    duration_s: float | None = None,
) -> VentStudyRow:
    """Run the pulsed-profile ``GP`` row against the synthesized square wave.

    The case is built by the same :func:`simulation_config_for` path as every
    constant row, at ``gravity_g = peak_g * duty_cycle`` — the true time mean
    of the train. That value is bookkeeping and reporting only: the profile
    handed to ``run_single_case`` is what the solver integrates.
    """

    spec = config.spaceflight
    if spec is None or not spec.pulsed_row:
        raise VentStudyConfigError(
            "the pulsed row requires gravity.spaceflight.pulsed_row true"
        )
    timestep = config.timestep_s if timestep_s is None else float(timestep_s)
    duration = config.duration_s if duration_s is None else float(duration_s)
    check_pulse_resolution(spec, timestep)
    prepared, descriptor = prepare_pulsed_gravity(config, duration_s=duration)
    mean_g = descriptor["mean_g"]
    simulation = simulation_config_for(
        config, mean_g, timestep_s=timestep, duration_s=duration
    )
    result = run_single_case(simulation, prepared_gravity=prepared)
    return _summarize_row(
        config,
        PULSED_ROW,
        mean_g,
        result,
        timestep,
        gravity_mode=PULSED_GRAVITY_MODE,
        pulse=descriptor,
        notes=pulse_timing_notes(spec, timestep),
    )


def run_dt_plateau(
    config: VentStudyConfig,
    *,
    row: str | None = None,
    base_row: VentStudyRow | None = None,
) -> DtPlateauReport:
    """Re-run the nominal-maximum row at dt/2 and dt/4 and report plateau %.

    ``delta_percent`` is the signed change in rise relative to the primary
    timestep, i.e. ``100 * (rise(dt_i) - rise(dt)) / rise(dt)``.

    ``GP`` is accepted like any other row: the synthesized table is a function
    of the pulse timing alone, so halving the timestep only resamples it.
    """

    name = config.dt_plateau_row if row is None else row
    levels = dict(config.gravity_levels)
    if name == PULSED_ROW and config.has_pulsed_row:
        gravity_g = config.spaceflight.duty_average_g

        def _run(timestep_s: float | None = None) -> VentStudyRow:
            return run_pulsed_row(config, timestep_s=timestep_s)

    elif name in levels:
        gravity_g = levels[name]

        def _run(timestep_s: float | None = None) -> VentStudyRow:
            return run_row(config, name, gravity_g, timestep_s=timestep_s)

    else:
        raise VentStudyConfigError(
            f"dt-plateau row {name!r} is not a configured gravity level"
        )
    base = base_row if base_row is not None else _run()
    base_rise = base.rise_mm
    points = [DtPlateauPoint(config.timestep_s, base_rise, 0.0)]
    for divisor in DT_PLATEAU_DIVISORS:
        timestep = config.timestep_s / divisor
        refined = _run(timestep_s=timestep)
        delta = (
            100.0 * (refined.rise_mm - base_rise) / base_rise
            if base_rise not in (0.0,) and np.isfinite(base_rise)
            else float("nan")
        )
        points.append(DtPlateauPoint(timestep, refined.rise_mm, delta))
    finite = [abs(point.delta_percent) for point in points if np.isfinite(point.delta_percent)]
    return DtPlateauReport(
        row=name,
        gravity_g=gravity_g,
        base_timestep_s=config.timestep_s,
        base_rise_mm=base_rise,
        points=tuple(points),
        max_abs_delta_percent=max(finite) if finite else float("nan"),
    )


def read_assumptions_block(path: str | Path | None = None) -> str:
    """Return §5 of the LOX vent test definition verbatim.

    The eight assumptions are quoted from the committed document rather than
    restated here, so the report cannot drift from the definition of record.
    """

    target = Path(path) if path is not None else ASSUMPTIONS_DOC_PATH
    try:
        text = target.read_text(encoding="utf-8")
    except OSError as exc:  # pragma: no cover - committed doc is always present
        raise RuntimeError(f"could not read assumptions document: {target}") from exc
    lines = text.splitlines()
    start = None
    for index, line in enumerate(lines):
        if line.startswith(f"## {ASSUMPTIONS_SECTION}."):
            start = index
            break
    if start is None:  # pragma: no cover - guarded by test
        raise RuntimeError(
            f"section {ASSUMPTIONS_SECTION} not found in {target}"
        )
    end = len(lines)
    for index in range(start + 1, len(lines)):
        if lines[index].startswith("## "):
            end = index
            break
    return "\n".join(lines[start:end]).strip()


def run_study(
    config: VentStudyConfig, *, stream: TextIO | None = None
) -> VentStudyResult:
    """Run every gravity row plus the automatic dt-plateau check."""

    rows: list[VentStudyRow] = []
    for name, gravity_g in config.gravity_levels:
        if stream is not None:
            print(
                f"[vent-study] running {name} (g = {gravity_g:.6e}) "
                f"dt = {config.timestep_s} s, duration = {config.duration_s} s",
                file=stream,
                flush=True,
            )
        rows.append(run_row(config, name, gravity_g))

    if config.has_pulsed_row:
        spec = config.spaceflight
        if stream is not None:
            print(
                f"[vent-study] running {PULSED_ROW} (pulsed square wave "
                f"{spec.on_ms:g}/{spec.period_ms:g} ms, phase "
                f"{spec.phase_ms:g} ms, mean g = {spec.duty_average_g:.6e}) "
                f"dt = {config.timestep_s} s, duration = {config.duration_s} s",
                file=stream,
                flush=True,
            )
        rows.append(run_pulsed_row(config))

    plateau: DtPlateauReport | None = None
    if config.dt_plateau_enabled:
        base = next(
            (row for row in rows if row.name == config.dt_plateau_row), None
        )
        if stream is not None:
            print(
                f"[vent-study] dt-plateau check on {config.dt_plateau_row} "
                f"at dt/{DT_PLATEAU_DIVISORS[0]} and dt/{DT_PLATEAU_DIVISORS[1]}",
                file=stream,
                flush=True,
            )
        plateau = run_dt_plateau(config, base_row=base)

    return VentStudyResult(
        config=config,
        rows=tuple(rows),
        dt_plateau=plateau,
        assumptions=read_assumptions_block(),
    )


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------


def _closure_verdict(summary: lox.LoxVentSummary) -> str:
    if not np.isfinite(summary.ullage_closure_max_relative):
        return "n/a"
    return "under" if summary.ullage_closure_within_tolerance else "EXCEEDS"


def format_report(result: VentStudyResult) -> str:
    """Render the results table, the dt-plateau block and §5 assumptions."""

    config = result.config
    fill = f"{config.fill_fraction:.6f} fraction"
    if config.fill_liters is not None:
        fill += f" ({config.fill_liters:g} L)"
    lines = [
        f"LIQLEV vent study — {config.name}",
        "=" * 78,
    ]
    if config.description:
        lines.append(config.description)
        lines.append("")
    if config.source_path is not None:
        lines.append(f"Config          : {config.source_path}")
    lines += [
        f"Fluid           : {config.fluid}",
        f"Geometry        : {_repo_relative(config.geometry_path)}",
        f"Fill            : {fill}",
        f"Pressure        : {config.initial_pressure_psia:g} -> "
        f"{config.final_pressure_psia:g} psia",
        f"Vent rate       : {config.vent_rate_lbm_s:.9g} lbm/s "
        f"({config.vent_rate_g_s:.6g} g/s)",
        f"Duration        : {config.duration_s:g} s   "
        f"Timestep: {config.timestep_s:g} s   Epsilon: height_dep",
        f"Tank            : D = {config.tank_diameter_ft:g} ft, "
        f"H = {config.tank_height_ft:g} ft",
        f"Gravity source  : {config.gravity_source}",
    ]
    if config.spaceflight is not None:
        spec = config.spaceflight
        pulse = ""
        if spec.on_ms is not None and spec.period_ms is not None:
            pulse = f", {spec.on_ms:g}/{spec.period_ms:g} ms"
        lines.append(
            f"Spaceflight     : {spec.thrust_n:g} N / {spec.vehicle_mass_kg:g} kg, "
            f"duty {spec.duty_cycle * 100.0:.4f}%{pulse}"
        )
    pulsed_rows = [
        row for row in result.rows if row.gravity_mode == PULSED_GRAVITY_MODE
    ]
    # The mode column only appears when a pulsed row is present, so a
    # constant-g report keeps its pre-feature layout byte for byte.
    mode_width = 9 if pulsed_rows else 0
    mode_header = f"{'mode':<{mode_width}} " if pulsed_rows else ""
    lines += [
        f"Versions        : python {platform.python_version()} / "
        f"numpy {np.__version__} / numba {numba.__version__}",
        "",
        f"{'Row':<6} {mode_header}{'g (std)':>13} {'t_end (s)':>10} "
        f"{'rise (mm)':>11} "
        f"{'dh/h0':>8} {'max AK3':>9} {'ullage':>9} {'vs 5%':>8} "
        f"{'ConvFail':>9}  physical",
        "-" * (108 + mode_width + (1 if pulsed_rows else 0)),
    ]
    for row in result.rows:
        summary = row.summary
        physical = "yes" if summary.physical else (
            "NO (" + ", ".join(summary.failure_classifications) + ")"
        )
        mode = f"{row.gravity_mode:<{mode_width}} " if pulsed_rows else ""
        lines.append(
            f"{row.name:<6} {mode}{row.gravity_g:>13.6e} {summary.t_end_s:>10.2f} "
            f"{row.rise_mm:>11.3f} {summary.dh_over_h0:>8.4f} "
            f"{summary.max_ak3:>9.4f} "
            f"{summary.ullage_closure_max_relative * 100.0:>8.3f}% "
            f"{_closure_verdict(summary):>8} {summary.conv_failed_total:>9d}  "
            f"{physical}"
        )
    for row in pulsed_rows:
        descriptor = row.pulse or {}
        lines.append("")
        lines.append(
            f"Row {row.name} pulsed profile: "
            f"{descriptor.get('on_ms', 0.0):g}/"
            f"{descriptor.get('period_ms', 0.0):g} ms ON/period, phase "
            f"{descriptor.get('phase_ms', 0.0):g} ms"
        )
        lines.append(
            f"  peak g = {descriptor.get('peak_g', float('nan')):.6e}   "
            f"mean g = {descriptor.get('mean_g', float('nan')):.6e} (reported "
            "in the table above)"
        )
        lines.append(
            f"  table = {descriptor.get('n_table_points', 0)} points, "
            f"sha256 = {descriptor.get('table_sha256', '')}"
        )
        for note in row.notes:
            lines.append(f"  {note}")
    lines.append("")
    lines.append(
        f"Ullage-closure acceptance gate: "
        f"{lox.ULLAGE_CLOSURE_RELATIVE_TOLERANCE * 100.0:g}% "
        "(reported, never relaxed)."
    )
    if pulsed_rows:
        lines.append("")
        lines.append(PULSED_ROW_CAVEAT.rstrip("\n"))

    plateau = result.dt_plateau
    lines.append("")
    if plateau is None:
        lines.append("dt-plateau check: disabled for this run.")
    else:
        lines.append(
            f"dt-plateau check — row {plateau.row} "
            f"(g = {plateau.gravity_g:.6e}), base dt = "
            f"{plateau.base_timestep_s:g} s"
        )
        for point in plateau.points:
            marker = (
                "(base)"
                if point.timestep_s == plateau.base_timestep_s
                else f"{point.delta_percent:+.3f}%"
            )
            lines.append(
                f"  dt = {point.timestep_s:<9g} rise = "
                f"{point.rise_mm:>10.3f} mm   {marker}"
            )
        lines.append(
            f"  max |delta| vs base = {plateau.max_abs_delta_percent:.3f}%"
        )

    lines += [
        "",
        f"Assumptions and limitations — quoted verbatim from "
        f"{_repo_relative(ASSUMPTIONS_DOC_PATH)} §{ASSUMPTIONS_SECTION}",
        "-" * 78,
        result.assumptions,
    ]
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------
# Manifest (F8 provenance discipline, reused from the validation harness)
# --------------------------------------------------------------------------


def _git_worktree_is_dirty() -> bool:
    """Dirty-tree probe — delegates to the validation harness guard (F8)."""

    return lox._git_worktree_is_dirty()


def _git_describe() -> str:
    """``git describe --dirty --always`` — same stamp as the LOX manifest."""

    return lox._git_describe_dirty()


def _repo_relative(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT)).replace("\\", "/")
    except ValueError:
        return path.as_posix()


def _row_manifest(row: VentStudyRow) -> dict[str, Any]:
    payload: dict[str, Any] = {"gravity_row": row.name}
    # Reuse the production summary serializer so both manifests agree field
    # for field; only the reported rise in mm is added here.
    payload.update(lox._summary_manifest(row.summary))
    payload["rise_mm"] = lox._finite_or_none(row.rise_mm)
    payload["gravity_mode"] = row.gravity_mode
    if row.pulse is not None:
        # sha256 over the raw table bytes binds the profile the solver saw.
        payload["pulse"] = dict(row.pulse)
    if row.notes:
        payload["notes"] = list(row.notes)
    return payload


def _plateau_manifest(plateau: DtPlateauReport | None) -> dict[str, Any] | None:
    if plateau is None:
        return None
    return {
        "row": plateau.row,
        "gravity_g": plateau.gravity_g,
        "base_timestep_s": plateau.base_timestep_s,
        "base_rise_mm": lox._finite_or_none(plateau.base_rise_mm),
        "max_abs_delta_percent": lox._finite_or_none(plateau.max_abs_delta_percent),
        "points": [
            {
                "timestep_s": point.timestep_s,
                "rise_mm": lox._finite_or_none(point.rise_mm),
                "delta_percent": lox._finite_or_none(point.delta_percent),
            }
            for point in plateau.points
        ],
    }


def build_manifest(
    result: VentStudyResult, *, solver_describe: str | None = None
) -> dict[str, Any]:
    """Build the hash-bound study manifest payload."""

    config = result.config
    if solver_describe is None:
        solver_describe = _git_describe()
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "version": SCHEMA_VERSION,
        "study_name": config.name,
        "geometry_npz": _repo_relative(config.geometry_path),
        "geometry_npz_sha256": lox.sha256_file(config.geometry_path),
        "harness_module": _repo_relative(lox.LOX_VENT_MODULE_PATH),
        "harness_module_sha256": lox.sha256_file(lox.LOX_VENT_MODULE_PATH),
        "analysis_module": _repo_relative(ANALYSIS_MODULE_PATH),
        "analysis_module_sha256": lox.sha256_file(ANALYSIS_MODULE_PATH),
        "assumptions_doc": _repo_relative(ASSUMPTIONS_DOC_PATH),
        "assumptions_doc_sha256": lox.sha256_file(ASSUMPTIONS_DOC_PATH),
        "solver_describe": solver_describe,
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "numba": numba.__version__,
        },
        "case_definition": {
            "source": _repo_relative(ASSUMPTIONS_DOC_PATH),
            "description": config.description,
            "fluid": config.fluid,
            "fill_fraction": config.fill_fraction,
            "fill_liters": config.fill_liters,
            "initial_pressure_psia": config.initial_pressure_psia,
            "final_pressure_psia": config.final_pressure_psia,
            "vent_rate_lbm_s": config.vent_rate_lbm_s,
            "vent_rate_g_s": config.vent_rate_g_s,
            "duration_s": config.duration_s,
            "timestep_s_primary": config.timestep_s,
            "epsilon": "height_dep",
            "tank_diameter_ft": config.tank_diameter_ft,
            "tank_height_ft": config.tank_height_ft,
            "gravity_source": config.gravity_source,
            "spaceflight": (
                config.spaceflight.as_dict()
                if config.spaceflight is not None
                else None
            ),
            "gravity_matrix_g": config.gravity_matrix,
        },
        "dt_plateau": _plateau_manifest(result.dt_plateau),
        "results": {row.name: _row_manifest(row) for row in result.rows},
        "assumptions_section": f"§{ASSUMPTIONS_SECTION}",
    }
    if lox.FLUID_STEP_PATH.is_file():
        payload["fluid_step"] = _repo_relative(lox.FLUID_STEP_PATH)
        payload["fluid_step_sha256"] = lox.sha256_file(lox.FLUID_STEP_PATH)
    if config.source_path is not None and config.source_path.is_file():
        payload["config_source"] = _repo_relative(config.source_path)
        payload["config_sha256"] = lox.sha256_file(config.source_path)
    return payload


def write_manifest(
    result: VentStudyResult,
    path: str | Path,
    *,
    solver_describe: str | None = None,
) -> dict[str, Any]:
    """Write the study manifest, refusing a dirty git worktree (F8)."""

    if _git_worktree_is_dirty():
        raise RuntimeError(
            "Refusing to write vent-study manifest from a dirty git worktree "
            "(F8 provenance guard). Commit or stash changes first."
        )
    target = Path(path)
    payload = build_manifest(result, solver_describe=solver_describe)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return payload


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m liqlev.analysis.vent_study",
        description=(
            "Re-run a LIQLEV vent study from a JSON config, print the results "
            "table with the §5 assumptions block, and write a hash-bound "
            "manifest."
        ),
    )
    parser.add_argument("--config", required=True, help="study config JSON path")
    parser.add_argument(
        "--output-dir",
        default=None,
        help=f"manifest/report directory (default: {_repo_relative(DEFAULT_OUTPUT_DIR)})",
    )
    parser.add_argument(
        "--gravity",
        action="append",
        default=None,
        metavar="ROW",
        help="only run these gravity rows (repeatable or comma-separated)",
    )
    parser.add_argument(
        "--duration-s", type=float, default=None, help="override duration_s"
    )
    parser.add_argument(
        "--timestep-s", type=float, default=None, help="override timestep_s"
    )
    parser.add_argument(
        "--skip-dt-plateau",
        action="store_true",
        help="skip the automatic dt/2 and dt/4 plateau re-runs",
    )
    parser.add_argument(
        "--no-manifest",
        action="store_true",
        help="print the report without writing manifest/report files",
    )
    return parser


def _apply_overrides(
    config: VentStudyConfig, args: argparse.Namespace
) -> VentStudyConfig:
    if args.gravity:
        requested: list[str] = []
        for value in args.gravity:
            requested.extend(part.strip() for part in value.split(",") if part.strip())
        config = select_rows(config, requested)
    updates: dict[str, Any] = {}
    if args.duration_s is not None:
        updates["duration_s"] = _positive_float(args.duration_s, "--duration-s")
    if args.timestep_s is not None:
        updates["timestep_s"] = _positive_float(args.timestep_s, "--timestep-s")
    if args.skip_dt_plateau:
        updates["dt_plateau_enabled"] = False
    if updates:
        config = replace(config, **updates)
    if config.timestep_s >= config.duration_s:
        raise VentStudyConfigError(
            f"timestep_s ({config.timestep_s}) must be smaller than duration_s "
            f"({config.duration_s})"
        )
    if config.has_pulsed_row:
        check_pulse_resolution(config.spaceflight, config.timestep_s)
    return config


def _console_safe(text: str, stream: TextIO) -> str:
    """Down-convert characters the console encoding cannot represent."""

    encoding = getattr(stream, "encoding", None)
    if not encoding:
        return text
    try:
        text.encode(encoding)
    except (UnicodeEncodeError, LookupError):
        return text.encode(encoding, errors="replace").decode(encoding, "replace")
    return text


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point. Returns a process exit code."""

    args = build_parser().parse_args(argv)
    try:
        config = _apply_overrides(load_study_config(args.config), args)
    except VentStudyConfigError as exc:
        print(f"config error: {exc}", file=sys.stderr)
        return 2

    result = run_study(config, stream=sys.stderr)
    report = format_report(result)
    print(_console_safe(report, sys.stdout), end="")

    if args.no_manifest:
        return 0

    output_dir = (
        Path(args.output_dir).resolve()
        if args.output_dir is not None
        else DEFAULT_OUTPUT_DIR
    )
    manifest_path = output_dir / f"{config.name}_manifest.json"
    try:
        write_manifest(result, manifest_path)
    except RuntimeError as exc:
        print(f"manifest refused: {exc}", file=sys.stderr)
        return 3
    (output_dir / f"{config.name}_report.txt").write_text(report, encoding="utf-8")
    print(f"[vent-study] manifest written: {manifest_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry
    raise SystemExit(main())
