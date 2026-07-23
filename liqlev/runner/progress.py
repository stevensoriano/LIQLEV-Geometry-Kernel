"""Structured progress events for headless and GUI runners."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from typing import Literal


ProgressKind = Literal[
    "log",
    "run_start",
    "solver_progress",
    "run_complete",
    "warning",
    "error",
    "complete",
]


@dataclass(frozen=True)
class ProgressEvent:
    kind: ProgressKind
    message: str = ""
    fraction: float | None = None
    scenario_key: str | None = None
    run_index: int | None = None
    total_runs: int | None = None
    stats: Mapping[str, float] = field(default_factory=dict)


ProgressCallback = Callable[[ProgressEvent], None]


def emit_progress(callback: ProgressCallback | None, event: ProgressEvent) -> None:
    """Emit an event only when a callback is present."""
    if callback is not None:
        callback(event)
