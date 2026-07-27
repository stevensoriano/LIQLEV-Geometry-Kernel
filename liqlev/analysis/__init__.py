"""Config-driven analysis entry points built on the validation machinery.

Modules here re-run committed acceptance cases from JSON configs. They own no
physics: the solver, case builder, ullage-closure metric, physicality
assessment and provenance guards all live in ``validation.lox_vent_cases`` and
``liqlev.runner`` and are reused unchanged.
"""

from __future__ import annotations

__all__ = ["vent_study"]
