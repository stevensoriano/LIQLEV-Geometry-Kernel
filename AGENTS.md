# AGENTS.md

## Repository Rules

- Do not include Codex or OpenAI as authors, co-authors, commit authors, or PR authors during GitHub pushes.
- Preserve the LIQLEV solver physics and numeric outputs unless the user explicitly approves a physics change.
- Treat `core.py` and `thermo_utils.py` as physics-critical. Refactor them only after the baseline validation scripts pass and a before/after comparison is saved.
- Keep generated outputs, local environments, build products, and cluster scratch files out of commits.
- Prefer small, reviewable changes with clear validation notes.

## Cluster Environment

- Install project dependencies under `/gen/home/sasorian`, not in system Python locations.
- Prefer a project virtual environment such as `/gen/home/sasorian/liqlev-modern-venv`.
- Keep package caches under `/gen/home/sasorian/.cache` when possible.
- Assume the cluster may be headless. GUI work should be testable without a display, and desktop launch should be documented separately for X11/Wayland forwarding or local execution.

## Modernization Target

The final application should be a modern, fast, engineer-facing cryogenic vent analysis tool. The look should be high-contrast, precise, and operations-console inspired: matte dark surfaces, thin technical lines, restrained accent colors, dense readable controls, and strong visual hierarchy. Do not name any specific aerospace company as the design source.

## Physics Preservation Gate

Before UI or architecture work changes behavior:

1. Run `python scripts/write_physics_baseline.py`.
2. Commit or otherwise preserve the generated baseline outside the working refactor.
3. After each phase, run `python scripts/check_physics_baseline.py`.
4. If the check fails, stop and classify the change as either a bug, environment drift, or an explicitly approved physics change.

## Recommended Work Order

1. Fix known GUI reliability bugs without changing solver output.
2. Add tests and baseline checks.
3. Extract typed configuration, validation, runner, and reporting modules.
4. Keep the existing GUI working while the new GUI is built.
5. Build the new PySide6 interface behind a separate entry point.
6. Migrate users only after result parity, packaging, and engineer workflow checks pass.

