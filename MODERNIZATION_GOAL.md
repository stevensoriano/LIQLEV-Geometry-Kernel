# LIQLEV Modernization Goal

## Objective

Modernize the LIQLEV cryogenic liquid level rise application into a fast, robust, engineer-ready cryovent analysis tool with a cutting-edge desktop GUI, richer physics visualizations, reliable exports, and strict preservation of existing solver results.

## Non-Negotiable Requirement

Core physics and computed results must not change during the application modernization. Any change that affects `core.py`, `thermo_utils.py`, solver inputs, result columns, units, or numeric outputs must be guarded by baseline comparison and explicitly called out.

## Final Definition of Done

- A modern desktop GUI exists as the primary user experience.
- The solver can be run headlessly from scripts/tests and interactively from the GUI.
- The application validates inputs before a run and gives engineers precise correction messages.
- Engineers can run sweeps, Monte Carlo analysis, CSV gravity profiles, vent-rate profiles, and validation presets.
- Result views include time histories, scenario comparison, convergence health, tank-fill risk, vapor generation, pressure evolution, gravity history, and a visual summary of how the vent event evolves.
- PDF/CSV/config exports work from the modern app.
- Existing physics baselines pass after the refactor.
- The legacy GUI remains runnable until the replacement is proven.
- Cluster setup and dependency installation are documented under `/gen/home/sasorian`.

## First Commands On The Cluster

Run these from the uploaded repository root:

```bash
python -m venv /gen/home/sasorian/liqlev-modern-venv
source /gen/home/sasorian/liqlev-modern-venv/bin/activate
python -m pip install --upgrade pip
PIP_CACHE_DIR=/gen/home/sasorian/.cache/pip python -m pip install -r requirements.txt -r requirements-modernization.txt
python scripts/write_physics_baseline.py
python scripts/check_physics_baseline.py
```

If GUI packages are not available on the headless cluster, still complete the engine extraction and validation work there. GUI launch can be verified later on a desktop or through X11/Wayland forwarding.

## Progressive Phases

1. Baseline and safety gates
   - Generate numeric baselines from the current code.
   - Add tests around parsing, unit conversion, input validation, and result schema.
   - Fix the known theme-toggle bug.

2. Architecture extraction
   - Move non-GUI logic out of `gui.py`.
   - Create typed simulation configuration models.
   - Create a headless runner for single cases, sweeps, Monte Carlo, and export.
   - Preserve the old GUI entry point.

3. Modern GUI foundation
   - Add a new PySide6-based app entry point.
   - Implement a workflow layout: Setup, Profiles, Run, Results, Export.
   - Keep solver execution off the UI thread.

4. Engineer visualizations
   - Add comparison plots for pressure, level rise, vent rate, epsilon, gravity, vapor generation, and convergence.
   - Add an event-evolution panel that shows tank state over time: liquid height, ullage, vent timeline, pressure drop, and warning thresholds.
   - Add scenario filters and linked cursors for comparing sweeps.

5. Packaging and release readiness
   - Add a reproducible build path.
   - Validate on the cluster and on a desktop.
   - Document launch, dependency, and troubleshooting steps.

## Stop Conditions

Stop and investigate before continuing if:

- Baseline comparisons fail.
- A solver result column is renamed, removed, or unit-shifted.
- GUI refactors require changes to physics formulas.
- Cluster dependency installation writes outside `/gen/home/sasorian`.
- A new visualization disagrees with exported CSV data.

