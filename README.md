# LIQLEV-Python: Cryogenic Liquid Level Rise Simulation

A modernized Python framework of the legacy LIQLEV boundary-layer model with a full-featured desktop GUI. This tool predicts transient liquid level rise (LLR) and two-phase fluid behavior during zero-gravity venting of cryogenic propellants.

## Overview

The management of cryogenic propellants under microgravity conditions presents significant challenges, particularly the safe venting of vapor without entraining liquid. When a cryogenic tank is vented to relieve internal pressure, the rapid pressure reduction initiates thermodynamic flashing and bulk boiling. Vapor generated within the superheated liquid bulk and along the wetted tank walls forms bubbles that displace the surrounding liquid, causing the liquid-vapor interface to rise. If this mixture reaches the vent port, liquid entrainment can occur, resulting in propellant loss.

This project is a complete Python reimplementation of the original Fortran LIQLEV model developed following the Saturn S-IVB AS-203 flight experiment[^1]. It builds upon subsequent adaptations for the Human Landing System (HLS) program[^2], upgrading the sequential Excel VBA environment into a robust, parallelized Python framework with an interactive desktop GUI.

## Features

### Simulation & Solver
- **Numba JIT Acceleration** — Core solver loop compiled with Numba (`@njit`) for ~10–30× speedup over pure Python
- **Pre-Computed Property Tables** — CoolProp saturation lookups cached into interpolation tables for ~100× faster thermodynamic evaluation per time step
- **Solver Convergence Tracking** — Flags time steps where the iterative boundary-layer solver did not converge (column 29 in result array)
- **Parametric Sweeps** — Sweep across multiple vent rates, fill fractions, and epsilon values using comma-separated lists, range notation (`0.1:0.1:0.5`), or `linspace(0.1, 0.5, 5)`
- **Monte Carlo Analysis** — Randomized sampling across vent rate, fill fraction, and gravity ranges with histogram output and 1σ/2σ confidence bands

### GUI & Visualization
- **Interactive GUI** — Dark theme built with CustomTkinter
- **Embedded Matplotlib Plots** — Interactive dark-themed plots with crosshair cursor, copy-to-clipboard, and adjustable font scaling (A+/A−)
- **Plot Overlay** — Superimpose results from multiple parameter combinations on a single plot for direct comparison
- **7 Plot Types** — Dimensionless liquid level rise (dh/h₀), liquid level rise (inches), tank pressure, vent rate profile, epsilon, diagnostics (superheat + dP/dt), and vapor generation (BL volume + gravity)
- **Progress Tracking** — Live progress bar with scenario count and percentage
- **Keyboard Shortcuts** — Ctrl+Enter / F5 to run, Esc for pointer mode

### Unit System
- **Dual Unit System** — Toggle between Imperial (psia, ft, lbm/s) and SI (bar, m, kg/s)
- **Live Conversion Display** — Alternate-unit value shown in real time next to each input field
- **Automatic Legend Conversion** — Plot legends display converted vent rate values in the active unit system
- **Solver Transparency** — Solver always operates in British units internally; all conversions handled transparently

### Tank Geometry & Presets
- **Tank Presets** — Quick-load geometry for common tanks:

| Preset | Diameter | Height | Fluid |
|---|---|---|---|
| S-IVB AS-203 LH2 | 6.604 m (21.67 ft) | 8.59 m (28.18 ft) | Hydrogen |
| Centaur LH2 | 3.05 m (10.0 ft) | 9.14 m (30.0 ft) | Hydrogen |
| Centaur LOX | 3.05 m (10.0 ft) | 3.66 m (12.0 ft) | Oxygen |
| SLS Core LH2 | 8.38 m (27.5 ft) | 39.7 m (130.2 ft) | Hydrogen |
| SLS Core LOX | 8.38 m (27.5 ft) | 17.0 m (55.8 ft) | Oxygen |
| MHTB | 3.05 m (10.0 ft) | 3.05 m (10.0 ft) | Hydrogen |

- **Measured Initial Conditions (Validation Overrides)** — Optional override fields for initial liquid mass (lbm) and initial temperature (°R). The solver normally computes these from cylindrical geometry and fill fraction, but real tanks have curved domes that change the usable volume. For flight validation (e.g. AS-203: 16,300 lbm, 38.3 °R), enter the measured values directly. An info popup explains when and why these overrides are needed.

### Epsilon (Boil-off Partitioning)
- **height_dep** — Geometric wetted-area ratio computed dynamically from liquid height
- **bulk_fake** — Large epsilon (~50) to emulate distributed volumetric boiling in large tanks
- **AS-203 Schedule** — Time-varying epsilon schedule from the Saturn S-IVB AS-203 validation case (11-point interpolation table)
- **Custom** — User-specified constant epsilon value

### Gravity Profiles
- **Constant** — Fixed g-level (ft/s²)
- **Custom Expression** — Python math expression as a function of time `t` (supports `sin`, `cos`, `exp`, `log`, `sqrt`, `pi`, etc.)
- **CSV Profile** — Load transient acceleration data from a CSV file (e.g. drop tower, parabolic flight). If the simulation exceeds the CSV duration, the final value is held constant.

### Vent Rate Profiles
- **Constant / Swept** — Single value or comma-separated array of constant vent rates
- **Ramp** — Linear ramp from initial rate to a target factor over a configurable duration
- **CSV Profile** — Load time-varying vent rate schedules from a CSV file

### Export & Reporting
- **PDF Report Generation** — Multi-page PDF with input parameters page and all diagnostic plots
- **CSV Export** — Full time-series results for each scenario
- **Summary CSV** — One-row-per-scenario summary with peak values and convergence failure counts
- **Clipboard Copy** — Copy any plot directly to the Windows clipboard (PNG→BMP conversion)
- **Configuration Save/Load** — Save and restore full simulation configurations as JSON files (includes unit mode, overrides, and all parameters)

## Supported Fluids

- Nitrogen (N₂)
- Hydrogen (H₂) — includes legacy AS-203 polynomial correlations
- Oxygen (O₂)
- Methane (CH₄)

## Technical Details & Mathematical Model

The LIQLEV model explicitly accounts for wall-driven boil-off and its contribution to interface rise by treating the sidewall as a surface on which a growing vapor film forms. Bubbles nucleate, grow, and detach within this film, displacing bulk liquid.

### Piston-Like Displacement
The framework retains the original one-dimensional piston-like displacement assumption, where all vapor generated is treated as uniformly distributed across the tank cross-section. The instantaneous interface location is tracked via the discrete relation:

$$Z_{\rm int}(t) = Z_{\rm int}(t-\Delta t) + \frac{\Delta V_{\rm BL} + \Delta m_{\rm liq}/\rho_{\rm liq}}{A_c}\Delta t$$

where $\Delta V_{\rm BL}$ is the incremental vapor volume generated in the wall thermal boundary layer, $\Delta m_{\rm liq}$ is the mass of liquid flashed to vapor, $\rho_{\rm liq}$ is the saturated liquid density, and $A_c$ is the tank cross-sectional area.

### Boil-Off Partitioning Ratio ($\epsilon$)
The parameter $\epsilon$ represents the instantaneous fraction of total vapor mass generated within the wall thermal boundary layer relative to the vapor generated at the liquid-vapor free surface. In the height-dependent mode, it is computed from the geometric wetted-area ratio:

$$\epsilon(t) = \frac{\pi D h(t)}{\pi D h(t) + \frac{\pi D^{2}}{4}}$$

*Note for Large-Scale Tanks:* For large-volume tanks at high fill fractions, distributed bulk flashing dominates. The model includes a "bulk_fake" mode (setting $\epsilon \approx 50$) which scales the effective vapor-generating surface area to emulate distributed volumetric boiling.

### Thermodynamic Integration
Outdated localized property lookups and curve fits have been replaced with the open-source `CoolProp` library[^3]. This provides consistent, high-accuracy real-fluid properties for Nitrogen, Hydrogen, Oxygen, and Methane across any pressure regime. The critical saturation-curve slope $(dP/dT)_{\rm sat}$ is computed dynamically via a central finite-difference scheme. For Hydrogen, legacy AS-203 polynomial correlations are retained for direct validation comparison.

## Project Structure

```text
LIQLEV-Python-Simulation/
├── data/                 # Input gravity and vent rate profiles (CSV)
├── results/              # Output directory for simulation results
├── gui.py                # Main application — desktop GUI and visualization
├── core.py               # Numba JIT-compiled mathematical solver loop
├── thermo_utils.py       # CoolProp thermodynamic property functions
├── requirements.txt      # Python dependencies
└── README.md             # Documentation
```

## Installation

1. Clone the repository:
```bash
git clone https://github.com/stevensoriano/LIQLEV-Python-Simulation.git
cd LIQLEV-Python-Simulation
```

2. Create a Conda environment (recommended):
```bash
conda create -n liqlev python=3.11
conda activate liqlev
```

3. Install the required dependencies:
```bash
pip install -r requirements.txt
```

> **Note:** CoolProp and Numba require a compatible Python version (3.9–3.12 recommended). A Conda/Miniconda environment is recommended for easiest dependency resolution.

## Usage

Launch the GUI:
```bash
python gui.py
```

### Quick Start (AS-203 Validation)

The program defaults to the Saturn S-IVB AS-203 LH2 validation case on startup:
1. Fluid is set to **Hydrogen** with initial/final pressures of **19.5 / 13.8 psia**
2. Tank geometry loads the **S-IVB AS-203 LH2** preset (D = 6.604 m, H = 8.59 m)
3. Measured initial conditions are pre-filled: **16,300 lbm** liquid mass, **38.3 °R**
4. Epsilon mode is set to **AS-203 Schedule** (time-varying 11-point table)
5. Three vent rates are pre-loaded: **3.3069, 2.2046, 1.1023 lbm/s**
6. Click **RUN SIMULATION** (or Ctrl+Enter) and compare results to published AS-203 data

### General Workflow

1. **Select a Fluid** — Choose from Nitrogen, Hydrogen, Oxygen, or Methane.
2. **Set Tank Geometry** — Enter diameter and height manually, or select a tank preset. Optionally set measured initial conditions for validation.
3. **Configure Venting** — Set initial/final pressure, vent rate(s), and fill fraction(s). Supports arrays for parametric sweeps.
4. **Choose Gravity Profile** — Use a constant g-level, enter a custom math expression, or load a transient CSV file.
5. **Select Epsilon Mode** — Choose height-dependent, bulk_fake, AS-203 Schedule, or custom constant.
6. **Run** — Click **RUN SIMULATION** (or press Ctrl+Enter / F5). Progress is shown in the sidebar.
7. **Analyze Results** — Browse interactive plots (liquid level rise, pressure, diagnostics, vapor generation), view the data table, or read the summary/log tabs.
8. **Export** — Generate a PDF report, export detailed CSV results, or export a one-row-per-scenario summary CSV. Copy individual plots to the clipboard with one click.

### Unit Toggle

Click the **SI / IMPERIAL** toggle to switch between unit systems. All input fields convert automatically, and a live alternate-unit display shows the equivalent value in the other system next to each field. Plot legends also update with converted vent rate values. The solver always operates in British units internally — conversions are handled transparently.

### Monte Carlo Mode

Enable the Monte Carlo tab to run randomized simulations across user-defined ranges for vent rate, fill fraction, and gravity level. Results are displayed as histograms with 1σ and 2σ confidence bands, plus 95th/99th percentile markers.

### Overlay Mode

Select **Overlay** in the plot controls to superimpose results from multiple parameter combinations on a single plot, enabling direct visual comparison across sweep cases.

---

## Contact

Steven Soriano — steven.a.soriano@nasa.gov

---

## References

[^1]: Bradshaw, R. D., "Evaluation and Application of Data from Low-Gravity Orbital Experiment, Phase I Final Report", General Dynamics, Convair Division, NASA CR-109847 (Report No. GDC-DDB-70-003), 1970.
[^2]: Moran, Matt, "Boundary Layer Model: Adaptation for HLS", Internal NASA Document, n.d.
[^3]: Bell, Ian H. et al., "Pure and Pseudo-pure Fluid Thermophysical Property Evaluation and the Open-Source Thermophysical Property Library CoolProp", Industrial & Engineering Chemistry Research, vol. 53, no. 6, pp. 2498-2508, 2014.
