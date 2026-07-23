"""Modern PySide6 entry point for LIQLEV cryovent analysis."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pyqtgraph as pg
from PySide6.QtCore import QObject, Qt, QThread, Signal, Slot
from PySide6.QtGui import QAction
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from PySide6.QtWidgets import (
    QApplication,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QSplitter,
    QStackedWidget,
    QStyle,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from liqlev.io.config_json import save_simulation_config
from liqlev.model.config import (
    EpsilonConfig,
    FluidConfig,
    GravityProfileConfig,
    RunControls,
    SimulationConfig,
    TankConfig,
    VentProfileConfig,
)
from liqlev.model.parsing import parse_numeric_array
from liqlev.model.validation import InputValidationError, validate_simulation_config
from liqlev.runner.monte_carlo import (
    MonteCarloRequest,
    MonteCarloResult,
    run_monte_carlo,
)
from liqlev.runner.progress import ProgressEvent
from liqlev.runner.sweep import SweepResult, run_sweep
from liqlev.viz.datasets import (
    boundary_layer_traces,
    convergence_traces,
    event_evolution_traces,
    pressure_level_trace,
    time_trace,
)
from liqlev.viz.summaries import summary_rows


APP_QSS = """
QMainWindow, QWidget {
    background: #06080b;
    color: #d8dee9;
    font-family: "Inter", "Segoe UI", "Arial";
    font-size: 12px;
}
QFrame#TopBar {
    background: #090d12;
    border-bottom: 1px solid #263140;
}
QLabel#AppTitle {
    color: #7dd3fc;
    font-size: 18px;
    font-weight: 700;
    letter-spacing: 0px;
}
QLabel#StatusPill {
    background: #102016;
    color: #8ff0a4;
    border: 1px solid #245c39;
    border-radius: 2px;
    padding: 4px 10px;
    font-weight: 700;
}
QListWidget {
    background: #090d12;
    border: 0;
    outline: 0;
}
QListWidget::item {
    border-left: 2px solid transparent;
    padding: 10px 12px;
    min-height: 28px;
}
QListWidget::item:selected {
    background: #111923;
    border-left: 2px solid #38bdf8;
    color: #f8fafc;
}
QGroupBox {
    border: 1px solid #263140;
    border-radius: 4px;
    margin-top: 18px;
    padding: 12px 10px 10px 10px;
    background: #0b1016;
    font-weight: 700;
    color: #cbd5e1;
}
QGroupBox::title {
    subcontrol-origin: margin;
    left: 10px;
    padding: 0 4px;
}
QLineEdit, QComboBox, QTextEdit, QTableWidget {
    background: #05070a;
    color: #e5eef7;
    border: 1px solid #263140;
    border-radius: 3px;
    selection-background-color: #155e75;
}
QLineEdit, QComboBox {
    min-height: 26px;
    padding: 2px 6px;
}
QPushButton {
    background: #101823;
    color: #e5eef7;
    border: 1px solid #2b394a;
    border-radius: 3px;
    min-height: 28px;
    padding: 4px 10px;
    font-weight: 600;
}
QPushButton:hover {
    background: #172233;
    border-color: #38bdf8;
}
QPushButton#PrimaryButton {
    background: #14532d;
    border-color: #22c55e;
}
QPushButton#DangerButton {
    background: #4a1018;
    border-color: #ef4444;
}
QTabWidget::pane {
    border: 1px solid #263140;
}
QTabBar::tab {
    background: #090d12;
    color: #94a3b8;
    padding: 7px 12px;
    border: 1px solid #263140;
}
QTabBar::tab:selected {
    color: #e5eef7;
    background: #101823;
    border-bottom-color: #38bdf8;
}
QHeaderView::section {
    background: #111923;
    color: #cbd5e1;
    border: 0;
    border-right: 1px solid #263140;
    padding: 5px;
}
"""


class RunnerWorker(QObject):
    progress = Signal(object)
    finished = Signal(object)
    failed = Signal(str)

    def __init__(
        self,
        config: SimulationConfig,
        *,
        mode: str = "sweep",
        mc_request: MonteCarloRequest | None = None,
    ):
        super().__init__()
        self._config = config
        self._mode = mode
        self._mc_request = mc_request
        self._cancel_requested = False

    @Slot()
    def cancel(self) -> None:
        self._cancel_requested = True

    def _emit_progress(self, event: ProgressEvent) -> None:
        if self._cancel_requested:
            raise RuntimeError("Run cancelled.")
        self.progress.emit(event)

    @Slot()
    def run(self) -> None:
        try:
            if self._mode == "monte_carlo":
                if self._mc_request is None:
                    raise ValueError("Monte Carlo request is missing.")
                result = run_monte_carlo(
                    self._config,
                    self._mc_request,
                    progress_cb=self._emit_progress,
                )
            else:
                result = run_sweep(self._config, progress_cb=self._emit_progress)
        except Exception as exc:
            self.failed.emit(str(exc))
            return
        self.finished.emit(result)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("LIQLEV Cryovent Analysis Console")
        self.resize(1500, 920)
        self.setMinimumSize(1160, 720)
        self._thread: QThread | None = None
        self._worker: RunnerWorker | None = None
        self._result: SweepResult | None = None
        self._mc_result: MonteCarloResult | None = None

        pg.setConfigOptions(antialias=True, background="#05070a", foreground="#d8dee9")

        central = QWidget()
        root = QVBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)
        root.addWidget(self._build_top_bar())
        root.addWidget(self._build_body(), 1)
        self.setCentralWidget(central)
        self._connect_actions()
        self._update_config_preview()

    def _build_top_bar(self) -> QWidget:
        bar = QFrame()
        bar.setObjectName("TopBar")
        layout = QHBoxLayout(bar)
        layout.setContentsMargins(14, 8, 14, 8)

        title = QLabel("LIQLEV")
        title.setObjectName("AppTitle")
        subtitle = QLabel("CRYOVENT ANALYSIS CONSOLE")
        subtitle.setStyleSheet("color: #94a3b8; font-weight: 600;")
        self.status = QLabel("READY")
        self.status.setObjectName("StatusPill")

        run_action = QAction("Run", self)
        run_action.triggered.connect(self._start_run)
        self.addAction(run_action)

        self.run_button = QPushButton("Run")
        self.run_button.setObjectName("PrimaryButton")
        self.run_button.setIcon(
            self.style().standardIcon(QStyle.StandardPixmap.SP_MediaPlay)
        )
        self.run_button.clicked.connect(self._start_run)

        layout.addWidget(title)
        layout.addWidget(subtitle)
        layout.addStretch(1)
        layout.addWidget(self.status)
        layout.addWidget(self.run_button)
        return bar

    def _build_body(self) -> QWidget:
        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setChildrenCollapsible(False)
        self.nav = QListWidget()
        self.nav.setFixedWidth(172)
        self.stack = QStackedWidget()
        for name in ["Setup", "Profiles", "Run", "Results", "Export"]:
            QListWidgetItem(name, self.nav)

        self.setup_page = self._build_setup_page()
        self.profiles_page = self._build_profiles_page()
        self.run_page = self._build_run_page()
        self.results_page = self._build_results_page()
        self.export_page = self._build_export_page()
        for page in [
            self.setup_page,
            self.profiles_page,
            self.run_page,
            self.results_page,
            self.export_page,
        ]:
            self.stack.addWidget(page)

        self.nav.currentRowChanged.connect(self.stack.setCurrentIndex)
        self.nav.setCurrentRow(0)
        splitter.addWidget(self.nav)
        splitter.addWidget(self.stack)
        splitter.setStretchFactor(1, 1)
        return splitter

    def _group(self, title: str) -> tuple[QGroupBox, QFormLayout]:
        group = QGroupBox(title)
        form = QFormLayout(group)
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
        form.setHorizontalSpacing(10)
        form.setVerticalSpacing(8)
        return group, form

    def _line(self, value: str) -> QLineEdit:
        field = QLineEdit(value)
        field.textChanged.connect(self._update_config_preview)
        return field

    def _combo(self, values: list[str], value: str) -> QComboBox:
        combo = QComboBox()
        combo.addItems(values)
        combo.setCurrentText(value)
        combo.currentTextChanged.connect(self._update_config_preview)
        return combo

    def _build_setup_page(self) -> QWidget:
        page = QWidget()
        layout = QGridLayout(page)
        layout.setContentsMargins(14, 14, 14, 14)
        layout.setHorizontalSpacing(12)
        layout.setVerticalSpacing(12)

        fluid_group, fluid_form = self._group("Fluid")
        self.preset = self._combo(
            ["AS-203 Default", "Hydrogen Mid Fill", "Nitrogen Custom"], "AS-203 Default"
        )
        self.fluid = self._combo(
            ["Hydrogen", "Nitrogen", "Oxygen", "Methane"], "Hydrogen"
        )
        self.pinit = self._line("19.5")
        self.pfinal = self._line("13.8")
        self.mass_override = self._line("16300.0")
        self.temp_override = self._line("38.3")
        fluid_form.addRow("Preset", self.preset)
        fluid_form.addRow("Fluid", self.fluid)
        fluid_form.addRow("Initial pressure, psia", self.pinit)
        fluid_form.addRow("Final pressure, psia", self.pfinal)
        fluid_form.addRow("Initial mass, lbm", self.mass_override)
        fluid_form.addRow("Initial temperature, R", self.temp_override)

        tank_group, tank_form = self._group("Tank")
        self.dtank = self._line("21.670")
        self.htank = self._line("28.18")
        self.fills = self._line("0.5116")
        tank_form.addRow("Diameter, ft", self.dtank)
        tank_form.addRow("Height, ft", self.htank)
        tank_form.addRow("Fill fractions", self.fills)

        run_group, run_form = self._group("Run Controls")
        self.duration = self._line("400.0")
        self.timestep = self._line("10.0")
        self.threshold = self._line("")
        run_form.addRow("Duration, s", self.duration)
        run_form.addRow("Time step, s", self.timestep)
        run_form.addRow("Risk threshold, dh/h0", self.threshold)

        self.preview = QTextEdit()
        self.preview.setReadOnly(True)
        self.preview.setMaximumHeight(150)

        layout.addWidget(fluid_group, 0, 0)
        layout.addWidget(tank_group, 0, 1)
        layout.addWidget(run_group, 1, 0)
        layout.addWidget(self.preview, 1, 1)
        layout.setRowStretch(2, 1)
        return page

    def _build_profiles_page(self) -> QWidget:
        page = QWidget()
        layout = QGridLayout(page)
        layout.setContentsMargins(14, 14, 14, 14)
        layout.setHorizontalSpacing(12)

        vent_group, vent_form = self._group("Vent Profile")
        self.vent_rates = self._line("3.3069")
        self.ramp_duration = self._line("400.0")
        self.ramp_factor = self._line("1.0")
        self.vent_csv = self._line("")
        vent_browse = QPushButton("Browse")
        vent_browse.setIcon(
            self.style().standardIcon(QStyle.StandardPixmap.SP_DirOpenIcon)
        )
        vent_browse.clicked.connect(lambda: self._browse_file(self.vent_csv))
        vent_csv_row = self._field_with_button(self.vent_csv, vent_browse)
        vent_form.addRow("Rates, lbm/s", self.vent_rates)
        vent_form.addRow("Ramp duration, s", self.ramp_duration)
        vent_form.addRow("Ramp factor", self.ramp_factor)
        vent_form.addRow("CSV", vent_csv_row)

        eps_group, eps_form = self._group("Epsilon")
        self.eps_mode = self._combo(
            ["height_dep", "bulk_fake", "AS-203 Schedule", "Custom"], "AS-203 Schedule"
        )
        self.eps_values = self._line("0.4")
        eps_form.addRow("Mode", self.eps_mode)
        eps_form.addRow("Custom values", self.eps_values)

        gravity_group, gravity_form = self._group("Gravity")
        self.gravity_mode = self._combo(
            ["Constant", "Function of Time", "CSV Profile"], "Constant"
        )
        self.constant_g = self._line("0.00000963")
        self.gravity_expr = self._line("0.001")
        self.gravity_csv = self._line("")
        self.hold_g = self._line("0.0014")
        gravity_browse = QPushButton("Browse")
        gravity_browse.setIcon(
            self.style().standardIcon(QStyle.StandardPixmap.SP_DirOpenIcon)
        )
        gravity_browse.clicked.connect(lambda: self._browse_file(self.gravity_csv))
        gravity_form.addRow("Mode", self.gravity_mode)
        gravity_form.addRow("Constant g", self.constant_g)
        gravity_form.addRow("g(t)", self.gravity_expr)
        gravity_form.addRow(
            "CSV", self._field_with_button(self.gravity_csv, gravity_browse)
        )
        gravity_form.addRow("Hold g", self.hold_g)

        layout.addWidget(vent_group, 0, 0)
        layout.addWidget(eps_group, 0, 1)
        layout.addWidget(gravity_group, 1, 0, 1, 2)
        layout.setRowStretch(2, 1)
        return page

    def _field_with_button(self, field: QLineEdit, button: QPushButton) -> QWidget:
        widget = QWidget()
        layout = QHBoxLayout(widget)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)
        layout.addWidget(field, 1)
        layout.addWidget(button)
        return widget

    def _build_run_page(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(14, 14, 14, 14)
        control_row = QHBoxLayout()
        self.run_page_button = QPushButton("Run Simulation")
        self.run_page_button.setObjectName("PrimaryButton")
        self.run_page_button.setIcon(
            self.style().standardIcon(QStyle.StandardPixmap.SP_MediaPlay)
        )
        self.run_page_button.clicked.connect(self._start_run)
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.setObjectName("DangerButton")
        self.cancel_button.setEnabled(False)
        self.cancel_button.setIcon(
            self.style().standardIcon(QStyle.StandardPixmap.SP_BrowserStop)
        )
        self.cancel_button.clicked.connect(self._cancel_run)
        self.progress_label = QLabel("0%")
        self.progress_label.setAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        control_row.addWidget(self.run_page_button)
        control_row.addWidget(self.cancel_button)
        control_row.addStretch(1)
        control_row.addWidget(self.progress_label)

        mc_group, mc_form = self._group("Monte Carlo")
        self.mc_n = self._line("50")
        self.mc_vent_min = self._line("0.001")
        self.mc_vent_max = self._line("0.005")
        self.mc_fill_min = self._line("0.3")
        self.mc_fill_max = self._line("0.7")
        self.mc_grav_min = self._line("0.0005")
        self.mc_grav_max = self._line("0.005")
        self.mc_seed = self._line("")
        self.mc_button = QPushButton("Run Monte Carlo")
        self.mc_button.setIcon(
            self.style().standardIcon(QStyle.StandardPixmap.SP_MediaPlay)
        )
        self.mc_button.clicked.connect(self._start_monte_carlo)
        mc_form.addRow("Samples", self.mc_n)
        mc_form.addRow("Vent min, lbm/s", self.mc_vent_min)
        mc_form.addRow("Vent max, lbm/s", self.mc_vent_max)
        mc_form.addRow("Fill min", self.mc_fill_min)
        mc_form.addRow("Fill max", self.mc_fill_max)
        mc_form.addRow("Gravity min, g", self.mc_grav_min)
        mc_form.addRow("Gravity max, g", self.mc_grav_max)
        mc_form.addRow("Seed", self.mc_seed)
        mc_form.addRow("", self.mc_button)
        self.log = QTextEdit()
        self.log.setReadOnly(True)
        layout.addLayout(control_row)
        layout.addWidget(mc_group)
        layout.addWidget(self.log, 1)
        return page

    def _build_results_page(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(14, 14, 14, 14)
        top = QHBoxLayout()
        self.scenario_combo = QComboBox()
        self.scenario_combo.currentTextChanged.connect(self._render_results)
        self.plot_combo = QComboBox()
        self.plot_combo.addItems(
            [
                "Pressure",
                "Liquid Level",
                "Event Evolution",
                "Pressure vs Level",
                "Ullage",
                "Vent Rate",
                "Epsilon",
                "Vapor Generation",
                "Gravity",
                "Convergence",
                "Boundary Layer",
            ]
        )
        self.plot_combo.currentTextChanged.connect(self._render_results)
        top.addWidget(QLabel("Scenario"))
        top.addWidget(self.scenario_combo, 1)
        top.addWidget(QLabel("View"))
        top.addWidget(self.plot_combo)

        tabs = QTabWidget()
        plot_tab = QWidget()
        plot_layout = QVBoxLayout(plot_tab)
        self.plot = pg.PlotWidget()
        self.plot.showGrid(x=True, y=True, alpha=0.22)
        plot_layout.addWidget(self.plot)
        self.summary = QTextEdit()
        self.summary.setReadOnly(True)
        self.table = QTableWidget()
        self.table.setAlternatingRowColors(True)
        self.table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.ResizeToContents
        )
        tabs.addTab(plot_tab, "Plots")
        tabs.addTab(self.table, "Table")
        tabs.addTab(self.summary, "Summary")

        layout.addLayout(top)
        layout.addWidget(tabs, 1)
        return page

    def _build_export_page(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(14, 14, 14, 14)
        group, form = self._group("Export")
        self.export_dir = self._line(str(Path.cwd() / "results"))
        browse = QPushButton("Browse")
        browse.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_DirOpenIcon))
        browse.clicked.connect(self._browse_export_dir)
        form.addRow("Folder", self._field_with_button(self.export_dir, browse))

        row = QHBoxLayout()
        csv_button = QPushButton("Export CSV")
        csv_button.setIcon(
            self.style().standardIcon(QStyle.StandardPixmap.SP_DialogSaveButton)
        )
        csv_button.clicked.connect(self._export_csv)
        summary_button = QPushButton("Summary CSV")
        summary_button.setIcon(
            self.style().standardIcon(QStyle.StandardPixmap.SP_FileDialogDetailedView)
        )
        summary_button.clicked.connect(self._export_summary_csv)
        pdf_button = QPushButton("PDF Report")
        pdf_button.setIcon(
            self.style().standardIcon(QStyle.StandardPixmap.SP_FileDialogContentsView)
        )
        pdf_button.clicked.connect(self._export_pdf)
        image_button = QPushButton("Plot PNG")
        image_button.setIcon(
            self.style().standardIcon(QStyle.StandardPixmap.SP_DialogSaveButton)
        )
        image_button.clicked.connect(self._export_plot_image)
        config_button = QPushButton("Save Config")
        config_button.setIcon(
            self.style().standardIcon(QStyle.StandardPixmap.SP_DialogSaveButton)
        )
        config_button.clicked.connect(self._save_config)
        row.addWidget(csv_button)
        row.addWidget(summary_button)
        row.addWidget(pdf_button)
        row.addWidget(image_button)
        row.addWidget(config_button)
        row.addStretch(1)
        layout.addWidget(group)
        layout.addLayout(row)
        layout.addStretch(1)
        return page

    def _connect_actions(self) -> None:
        for field in [
            self.pinit,
            self.pfinal,
            self.mass_override,
            self.temp_override,
            self.dtank,
            self.htank,
            self.fills,
            self.duration,
            self.timestep,
            self.threshold,
            self.vent_rates,
            self.ramp_duration,
            self.ramp_factor,
            self.vent_csv,
            self.eps_values,
            self.constant_g,
            self.gravity_expr,
            self.gravity_csv,
            self.hold_g,
            self.mc_n,
            self.mc_vent_min,
            self.mc_vent_max,
            self.mc_fill_min,
            self.mc_fill_max,
            self.mc_grav_min,
            self.mc_grav_max,
            self.mc_seed,
        ]:
            field.textChanged.connect(self._update_config_preview)
        self.preset.currentTextChanged.connect(self._apply_preset)

    def _apply_preset(self, name: str) -> None:
        presets = {
            "AS-203 Default": {
                "fluid": "Hydrogen",
                "pinit": "19.5",
                "pfinal": "13.8",
                "mass": "16300.0",
                "temp": "38.3",
                "dtank": "21.670",
                "htank": "28.18",
                "fills": "0.5116",
                "duration": "400.0",
                "dt": "10.0",
                "vent": "3.3069",
                "eps_mode": "AS-203 Schedule",
                "gravity": "0.00000963",
            },
            "Hydrogen Mid Fill": {
                "fluid": "Hydrogen",
                "pinit": "19.5",
                "pfinal": "13.8",
                "mass": "",
                "temp": "",
                "dtank": "21.670",
                "htank": "28.18",
                "fills": "0.50",
                "duration": "300.0",
                "dt": "10.0",
                "vent": "0.0015",
                "eps_mode": "height_dep",
                "gravity": "0.001",
            },
            "Nitrogen Custom": {
                "fluid": "Nitrogen",
                "pinit": "30.0",
                "pfinal": "20.0",
                "mass": "",
                "temp": "",
                "dtank": "5.0",
                "htank": "10.0",
                "fills": "0.45",
                "duration": "160.0",
                "dt": "5.0",
                "vent": "0.02",
                "eps_mode": "Custom",
                "gravity": "0.002",
            },
        }
        preset = presets[name]
        self.fluid.setCurrentText(preset["fluid"])
        self.pinit.setText(preset["pinit"])
        self.pfinal.setText(preset["pfinal"])
        self.mass_override.setText(preset["mass"])
        self.temp_override.setText(preset["temp"])
        self.dtank.setText(preset["dtank"])
        self.htank.setText(preset["htank"])
        self.fills.setText(preset["fills"])
        self.duration.setText(preset["duration"])
        self.timestep.setText(preset["dt"])
        self.vent_rates.setText(preset["vent"])
        self.ramp_duration.setText(preset["duration"])
        self.ramp_factor.setText("1.0")
        self.eps_mode.setCurrentText(preset["eps_mode"])
        self.eps_values.setText("0.4")
        self.gravity_mode.setCurrentText("Constant")
        self.constant_g.setText(preset["gravity"])
        self._update_config_preview()

    def _browse_file(self, target: QLineEdit) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Select CSV", "", "CSV Files (*.csv);;All Files (*)"
        )
        if path:
            target.setText(path)

    def _browse_export_dir(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self, "Select Export Folder", self.export_dir.text()
        )
        if path:
            self.export_dir.setText(path)

    def _optional_float(self, field: QLineEdit) -> float | None:
        text = field.text().strip()
        return float(text) if text else None

    def _build_config(self) -> SimulationConfig:
        eps_mode = self.eps_mode.currentText()
        eps_values = (
            tuple(parse_numeric_array(self.eps_values.text()))
            if eps_mode == "Custom"
            else (0.4,)
        )
        threshold = self._optional_float(self.threshold)

        config = SimulationConfig(
            fluid=FluidConfig(
                name=self.fluid.currentText(),
                initial_pressure_psia=float(self.pinit.text()),
                final_pressure_psia=float(self.pfinal.text()),
                initial_mass_lbm=self._optional_float(self.mass_override),
                initial_temperature_r=self._optional_float(self.temp_override),
            ),
            tank=TankConfig(
                diameter_ft=float(self.dtank.text()),
                height_ft=float(self.htank.text()),
                fill_fractions=tuple(parse_numeric_array(self.fills.text())),
            ),
            vent=VentProfileConfig(
                rates_lbm_s=tuple(parse_numeric_array(self.vent_rates.text())),
                ramp_duration_s=float(self.ramp_duration.text()),
                ramp_target_factor=float(self.ramp_factor.text()),
                csv_path=self.vent_csv.text().strip(),
            ),
            gravity=GravityProfileConfig(
                mode=self.gravity_mode.currentText(),
                constant_g=float(self.constant_g.text()),
                expression=self.gravity_expr.text().strip(),
                csv_path=self.gravity_csv.text().strip(),
                hold_g=float(self.hold_g.text()),
            ),
            epsilon=EpsilonConfig(mode=eps_mode, values=eps_values),
            run=RunControls(
                duration_s=float(self.duration.text()),
                timestep_s=float(self.timestep.text()),
                threshold_dh_h0=threshold,
            ),
        )
        validate_simulation_config(config)
        return config

    def _update_config_preview(self) -> None:
        try:
            fills = parse_numeric_array(self.fills.text())
            eps_count = (
                len(parse_numeric_array(self.eps_values.text()))
                if self.eps_mode.currentText() == "Custom"
                else 1
            )
            vent_count = (
                1
                if self.vent_csv.text().strip()
                else len(parse_numeric_array(self.vent_rates.text()))
            )
            runs = len(fills) * eps_count * vent_count
            self.preview.setPlainText(
                "\n".join(
                    [
                        f"Run count: {runs}",
                        f"Fluid: {self.fluid.currentText()}",
                        f"Pressure: {self.pinit.text()} -> {self.pfinal.text()} psia",
                        f"Tank: D={self.dtank.text()} ft, H={self.htank.text()} ft",
                        f"Duration: {self.duration.text()} s, dt={self.timestep.text()} s",
                    ]
                )
            )
        except Exception as exc:
            self.preview.setPlainText(str(exc))

    def _start_run(self) -> None:
        if self._thread is not None:
            return
        try:
            config = self._build_config()
        except (ValueError, InputValidationError) as exc:
            QMessageBox.critical(self, "Input Error", str(exc))
            return

        self._begin_worker(config, mode="sweep")

    def _start_monte_carlo(self) -> None:
        if self._thread is not None:
            return
        try:
            config = self._build_config()
            request = MonteCarloRequest(
                n=int(self.mc_n.text()),
                vent_min_lbm_s=float(self.mc_vent_min.text()),
                vent_max_lbm_s=float(self.mc_vent_max.text()),
                fill_min=float(self.mc_fill_min.text()),
                fill_max=float(self.mc_fill_max.text()),
                gravity_min_g=float(self.mc_grav_min.text()),
                gravity_max_g=float(self.mc_grav_max.text()),
                seed=int(self.mc_seed.text()) if self.mc_seed.text().strip() else None,
            )
        except (ValueError, InputValidationError) as exc:
            QMessageBox.critical(self, "Input Error", str(exc))
            return

        self._begin_worker(config, mode="monte_carlo", mc_request=request)

    def _begin_worker(
        self,
        config: SimulationConfig,
        *,
        mode: str,
        mc_request: MonteCarloRequest | None = None,
    ) -> None:
        self.nav.setCurrentRow(2)
        self.log.clear()
        self.progress_label.setText("0%")
        self.status.setText("RUNNING")
        self.run_button.setEnabled(False)
        self.run_page_button.setEnabled(False)
        self.mc_button.setEnabled(False)
        self.cancel_button.setEnabled(True)

        self._thread = QThread(self)
        self._worker = RunnerWorker(config, mode=mode, mc_request=mc_request)
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run)
        self._worker.progress.connect(self._on_progress)
        self._worker.finished.connect(self._on_finished)
        self._worker.failed.connect(self._on_failed)
        self._worker.finished.connect(self._thread.quit)
        self._worker.failed.connect(self._thread.quit)
        self._thread.finished.connect(self._thread.deleteLater)
        self._thread.finished.connect(self._clear_worker)
        self._thread.start()

    def _cancel_run(self) -> None:
        if self._worker is not None:
            self.status.setText("CANCELLING")
            self._worker.cancel()

    @Slot(object)
    def _on_progress(self, event: ProgressEvent) -> None:
        if event.message:
            self.log.append(event.message)
        if event.kind == "solver_progress" and event.stats:
            stats = event.stats
            if "sim_time" in stats:
                self.log.append(
                    "t={sim_time:7.2f}s  P={pressure:8.3f} psia  dh/h0={dh_h0:+.5f}".format(
                        **stats
                    )
                )
        if event.fraction is not None:
            self.progress_label.setText(f"{int(event.fraction * 100)}%")

    @Slot(object)
    def _on_finished(self, result: object) -> None:
        self.status.setText("READY")
        self.run_button.setEnabled(True)
        self.run_page_button.setEnabled(True)
        self.mc_button.setEnabled(True)
        self.cancel_button.setEnabled(False)
        self.progress_label.setText("100%")
        if isinstance(result, MonteCarloResult):
            self._mc_result = result
            self.log.append(
                f"Monte Carlo complete ({result.n} samples, {result.elapsed_s:.2f}s)"
            )
            self._load_monte_carlo_results()
        else:
            self._result = result
            self.log.append(
                f"Completed {result.run_count} run(s) in {result.elapsed_s:.2f}s"
            )
            self._load_results()
        self.nav.setCurrentRow(3)

    @Slot(str)
    def _on_failed(self, message: str) -> None:
        self.status.setText("CANCELLED" if "cancelled" in message.lower() else "ERROR")
        self.run_button.setEnabled(True)
        self.run_page_button.setEnabled(True)
        self.mc_button.setEnabled(True)
        self.cancel_button.setEnabled(False)
        self.log.append(message)
        if "cancelled" not in message.lower():
            QMessageBox.critical(self, "Run Error", message)

    @Slot()
    def _clear_worker(self) -> None:
        self._thread = None
        self._worker = None

    def _load_results(self) -> None:
        if self._result is None:
            return
        self.scenario_combo.blockSignals(True)
        self.scenario_combo.clear()
        self.scenario_combo.addItems(list(self._result.scenarios.keys()))
        self.scenario_combo.blockSignals(False)
        self._render_results()
        self._render_summary()

    def _load_monte_carlo_results(self) -> None:
        if self._mc_result is None:
            return
        self._result = None
        self.scenario_combo.blockSignals(True)
        self.scenario_combo.clear()
        self.scenario_combo.addItem("Monte Carlo")
        self.scenario_combo.blockSignals(False)
        self.plot.clear()
        categories, edges = pd.cut(
            self._mc_result.all_dh,
            bins=min(30, max(10, self._mc_result.n // 5)),
            retbins=True,
        )
        histogram = pd.Series(categories).value_counts(sort=False)
        centers = [
            (float(edges[i]) + float(edges[i + 1])) / 2 for i in range(len(edges) - 1)
        ]
        widths = [float(edges[i + 1]) - float(edges[i]) for i in range(len(edges) - 1)]
        bar = pg.BarGraphItem(
            x=centers, height=histogram.to_numpy(), width=widths, brush="#38bdf8"
        )
        self.plot.addItem(bar)
        self.plot.setLabel("bottom", "Max dh/h0")
        self.plot.setLabel("left", "Count")
        self.summary.setPlainText(
            "\n".join(
                [
                    f"Samples: {self._mc_result.n}",
                    f"Max dh/h0: {self._mc_result.max_dh:.6f}",
                    f"Mean dh/h0: {self._mc_result.mean_dh:.6f}",
                    f"Std dh/h0: {self._mc_result.std_dh:.6f}",
                    f"95th percentile: {self._mc_result.p95:.6f}",
                    f"99th percentile: {self._mc_result.p99:.6f}",
                    "Worst case: "
                    f"vent={self._mc_result.worst['vent']:.6f} lbm/s, "
                    f"fill={self._mc_result.worst['fill']:.4f}, "
                    f"g={self._mc_result.worst['grav']:.6f}",
                ]
            )
        )

    def _current_scenario(self) -> dict[str, Any] | None:
        if self._result is None:
            return None
        return self._result.scenarios.get(self.scenario_combo.currentText())

    def _render_results(self) -> None:
        scenario = self._current_scenario()
        if not scenario:
            return
        self.plot.clear()
        view = self.plot_combo.currentText()
        colors = ["#38bdf8", "#22c55e", "#f59e0b", "#f472b6", "#a78bfa"]
        threshold = self._optional_float(self.threshold)
        self.plot.setLabel("bottom", "Time", units="s")
        for index, (dataframe, rate) in enumerate(
            zip(scenario["dfs"], scenario["vent_rates"])
        ):
            if dataframe.empty:
                continue
            color = colors[index % len(colors)]
            pen = pg.mkPen(color=color, width=2)
            name = f"{rate:g} lbm/s"
            x = dataframe["Time"].to_numpy()
            if view == "Pressure":
                self.plot.plot(x, dataframe["Press"].to_numpy(), pen=pen, name=name)
                self.plot.setLabel("left", "Pressure", units="psia")
            elif view == "Liquid Level":
                self.plot.plot(x, dataframe["Hratio"].to_numpy(), pen=pen, name=name)
                self.plot.setLabel("left", "dh/h0")
                if threshold is not None:
                    self.plot.addLine(
                        y=threshold,
                        pen=pg.mkPen(color="#ef4444", style=Qt.PenStyle.DashLine),
                    )
            elif view == "Event Evolution":
                for trace_index, trace in enumerate(event_evolution_traces(dataframe)):
                    trace_pen = pg.mkPen(
                        color=colors[(index + trace_index) % len(colors)], width=2
                    )
                    self.plot.plot(trace.x, trace.y, pen=trace_pen, name=trace.name)
                self.plot.setLabel("left", "Normalized State")
            elif view == "Pressure vs Level":
                trace = pressure_level_trace(dataframe)
                self.plot.plot(trace.x, trace.y, pen=pen, name=name)
                self.plot.setLabel("bottom", "dh/h0")
                self.plot.setLabel("left", trace.y_label, units=trace.unit)
            elif view == "Ullage":
                trace = time_trace(
                    dataframe, "Ullage Mass", "Ullage mass", "Mass", "lbm"
                )
                self.plot.plot(trace.x, trace.y, pen=pen, name=name)
                self.plot.setLabel("left", trace.y_label, units=trace.unit)
            elif view == "Vent Rate":
                self.plot.plot(x, dataframe["Vent Rate"].to_numpy(), pen=pen, name=name)
                self.plot.setLabel("left", "Vent Rate", units="lbm/s")
            elif view == "Epsilon":
                self.plot.plot(x, dataframe["eps"].to_numpy(), pen=pen, name=name)
                self.plot.setLabel("left", "Epsilon")
            elif view == "Vapor Generation":
                self.plot.plot(
                    x, dataframe["Vap Gen Rate (kg/s)"].to_numpy(), pen=pen, name=name
                )
                self.plot.setLabel("left", "Vapor Generation", units="kg/s")
            elif view == "Gravity":
                self.plot.plot(x, dataframe["Gravity_g"].to_numpy(), pen=pen, name=name)
                self.plot.setLabel("left", "Gravity", units="g")
            elif view == "Convergence":
                for trace_index, trace in enumerate(convergence_traces(dataframe)):
                    trace_pen = pg.mkPen(
                        color=colors[(index + trace_index) % len(colors)], width=2
                    )
                    self.plot.plot(trace.x, trace.y, pen=trace_pen, name=trace.name)
                self.plot.setLabel("left", "Iterations")
            elif view == "Boundary Layer":
                for trace_index, trace in enumerate(boundary_layer_traces(dataframe)):
                    trace_pen = pg.mkPen(
                        color=colors[(index + trace_index) % len(colors)], width=2
                    )
                    self.plot.plot(trace.x, trace.y, pen=trace_pen, name=trace.name)
                self.plot.setLabel("left", "Boundary Layer Diagnostics")
        self._render_table(scenario["dfs"][0])

    def _render_table(self, dataframe: pd.DataFrame) -> None:
        if dataframe.empty:
            self.table.clear()
            return
        rows = min(len(dataframe), 500)
        self.table.setRowCount(rows)
        self.table.setColumnCount(len(dataframe.columns))
        self.table.setHorizontalHeaderLabels(list(dataframe.columns))
        for row in range(rows):
            for col, column in enumerate(dataframe.columns):
                value = dataframe.iat[row, col]
                text = f"{value:.6g}" if isinstance(value, float) else str(value)
                self.table.setItem(row, col, QTableWidgetItem(text))

    def _render_summary(self) -> None:
        if self._result is None:
            return
        lines = []
        threshold = self._optional_float(self.threshold)
        for row in summary_rows(self._result.scenarios, threshold):
            threshold_text = ""
            if threshold is not None:
                threshold_time = (
                    f"{row['Threshold Time (s)']:.2f}s"
                    if row["Threshold Time (s)"] is not None
                    else "n/a"
                )
                threshold_text = (
                    f"  threshold={row['Threshold Crossed']} at {threshold_time}"
                )
            tank_time = (
                f"{row['Tank Exceeded Time (s)']:.2f}s"
                if row["Tank Exceeded Time (s)"] is not None
                else "n/a"
            )
            lines.append(row["Scenario"])
            lines.append(
                f"  vent={row['Vent Rate (lbm/s)']:g} lbm/s"
                f"  max dh/h0={row['Max dh/h0']:.6f}"
                f"  t={row['Time to Peak (s)']:.2f}s"
                f"  final P={row['Final Pressure (psia)']:.4f} psia"
                f"{threshold_text}"
                f"  tank={row['Tank Exceeded']} at {tank_time}"
                f"  conv_fail={row['Convergence Failures']}"
            )
            lines.append("")
        self.summary.setPlainText("\n".join(lines))

    def _export_csv(self) -> None:
        if self._result is None:
            QMessageBox.information(self, "Export", "No results available.")
            return
        folder = Path(self.export_dir.text())
        folder.mkdir(parents=True, exist_ok=True)
        count = 0
        for key, scenario in self._result.scenarios.items():
            safe = key.replace(" ", "_").replace(",", "").replace("%", "pct")
            for dataframe, rate in zip(scenario["dfs"], scenario["vent_rates"]):
                dataframe.to_csv(folder / f"liqlev_{safe}_vent{rate}.csv", index=False)
                count += 1
        QMessageBox.information(self, "Export", f"Exported {count} CSV file(s).")

    def _export_summary_csv(self) -> None:
        if self._result is None:
            QMessageBox.information(self, "Export", "No results available.")
            return
        folder = Path(self.export_dir.text())
        folder.mkdir(parents=True, exist_ok=True)
        threshold = self._optional_float(self.threshold)
        rows = summary_rows(self._result.scenarios, threshold)
        pd.DataFrame(rows).to_csv(folder / "liqlev_summary.csv", index=False)
        QMessageBox.information(self, "Export", "Summary CSV exported.")

    def _export_plot_image(self) -> None:
        folder = Path(self.export_dir.text())
        folder.mkdir(parents=True, exist_ok=True)
        image = self.plot.grab()
        image.save(str(folder / "liqlev_current_plot.png"))
        QMessageBox.information(self, "Export", "Plot image exported.")

    def _export_pdf(self) -> None:
        if self._result is None and self._mc_result is None:
            QMessageBox.information(self, "Export", "No results available.")
            return
        folder = Path(self.export_dir.text())
        folder.mkdir(parents=True, exist_ok=True)
        path = folder / "liqlev_report.pdf"
        with PdfPages(path) as pdf:
            fig = Figure(figsize=(11, 8.5), dpi=100, facecolor="#06080b")
            ax = fig.add_subplot(111)
            ax.set_facecolor("#06080b")
            ax.axis("off")
            ax.text(
                0.03,
                0.96,
                "LIQLEV CRYOVENT ANALYSIS REPORT",
                color="#7dd3fc",
                fontsize=18,
                fontweight="bold",
                transform=ax.transAxes,
                va="top",
            )
            ax.text(
                0.03,
                0.88,
                self.summary.toPlainText(),
                color="#d8dee9",
                fontsize=8,
                family="monospace",
                transform=ax.transAxes,
                va="top",
            )
            pdf.savefig(fig, facecolor=fig.get_facecolor())

            if self._result is not None:
                for scenario_key, scenario in self._result.scenarios.items():
                    fig = Figure(figsize=(11, 8.5), dpi=100, facecolor="#06080b")
                    ax = fig.add_subplot(111)
                    ax.set_facecolor("#05070a")
                    for dataframe, rate in zip(scenario["dfs"], scenario["vent_rates"]):
                        if dataframe.empty:
                            continue
                        ax.plot(
                            dataframe["Time"],
                            dataframe["Press"],
                            label=f"{rate:g} lbm/s",
                        )
                    ax.set_title(scenario_key, color="#d8dee9")
                    ax.set_xlabel("Time (s)", color="#d8dee9")
                    ax.set_ylabel("Pressure (psia)", color="#d8dee9")
                    ax.tick_params(colors="#d8dee9")
                    ax.grid(True, alpha=0.25)
                    ax.legend()
                    pdf.savefig(fig, facecolor=fig.get_facecolor())
        QMessageBox.information(self, "Export", "PDF report exported.")

    def _save_config(self) -> None:
        try:
            config = self._build_config()
        except (ValueError, InputValidationError) as exc:
            QMessageBox.critical(self, "Input Error", str(exc))
            return
        folder = Path(self.export_dir.text())
        folder.mkdir(parents=True, exist_ok=True)
        save_simulation_config(config, folder / "liqlev_config.json")
        QMessageBox.information(self, "Export", "Configuration saved.")


def build_app(argv: list[str] | None = None) -> QApplication:
    app = QApplication(argv or [])
    app.setStyleSheet(APP_QSS)
    return app


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--smoke-test",
        action="store_true",
        help="Create the window and exit without showing it.",
    )
    args = parser.parse_args(argv)

    app = build_app(sys.argv[:1])
    window = MainWindow()
    if args.smoke_test:
        app.processEvents()
        return 0
    window.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
