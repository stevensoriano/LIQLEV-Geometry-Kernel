# -*- coding: utf-8 -*-
"""
gui.py  –  LIQLEV Cryogenic Liquid Level Rise Simulator (v2)

GitHub-styled dark-theme engineering GUI.
Features:
  - Array / sweep inputs for vent rates, fill fractions, epsilon values
  - Pre-computed property tables for ~100x faster simulation
  - Live run-count preview and progress tracking
  - Embedded matplotlib with matching dark theme
  - Data table, summary, and log tabs
"""

import os, sys, json, time, threading, queue, re, io, ctypes, datetime
import numpy as np
import pandas as pd
import customtkinter as ctk
from tkinter import ttk, filedialog, messagebox
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
import matplotlib.pyplot as _plt
matplotlib.use("TkAgg")

from core import liqlev_simulation
from thermo_utils import Tsat, Psat, DensitySat, build_property_table

# ═════════════════════════════════════════════════════════════════════════════
#  THEME — Dark / Light
# ═════════════════════════════════════════════════════════════════════════════
_THEME_DARK = dict(
    bg          = "#000000",
    surface     = "#0a0a0f",
    overlay     = "#111118",
    border      = "#1e1e2a",
    text        = "#e8eaed",
    text2       = "#9aa0a6",
    muted       = "#5f6368",
    blue        = "#4da6ff",
    green       = "#00c853",
    green_hover = "#00e676",
    red         = "#ff1744",
    orange      = "#ffab00",
    purple      = "#7c4dff",
    header      = "#000000",
    input_bg    = "#05050a",
    btn_secondary     = "#141420",
    btn_secondary_hov = "#1a1a2e",
    accent      = "#005288",
    cyan        = "#18ffff",
)

_THEME_LIGHT = dict(
    bg          = "#f0f2f5",
    surface     = "#ffffff",
    overlay     = "#e8ecf0",
    border      = "#c0c8d0",
    text        = "#1a1a2e",
    text2       = "#3c4450",
    muted       = "#6b7685",
    blue        = "#0066cc",
    green       = "#00873e",
    green_hover = "#00a64a",
    red         = "#cc1430",
    orange      = "#cc8800",
    purple      = "#5c33cc",
    header      = "#ffffff",
    input_bg    = "#f7f8fa",
    btn_secondary     = "#dde2e8",
    btn_secondary_hov = "#cdd3da",
    accent      = "#005288",
    cyan        = "#007a8c",
)

_PLOT_DARK = dict(bg="#000000", face="#05050a", grid="#1a1a2e", text="#e8eaed")
_PLOT_LIGHT = dict(bg="#f0f2f5", face="#ffffff", grid="#c0c8d0", text="#1a1a2e")

# Active theme dicts — mutable, swapped at runtime
GH = dict(_THEME_DARK)

PLOT_BG   = _PLOT_DARK["bg"]
PLOT_FACE = _PLOT_DARK["face"]
PLOT_GRID = _PLOT_DARK["grid"]
PLOT_TEXT  = _PLOT_DARK["text"]
PLOT_COLORS = ["#4da6ff", "#00c853", "#ffab00", "#ff4081",
               "#7c4dff", "#ff6e40", "#18ffff", "#69f0ae"]

FLUIDS = ["Nitrogen", "Oxygen", "Hydrogen", "Methane"]
PSI_TO_KPA = 6.89475729
G_TO_FT_S2 = 32.174
EPSILON_MODES = ["height_dep", "bulk_fake", "AS-203 Schedule", "Custom"]

# ── Tank geometry presets ──
TANK_PRESETS = {
    "Custom": {},
    "S-IVB AS-203 LH2":  {"dtank": "6.604", "htank": "8.59",  "fluid": "Hydrogen",
                          "xmlzro_override": "16300.0", "tinit_override": "38.3"},
    "Centaur LH2":       {"dtank": "3.05",  "htank": "9.14",  "fluid": "Hydrogen"},
    "Centaur LOX":       {"dtank": "3.05",  "htank": "3.66",  "fluid": "Oxygen"},
    "SLS Core LH2":      {"dtank": "8.38",  "htank": "39.7",  "fluid": "Hydrogen"},
    "SLS Core LOX":      {"dtank": "8.38",  "htank": "17.0",  "fluid": "Oxygen"},
    "MHTB":              {"dtank": "3.05",  "htank": "3.05",  "fluid": "Hydrogen"},
}

# ── Font scaling base ──
BASE_FONT = 12

# ── Default field values (for Reset) ──
_DEFAULTS = {
    "fluid": "Hydrogen", "pinit": "19.5", "pfinal": "13.8",
    "dtank": "21.670", "htank": "28.18", "duration": "400.0", "dt": "10.0",
    "vent": "3.3069, 2.2046, 1.1023", "fill": "0.5116",
    "eps_mode": "AS-203 Schedule",
    "eps_custom": "0.4", "ramp_dur": "400.0", "ramp_factor": "1.0",
    "grav_mode": "Constant", "const_g": "0.00000963", "grav_func": "0.001",
    "grav_csv": "./data/5s_drop_tower_extracted_az_positive_data.csv",
    "hold_g": "0.0014", "mc_n": "50", "mc_vent_min": "0.001",
    "mc_vent_max": "0.005", "mc_fill_min": "0.3", "mc_fill_max": "0.7",
    "mc_grav_min": "0.0005", "mc_grav_max": "0.005", "threshold": "",
    "tank_preset": "S-IVB AS-203 LH2",
    "xmlzro_override": "16300.0", "tinit_override": "38.3",
}

# ═════════════════════════════════════════════════════════════════════════════
#  UNIT CONVERSION SYSTEM
# ═════════════════════════════════════════════════════════════════════════════
# Simulation always runs in British units internally.
# SI conversions applied on display only.
UNIT_CONV = {
    # key: (british_unit, si_unit, multiply_by_to_get_si)
    "length_ft":     ("ft",    "m",     0.3048),
    "length_in":     ("in",    "cm",    2.54),
    "pressure":      ("psia",  "bar",   0.0689476),
    "mass_flow":     ("lbm/s", "kg/s",  0.453592),
    "temperature":   ("R",     "K",     1 / 1.8),
    "area":          ("ft\u00b2",   "m\u00b2",    0.092903),
    "volume":        ("ft\u00b3",   "m\u00b3",    0.028317),
    "density":       ("lbm/ft\u00b3", "kg/m\u00b3", 16.0185),
}


def _u(key, si_mode):
    """Return the display unit string for a given key."""
    brit, si, _ = UNIT_CONV[key]
    return si if si_mode else brit


def _cv(value, key, si_mode):
    """Convert a value from British to SI if si_mode is True."""
    if not si_mode:
        return value
    _, _, factor = UNIT_CONV[key]
    return value * factor

PLOT_TYPES = [
    "Liquid Level Rise (dh/h0)",
    "Liquid Level Rise (inches)",
    "Tank Pressure",
    "Vent Rate Profile",
    "Epsilon",
    "Diagnostics",
    "Vapor Generation",
]


# ═════════════════════════════════════════════════════════════════════════════
#  HELPERS
# ═════════════════════════════════════════════════════════════════════════════
def parse_array(text):
    """Parse user input into a list of floats.

    Supports:
      - Comma separated: "0.1, 0.3, 0.5"
      - Range notation:  "0.1:0.1:0.5"  (start:step:stop inclusive)
      - Linspace:        "linspace(0.1, 0.5, 5)"
    """
    text = text.strip()
    if not text:
        return []

    # linspace(start, stop, n)
    m = re.match(r"linspace\(\s*([^,]+),\s*([^,]+),\s*([^)]+)\)", text, re.I)
    if m:
        return np.linspace(float(m.group(1)), float(m.group(2)),
                           int(m.group(3))).tolist()

    # start:step:stop
    parts = text.split(":")
    if len(parts) == 3:
        start, step, stop = float(parts[0]), float(parts[1]), float(parts[2])
        arr = np.arange(start, stop + step * 0.5, step)
        return arr.tolist()

    # Comma separated
    return [float(x.strip()) for x in text.split(",") if x.strip()]


def build_inputs(fluid, pinit_psia, pfinal_psia, dtank, htank, fill_fraction,
                 duration, delta_t, vent_rate, neps, teps, xeps,
                 ramp_duration, ramp_target_factor, nggo, tggo, xggo,
                 gravity_function=None,
                 xmlzro_override=None, tinit_override=None):
    """Build the inputs dict expected by core.liqlev_simulation."""
    volt = (np.pi / 4) * (dtank ** 2) * htank
    ac = 0.7854 * (dtank ** 2)
    htzero = fill_fraction * htank

    press_kpa = pinit_psia * PSI_TO_KPA
    tinit = Tsat(fluid, press_kpa) * 1.8  # K -> Rankine

    # Override tinit with measured value (e.g. AS-203 flight data)
    if tinit_override is not None:
        tinit = tinit_override

    if fluid == "Hydrogen":
        rhol = (0.1709 + 0.7454 * tinit - 0.04421 * tinit ** 2
                + 0.001248 * tinit ** 3 - 1.738e-5 * tinit ** 4
                + 9.424e-8 * tinit ** 5)
    else:
        t_k = tinit / 1.8
        ps_kpa = Psat(fluid, t_k)
        rhol = DensitySat(fluid, "liquid", ps_kpa) * 0.0624279606

    xmlzro = rhol * (htzero * ac)

    # Override xmlzro with measured value (e.g. AS-203 flight data)
    if xmlzro_override is not None:
        xmlzro = xmlzro_override

    if duration > ramp_duration:
        tvmdot = np.array([0.0, ramp_duration, duration])
        xvmdot = np.array([vent_rate,
                           vent_rate * ramp_target_factor,
                           vent_rate * ramp_target_factor])
    else:
        slope = (vent_rate * ramp_target_factor - vent_rate) / ramp_duration
        end_rate = vent_rate + slope * duration
        tvmdot = np.array([0.0, duration])
        xvmdot = np.array([vent_rate, end_rate])

    if neps is None:
        neps = 11
        teps = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, duration])
        xeps = np.array([0.0, 0.0513, 0.178, 0.28, 0.362,
                         0.422, 0.47, 0.52, 0.56, 0.6, 0.6])

    return {
        "Title": "LIQLEV GUI Simulation",
        "Liquid": fluid, "Units": "British", "Delta": delta_t,
        "Dtank": dtank, "Htzero": htzero, "Volt": volt,
        "Xmlzro": xmlzro, "Pinit": pinit_psia, "Pfinal": pfinal_psia,
        "Tinit": tinit, "Thetin": 0.0,
        "Nvmd": len(tvmdot), "Neps": neps,
        "Nlattm": 2, "Nvertm": 2, "Nggo": nggo,
        "Tvmdot": tvmdot, "Xvmdot": xvmdot,
        "Teps": teps, "Xeps": xeps,
        "Tspal": np.array([0.0, duration]), "Xspacl": np.array([1.0, 1.0]),
        "Tspav": np.array([0.0, duration]), "Xspacv": np.array([1.0, 1.0]),
        "Tggo": tggo, "Xggo": xggo,
        "gravity_function": gravity_function,
    }


def calculate_epsilon(h, dtank):
    perim = np.pi * dtank
    a_wall = perim * h
    ac = 0.7854 * (dtank ** 2)
    return a_wall / (a_wall + ac) if (a_wall + ac) != 0 else 0


# ── Safe math evaluator for user-defined gravity functions ───────────────────
import math as _math

_GRAV_SAFE_DICT = {
    "sin": _math.sin, "cos": _math.cos, "tan": _math.tan,
    "asin": _math.asin, "acos": _math.acos, "atan": _math.atan,
    "exp": _math.exp, "log": _math.log, "log10": _math.log10,
    "sqrt": _math.sqrt, "abs": abs,
    "pi": _math.pi, "e": _math.e,
    "min": min, "max": max,
    "pow": pow,
}


def _safe_eval_gravity(expr, t):
    """Evaluate a user gravity expression safely. Only math functions + t allowed."""
    safe = dict(_GRAV_SAFE_DICT)
    safe["t"] = t
    safe["__builtins__"] = {}  # block all builtins for safety
    return float(eval(expr, safe))


def _make_gravity_function(expr):
    """Create a gravity_function(t) callable from a user expression.

    The returned function takes time in seconds and returns gravity in ft/s^2
    (the expression is in g's, so we convert internally).
    """
    def gravity_func(t_sec):
        g_level = _safe_eval_gravity(expr, t_sec)
        return g_level * G_TO_FT_S2  # convert g's to ft/s^2 for the solver
    return gravity_func


# ═════════════════════════════════════════════════════════════════════════════
#  MAIN APPLICATION
# ═════════════════════════════════════════════════════════════════════════════
class LIQLEVApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("LIQLEV / Cryogenic Liquid Level Rise Simulator")
        self.geometry("1720x980")
        self.minsize(1280, 720)

        # Set window icon (taskbar + title bar)
        icon_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "liqlev.ico")
        if os.path.exists(icon_path):
            self.iconbitmap(icon_path)

        ctk.set_appearance_mode("dark")
        self.configure(fg_color=GH["bg"])

        # Data stores
        self.sweep_results = {}
        self.scenario_keys = []
        self.log_queue = queue.Queue()
        self.is_running = False
        self._mc_results = None
        self._crosshair_lines = []  # for crosshair cursor

        # Unit-aware field registry: list of (var, unit_label, unit_key)
        # Populated by _field() when unit_key is provided.
        # Used by _toggle_units() to convert values + relabel.
        self._unit_fields = []

        # Layout: header row 0, body row 1
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(1, weight=1)

        self._build_header()
        self._build_sidebar()
        self._build_results_panel()
        self._poll_log()

        # Keyboard shortcuts
        self.bind("<Control-Return>", lambda e: self._run_simulation())
        self.bind("<F5>", lambda e: self._run_simulation())
        self.bind("<Escape>", lambda e: self._tool_pointer())

    # ─────────────────────────────────────────────────────────────────────
    #  HEADER BAR
    # ─────────────────────────────────────────────────────────────────────
    def _build_header(self):
        hdr = ctk.CTkFrame(self, fg_color=GH["header"], height=52,
                           corner_radius=0)
        hdr.grid(row=0, column=0, columnspan=2, sticky="ew")
        hdr.grid_propagate(False)

        # Accent stripe (thin blue line at top — SpaceX webcast style)
        stripe = ctk.CTkFrame(hdr, fg_color=GH["blue"], height=2,
                              corner_radius=0)
        stripe.pack(side="top", fill="x")

        # App name — bold uppercase aerospace typography
        ctk.CTkLabel(
            hdr, text="[ LIQLEV ]",
            font=ctk.CTkFont(family="Consolas", size=22, weight="bold"),
            text_color=GH["blue"]
        ).pack(side="left", padx=(16, 0), pady=8)

        ctk.CTkLabel(
            hdr, text="/",
            font=ctk.CTkFont(family="Consolas", size=18),
            text_color=GH["muted"]
        ).pack(side="left", padx=(4, 4), pady=8)

        ctk.CTkLabel(
            hdr, text="CRYOGENIC LIQUID LEVEL RISE SIMULATOR",
            font=ctk.CTkFont(family="Consolas", size=11),
            text_color=GH["text2"]
        ).pack(side="left", padx=(0, 16), pady=8)

        # Theme toggle button
        self._dark_mode = True
        self._theme_btn = ctk.CTkButton(
            hdr, text="\u263e DARK", width=80, height=24,
            corner_radius=12,
            font=ctk.CTkFont(family="Consolas", size=10, weight="bold"),
            fg_color=GH["btn_secondary"],
            hover_color=GH["btn_secondary_hov"],
            text_color=GH["text2"],
            border_color=GH["border"], border_width=1,
            command=self._toggle_theme)
        self._theme_btn.pack(side="right", padx=(0, 12), pady=8)

        # Unit toggle button — pill style
        self._si_mode = False
        self._unit_btn = ctk.CTkButton(
            hdr, text="IMPERIAL", width=100, height=24,
            corner_radius=12,
            font=ctk.CTkFont(family="Consolas", size=10, weight="bold"),
            fg_color=GH["btn_secondary"],
            hover_color=GH["btn_secondary_hov"],
            text_color=GH["text2"],
            border_color=GH["border"], border_width=1,
            command=self._toggle_units)
        self._unit_btn.pack(side="right", padx=(0, 4), pady=8)

        ctk.CTkLabel(
            hdr, text="UNITS",
            font=ctk.CTkFont(family="Consolas", size=9),
            text_color=GH["muted"]
        ).pack(side="right", padx=(0, 4), pady=8)

        # Font scale buttons
        ctk.CTkButton(
            hdr, text="A+", width=32, height=24, corner_radius=3,
            font=ctk.CTkFont(family="Consolas", size=10, weight="bold"),
            fg_color=GH["btn_secondary"],
            hover_color=GH["btn_secondary_hov"],
            text_color=GH["text2"],
            border_color=GH["border"], border_width=1,
            command=lambda: self._scale_fonts(1)
        ).pack(side="right", padx=(0, 2), pady=8)
        ctk.CTkButton(
            hdr, text="A-", width=32, height=24, corner_radius=3,
            font=ctk.CTkFont(family="Consolas", size=10, weight="bold"),
            fg_color=GH["btn_secondary"],
            hover_color=GH["btn_secondary_hov"],
            text_color=GH["text2"],
            border_color=GH["border"], border_width=1,
            command=lambda: self._scale_fonts(-1)
        ).pack(side="right", padx=(0, 4), pady=8)

        # Status pill
        self.status_label = ctk.CTkLabel(
            hdr, text="  READY  ",
            font=ctk.CTkFont(family="Consolas", size=10, weight="bold"),
            fg_color=GH["green"], corner_radius=12,
            text_color="#000000", height=24)
        self.status_label.pack(side="right", padx=(0, 12), pady=8)

    # ─────────────────────────────────────────────────────────────────────
    #  SIDEBAR (Configuration)
    # ─────────────────────────────────────────────────────────────────────
    def _build_sidebar(self):
        sb = ctk.CTkScrollableFrame(
            self, width=360, fg_color=GH["surface"],
            border_color=GH["border"], border_width=1,
            corner_radius=8,
            scrollbar_button_color=GH["border"],
            scrollbar_button_hover_color=GH["muted"])
        sb.grid(row=1, column=0, sticky="nsew", padx=(8, 4), pady=(4, 8))

        # ── Fluid & Thermodynamic ──
        self._section(sb, "Fluid & Thermodynamic")

        self.fluid_var = ctk.StringVar(value=_DEFAULTS["fluid"])
        self._dropdown(sb, "Fluid", self.fluid_var, FLUIDS)
        self.pinit_var  = self._field(sb, "Initial Pressure", _DEFAULTS["pinit"], "psia",
                                      unit_key="pressure")
        self.pfinal_var = self._field(sb, "Final Pressure",   _DEFAULTS["pfinal"],  "psia",
                                      unit_key="pressure")

        # ── Tank Geometry ──
        self._section(sb, "Tank Geometry")

        # Preset dropdown
        self.tank_preset_var = ctk.StringVar(value=_DEFAULTS.get("tank_preset", "Custom"))
        self._dropdown(sb, "Preset", self.tank_preset_var,
                       list(TANK_PRESETS.keys()),
                       command=self._on_tank_preset)

        self.dtank_var = self._field(sb, "Diameter", _DEFAULTS["dtank"], "ft",
                                     unit_key="length_ft")
        self.htank_var = self._field(sb, "Height",   _DEFAULTS["htank"], "ft",
                                     unit_key="length_ft")

        self.vol_label = ctk.CTkLabel(
            sb, text="VOL: -- ft\u00b3   |   Ac: -- ft\u00b2",
            font=ctk.CTkFont(family="Consolas", size=10),
            text_color=GH["cyan"])
        self.vol_label.pack(anchor="w", padx=8, pady=(0, 4))

        for v in (self.dtank_var, self.htank_var):
            v.trace_add("write", lambda *_: self._update_computed())
        self._update_computed()

        # ── Measured Initial Conditions (inside Tank Geometry) ──
        # Separator label
        ctk.CTkLabel(
            sb, text="\u2500\u2500  Measured Initial Conditions  \u2500\u2500",
            font=ctk.CTkFont(family="Consolas", size=11, weight="bold"),
            text_color=GH["orange"]).pack(anchor="w", padx=8, pady=(6, 2))

        self.xmlzro_override_var = self._field(
            sb, "Initial Liquid Mass*", _DEFAULTS["xmlzro_override"], "lbm")
        self.tinit_override_var = self._field(
            sb, "Initial Temperature", _DEFAULTS["tinit_override"], "\u00b0R")

        # Info button — opens a popup dialog
        override_info_row = ctk.CTkFrame(sb, fg_color="transparent")
        override_info_row.pack(fill="x", padx=8, pady=(2, 2))
        self._override_info_btn = ctk.CTkButton(
            override_info_row, text="\u24d8  Why are these needed?", width=200, height=22,
            font=ctk.CTkFont(family="Consolas", size=10),
            fg_color=GH["btn_secondary"], hover_color=GH["btn_secondary_hov"],
            text_color=GH["cyan"], border_color=GH["border"], border_width=1,
            corner_radius=3, command=self._show_override_info)
        self._override_info_btn.pack(side="left")

        self._hint(sb, "*Leave blank \u2192 auto-computed from geometry and fill fraction")

        # ── Simulation Control ──
        self._section(sb, "Simulation Control")
        self.duration_var = self._field(sb, "Duration",       _DEFAULTS["duration"], "s")
        self.dt_var       = self._field(sb, "Time Step (dt)", _DEFAULTS["dt"], "s")

        # ── Sweep Parameters ──
        self._section(sb, "Sweep Parameters")

        self.vent_var = self._field(sb, "Vent Rate(s)", _DEFAULTS["vent"], "lbm/s",
                                    unit_key="mass_flow")
        self._hint(sb, "Comma, range (0.001:0.001:0.005), or linspace(0.001,0.005,5)")

        # Vent rate CSV profile import
        vent_csv_row = ctk.CTkFrame(sb, fg_color="transparent")
        vent_csv_row.pack(fill="x", padx=8, pady=(2, 0))
        self.vent_csv_var = ctk.StringVar(value="")
        ctk.CTkLabel(vent_csv_row, text="or CSV Profile:", width=100,
                     anchor="w", text_color=GH["muted"],
                     font=ctk.CTkFont(size=10)).pack(side="left")
        ctk.CTkEntry(vent_csv_row, textvariable=self.vent_csv_var, width=140,
                     fg_color=GH["input_bg"], border_color=GH["border"],
                     text_color=GH["text"],
                     placeholder_text="time,vent_rate CSV",
                     font=ctk.CTkFont(family="Consolas", size=10)
                     ).pack(side="left", padx=2)
        ctk.CTkButton(vent_csv_row, text="...", width=28, height=24,
                      fg_color=GH["btn_secondary"],
                      hover_color=GH["btn_secondary_hov"],
                      text_color=GH["text"],
                      command=self._browse_vent_csv).pack(side="left")

        self.fill_var = self._field(sb, "Fill Fraction(s)", _DEFAULTS["fill"], "")
        self._hint(sb, "Range 0-1.  e.g.  0.25, 0.50, 0.75")

        # Sweep preview
        self.sweep_label = ctk.CTkLabel(
            sb, text="",
            font=ctk.CTkFont(family="Consolas", size=10, weight="bold"),
            text_color=GH["orange"])
        self.sweep_label.pack(anchor="w", padx=8, pady=(4, 2))

        for v in (self.vent_var, self.fill_var):
            v.trace_add("write", lambda *_: self._update_sweep_count())

        # ── Epsilon ──
        self._section(sb, "Epsilon (Boil-off Partitioning)")
        self.eps_mode_var = ctk.StringVar(value=_DEFAULTS["eps_mode"])
        self._dropdown(sb, "Mode", self.eps_mode_var, EPSILON_MODES,
                       command=self._on_eps_mode_change)

        # Description label that updates per mode — this is the anchor widget
        self.eps_desc_label = ctk.CTkLabel(
            sb, text="",
            font=ctk.CTkFont(size=10), text_color=GH["muted"],
            wraplength=320, justify="left")
        self.eps_desc_label.pack(anchor="w", padx=12, pady=(0, 4))

        # Custom value entry — packed/forgotten AFTER eps_desc_label
        self.eps_custom_frame = ctk.CTkFrame(sb, fg_color="transparent")
        self.eps_custom_var = self._field(self.eps_custom_frame,
                                          "Custom Value(s)", _DEFAULTS["eps_custom"], "")
        self._hint(self.eps_custom_frame,
                   "Comma-separated, range (0.2:0.2:0.8), or linspace(0.2,0.8,4)")
        # Store the anchor so we can re-pack after it
        self._eps_anchor = self.eps_desc_label

        # Update sweep count when custom epsilon values change
        self.eps_custom_var.trace_add("write", lambda *_: self._update_sweep_count())

        # Initialize visibility
        self._on_eps_mode_change(self.eps_mode_var.get())

        # ── Vent Valve Ramping ──
        self._section(sb, "Vent Valve Ramping")
        self.ramp_dur_var    = self._field(sb, "Ramp Duration",  _DEFAULTS["ramp_dur"], "s")
        self.ramp_factor_var = self._field(sb, "Target Factor",  _DEFAULTS["ramp_factor"], "")
        self._hint(sb, "1.0 = constant rate, 0.8 = 20% linear reduction")

        # ── Gravity ──
        self._section(sb, "Gravity Settings")

        GRAVITY_MODES = ["Constant", "Function of Time", "CSV Profile"]
        self.grav_mode_var = ctk.StringVar(value="Constant")
        self._dropdown(sb, "Mode", self.grav_mode_var, GRAVITY_MODES,
                       command=self._on_gravity_mode_change)

        # Description label
        self.grav_desc_label = ctk.CTkLabel(
            sb, text="",
            font=ctk.CTkFont(size=10), text_color=GH["muted"],
            wraplength=320, justify="left")
        self.grav_desc_label.pack(anchor="w", padx=12, pady=(0, 4))
        self._grav_anchor = self.grav_desc_label

        # ── Mode 1: Constant ──
        self.grav_const_frame = ctk.CTkFrame(sb, fg_color="transparent")
        self.grav_const_val_var = self._field(self.grav_const_frame,
                                              "Gravity Level", _DEFAULTS["const_g"], "g's")

        # ── Mode 2: Function of Time ──
        self.grav_func_frame = ctk.CTkFrame(sb, fg_color="transparent")
        self.grav_func_var = ctk.StringVar(value=_DEFAULTS["grav_func"])
        ff = ctk.CTkFrame(self.grav_func_frame, fg_color="transparent")
        ff.pack(fill="x", padx=8, pady=2)
        ctk.CTkLabel(ff, text="g(t) =", width=50, anchor="w",
                     text_color=GH["blue"],
                     font=ctk.CTkFont(family="Consolas", size=13,
                                      weight="bold")).pack(side="left")
        ctk.CTkEntry(ff, textvariable=self.grav_func_var, width=220,
                     fg_color=GH["input_bg"], border_color=GH["border"],
                     text_color=GH["text"],
                     font=ctk.CTkFont(family="Consolas", size=12)
                     ).pack(side="left", padx=4)
        ctk.CTkLabel(self.grav_func_frame, text="  g's",
                     font=ctk.CTkFont(size=11), text_color=GH["muted"]
                     ).pack(anchor="w", padx=60)
        # Collapsible hints behind ? button
        hint_row = ctk.CTkFrame(self.grav_func_frame, fg_color="transparent")
        hint_row.pack(fill="x", padx=8, pady=2)
        self._grav_hints_visible = False
        self._grav_hint_btn = ctk.CTkButton(
            hint_row, text="?  Show Examples", width=130, height=22,
            font=ctk.CTkFont(family="Consolas", size=10),
            fg_color=GH["btn_secondary"], hover_color=GH["btn_secondary_hov"],
            text_color=GH["text2"], border_color=GH["border"], border_width=1,
            corner_radius=3, command=self._toggle_grav_hints)
        self._grav_hint_btn.pack(side="left")

        self._grav_hints_frame = ctk.CTkFrame(self.grav_func_frame,
                                               fg_color=GH["overlay"],
                                               corner_radius=4)
        hints = [
            "Variable: t (time in seconds)",
            "Math: sin, cos, exp, log, sqrt, abs, pi, e",
            "",
            "Examples:",
            "  Constant:     0.001",
            "  Step at 5s:   0.01 if t < 5 else 0.001",
            "  Linear ramp:  0.01 - 0.009 * t / 20",
            "  Oscillating:  0.001 + 0.0005*sin(2*pi*t/10)",
            "  Drop tower:   0.01*exp(-t/2) + 1e-4",
        ]
        for h in hints:
            ctk.CTkLabel(self._grav_hints_frame, text="  " + h,
                         font=ctk.CTkFont(family="Consolas", size=10),
                         text_color=GH["muted"]).pack(anchor="w", padx=4)

        # ── Mode 3: CSV Profile ──
        self.grav_csv_frame = ctk.CTkFrame(sb, fg_color="transparent")
        cf = ctk.CTkFrame(self.grav_csv_frame, fg_color="transparent")
        cf.pack(fill="x", padx=8, pady=2)
        ctk.CTkLabel(cf, text="CSV File", width=120, anchor="w",
                     text_color=GH["text2"],
                     font=ctk.CTkFont(size=12)).pack(side="left")
        self.grav_csv_var = ctk.StringVar(value=_DEFAULTS["grav_csv"])
        ctk.CTkEntry(cf, textvariable=self.grav_csv_var, width=130,
                     fg_color=GH["input_bg"], border_color=GH["border"],
                     text_color=GH["text"]).pack(side="left", padx=4)
        ctk.CTkButton(cf, text="...", width=32, height=28,
                      fg_color=GH["btn_secondary"],
                      hover_color=GH["btn_secondary_hov"],
                      text_color=GH["text"],
                      command=self._browse_gravity).pack(side="left")
        self.hold_g_var = self._field(self.grav_csv_frame, "Hold G Value",
                                      _DEFAULTS["hold_g"], "g's")
        self._hint(self.grav_csv_frame,
                   "CSV must have columns: normalized_time, az_positive (in g's)")
        self._hint(self.grav_csv_frame,
                   "Hold G: if sim runs longer than CSV data, hold this value")

        # Initialize
        self._on_gravity_mode_change(self.grav_mode_var.get())

        # ── Monte Carlo / Sensitivity ──
        self._section(sb, "Monte Carlo / Sensitivity")
        self._hint(sb, "Randomize parameters within ranges to find worst-case")

        self.mc_n_var = self._field(sb, "N Samples", _DEFAULTS["mc_n"], "")
        self.mc_vent_min_var = self._field(sb, "Vent Rate Min", _DEFAULTS["mc_vent_min"], "lbm/s",
                                          unit_key="mass_flow")
        self.mc_vent_max_var = self._field(sb, "Vent Rate Max", _DEFAULTS["mc_vent_max"], "lbm/s",
                                           unit_key="mass_flow")
        self.mc_fill_min_var = self._field(sb, "Fill Frac Min",  _DEFAULTS["mc_fill_min"], "")
        self.mc_fill_max_var = self._field(sb, "Fill Frac Max",  _DEFAULTS["mc_fill_max"], "")
        self.mc_grav_min_var = self._field(sb, "Gravity Min", _DEFAULTS["mc_grav_min"], "g's")
        self.mc_grav_max_var = self._field(sb, "Gravity Max", _DEFAULTS["mc_grav_max"],  "g's")

        self.mc_btn = ctk.CTkButton(
            sb, text="RUN MONTE CARLO", height=36,
            font=ctk.CTkFont(family="Consolas", size=12, weight="bold"),
            fg_color=GH["purple"], hover_color="#7c3aed",
            text_color="#ffffff", corner_radius=4,
            command=self._run_monte_carlo)
        self.mc_btn.pack(fill="x", padx=8, pady=(6, 4))

        # ── Action Buttons ──
        ctk.CTkFrame(sb, fg_color="transparent", height=12).pack()

        self.run_btn = ctk.CTkButton(
            sb, text="\u25B6  RUN SIMULATION", height=42,
            font=ctk.CTkFont(family="Consolas", size=13, weight="bold"),
            fg_color=GH["green"], hover_color=GH["green_hover"],
            text_color="#000000", corner_radius=4,
            command=self._run_simulation)
        self.run_btn.pack(fill="x", padx=8, pady=(4, 4))

        prog_row = ctk.CTkFrame(sb, fg_color="transparent")
        prog_row.pack(fill="x", padx=8, pady=(0, 6))
        self.progress = ctk.CTkProgressBar(
            prog_row, mode="determinate", height=8,
            progress_color=GH["blue"], fg_color=GH["border"])
        self.progress.pack(side="left", fill="x", expand=True)
        self.progress.set(0)
        self._progress_label = ctk.CTkLabel(
            prog_row, text="0%", width=40,
            font=ctk.CTkFont(family="Consolas", size=10, weight="bold"),
            text_color=GH["text2"])
        self._progress_label.pack(side="left", padx=(4, 0))

        btn_row = ctk.CTkFrame(sb, fg_color="transparent")
        btn_row.pack(fill="x", padx=8, pady=2)
        for txt, cmd, col in [
            ("EXPORT CSV",    self._export_results,  GH["btn_secondary"]),
            ("EXPORT SUMMARY", self._export_summary_csv, GH["btn_secondary"]),
            ("CLEAR",         self._clear_results,   GH["btn_secondary"]),
        ]:
            ctk.CTkButton(
                btn_row, text=txt, height=30, fg_color=col,
                hover_color=GH["btn_secondary_hov"],
                text_color=GH["text2"], corner_radius=4,
                font=ctk.CTkFont(family="Consolas", size=10),
                border_color=GH["border"], border_width=1,
                command=cmd
            ).pack(side="left", expand=True, fill="x", padx=2)

        btn_row2 = ctk.CTkFrame(sb, fg_color="transparent")
        btn_row2.pack(fill="x", padx=8, pady=2)
        for txt, cmd in [("SAVE CONFIG", self._save_config),
                         ("LOAD CONFIG", self._load_config)]:
            ctk.CTkButton(
                btn_row2, text=txt, height=28,
                fg_color=GH["btn_secondary"],
                hover_color=GH["btn_secondary_hov"],
                text_color=GH["text2"], corner_radius=4,
                font=ctk.CTkFont(family="Consolas", size=10),
                border_color=GH["border"], border_width=1,
                command=cmd
            ).pack(side="left", expand=True, fill="x", padx=2)

        btn_row3 = ctk.CTkFrame(sb, fg_color="transparent")
        btn_row3.pack(fill="x", padx=8, pady=2)
        ctk.CTkButton(
            btn_row3, text="GENERATE PDF REPORT", height=30,
            fg_color=GH["blue"], hover_color="#4493f8",
            text_color="#ffffff", corner_radius=4,
            font=ctk.CTkFont(family="Consolas", size=10, weight="bold"),
            command=self._generate_report
        ).pack(fill="x", padx=2)

        btn_row4 = ctk.CTkFrame(sb, fg_color="transparent")
        btn_row4.pack(fill="x", padx=8, pady=2)
        ctk.CTkButton(
            btn_row4, text="\u21BA  RESET DEFAULTS", height=28,
            fg_color=GH["btn_secondary"],
            hover_color=GH["btn_secondary_hov"],
            text_color=GH["text2"], corner_radius=4,
            font=ctk.CTkFont(family="Consolas", size=10),
            border_color=GH["border"], border_width=1,
            command=self._reset_defaults
        ).pack(fill="x", padx=2)

        self._update_sweep_count()

    # ── Sidebar widget helpers ──
    def _section(self, parent, text):
        f = ctk.CTkFrame(parent, fg_color="transparent", height=28)
        f.pack(fill="x", padx=4, pady=(14, 2))
        # Accent bar
        bar = ctk.CTkFrame(f, fg_color=GH["blue"], width=3, height=16,
                           corner_radius=2)
        bar.pack(side="left", padx=(4, 8))
        ctk.CTkLabel(f, text=f"[ {text.upper()} ]",
                     font=ctk.CTkFont(family="Consolas", size=11,
                                      weight="bold"),
                     text_color=GH["text2"]).pack(side="left")

    def _field(self, parent, label, default, unit, numeric=True,
               unit_key=None):
        """Create a labeled input field.

        Parameters
        ----------
        unit_key : str or None
            Key into UNIT_CONV (e.g. "pressure", "length_ft", "mass_flow").
            When provided the field is registered for automatic conversion
            when the user toggles Imperial / SI.  A small secondary label
            shows the live-converted value in the alternate unit system.
        """
        f = ctk.CTkFrame(parent, fg_color="transparent")
        f.pack(fill="x", padx=8, pady=2)
        ctk.CTkLabel(f, text=label, width=120, anchor="w",
                     text_color=GH["muted"],
                     font=ctk.CTkFont(family="Consolas", size=11)).pack(
                         side="left")
        var = ctk.StringVar(value=default)
        entry = ctk.CTkEntry(f, textvariable=var, width=120,
                     fg_color=GH["input_bg"],
                     border_color=GH["border"],
                     text_color=GH["cyan"],
                     font=ctk.CTkFont(family="Consolas", size=12))
        entry.pack(side="left", padx=4)
        unit_lbl = None
        alt_lbl = None  # secondary label for live-converted value
        if unit:
            unit_lbl = ctk.CTkLabel(f, text=unit, width=34, anchor="w",
                         text_color=GH["muted"],
                         font=ctk.CTkFont(family="Consolas", size=10))
            unit_lbl.pack(side="left")
        # Register for unit conversion if unit_key provided
        if unit_key and unit_lbl:
            alt_lbl = ctk.CTkLabel(f, text="", width=90, anchor="w",
                                   text_color=GH["muted"],
                                   font=ctk.CTkFont(family="Consolas",
                                                     size=9))
            alt_lbl.pack(side="left", padx=(0, 0))
            self._unit_fields.append((var, unit_lbl, unit_key, alt_lbl))
            # Live trace — update alt label whenever value changes
            var.trace_add("write",
                          lambda *_, v=var, uk=unit_key, al=alt_lbl:
                          self._update_alt_unit(v, uk, al))
            # Initialize the alt label
            self._update_alt_unit(var, unit_key, alt_lbl)
        # Live validation for numeric fields
        if numeric:
            var.trace_add("write", lambda *_, e=entry, v=var:
                          self._validate_field(e, v))
        return var

    def _update_alt_unit(self, var, unit_key, alt_lbl):
        """Update the secondary label with the value converted to the
        alternate unit system."""
        brit_u, si_u, factor = UNIT_CONV[unit_key]
        raw = var.get().strip()
        if not raw:
            alt_lbl.configure(text="")
            return
        try:
            val = float(raw)
            if self._si_mode:
                # Currently SI → show British equivalent
                converted = val / factor
                alt_lbl.configure(text=f"({converted:.4g} {brit_u})")
            else:
                # Currently British → show SI equivalent
                converted = val * factor
                alt_lbl.configure(text=f"({converted:.4g} {si_u})")
        except ValueError:
            # Array input — try first value
            try:
                vals = parse_array(raw)
                if vals:
                    v0 = vals[0]
                    if self._si_mode:
                        c = v0 / factor
                        alt_lbl.configure(
                            text=f"({c:.4g} {brit_u}"
                                 f"{', ...' if len(vals) > 1 else ''})")
                    else:
                        c = v0 * factor
                        alt_lbl.configure(
                            text=f"({c:.4g} {si_u}"
                                 f"{', ...' if len(vals) > 1 else ''})")
                else:
                    alt_lbl.configure(text="")
            except Exception:
                alt_lbl.configure(text="")

    def _validate_field(self, entry, var):
        """Color entry border red if value is not a valid number/array."""
        val = var.get().strip()
        if not val:
            entry.configure(border_color=GH["border"])
            return
        try:
            # Try as single float first
            float(val)
            entry.configure(border_color=GH["border"])
        except ValueError:
            # Try as array
            try:
                result = parse_array(val)
                if result:
                    entry.configure(border_color=GH["border"])
                else:
                    entry.configure(border_color=GH["red"])
            except Exception:
                entry.configure(border_color=GH["red"])

    def _dropdown(self, parent, label, var, values, command=None):
        f = ctk.CTkFrame(parent, fg_color="transparent")
        f.pack(fill="x", padx=8, pady=2)
        ctk.CTkLabel(f, text=label, width=120, anchor="w",
                     text_color=GH["muted"],
                     font=ctk.CTkFont(family="Consolas", size=11)).pack(
                         side="left")
        ctk.CTkOptionMenu(
            f, variable=var, values=values, width=160,
            fg_color=GH["btn_secondary"], button_color=GH["border"],
            button_hover_color=GH["muted"],
            dropdown_fg_color=GH["overlay"],
            dropdown_hover_color=GH["border"],
            text_color=GH["text"],
            font=ctk.CTkFont(family="Consolas", size=11),
            command=command
        ).pack(side="left", padx=4)

    def _hint(self, parent, text):
        ctk.CTkLabel(parent, text="  " + text,
                     font=ctk.CTkFont(size=10),
                     text_color=GH["muted"]).pack(anchor="w", padx=8)

    def _scale_fonts(self, direction):
        """Scale all widget fonts up (+1) or down (-1)."""
        global BASE_FONT
        BASE_FONT = max(8, min(18, BASE_FONT + direction))
        # Update plot fonts on next render
        if self.sweep_results:
            self._render_current_plot()

    def _on_tank_preset(self, choice):
        """Apply a tank geometry preset."""
        preset = TANK_PRESETS.get(choice, {})
        if not preset:
            return
        if "dtank" in preset:
            self.dtank_var.set(preset["dtank"])
        if "htank" in preset:
            self.htank_var.set(preset["htank"])
        if "fluid" in preset:
            self.fluid_var.set(preset["fluid"])
        self.xmlzro_override_var.set(preset.get("xmlzro_override", ""))
        self.tinit_override_var.set(preset.get("tinit_override", ""))

    def _browse_vent_csv(self):
        path = filedialog.askopenfilename(
            title="Select Vent Rate CSV",
            filetypes=[("CSV", "*.csv"), ("All", "*.*")])
        if path:
            self.vent_csv_var.set(path)

    def _toggle_grav_hints(self):
        """Show/hide the gravity function examples."""
        self._grav_hints_visible = not self._grav_hints_visible
        if self._grav_hints_visible:
            self._grav_hints_frame.pack(fill="x", padx=8, pady=(2, 4))
            self._grav_hint_btn.configure(text="\u25B2  Hide Examples")
        else:
            self._grav_hints_frame.pack_forget()
            self._grav_hint_btn.configure(text="?  Show Examples")

    def _show_override_info(self):
        """Pop up a dialog explaining the validation override fields."""
        sep = "\u2500" * 48
        msg = (
            "HOW THE SOLVER COMPUTES INITIAL CONDITIONS\n"
            f"{sep}\n"
            "The solver models the tank as a simple cylinder:\n\n"
            "    Vol  = (\u03c0/4) \u00d7 D\u00b2 \u00d7 H\n"
            "    T\u2080    = Tsat(fluid, P\u2080)   via CoolProp\n"
            "    m\u2080    = \u03c1_liq(T\u2080) \u00d7 fill fraction \u00d7 Vol\n\n\n"
            "WHY OVERRIDES MAY BE NEEDED\n"
            f"{sep}\n"
            "Real tanks have curved domes (ellipsoidal, hemispherical,\n"
            "or torispherical heads) that reduce the usable volume\n"
            "compared to a cylinder of the same diameter and height.\n\n"
            "For flight validation (e.g. Saturn S-IVB AS-203), the\n"
            "actual propellant mass and temperature were measured\n"
            "\u2014 not derived from geometry.  These measured values\n"
            "differ from what the cylindrical formula computes.\n\n\n"
            "HOW TO USE\n"
            f"{sep}\n"
            "\u25b6  Leave blank  \u2192  solver computes from geometry\n"
            "\u25b6  Enter value  \u2192  overrides solver computation\n\n\n"
            "AS-203 FLIGHT DATA\n"
            f"{sep}\n"
            "  Liquid mass:      16,300 lbm\n"
            "  Temperature:      38.3 \u00b0R   (21.3 K)\n"
            "  Tank diameter:    21.67 ft  (6.604 m)\n"
            "  Tank volume:      10,392 ft\u00b3\n"
            "  Fill level:       ~51%"
        )
        messagebox.showinfo("Measured Initial Conditions", msg)

    def _reset_defaults(self):
        """Reset all input fields to default values (British).

        If currently in SI mode, switch to Imperial first, set defaults,
        then switch back so the user sees correct SI values.
        """
        was_si = self._si_mode
        # Temporarily force Imperial so defaults (which are British) set cleanly
        if was_si:
            self._si_mode = True  # _toggle_units will flip to False
            self._toggle_units()

        self.fluid_var.set(_DEFAULTS["fluid"])
        self.pinit_var.set(_DEFAULTS["pinit"])
        self.pfinal_var.set(_DEFAULTS["pfinal"])
        self.dtank_var.set(_DEFAULTS["dtank"])
        self.htank_var.set(_DEFAULTS["htank"])
        self.duration_var.set(_DEFAULTS["duration"])
        self.dt_var.set(_DEFAULTS["dt"])
        self.vent_var.set(_DEFAULTS["vent"])
        self.fill_var.set(_DEFAULTS["fill"])
        self.eps_mode_var.set(_DEFAULTS["eps_mode"])
        self.eps_custom_var.set(_DEFAULTS["eps_custom"])
        self.ramp_dur_var.set(_DEFAULTS["ramp_dur"])
        self.ramp_factor_var.set(_DEFAULTS["ramp_factor"])
        self.grav_mode_var.set(_DEFAULTS["grav_mode"])
        self.grav_const_val_var.set(_DEFAULTS["const_g"])
        self.grav_func_var.set(_DEFAULTS["grav_func"])
        self.grav_csv_var.set(_DEFAULTS["grav_csv"])
        self.hold_g_var.set(_DEFAULTS["hold_g"])
        self.mc_n_var.set(_DEFAULTS["mc_n"])
        self.mc_vent_min_var.set(_DEFAULTS["mc_vent_min"])
        self.mc_vent_max_var.set(_DEFAULTS["mc_vent_max"])
        self.mc_fill_min_var.set(_DEFAULTS["mc_fill_min"])
        self.mc_fill_max_var.set(_DEFAULTS["mc_fill_max"])
        self.mc_grav_min_var.set(_DEFAULTS["mc_grav_min"])
        self.mc_grav_max_var.set(_DEFAULTS["mc_grav_max"])
        self._threshold_var.set(_DEFAULTS["threshold"])
        self.vent_csv_var.set("")
        self.tank_preset_var.set(_DEFAULTS.get("tank_preset", "Custom"))
        self.xmlzro_override_var.set(_DEFAULTS["xmlzro_override"])
        self.tinit_override_var.set(_DEFAULTS["tinit_override"])
        self._on_gravity_mode_change(_DEFAULTS["grav_mode"])
        self._on_eps_mode_change(_DEFAULTS["eps_mode"])
        self._update_computed()
        self._update_sweep_count()

        # Restore SI mode if that's where the user was
        if was_si:
            self._toggle_units()  # flips False → True, converts values

    def _update_computed(self):
        try:
            d = float(self.dtank_var.get())
            h = float(self.htank_var.get())
            # d and h are in whatever units the user currently sees
            # compute vol/ac in those units
            vol = (np.pi / 4) * d ** 2 * h
            ac = 0.7854 * d ** 2
            if self._si_mode:
                vol_u, area_u = "m\u00b3", "m\u00b2"
            else:
                vol_u, area_u = "ft\u00b3", "ft\u00b2"
            self.vol_label.configure(
                text=f"VOL: {vol:.5f} {vol_u}  |  Ac: {ac:.5f} {area_u}")
        except ValueError:
            self.vol_label.configure(
                text="VOL: --  |  Ac: --")

    def _update_sweep_count(self):
        try:
            nv = len(parse_array(self.vent_var.get()))
            nf = len(parse_array(self.fill_var.get()))
        except Exception:
            nv = nf = 0
        ne = 1  # epsilon mode counted as 1 unless custom has multiple
        if self.eps_mode_var.get() == "Custom":
            try:
                ne = len(parse_array(self.eps_custom_var.get()))
            except Exception:
                ne = 1
        total = nv * nf * ne
        if total > 0:
            self.sweep_label.configure(
                text=f"{total} run{'s' if total != 1 else ''}  "
                     f"({nf} fill x {nv} vent x {ne} eps)")
        else:
            self.sweep_label.configure(text="")

    def _on_eps_mode_change(self, choice):
        EPS_DESCRIPTIONS = {
            "height_dep": "Dynamic: eps = wetted wall area / (wall + free surface). "
                          "Recalculated every timestep based on current liquid height.",
            "bulk_fake":  "Bulk boiling: eps = 50. Mimics distributed volumetric "
                          "boiling in large tanks (e.g. JAXA 30m\u00b3 experiment).",
            "AS-203 Schedule": "Time-varying schedule from the Saturn S-IVB AS-203 "
                          "validation study. Ramps 0\u21920.6 over 180s (11 points).",
            "Custom":     "Fixed constant value(s). Enter one or more epsilon values "
                          "to sweep. Each value runs as a separate scenario.",
        }
        self.eps_desc_label.configure(text=EPS_DESCRIPTIONS.get(choice, ""))

        # Show/hide custom entry — re-pack after the description label anchor
        self.eps_custom_frame.pack_forget()
        if choice == "Custom":
            self.eps_custom_frame.pack(fill="x", pady=2,
                                       after=self._eps_anchor)

        self._update_sweep_count()

    def _on_gravity_mode_change(self, choice):
        GRAV_DESCRIPTIONS = {
            "Constant":         "Fixed gravity level for the entire simulation duration.",
            "Function of Time": "Define gravity as a Python math expression of time (t). "
                                "The expression is evaluated at every solver timestep. "
                                "Output must be in g's (1 g = 32.174 ft/s\u00b2).",
            "CSV Profile":      "Import a time-vs-gravity profile from a CSV file. "
                                "Typically used for drop tower or parabolic flight data. "
                                "The solver linearly interpolates between data points.",
        }
        self.grav_desc_label.configure(text=GRAV_DESCRIPTIONS.get(choice, ""))

        # Hide all mode frames
        for frame in (self.grav_const_frame, self.grav_func_frame,
                      self.grav_csv_frame):
            frame.pack_forget()

        # Show the selected one
        if choice == "Constant":
            self.grav_const_frame.pack(fill="x", pady=2,
                                       after=self._grav_anchor)
        elif choice == "Function of Time":
            self.grav_func_frame.pack(fill="x", pady=2,
                                      after=self._grav_anchor)
        elif choice == "CSV Profile":
            self.grav_csv_frame.pack(fill="x", pady=2,
                                     after=self._grav_anchor)

    def _browse_gravity(self):
        path = filedialog.askopenfilename(
            title="Select Gravity CSV",
            filetypes=[("CSV", "*.csv"), ("All", "*.*")])
        if path:
            self.grav_csv_var.set(path)

    # ─────────────────────────────────────────────────────────────────────
    #  RESULTS PANEL
    # ─────────────────────────────────────────────────────────────────────
    def _build_results_panel(self):
        rp = ctk.CTkFrame(self, fg_color=GH["surface"],
                          border_color=GH["border"], border_width=1,
                          corner_radius=8)
        rp.grid(row=1, column=1, sticky="nsew", padx=(4, 8), pady=(4, 8))
        rp.grid_columnconfigure(0, weight=1)
        rp.grid_rowconfigure(1, weight=1)

        # ── Top control bar ──
        top = ctk.CTkFrame(rp, fg_color=GH["overlay"], height=44,
                           corner_radius=0)
        top.grid(row=0, column=0, sticky="ew", padx=0, pady=0)

        ctk.CTkLabel(top, text="SCENARIO",
                     font=ctk.CTkFont(family="Consolas", size=10,
                                      weight="bold"),
                     text_color=GH["muted"]).pack(side="left", padx=(12, 4))
        self.scenario_var = ctk.StringVar(value="(run a simulation)")
        self.scenario_menu = ctk.CTkOptionMenu(
            top, variable=self.scenario_var,
            values=["(run a simulation)"], width=340,
            fg_color=GH["btn_secondary"], button_color=GH["border"],
            button_hover_color=GH["muted"],
            dropdown_fg_color=GH["overlay"],
            dropdown_hover_color=GH["border"],
            text_color=GH["text"],
            font=ctk.CTkFont(family="Consolas", size=11),
            command=self._on_scenario_change)
        self.scenario_menu.pack(side="left", padx=(0, 20))

        ctk.CTkLabel(top, text="PLOT",
                     font=ctk.CTkFont(family="Consolas", size=10,
                                      weight="bold"),
                     text_color=GH["muted"]).pack(side="left", padx=(0, 4))
        self.plot_var = ctk.StringVar(value=PLOT_TYPES[0])
        ctk.CTkOptionMenu(
            top, variable=self.plot_var, values=PLOT_TYPES, width=260,
            fg_color=GH["btn_secondary"], button_color=GH["border"],
            button_hover_color=GH["muted"],
            dropdown_fg_color=GH["overlay"],
            dropdown_hover_color=GH["border"],
            text_color=GH["text"],
            font=ctk.CTkFont(family="Consolas", size=11),
            command=self._on_plot_change
        ).pack(side="left")

        # Overlay toggle
        self._overlay_var = ctk.BooleanVar(value=False)
        ctk.CTkCheckBox(
            top, text="OVERLAY",
            variable=self._overlay_var,
            font=ctk.CTkFont(family="Consolas", size=10),
            text_color=GH["text2"],
            fg_color=GH["blue"], hover_color=GH["blue"],
            border_color=GH["border"],
            command=self._on_plot_change
        ).pack(side="left", padx=(16, 4))

        # Threshold
        ctk.CTkLabel(top, text="THR dh/h\u2080:",
                     font=ctk.CTkFont(family="Consolas", size=10),
                     text_color=GH["muted"]).pack(side="left", padx=(12, 2))
        self._threshold_var = ctk.StringVar(value="")
        self._thresh_entry_ref = ctk.CTkEntry(
            top, textvariable=self._threshold_var, width=60, height=26,
            fg_color=GH["input_bg"], border_color=GH["border"],
            text_color=GH["orange"],
            placeholder_text="e.g. 0.3",
            font=ctk.CTkFont(family="Consolas", size=11))
        self._thresh_entry_ref.pack(side="left", padx=2)
        self._threshold_var.trace_add("write",
                                       lambda *_: self._on_plot_change())

        # ── Tabview ──
        self.tabs = ctk.CTkTabview(
            rp, fg_color=GH["surface"],
            segmented_button_fg_color=GH["overlay"],
            segmented_button_selected_color=GH["btn_secondary"],
            segmented_button_selected_hover_color=GH["border"],
            segmented_button_unselected_color=GH["overlay"],
            segmented_button_unselected_hover_color=GH["border"],
            text_color=GH["text"],
            corner_radius=6)
        self.tabs.grid(row=1, column=0, sticky="nsew", padx=4, pady=(0, 4))

        # -- Plot tab --
        self.tab_plot = self.tabs.add("Plots")
        self.tab_plot.grid_columnconfigure(0, weight=1)
        self.tab_plot.grid_rowconfigure(0, weight=1)

        self.fig = Figure(figsize=(10, 6), dpi=100, facecolor=PLOT_BG)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.tab_plot)
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")

        # ── Custom themed toolbar ──
        toolbar_frame = ctk.CTkFrame(self.tab_plot, fg_color=GH["overlay"],
                                     height=36, corner_radius=0)
        toolbar_frame.grid(row=1, column=0, sticky="ew")

        # Hidden NavigationToolbar (we use it for its backend methods)
        self._mpl_toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        self._mpl_toolbar.pack_forget()  # hide the ugly default toolbar

        self._plot_mode = ctk.StringVar(value="pointer")  # pointer|pan|zoom

        btn_style = dict(
            height=26, width=70, corner_radius=3,
            font=ctk.CTkFont(family="Consolas", size=10),
            text_color=GH["text2"],
            fg_color=GH["btn_secondary"],
            hover_color=GH["btn_secondary_hov"],
            border_color=GH["border"], border_width=1,
        )

        # Tool buttons
        self._btn_pointer = ctk.CTkButton(
            toolbar_frame, text="\u25C7 Pointer", command=self._tool_pointer,
            **btn_style)
        self._btn_pointer.pack(side="left", padx=(8, 2), pady=4)

        self._btn_pan = ctk.CTkButton(
            toolbar_frame, text="\u2725 Pan", command=self._tool_pan,
            **btn_style)
        self._btn_pan.pack(side="left", padx=2, pady=4)

        self._btn_zoom = ctk.CTkButton(
            toolbar_frame, text="\u2922 Zoom", command=self._tool_zoom,
            **btn_style)
        self._btn_zoom.pack(side="left", padx=2, pady=4)

        # Separator
        ctk.CTkLabel(toolbar_frame, text="|", text_color=GH["muted"],
                     font=ctk.CTkFont(size=14)).pack(side="left", padx=4)

        ctk.CTkButton(
            toolbar_frame, text="\u21BA Reset View",
            command=self._reset_view, width=90,
            **{k: v for k, v in btn_style.items() if k != "width"}
        ).pack(side="left", padx=2, pady=4)

        ctk.CTkButton(
            toolbar_frame, text="\U0001F4BE Save Plot",
            command=self._save_plot, width=90,
            **{k: v for k, v in btn_style.items() if k != "width"}
        ).pack(side="left", padx=2, pady=4)

        ctk.CTkButton(
            toolbar_frame, text="\U0001F4CB Copy",
            command=self._copy_plot_to_clipboard, width=70,
            **{k: v for k, v in btn_style.items() if k != "width"}
        ).pack(side="left", padx=2, pady=4)

        # Coordinate readout label (right side)
        self._coord_label = ctk.CTkLabel(
            toolbar_frame, text="",
            font=ctk.CTkFont(family="Consolas", size=10),
            text_color=GH["cyan"])
        self._coord_label.pack(side="right", padx=(0, 12))

        ctk.CTkLabel(toolbar_frame, text="XY:",
                     text_color=GH["muted"],
                     font=ctk.CTkFont(family="Consolas", size=9)).pack(
                         side="right", padx=(8, 2))

        # Connect matplotlib canvas events for interactivity
        self.canvas.mpl_connect("scroll_event", self._on_scroll_zoom)
        self.canvas.mpl_connect("motion_notify_event", self._on_mouse_move)

        # Store default axis limits for reset
        self._default_view_limits = {}
        self._update_toolbar_highlight()

        # -- Data table tab --
        self.tab_data = self.tabs.add("Data Table")
        self.tab_data.grid_columnconfigure(0, weight=1)
        self.tab_data.grid_rowconfigure(0, weight=1)

        table_frame = ctk.CTkFrame(self.tab_data, fg_color="transparent")
        table_frame.grid(row=0, column=0, sticky="nsew")
        table_frame.grid_columnconfigure(0, weight=1)
        table_frame.grid_rowconfigure(0, weight=1)

        style = ttk.Style()
        style.theme_use("clam")
        style.configure("GH.Treeview",
                        background=GH["bg"],
                        foreground=GH["text"],
                        fieldbackground=GH["bg"],
                        rowheight=22,
                        borderwidth=0,
                        font=("Consolas", 10))
        style.configure("GH.Treeview.Heading",
                        background=GH["overlay"],
                        foreground=GH["text2"],
                        borderwidth=0,
                        font=("Consolas", 9, "bold"))
        style.map("GH.Treeview",
                  background=[("selected", GH["blue"])])

        self.tree = ttk.Treeview(table_frame, style="GH.Treeview",
                                 show="headings")
        vsb = ttk.Scrollbar(table_frame, orient="vertical",
                            command=self.tree.yview)
        hsb = ttk.Scrollbar(table_frame, orient="horizontal",
                            command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")

        # Data table controls
        dt_ctrl = ctk.CTkFrame(self.tab_data, fg_color="transparent", height=36)
        dt_ctrl.grid(row=2, column=0, sticky="ew", pady=4)
        ctk.CTkLabel(dt_ctrl, text="VENT RATE:",
                     text_color=GH["muted"],
                     font=ctk.CTkFont(family="Consolas", size=10)).pack(
                         side="left", padx=(8, 4))
        self.table_rate_var = ctk.StringVar(value="--")
        self.table_rate_menu = ctk.CTkOptionMenu(
            dt_ctrl, variable=self.table_rate_var, values=["--"], width=160,
            fg_color=GH["btn_secondary"], button_color=GH["border"],
            text_color=GH["text"],
            command=lambda _: self._update_data_table())
        self.table_rate_menu.pack(side="left", padx=4)

        # -- Summary tab --
        self.tab_summary = self.tabs.add("Summary")
        self.summary_text = ctk.CTkTextbox(
            self.tab_summary,
            font=ctk.CTkFont(family="Consolas", size=12),
            fg_color=GH["bg"], text_color=GH["text"],
            state="disabled")
        self.summary_text.pack(fill="both", expand=True, padx=4, pady=4)

        # -- Log tab --
        self.tab_log = self.tabs.add("Log")
        self.log_text = ctk.CTkTextbox(
            self.tab_log,
            font=ctk.CTkFont(family="Consolas", size=11),
            fg_color=GH["bg"], text_color=GH["text2"],
            state="disabled")
        self.log_text.pack(fill="both", expand=True, padx=4, pady=4)

        self._show_welcome_plot()

    def _show_welcome_plot(self):
        self.fig.clear()
        self.fig.set_facecolor(PLOT_BG)
        ax = self.fig.add_subplot(111)
        ax.set_facecolor(PLOT_BG)

        # Thin accent line
        ax.axhline(y=0.50, xmin=0.3, xmax=0.7,
                    color=GH["blue"], linewidth=1.5, alpha=0.6)

        ax.text(0.5, 0.60, "LIQLEV", transform=ax.transAxes, ha="center",
                va="center", fontsize=52, fontweight="bold", color=GH["blue"],
                family="Consolas", alpha=0.95)
        ax.text(0.5, 0.44, "CRYOGENIC LIQUID LEVEL RISE SIMULATOR",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=12, color=GH["text2"], family="Consolas")
        ax.text(0.5, 0.33,
                "Configure parameters  \u2192  Run Simulation",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=11, color=GH["muted"], family="Consolas")
        ax.text(0.5, 0.22, "v2.0  |  Numba JIT  |  CoolProp",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=9, color=GH["muted"], family="Consolas", alpha=0.5)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        self.canvas.draw_idle()

    # ─────────────────────────────────────────────────────────────────────
    #  LOGGING (thread-safe)
    # ─────────────────────────────────────────────────────────────────────
    def _log(self, msg):
        self.log_queue.put(msg)

    def _poll_log(self):
        batch = []
        try:
            while True:
                batch.append(self.log_queue.get_nowait())
        except queue.Empty:
            pass
        if batch:
            self.log_text.configure(state="normal")
            self.log_text.insert("end", "\n".join(batch) + "\n")
            self.log_text.see("end")
            self.log_text.configure(state="disabled")
        self.after(80, self._poll_log)

    # ─────────────────────────────────────────────────────────────────────
    #  SIMULATION EXECUTION
    # ─────────────────────────────────────────────────────────────────────
    def _run_simulation(self):
        if self.is_running:
            return

        try:
            params = self._collect_params()
        except ValueError as e:
            messagebox.showerror("Input Error", str(e))
            return

        self.is_running = True
        self.run_btn.configure(state="disabled", text="RUNNING...",
                               fg_color=GH["orange"])
        self.status_label.configure(text="  RUNNING  ", fg_color=GH["orange"])
        self.progress.set(0)
        self.tabs.set("Log")

        self._log("=" * 60)
        self._log("  LIQLEV SIMULATION STARTED")
        self._log("=" * 60)

        threading.Thread(target=self._simulation_worker,
                         args=(params,), daemon=True).start()

    def _to_british(self, value, unit_key):
        """Convert a value from current display units to British.

        If already in British (si_mode is False), returns value unchanged.
        If in SI mode, divides by the SI conversion factor.
        """
        if not self._si_mode:
            return value
        _, _, factor = UNIT_CONV[unit_key]
        return value / factor

    def _array_to_british(self, values, unit_key):
        """Convert a list of values from current display units to British."""
        if not self._si_mode:
            return values
        _, _, factor = UNIT_CONV[unit_key]
        return [v / factor for v in values]

    def _collect_params(self):
        p = {}
        try:
            p["fluid"]       = self.fluid_var.get()
            p["pinit"]       = self._to_british(float(self.pinit_var.get()),
                                                "pressure")
            p["pfinal"]      = self._to_british(float(self.pfinal_var.get()),
                                                "pressure")
            p["dtank"]       = self._to_british(float(self.dtank_var.get()),
                                                "length_ft")
            p["htank"]       = self._to_british(float(self.htank_var.get()),
                                                "length_ft")
            p["duration"]    = float(self.duration_var.get())
            p["delta_t"]     = float(self.dt_var.get())
            p["ramp_dur"]    = float(self.ramp_dur_var.get())
            p["ramp_factor"] = float(self.ramp_factor_var.get())
        except ValueError:
            raise ValueError("All numeric fields must contain valid numbers.")

        # Validation overrides (optional — blank means use computed values)
        xmlzro_str = self.xmlzro_override_var.get().strip()
        tinit_str = self.tinit_override_var.get().strip()
        p["xmlzro_override"] = float(xmlzro_str) if xmlzro_str else None
        p["tinit_override"] = float(tinit_str) if tinit_str else None

        if p["pinit"] <= p["pfinal"]:
            raise ValueError("Initial pressure must exceed final pressure.")
        if p["dtank"] <= 0 or p["htank"] <= 0:
            raise ValueError("Tank dimensions must be positive.")
        if p["duration"] <= 0 or p["delta_t"] <= 0:
            raise ValueError("Duration and time step must be positive.")
        if p["delta_t"] >= p["duration"]:
            raise ValueError(
                f"Time step ({p['delta_t']}s) must be less than "
                f"duration ({p['duration']}s).")

        # Vent rate: either CSV profile or manual entry
        vent_csv = self.vent_csv_var.get().strip() if hasattr(self, 'vent_csv_var') else ""
        if vent_csv:
            if not os.path.exists(vent_csv):
                raise ValueError(f"Vent rate CSV not found: {vent_csv}")
            p["vent_csv"] = vent_csv
            p["vent_rates"] = [0.0]  # placeholder — CSV overrides
        else:
            p["vent_csv"] = ""
            try:
                raw_vents = parse_array(self.vent_var.get())
                if not raw_vents:
                    raise ValueError()
                p["vent_rates"] = self._array_to_british(raw_vents,
                                                         "mass_flow")
            except Exception:
                raise ValueError("Vent rates: enter comma-separated numbers, "
                                 "range (start:step:stop), or linspace(a,b,n).")

        try:
            p["fills"] = parse_array(self.fill_var.get())
            if not p["fills"]:
                raise ValueError()
        except Exception:
            raise ValueError("Fill fractions: enter comma-separated numbers "
                             "between 0 and 1.")

        for ff in p["fills"]:
            if ff <= 0 or ff > 1:
                raise ValueError(f"Fill fraction {ff} out of range (0, 1].")

        eps_mode = self.eps_mode_var.get()
        if eps_mode == "Custom":
            try:
                p["epsilons"] = parse_array(self.eps_custom_var.get())
                if not p["epsilons"]:
                    raise ValueError()
            except Exception:
                raise ValueError("Custom epsilon must be comma-separated numbers.")
        else:
            p["epsilons"] = [eps_mode]

        p["grav_mode"] = self.grav_mode_var.get()

        if p["grav_mode"] == "Constant":
            try:
                p["const_g"] = float(self.grav_const_val_var.get())
            except ValueError:
                raise ValueError("Gravity level must be a valid number.")

        elif p["grav_mode"] == "Function of Time":
            expr = self.grav_func_var.get().strip()
            if not expr:
                raise ValueError("Gravity function expression cannot be empty.")
            # Validate the expression by test-evaluating at t=0
            try:
                _safe_eval_gravity(expr, 0.0)
            except Exception as ex:
                raise ValueError(
                    f"Gravity function error at t=0:\n{ex}\n\n"
                    f"Expression: g(t) = {expr}\n\n"
                    f"Allowed variable: t\n"
                    f"Allowed functions: sin, cos, tan, exp, log, log10, "
                    f"sqrt, abs, pi, e, min, max")
            p["grav_func_expr"] = expr

        elif p["grav_mode"] == "CSV Profile":
            p["grav_csv"] = self.grav_csv_var.get()
            if not p["grav_csv"]:
                raise ValueError("Please select a gravity CSV file.")
            try:
                p["hold_g"] = float(self.hold_g_var.get())
            except ValueError:
                raise ValueError("Hold G value must be a valid number.")

        return p

    def _simulation_worker(self, p):
        """Run the full parameter sweep in a background thread."""
        t_start = time.perf_counter()
        old_stdout = sys.stdout
        sys.stdout = _LogRedirector(self._log)

        try:
            duration = p["duration"]

            # ── Gravity setup ──
            gravity_function = None  # default: use tggo/xggo arrays

            if p["grav_mode"] == "Constant":
                const_g_ft = p["const_g"] * G_TO_FT_S2
                nggo = 2
                tggo = np.array([0.0, duration])
                xggo = np.array([const_g_ft, const_g_ft])
                self._log(f"[*] Gravity: Constant {p['const_g']} g's")

            elif p["grav_mode"] == "Function of Time":
                expr = p["grav_func_expr"]
                self._log(f"[*] Gravity: g(t) = {expr}")
                gravity_function = _make_gravity_function(expr)
                # Test at a few points and log
                test_times = [0, duration / 4, duration / 2, duration]
                for tt in test_times:
                    gval = _safe_eval_gravity(expr, tt)
                    self._log(f"    g({tt:.1f}s) = {gval:.6f} g's")
                # Still need dummy tggo/xggo (core.py unpacks them even
                # though gravity_function overrides the interpolation)
                nggo = 2
                tggo = np.array([0.0, duration])
                xggo = np.array([0.0, 0.0])  # unused when gravity_function set

            elif p["grav_mode"] == "CSV Profile":
                self._log(f"[*] Gravity: CSV Profile from {p['grav_csv']}")
                g_df = pd.read_csv(p["grav_csv"])
                tggo_g = g_df["normalized_time"].to_numpy()
                xggo_g = g_df["az_positive"].to_numpy()
                last_t = tggo_g[-1]
                self._log(f"    CSV data: {len(tggo_g)} points, "
                          f"0 to {last_t:.2f}s, "
                          f"g range [{xggo_g.min():.6f}, {xggo_g.max():.6f}] g's")
                if duration > last_t:
                    self._log(f"    Extending with hold = "
                              f"{p['hold_g']} g's to {duration}s")
                    tggo_g = np.append(tggo_g, [last_t + 1e-9, duration])
                    xggo_g = np.append(xggo_g, [p["hold_g"], p["hold_g"]])
                nggo = len(tggo_g)
                tggo = tggo_g
                xggo = xggo_g * G_TO_FT_S2

            # ── Pre-compute property table (SPEED BOOST) ──
            fluid = p["fluid"]
            prop_table = None
            if fluid != "Hydrogen":
                self._log("[*] Pre-computing thermodynamic property table...")
                pt_start = time.perf_counter()
                prop_table = build_property_table(fluid, p["pfinal"], p["pinit"])
                pt_time = time.perf_counter() - pt_start
                self._log(f"[+] Property table built in {pt_time:.2f}s "
                          f"(400 points, eliminates per-step CoolProp calls)")

            htank = p["htank"]  # already the tank height in feet
            total_runs = (len(p["fills"]) * len(p["epsilons"])
                          * len(p["vent_rates"]))
            run_count = 0
            new_results = {}

            for fill_frac in p["fills"]:
                for eps_spec in p["epsilons"]:
                    # Determine epsilon schedule
                    if eps_spec == "height_dep":
                        neps, teps, xeps = 0, None, None
                        eps_label = "height_dep"
                    elif eps_spec == "bulk_fake":
                        neps = 2
                        teps = np.array([0.0, duration])
                        xeps = np.array([50.0, 50.0])
                        eps_label = "bulk_fake"
                    elif eps_spec == "AS-203 Schedule":
                        neps = 11
                        teps = np.array([0.0, 20.0, 40.0, 60.0, 80.0,
                                         100.0, 120.0, 140.0, 160.0,
                                         180.0, duration])
                        xeps = np.array([0.0000, 0.0513, 0.1780, 0.2800,
                                         0.3620, 0.4220, 0.4700, 0.5200,
                                         0.5600, 0.6000, 0.6000])
                        eps_label = "AS-203 Schedule"
                    else:
                        neps = 2
                        teps = np.array([0.0, duration])
                        xeps = np.array([float(eps_spec), float(eps_spec)])
                        eps_label = f"{float(eps_spec):.4f}"

                    scenario_key = (f"Fill {fill_frac*100:.0f}%, "
                                    f"eps={eps_label}")
                    result_dfs = []

                    # Vent rate CSV profile override
                    vent_csv_tvmdot = None
                    vent_csv_xvmdot = None
                    if p.get("vent_csv"):
                        vr_df = pd.read_csv(p["vent_csv"])
                        vent_csv_tvmdot = vr_df.iloc[:, 0].to_numpy()
                        vent_csv_xvmdot = vr_df.iloc[:, 1].to_numpy()
                        self._log(f"[*] Vent rate from CSV: "
                                  f"{len(vent_csv_tvmdot)} points")

                    for vent_rate in p["vent_rates"]:
                        run_count += 1
                        self._log(f"\n--- Run {run_count}/{total_runs}: "
                                  f"Fill={fill_frac*100:.0f}%, "
                                  f"eps={eps_label}, "
                                  f"Vent={vent_rate} lbm/s ---")

                        inputs = build_inputs(
                            fluid=fluid,
                            pinit_psia=p["pinit"],
                            pfinal_psia=p["pfinal"],
                            dtank=p["dtank"],
                            htank=p["htank"],
                            fill_fraction=fill_frac,
                            duration=duration,
                            delta_t=p["delta_t"],
                            vent_rate=vent_rate,
                            neps=neps, teps=teps, xeps=xeps,
                            ramp_duration=p["ramp_dur"],
                            ramp_target_factor=p["ramp_factor"],
                            nggo=nggo, tggo=tggo, xggo=xggo,
                            gravity_function=gravity_function,
                            xmlzro_override=p.get("xmlzro_override"),
                            tinit_override=p.get("tinit_override"),
                        )

                        # Override vent arrays if CSV provided
                        if vent_csv_tvmdot is not None:
                            inputs["Nvmd"] = len(vent_csv_tvmdot)
                            inputs["Tvmdot"] = vent_csv_tvmdot
                            inputs["Xvmdot"] = vent_csv_xvmdot

                        # Progress callback — logs live solver stats to GUI
                        def _on_progress(stats):
                            self._log(
                                f"    t={stats['sim_time']:6.2f}s  "
                                f"P={stats['pressure']:7.2f} psia  "
                                f"h={stats['height']:.4f} ft  "
                                f"dh/h0={stats['dh_h0']:+.4f}  "
                                f"[{stats['pct_complete']:4.0f}%  "
                                f"{stats['wall_time']:.1f}s elapsed]")

                        run_start = time.perf_counter()
                        df = liqlev_simulation(inputs, verbose=False,
                                               prop_table=prop_table,
                                               progress_cb=_on_progress)
                        run_time = time.perf_counter() - run_start

                        result_dfs.append(df)

                        if not df.empty:
                            max_row = df.loc[df["Hratio"].idxmax()]
                            exceeds = ("YES" if (df["Height"] >= htank).any()
                                       else "No")
                            conv_fails = (int(df["Conv Failed"].sum())
                                          if "Conv Failed" in df.columns
                                          else 0)
                            self._log(
                                f"  Max dh/h0={max_row['Hratio']:.4f}  "
                                f"t={max_row['Time']:.2f}s  "
                                f"Exceeds={exceeds}  "
                                f"({run_time:.2f}s)")
                            if conv_fails > 0:
                                self._log(
                                    f"  [WARNING] {conv_fails} timesteps "
                                    f"had solver convergence failures")

                        # Update progress bar + % label
                        frac = run_count / total_runs
                        self.after(0, lambda v=frac: (
                            self.progress.set(v),
                            self._progress_label.configure(
                                text=f"{int(v*100)}%")))

                    new_results[scenario_key] = {
                        "dfs": result_dfs,
                        "vent_rates": list(p["vent_rates"]),
                        "fill": fill_frac,
                        "eps_label": eps_label,
                        "htank": htank,
                    }

            self.sweep_results.update(new_results)
            self.scenario_keys = list(self.sweep_results.keys())

            elapsed = time.perf_counter() - t_start
            self._log(f"\n{'=' * 60}")
            self._log(f"  ALL {run_count} SIMULATIONS COMPLETE  "
                      f"({elapsed:.1f}s total)")
            self._log(f"{'=' * 60}")

        except Exception as e:
            self._log(f"\n[ERROR] {e}")
            import traceback
            self._log(traceback.format_exc())
        finally:
            sys.stdout = old_stdout
            self.after(0, self._on_simulation_done)

    def _on_simulation_done(self):
        self.is_running = False
        self.run_btn.configure(state="normal", text="\u25B6  RUN SIMULATION",
                               fg_color=GH["green"])
        self.status_label.configure(text="  READY  ", fg_color=GH["green"])
        self.progress.set(1.0)
        self._progress_label.configure(text="100%")

        if not self.scenario_keys:
            return

        self.scenario_menu.configure(values=self.scenario_keys)
        self.scenario_var.set(self.scenario_keys[-1])
        self._render_current_plot()
        self._update_data_table_controls()
        self._update_summary()
        self.tabs.set("Plots")

    # ─────────────────────────────────────────────────────────────────────
    #  PLOT RENDERING
    # ─────────────────────────────────────────────────────────────────────
    def _on_scenario_change(self, _=None):
        self._render_current_plot()
        self._update_data_table_controls()

    # ─────────────────────────────────────────────────────────────────────
    #  INTERACTIVE PLOT TOOLS
    # ─────────────────────────────────────────────────────────────────────
    def _update_toolbar_highlight(self):
        """Highlight the active tool button."""
        mode = self._plot_mode.get()
        for btn, name in [(self._btn_pointer, "pointer"),
                          (self._btn_pan, "pan"),
                          (self._btn_zoom, "zoom")]:
            if name == mode:
                btn.configure(fg_color=GH["blue"],
                              hover_color=GH["blue"])
            else:
                btn.configure(fg_color=GH["btn_secondary"],
                              hover_color=GH["btn_secondary_hov"])

    def _tool_pointer(self):
        """Switch to pointer mode (no pan/zoom, just coordinate readout)."""
        # Deactivate any active mpl tool
        if self._mpl_toolbar.mode == "pan/zoom":
            self._mpl_toolbar.pan()
        elif self._mpl_toolbar.mode == "zoom rect":
            self._mpl_toolbar.zoom()
        self._plot_mode.set("pointer")
        self._update_toolbar_highlight()

    def _tool_pan(self):
        """Toggle click-and-drag pan mode."""
        if self._plot_mode.get() == "pan":
            # Deactivate
            self._mpl_toolbar.pan()
            self._plot_mode.set("pointer")
        else:
            # Deactivate zoom if active
            if self._mpl_toolbar.mode == "zoom rect":
                self._mpl_toolbar.zoom()
            self._mpl_toolbar.pan()
            self._plot_mode.set("pan")
        self._update_toolbar_highlight()

    def _tool_zoom(self):
        """Toggle rectangle-zoom mode."""
        if self._plot_mode.get() == "zoom":
            self._mpl_toolbar.zoom()
            self._plot_mode.set("pointer")
        else:
            if self._mpl_toolbar.mode == "pan/zoom":
                self._mpl_toolbar.pan()
            self._mpl_toolbar.zoom()
            self._plot_mode.set("zoom")
        self._update_toolbar_highlight()

    def _reset_view(self):
        """Reset all axes to their original limits (after plot was drawn)."""
        for ax_id, (xl, yl) in self._default_view_limits.items():
            for ax in self.fig.get_axes():
                if id(ax) == ax_id:
                    ax.set_xlim(xl)
                    ax.set_ylim(yl)
        self.canvas.draw_idle()

    def _save_plot(self):
        """Save the current figure to a file."""
        path = filedialog.asksaveasfilename(
            title="Save Plot",
            defaultextension=".png",
            filetypes=[("PNG Image", "*.png"), ("PDF", "*.pdf"),
                       ("SVG", "*.svg"), ("All", "*.*")])
        if path:
            self.fig.savefig(path, dpi=200, facecolor=PLOT_BG,
                             edgecolor="none", bbox_inches="tight")
            messagebox.showinfo("Saved", f"Plot saved to:\n{path}")

    def _copy_plot_to_clipboard(self):
        """Copy the current plot to the Windows clipboard as an image."""
        try:
            from PIL import Image as _PILImage, BmpImagePlugin
            # Render plot to PNG, convert to RGB BMP via PIL
            png_buf = io.BytesIO()
            self.fig.savefig(png_buf, format="png", dpi=150, facecolor=PLOT_BG,
                             edgecolor="none", bbox_inches="tight")
            png_buf.seek(0)
            img = _PILImage.open(png_buf).convert("RGB")
            bmp_buf = io.BytesIO()
            img.save(bmp_buf, format="BMP")
            bmp_bytes = bmp_buf.getvalue()
            dib_data = bmp_bytes[14:]  # strip 14-byte BMP file header -> DIB
            png_buf.close()
            bmp_buf.close()

            CF_DIB = 8
            GMEM_MOVEABLE = 0x0002
            user32 = ctypes.windll.user32
            kernel32 = ctypes.windll.kernel32

            kernel32.GlobalAlloc.argtypes = [ctypes.c_uint, ctypes.c_size_t]
            kernel32.GlobalAlloc.restype = ctypes.c_void_p
            kernel32.GlobalLock.argtypes = [ctypes.c_void_p]
            kernel32.GlobalLock.restype = ctypes.c_void_p
            kernel32.GlobalUnlock.argtypes = [ctypes.c_void_p]
            user32.SetClipboardData.argtypes = [ctypes.c_uint, ctypes.c_void_p]

            user32.OpenClipboard(0)
            user32.EmptyClipboard()
            sz = len(dib_data)
            h = kernel32.GlobalAlloc(GMEM_MOVEABLE, sz)
            ptr = kernel32.GlobalLock(h)
            ctypes.memmove(ptr, dib_data, sz)
            kernel32.GlobalUnlock(h)
            user32.SetClipboardData(CF_DIB, h)
            user32.CloseClipboard()

            self._coord_label.configure(text="COPIED TO CLIPBOARD",
                                        text_color=GH["green"])
            self.after(2000, lambda: self._coord_label.configure(
                text="", text_color=GH["cyan"]))
        except Exception as e:
            messagebox.showerror("Clipboard Error", str(e))

    # ─────────────────────────────────────────────────────────────────────
    #  THEME TOGGLE
    # ─────────────────────────────────────────────────────────────────────
    def _toggle_theme(self):
        """Switch between dark and light mode."""
        global PLOT_BG, PLOT_FACE, PLOT_GRID, PLOT_TEXT

        self._dark_mode = not self._dark_mode
        src = _THEME_DARK if self._dark_mode else _THEME_LIGHT
        psrc = _PLOT_DARK if self._dark_mode else _PLOT_LIGHT

        # Swap global theme dicts
        GH.update(src)
        PLOT_BG   = psrc["bg"]
        PLOT_FACE = psrc["face"]
        PLOT_GRID = psrc["grid"]
        PLOT_TEXT  = psrc["text"]

        # Update button label
        if self._dark_mode:
            self._theme_btn.configure(text="\u263e DARK")
        else:
            self._theme_btn.configure(text="\u2600 LIGHT")

        # Apply to all widgets recursively
        self._apply_theme()

        # Redraw current plot if results exist
        if self.sweep_results:
            self._update_plot()

    def _apply_theme(self):
        """Recursively apply the current GH theme to all widgets."""
        # Root window
        self.configure(fg_color=GH["bg"])

        # Walk all widgets and re-color based on type
        self._recolor_widgets(self)

        # Re-style the figure canvas background
        self.fig.set_facecolor(PLOT_BG)
        for ax in self.fig.get_axes():
            ax.set_facecolor(PLOT_FACE)
        self.canvas.draw_idle()

    def _recolor_widgets(self, parent):
        """Recursively update widget colors for the active theme."""
        for child in parent.winfo_children():
            cls = type(child).__name__

            try:
                if cls in ("CTkFrame", "CTkScrollableFrame"):
                    fg = child.cget("fg_color")
                    # Map old colors to new theme
                    if fg in (_THEME_DARK["bg"], _THEME_LIGHT["bg"]):
                        child.configure(fg_color=GH["bg"])
                    elif fg in (_THEME_DARK["surface"], _THEME_LIGHT["surface"]):
                        child.configure(fg_color=GH["surface"])
                    elif fg in (_THEME_DARK["overlay"], _THEME_LIGHT["overlay"]):
                        child.configure(fg_color=GH["overlay"])
                    elif fg in (_THEME_DARK["header"], _THEME_LIGHT["header"]):
                        child.configure(fg_color=GH["header"])
                    elif fg in (_THEME_DARK["btn_secondary"],
                                _THEME_LIGHT["btn_secondary"]):
                        child.configure(fg_color=GH["btn_secondary"])
                    # Update border if present
                    try:
                        bc = child.cget("border_color")
                        if bc in (_THEME_DARK["border"], _THEME_LIGHT["border"]):
                            child.configure(border_color=GH["border"])
                    except Exception:
                        pass
                    # Scrollbar colors
                    try:
                        child.configure(
                            scrollbar_button_color=GH["border"],
                            scrollbar_button_hover_color=GH["muted"])
                    except Exception:
                        pass

                elif cls == "CTkLabel":
                    tc = child.cget("text_color")
                    if tc in (_THEME_DARK["text"], _THEME_LIGHT["text"]):
                        child.configure(text_color=GH["text"])
                    elif tc in (_THEME_DARK["text2"], _THEME_LIGHT["text2"]):
                        child.configure(text_color=GH["text2"])
                    elif tc in (_THEME_DARK["muted"], _THEME_LIGHT["muted"]):
                        child.configure(text_color=GH["muted"])
                    elif tc in (_THEME_DARK["cyan"], _THEME_LIGHT["cyan"]):
                        child.configure(text_color=GH["cyan"])
                    elif tc in (_THEME_DARK["blue"], _THEME_LIGHT["blue"]):
                        child.configure(text_color=GH["blue"])
                    elif tc in (_THEME_DARK["orange"], _THEME_LIGHT["orange"]):
                        child.configure(text_color=GH["orange"])
                    # Label backgrounds
                    try:
                        fg = child.cget("fg_color")
                        if fg in (_THEME_DARK["green"], _THEME_LIGHT["green"]):
                            child.configure(fg_color=GH["green"])
                    except Exception:
                        pass

                elif cls == "CTkButton":
                    fg = child.cget("fg_color")
                    if fg in (_THEME_DARK["btn_secondary"],
                              _THEME_LIGHT["btn_secondary"]):
                        child.configure(
                            fg_color=GH["btn_secondary"],
                            hover_color=GH["btn_secondary_hov"],
                            text_color=GH["text2"],
                            border_color=GH["border"])
                    elif fg in (_THEME_DARK["green"], _THEME_LIGHT["green"]):
                        child.configure(
                            fg_color=GH["green"],
                            hover_color=GH["green_hover"])
                    elif fg in (_THEME_DARK["blue"], _THEME_LIGHT["blue"]):
                        child.configure(fg_color=GH["blue"])
                    elif fg in (_THEME_DARK["red"], _THEME_LIGHT["red"]):
                        child.configure(fg_color=GH["red"])

                elif cls == "CTkEntry":
                    child.configure(
                        fg_color=GH["input_bg"],
                        text_color=GH["text"],
                        border_color=GH["border"])

                elif cls == "CTkOptionMenu":
                    child.configure(
                        fg_color=GH["input_bg"],
                        button_color=GH["border"],
                        button_hover_color=GH["muted"],
                        text_color=GH["text"])

                elif cls == "CTkTextbox":
                    child.configure(
                        fg_color=GH["surface"],
                        text_color=GH["text"])

                elif cls == "CTkProgressBar":
                    child.configure(
                        fg_color=GH["border"],
                        progress_color=GH["blue"])

            except Exception:
                pass  # Some widgets may not support all cget/configure calls

            # Recurse into children
            self._recolor_widgets(child)

    def _toggle_units(self):
        """Toggle between Imperial and SI display units.

        Converts every registered input field value and relabels units.
        The solver always receives British, so _collect_params converts
        back before building inputs.
        """
        going_to_si = not self._si_mode
        self._si_mode = going_to_si

        if going_to_si:
            self._unit_btn.configure(text="SI", fg_color=GH["blue"],
                                     text_color="#ffffff")
        else:
            self._unit_btn.configure(text="IMPERIAL",
                                     fg_color=GH["btn_secondary"],
                                     text_color=GH["text2"])

        # Convert every registered field
        for var, lbl, ukey, alt_lbl in self._unit_fields:
            brit_u, si_u, factor = UNIT_CONV[ukey]
            # Update the primary unit label
            if going_to_si:
                lbl.configure(text=si_u)
            else:
                lbl.configure(text=brit_u)
            # Convert the value(s) in the field
            raw = var.get().strip()
            if not raw:
                continue
            try:
                vals = parse_array(raw)
                if going_to_si:
                    converted = [v * factor for v in vals]
                else:
                    converted = [v / factor for v in vals]
                if len(converted) == 1:
                    var.set(f"{converted[0]:.6g}")
                else:
                    var.set(", ".join(f"{v:.6g}" for v in converted))
            except Exception:
                pass  # leave unparseable text alone
            # Alt label refreshes automatically via the var trace

        # Update computed labels (vol, Ac)
        self._update_computed()

        # Re-render results with new units
        if self.sweep_results:
            self._render_current_plot()
            self._update_summary()
            self._update_data_table()

    def _store_default_limits(self):
        """Capture axis limits right after drawing, for reset."""
        self._default_view_limits.clear()
        for ax in self.fig.get_axes():
            self._default_view_limits[id(ax)] = (
                ax.get_xlim(), ax.get_ylim())

    def _on_scroll_zoom(self, event):
        """Scroll-wheel zoom centered on the cursor position."""
        if event.inaxes is None:
            return
        ax = event.inaxes
        scale = 0.8 if event.button == "up" else 1.25  # zoom in / out

        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        # Zoom centered on cursor
        xdata, ydata = event.xdata, event.ydata
        new_w = (xlim[1] - xlim[0]) * scale
        new_h = (ylim[1] - ylim[0]) * scale
        relx = (xdata - xlim[0]) / (xlim[1] - xlim[0])
        rely = (ydata - ylim[0]) / (ylim[1] - ylim[0])

        ax.set_xlim(xdata - new_w * relx, xdata + new_w * (1 - relx))
        ax.set_ylim(ydata - new_h * rely, ydata + new_h * (1 - rely))
        self.canvas.draw_idle()

    def _on_mouse_move(self, event):
        """Update coordinate readout and crosshair on cursor movement."""
        # Remove old crosshair lines
        for line in self._crosshair_lines:
            try:
                line.remove()
            except Exception:
                pass
        self._crosshair_lines.clear()

        if event.inaxes is not None:
            self._coord_label.configure(
                text=f"x = {event.xdata:.3f}   y = {event.ydata:.5f}")
            # Draw crosshair
            ax = event.inaxes
            h_line = ax.axhline(y=event.ydata, color=GH["cyan"],
                                linewidth=0.5, alpha=0.4, linestyle=":")
            v_line = ax.axvline(x=event.xdata, color=GH["cyan"],
                                linewidth=0.5, alpha=0.4, linestyle=":")
            self._crosshair_lines = [h_line, v_line]
            self.canvas.draw_idle()
        else:
            self._coord_label.configure(text="")
            self.canvas.draw_idle()

    # ─────────────────────────────────────────────────────────────────────
    #  PLOT RENDERING
    # ─────────────────────────────────────────────────────────────────────
    def _on_plot_change(self, _=None):
        self._render_current_plot()

    def _get_current_scenario(self):
        return self.sweep_results.get(self.scenario_var.get())

    def _get_threshold(self):
        """Return the threshold value or None."""
        try:
            v = float(self._threshold_var.get())
            return v if v > 0 else None
        except (ValueError, TypeError):
            return None

    def _get_overlay_scenarios(self):
        """Return list of (label, scenario) for overlay mode."""
        if self._overlay_var.get() and len(self.sweep_results) > 1:
            return list(self.sweep_results.items())
        s = self._get_current_scenario()
        if s:
            return [(self.scenario_var.get(), s)]
        return []

    def _render_current_plot(self):
        scenarios = self._get_overlay_scenarios()
        if not scenarios:
            return

        self.fig.clear()

        dispatch = {
            "Liquid Level Rise (dh/h0)":    self._plot_level_rise,
            "Liquid Level Rise (inches)":   self._plot_level_inches,
            "Tank Pressure":                self._plot_pressure,
            "Vent Rate Profile":            self._plot_vent_rate,
            "Epsilon":                      self._plot_epsilon,
            "Diagnostics":                  self._plot_diagnostics,
            "Vapor Generation":             self._plot_vapor_gen,
        }

        renderer = dispatch.get(self.plot_var.get(), self._plot_level_rise)
        renderer(scenarios)
        self.fig.set_facecolor(PLOT_BG)
        self.canvas.draw_idle()

        # Check threshold warning for pulsing aesthetic
        thr = self._get_threshold()
        warning = False
        if thr is not None:
            for lbl, s in scenarios:
                for df in s["dfs"]:
                    if not df.empty and df["Hratio"].max() > thr:
                        warning = True
                        break
        self._set_warning_state(warning)

        # Capture limits for reset button
        self._store_default_limits()

    def _set_warning_state(self, is_warning):
        if is_warning:
            if not getattr(self, "_warning_active", False):
                self._warning_active = True
                self._pulse_warning()
        else:
            self._warning_active = False
            if hasattr(self, "_thresh_entry_ref"):
                self._thresh_entry_ref.configure(border_color=GH["border"])

    def _pulse_warning(self):
        if not getattr(self, "_warning_active", False):
            return
        if not hasattr(self, "_thresh_entry_ref"):
            return
        
        current = self._thresh_entry_ref.cget("border_color")
        next_color = GH["red"] if current == GH["border"] else GH["border"]
        self._thresh_entry_ref.configure(border_color=next_color)
        self.after(600, self._pulse_warning)

    def _colors(self, n):
        return PLOT_COLORS[:max(n, 1)] if n <= len(PLOT_COLORS) else \
            [matplotlib.cm.viridis(i) for i in np.linspace(0.15, 0.85, n)]

    def _style_ax(self, ax, title, xlabel, ylabel):
        ax.set_facecolor(PLOT_FACE)
        ax.set_title(title.upper(), fontsize=BASE_FONT + 1, fontweight="bold",
                     color=GH["text"], pad=10, family="Consolas")
        ax.set_xlabel(xlabel, fontsize=BASE_FONT, color=GH["text2"],
                      family="Consolas")
        ax.set_ylabel(ylabel, fontsize=BASE_FONT, color=GH["text2"],
                      family="Consolas")
        ax.tick_params(colors=GH["text"], labelsize=BASE_FONT - 2,
                       labelcolor=GH["text2"])
        ax.grid(True, alpha=0.15, color=PLOT_GRID, linewidth=0.5)
        for spine in ax.spines.values():
            spine.set_color(GH["border"])
        # Unit badge
        unit_str = "SI" if self._si_mode else "IMPERIAL"
        ax.text(0.99, 0.99, unit_str, transform=ax.transAxes,
                ha="right", va="top", fontsize=8, family="Consolas",
                color=GH["muted"], alpha=0.6,
                bbox=dict(boxstyle="round,pad=0.2",
                          facecolor=GH["overlay"], edgecolor=GH["border"],
                          alpha=0.5))

    def _legend(self, ax):
        leg = ax.legend(fontsize=BASE_FONT - 2, facecolor=GH["overlay"],
                        edgecolor=GH["border"],
                        labelcolor=GH["text"], framealpha=0.92,
                        prop={"family": "Consolas", "size": BASE_FONT - 2})

    def _mark_extremes(self, ax, df, htank, y_col, h0_offset=0, scale=1):
        """Add red X (tank exceeded) or green star (max level) markers
        with text annotation showing peak value and time."""
        if (df["Height"] >= htank).any():
            idx = df[df["Height"] >= htank].index[0]
            t_val = df.loc[idx, "Time"]
            y = (df.loc[idx, y_col] - h0_offset) * scale
            ax.plot(t_val, y, "x",
                    color=GH["red"], markersize=12, markeredgewidth=2.5,
                    zorder=5)
            ax.annotate(f"EXCEEDED @ {t_val:.1f}s",
                        xy=(t_val, y), xytext=(8, 8),
                        textcoords="offset points",
                        fontsize=9, fontweight="bold", color=GH["red"],
                        family="Consolas",
                        bbox=dict(boxstyle="round,pad=0.3",
                                  facecolor=GH["overlay"], edgecolor=GH["red"],
                                  alpha=0.9))
        else:
            idx = df["Hratio"].idxmax()
            t_val = df.loc[idx, "Time"]
            y = (df.loc[idx, y_col] - h0_offset) * scale
            ax.plot(t_val, y, marker="*",
                    color="#3fb950", markersize=14, markeredgecolor="black",
                    zorder=5)
            ax.axhline(y=y, color="#3fb950", linestyle="--",
                       linewidth=0.7, alpha=0.4)
            ax.annotate(f"Max: {y:.4f} @ {t_val:.1f}s",
                        xy=(t_val, y), xytext=(8, -14),
                        textcoords="offset points",
                        fontsize=9, fontweight="bold", color="#3fb950",
                        family="Consolas",
                        bbox=dict(boxstyle="round,pad=0.3",
                                  facecolor=GH["overlay"],
                                  edgecolor="#3fb950", alpha=0.9))

    def _iter_traces(self, scenarios):
        """Yield (label, df, rate, color_idx) for each trace across scenarios."""
        ci = 0
        overlay = len(scenarios) > 1
        unit_str = _u('mass_flow', self._si_mode)
        for key, s in scenarios:
            for df, rate in zip(s["dfs"], s["vent_rates"]):
                if df.empty:
                    continue
                disp_rate = _cv(rate, 'mass_flow', self._si_mode)
                disp_rate_str = f"{disp_rate:g}"
                if overlay:
                    lbl = f"{key} | {disp_rate_str} {unit_str}"
                else:
                    lbl = f"{disp_rate_str} {unit_str}"
                yield lbl, df, rate, s, ci
                ci += 1

    def _draw_threshold(self, ax):
        """Draw threshold line on level-rise plots."""
        thr = self._get_threshold()
        if thr is not None:
            ax.axhline(y=thr, color=GH["red"], linestyle="--",
                       linewidth=1.5, alpha=0.8, label=f"Threshold = {thr}")

    def _plot_level_rise(self, scenarios):
        ax = self.fig.add_subplot(111)
        traces = list(self._iter_traces(scenarios))
        colors = self._colors(len(traces))

        for lbl, df, rate, s, ci in traces:
            ax.plot(df["Time"], df["Hratio"], linewidth=1.8,
                    color=colors[ci], label=lbl)
            ax.fill_between(df["Time"], df["Hratio"], color=colors[ci], alpha=0.15)
            self._mark_extremes(ax, df, s["htank"], "Hratio")

        self._draw_threshold(ax)
        title = "Dimensionless Liquid Level Rise"
        if len(scenarios) == 1:
            s = scenarios[0][1]
            title += f"  (Fill: {s['fill']*100:.0f}%, eps: {s['eps_label']})"
        self._style_ax(ax, title, "Time (s)", "dh/h\u2080")
        self._legend(ax)
        self.fig.tight_layout()

    def _plot_level_inches(self, scenarios):
        ax = self.fig.add_subplot(111)
        si = self._si_mode
        traces = list(self._iter_traces(scenarios))
        colors = self._colors(len(traces))

        for lbl, df, rate, s, ci in traces:
            h0 = df["Height"].iloc[0]
            dh = (df["Height"] - h0) * 12  # inches
            if si:
                dh = dh * 2.54  # -> cm
            ax.plot(df["Time"], dh, linewidth=1.8,
                    color=colors[ci], label=lbl)
            ax.fill_between(df["Time"], dh, color=colors[ci], alpha=0.15)
            scale = 12 * 2.54 if si else 12
            self._mark_extremes(ax, df, s["htank"], "Height",
                                h0_offset=h0, scale=scale)

        self._draw_threshold(ax)
        unit = "cm" if si else "inches"
        title = "Liquid Level Increase"
        if len(scenarios) == 1:
            s = scenarios[0][1]
            title += f"  (Fill: {s['fill']*100:.0f}%, eps: {s['eps_label']})"
        self._style_ax(ax, title, "Time (s)", f"\u0394h ({unit})")
        self._legend(ax)
        self.fig.tight_layout()

    def _plot_pressure(self, scenarios):
        ax = self.fig.add_subplot(111)
        si = self._si_mode
        traces = list(self._iter_traces(scenarios))
        colors = self._colors(len(traces))

        for lbl, df, rate, s, ci in traces:
            p = _cv(df["Press"], "pressure", si) if si else df["Press"]
            ax.plot(df["Time"], p, linewidth=1.8,
                    color=colors[ci], label=lbl)
            ax.fill_between(df["Time"], p, color=colors[ci], alpha=0.15)
            if not df.empty:
                idx = df["Hratio"].idxmax()
                pv = _cv(df.loc[idx, "Press"], "pressure", si)
                ax.plot(df.loc[idx, "Time"], pv,
                        marker="*", color="#3fb950", markersize=14,
                        markeredgecolor="black", zorder=5)

        title = "Tank Pressure"
        if len(scenarios) == 1:
            s = scenarios[0][1]
            title += f"  (Fill: {s['fill']*100:.0f}%, eps: {s['eps_label']})"
        self._style_ax(ax, title, "Time (s)",
                       f"Pressure ({_u('pressure', si)})")
        self._legend(ax)
        self.fig.tight_layout()

    def _plot_vent_rate(self, scenarios):
        ax = self.fig.add_subplot(111)
        si = self._si_mode
        traces = list(self._iter_traces(scenarios))
        colors = self._colors(len(traces))

        for lbl, df, rate, s, ci in traces:
            vr = _cv(df["Vent Rate"], "mass_flow", si) if si \
                 else df["Vent Rate"]
            ax.plot(df["Time"], vr, linewidth=1.8,
                    color=colors[ci], label=lbl)
            ax.fill_between(df["Time"], vr, color=colors[ci], alpha=0.15)

        self._style_ax(ax, "Vent Rate Profile", "Time (s)",
                       f"Vent Rate ({_u('mass_flow', si)})")
        self._legend(ax)
        self.fig.tight_layout()

    def _plot_epsilon(self, scenarios):
        ax = self.fig.add_subplot(111)
        traces = list(self._iter_traces(scenarios))
        colors = self._colors(len(traces))

        for lbl, df, rate, s, ci in traces:
            ax.plot(df["Time"], df["eps"], linewidth=1.8,
                    color=colors[ci], label=lbl)
            ax.fill_between(df["Time"], df["eps"], color=colors[ci], alpha=0.15)

        title = "Epsilon vs Time"
        if len(scenarios) == 1:
            title += f"  (Mode: {scenarios[0][1]['eps_label']})"
        self._style_ax(ax, title, "Time (s)", "\u03b5")
        self._legend(ax)
        self.fig.tight_layout()

    def _plot_diagnostics(self, scenarios):
        ax1 = self.fig.add_subplot(211)
        ax2 = self.fig.add_subplot(212, sharex=ax1)
        si = self._si_mode
        traces = list(self._iter_traces(scenarios))
        colors = self._colors(len(traces))

        ax1t = ax1.twinx()
        first_df = None
        for lbl, df, rate, s, ci in traces:
            if first_df is None:
                first_df = df
            ax1.plot(df["Time"], df["Superheat"], linewidth=1.8,
                     color=colors[ci], label=f"SH: {lbl}")
            dpdt = _cv(df["dP/dtha"], "pressure", si) if si \
                   else df["dP/dtha"]
            ax1t.plot(df["Time"], dpdt, linewidth=1.2,
                      linestyle=":", color=colors[ci], alpha=0.7,
                      label=f"dP/dt: {lbl}")

        p_unit = _u("pressure", si)
        self._style_ax(ax1, "Superheat & Depressurization", "",
                       "Superheat (K)")
        ax1t.set_ylabel(f"dP/dt ({p_unit}/s)", fontsize=12,
                        color=GH["text2"])
        ax1t.tick_params(colors=GH["text"], labelsize=10,
                         labelcolor=GH["text2"])
        lines1, lab1 = ax1.get_legend_handles_labels()
        lines2, lab2 = ax1t.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, lab1 + lab2, fontsize=7,
                   facecolor=GH["overlay"], edgecolor=GH["border"],
                   labelcolor=GH["text"], ncol=2, loc="upper left")

        ax2t = ax2.twinx()
        for lbl, df, rate, s, ci in traces:
            vol_unit = _u("volume", si)
            vbl = _cv(df["VBL vol"], "volume", si) if si else df["VBL vol"]
            ax2.plot(df["Time"], vbl, linewidth=1.8,
                     color=colors[ci], label=f"VBL: {lbl}")

        if first_df is not None and not first_df.empty:
            ax2t.plot(first_df["Time"], first_df["Gravity_g"],
                      color=GH["muted"], linestyle="--", linewidth=1.5,
                      alpha=0.6, label="Gravity")

        vol_u = _u("volume", si)
        self._style_ax(ax2, "Boundary Layer Volume & Gravity",
                       "Time (s)", f"BL Volume ({vol_u})")
        ax2.set_yscale("log")
        ax2t.set_ylabel("Gravity (g's)", fontsize=12, color=GH["text2"])
        ax2t.tick_params(colors=GH["text"], labelsize=10,
                         labelcolor=GH["text2"])
        lines1, lab1 = ax2.get_legend_handles_labels()
        lines2, lab2 = ax2t.get_legend_handles_labels()
        ax2.legend(lines1 + lines2, lab1 + lab2, fontsize=8,
                   facecolor=GH["overlay"], edgecolor=GH["border"],
                   labelcolor=GH["text"])

        self.fig.tight_layout()

    def _plot_vapor_gen(self, scenarios):
        ax1 = self.fig.add_subplot(211)
        ax2 = self.fig.add_subplot(212, sharex=ax1)
        traces = list(self._iter_traces(scenarios))
        colors = self._colors(len(traces))

        for lbl, df, rate, s, ci in traces:
            ax1.plot(df["Time"], df["Vap Gen Rate (kg/s)"], linewidth=1.8,
                     color=colors[ci], label=lbl)
            ax2.plot(df["Time"], df["Total Vap Gen (kg)"], linewidth=1.8,
                     color=colors[ci], label=lbl)

        title = "Vapor Generation Rate"
        if len(scenarios) == 1:
            s = scenarios[0][1]
            title += f"  (Fill: {s['fill']*100:.0f}%, eps: {s['eps_label']})"
        self._style_ax(ax1, title, "", "Rate (kg/s)")
        self._style_ax(ax2, "Cumulative Vapor Mass",
                       "Time (s)", "Total (kg)")
        self._legend(ax1)
        self.fig.tight_layout()

    # ─────────────────────────────────────────────────────────────────────
    #  DATA TABLE
    # ─────────────────────────────────────────────────────────────────────
    def _update_data_table_controls(self):
        scenario = self._get_current_scenario()
        if not scenario:
            return
        rate_strs = [str(r) for r in scenario["vent_rates"]]
        self.table_rate_menu.configure(values=rate_strs)
        self.table_rate_var.set(rate_strs[0])
        self._update_data_table()

    def _update_data_table(self):
        scenario = self._get_current_scenario()
        if not scenario:
            return

        rate_str = self.table_rate_var.get()
        try:
            idx = scenario["vent_rates"].index(float(rate_str))
        except (ValueError, IndexError):
            idx = 0
        df = scenario["dfs"][idx]

        if df.empty:
            return

        # Clear
        self.tree.delete(*self.tree.get_children())

        # Columns
        cols = list(df.columns)
        self.tree["columns"] = cols
        for col in cols:
            short = col if len(col) <= 18 else col[:16] + ".."
            self.tree.heading(col, text=short,
                              command=lambda c=col: self._sort_table(c))
            self.tree.column(col, width=110, minwidth=70, anchor="e")

        # Batch insert (much faster than one-by-one)
        max_rows = min(len(df), 5000)
        data = df.iloc[:max_rows]

        # Pre-format all values at once
        rows = []
        for i in range(max_rows):
            vals = []
            for col in cols:
                v = data.iat[i, data.columns.get_loc(col)]
                if isinstance(v, float):
                    if abs(v) > 1e4 or (0 < abs(v) < 1e-3):
                        vals.append(f"{v:.4e}")
                    else:
                        vals.append(f"{v:.6f}")
                else:
                    vals.append(str(v))
            rows.append(vals)

        for i, vals in enumerate(rows):
            tag = "even" if i % 2 == 0 else "odd"
            self.tree.insert("", "end", values=vals, tags=(tag,))

        self.tree.tag_configure("even", background=GH["bg"])
        self.tree.tag_configure("odd", background=GH["overlay"])

    def _sort_table(self, col):
        items = [(self.tree.set(k, col), k) for k in self.tree.get_children()]
        try:
            items.sort(key=lambda t: float(t[0]))
        except ValueError:
            items.sort(key=lambda t: t[0])
        for index, (_, k) in enumerate(items):
            self.tree.move(k, "", index)

    # ─────────────────────────────────────────────────────────────────────
    #  SUMMARY
    # ─────────────────────────────────────────────────────────────────────
    def _update_summary(self):
        self.summary_text.configure(state="normal")
        self.summary_text.delete("1.0", "end")

        si = self._si_mode
        thr = self._get_threshold()

        L = []
        L.append("=" * 76)
        mode_str = "SI" if si else "Imperial"
        L.append(f"  LIQLEV SIMULATION RESULTS SUMMARY  [{mode_str}]")
        L.append("=" * 76)
        if thr is not None:
            L.append(f"  ** THRESHOLD: dh/h0 = {thr} **")
        L.append("")

        p_u = _u("pressure", si)
        mf_u = _u("mass_flow", si)
        dh_u = "cm" if si else "inches"

        for key, s in self.sweep_results.items():
            L.append(f"  SCENARIO: {key}")
            L.append("  " + "-" * 72)
            hdr_thr = "  Thresh?" if thr else ""
            L.append(f"  {'Vent Rate':<16} {'Max dh/h0':<14} "
                     f"{'Max dh':<14} {'Time (s)':<12} "
                     f"{'Final P':<12} {'Exceeds?':<10}{hdr_thr}")
            L.append(f"  {'(' + mf_u + ')':<16} {'(-)':<14} "
                     f"{'(' + dh_u + ')':<14} {'(s)':<12} "
                     f"{'(' + p_u + ')':<12} {'(Y/N)':<10}")
            L.append("  " + "-" * 72)

            htank = s["htank"]
            for df, rate in zip(s["dfs"], s["vent_rates"]):
                if df.empty:
                    L.append(f"  {rate:<16} {'N/A'}")
                    continue

                mr = df.loc[df["Hratio"].idxmax()]
                h0 = df["Height"].iloc[0]
                dh_raw = (mr["Height"] - h0) * 12  # inches
                dh_disp = dh_raw * 2.54 if si else dh_raw
                final_p = _cv(df["Press"].iloc[-1], "pressure", si)
                rate_disp = _cv(rate, "mass_flow", si)
                exceeds = "YES" if (df["Height"] >= htank).any() else "No"
                thr_flag = ""
                if thr is not None:
                    thr_flag = "  FAIL" if mr["Hratio"] > thr else "  PASS"

                L.append(f"  {rate_disp:<16.6f} {mr['Hratio']:<14.4f} "
                         f"{dh_disp:<14.4f} {mr['Time']:<12.2f} "
                         f"{final_p:<12.2f} {exceeds:<10}{thr_flag}")
            L.append("")

        # Reference epsilon
        if self.sweep_results:
            try:
                dtank = float(self.dtank_var.get())
                htank_ref = float(self.htank_var.get())
                h_u = _u("length_ft", si)
                L.append("  REFERENCE: Geometry-Based Epsilon")
                L.append("  " + "-" * 42)
                L.append(f"  {'Fill %':<10} {'Height':<14} {'Epsilon':<10}")
                L.append(f"  {'(%)':<10} {'(' + h_u + ')':<14} {'(-)':<10}")
                L.append("  " + "-" * 42)
                for ff in [0.25, 0.50, 0.75, 1.00]:
                    h_c = ff * htank_ref
                    h_disp = _cv(h_c, "length_ft", si)
                    ev = calculate_epsilon(h_c, dtank)
                    L.append(f"  {ff*100:<10.0f} {h_disp:<14.4f} {ev:<10.4f}")
                L.append("")
            except Exception:
                pass

        # Monte Carlo results summary
        if hasattr(self, '_mc_results') and self._mc_results is not None:
            mc = self._mc_results
            L.append("  MONTE CARLO SENSITIVITY ANALYSIS")
            L.append("  " + "-" * 50)
            L.append(f"  Samples:      {mc['n']}")
            L.append(f"  Max dh/h0:    {mc['max_dh']:.4f}  "
                     f"(vent={mc['worst']['vent']:.6f}, "
                     f"fill={mc['worst']['fill']:.2f}, "
                     f"g={mc['worst']['grav']:.6f})")
            L.append(f"  Mean dh/h0:   {mc['mean_dh']:.4f}")
            L.append(f"  Std dh/h0:    {mc['std_dh']:.4f}")
            L.append(f"  95th pctile:  {mc['p95']:.4f}")
            L.append(f"  99th pctile:  {mc['p99']:.4f}")
            if thr is not None:
                n_fail = sum(1 for v in mc['all_dh'] if v > thr)
                L.append(f"  Exceed threshold ({thr}):  "
                         f"{n_fail}/{mc['n']} "
                         f"({100*n_fail/mc['n']:.1f}%)")
            L.append("")

        L.append("=" * 76)

        self.summary_text.insert("1.0", "\n".join(L))
        self.summary_text.configure(state="disabled")

    # ─────────────────────────────────────────────────────────────────────
    #  MONTE CARLO / SENSITIVITY
    # ─────────────────────────────────────────────────────────────────────
    def _run_monte_carlo(self):
        if self.is_running:
            return
        try:
            mc_p = {
                "n": int(self.mc_n_var.get()),
                "vent_min": self._to_british(
                    float(self.mc_vent_min_var.get()), "mass_flow"),
                "vent_max": self._to_british(
                    float(self.mc_vent_max_var.get()), "mass_flow"),
                "fill_min": float(self.mc_fill_min_var.get()),
                "fill_max": float(self.mc_fill_max_var.get()),
                "grav_min": float(self.mc_grav_min_var.get()),
                "grav_max": float(self.mc_grav_max_var.get()),
            }
            # Validate ranges
            if mc_p["vent_min"] >= mc_p["vent_max"]:
                raise ValueError("Vent Rate Min must be less than Max.")
            if mc_p["fill_min"] >= mc_p["fill_max"]:
                raise ValueError("Fill Frac Min must be less than Max.")
            if mc_p["grav_min"] >= mc_p["grav_max"]:
                raise ValueError("Gravity Min must be less than Max.")
            if mc_p["n"] < 2:
                raise ValueError("N Samples must be at least 2.")
            # Also need base simulation params
            base = self._collect_params()
        except ValueError as e:
            messagebox.showerror("Input Error", str(e))
            return

        self.is_running = True
        self.mc_btn.configure(state="disabled", text="RUNNING MC...",
                              fg_color=GH["orange"])
        self.status_label.configure(text="  MONTE CARLO  ",
                                    fg_color=GH["purple"])
        self.progress.set(0)
        self.tabs.set("Log")

        self._log("=" * 60)
        self._log("  MONTE CARLO SENSITIVITY ANALYSIS")
        self._log(f"  N={mc_p['n']} samples")
        self._log(f"  Vent: [{mc_p['vent_min']}, {mc_p['vent_max']}] lbm/s")
        self._log(f"  Fill: [{mc_p['fill_min']}, {mc_p['fill_max']}]")
        self._log(f"  Grav: [{mc_p['grav_min']}, {mc_p['grav_max']}] g's")
        self._log("=" * 60)

        threading.Thread(target=self._mc_worker,
                         args=(mc_p, base), daemon=True).start()

    def _mc_worker(self, mc_p, base):
        old_stdout = sys.stdout
        sys.stdout = _LogRedirector(self._log)
        t_start = time.perf_counter()

        try:
            fluid = base["fluid"]
            duration = base["duration"]
            htank = base["htank"]

            prop_table = None
            if fluid != "Hydrogen":
                self._log("[*] Building property table...")
                prop_table = build_property_table(
                    fluid, base["pfinal"], base["pinit"])

            n = mc_p["n"]
            rng = np.random.default_rng()

            all_dh = []
            all_params = []
            worst = {"dh": 0, "vent": 0, "fill": 0, "grav": 0}

            for i in range(n):
                vent = rng.uniform(mc_p["vent_min"], mc_p["vent_max"])
                fill = rng.uniform(mc_p["fill_min"], mc_p["fill_max"])
                grav = rng.uniform(mc_p["grav_min"], mc_p["grav_max"])

                grav_ft = grav * G_TO_FT_S2
                nggo = 2
                tggo = np.array([0.0, duration])
                xggo = np.array([grav_ft, grav_ft])

                eps_mode = base["epsilons"][0]
                if eps_mode == "height_dep":
                    neps, teps, xeps = 0, None, None
                elif eps_mode == "bulk_fake":
                    neps = 2
                    teps = np.array([0.0, duration])
                    xeps = np.array([50.0, 50.0])
                elif eps_mode == "AS-203 Schedule":
                    neps = 11
                    teps = np.array([0.0, 20.0, 40.0, 60.0, 80.0,
                                     100.0, 120.0, 140.0, 160.0,
                                     180.0, duration])
                    xeps = np.array([0.0000, 0.0513, 0.1780, 0.2800,
                                     0.3620, 0.4220, 0.4700, 0.5200,
                                     0.5600, 0.6000, 0.6000])
                else:
                    neps = 2
                    teps = np.array([0.0, duration])
                    xeps = np.array([float(eps_mode), float(eps_mode)])

                inputs = build_inputs(
                    fluid=fluid,
                    pinit_psia=base["pinit"],
                    pfinal_psia=base["pfinal"],
                    dtank=base["dtank"],
                    htank=htank,
                    fill_fraction=fill,
                    duration=duration,
                    delta_t=base["delta_t"],
                    vent_rate=vent,
                    neps=neps, teps=teps, xeps=xeps,
                    ramp_duration=base["ramp_dur"],
                    ramp_target_factor=base["ramp_factor"],
                    nggo=nggo, tggo=tggo, xggo=xggo,
                )

                df = liqlev_simulation(inputs, verbose=False,
                                       prop_table=prop_table)
                if not df.empty:
                    max_dh = df["Hratio"].max()
                else:
                    max_dh = 0.0

                all_dh.append(max_dh)
                all_params.append({"vent": vent, "fill": fill, "grav": grav})

                if max_dh > worst["dh"]:
                    worst = {"dh": max_dh, "vent": vent,
                             "fill": fill, "grav": grav}

                if (i + 1) % max(1, n // 20) == 0 or i == n - 1:
                    self._log(f"  [{i+1}/{n}] max_so_far={worst['dh']:.4f}")
                    frac = (i + 1) / n
                    self.after(0, lambda v=frac: (
                        self.progress.set(v),
                        self._progress_label.configure(
                            text=f"{int(v*100)}%")))

            arr = np.array(all_dh)
            self._mc_results = {
                "n": n,
                "all_dh": all_dh,
                "all_params": all_params,
                "max_dh": arr.max(),
                "mean_dh": arr.mean(),
                "std_dh": arr.std(),
                "p95": np.percentile(arr, 95),
                "p99": np.percentile(arr, 99),
                "worst": worst,
            }

            elapsed = time.perf_counter() - t_start
            self._log(f"\n{'=' * 60}")
            self._log(f"  MONTE CARLO COMPLETE  ({n} samples, "
                      f"{elapsed:.1f}s)")
            self._log(f"  Worst case: dh/h0={worst['dh']:.4f}")
            self._log(f"  Mean: {arr.mean():.4f}  "
                      f"Std: {arr.std():.4f}  "
                      f"95th: {np.percentile(arr, 95):.4f}")
            self._log(f"{'=' * 60}")

        except Exception as e:
            self._log(f"\n[ERROR] {e}")
            import traceback
            self._log(traceback.format_exc())
        finally:
            sys.stdout = old_stdout
            self.after(0, self._on_mc_done)

    def _on_mc_done(self):
        self.is_running = False
        self.mc_btn.configure(state="normal", text="RUN MONTE CARLO",
                              fg_color=GH["purple"])
        self.status_label.configure(text="  READY  ", fg_color=GH["green"])
        self.progress.set(1.0)
        self._progress_label.configure(text="100%")

        if hasattr(self, '_mc_results') and self._mc_results:
            self._update_summary()
            # Show MC histogram
            self._plot_mc_histogram()
            self.tabs.set("Plots")

    def _plot_mc_histogram(self):
        """Show the Monte Carlo results as a histogram."""
        if not hasattr(self, '_mc_results') or not self._mc_results:
            return
        mc = self._mc_results
        self.fig.clear()

        ax = self.fig.add_subplot(111)
        ax.set_facecolor(PLOT_FACE)

        arr = np.array(mc["all_dh"])
        n_bins = min(30, max(10, mc["n"] // 5))
        ax.hist(arr, bins=n_bins, color=GH["blue"], alpha=0.8,
                edgecolor=GH["border"], linewidth=0.5)

        # 1σ and 2σ confidence bands
        mean = mc["mean_dh"]
        std = mc["std_dh"]
        ylims = ax.get_ylim()
        ax.axvspan(mean - 2 * std, mean + 2 * std,
                   alpha=0.08, color=GH["purple"],
                   label=f"2\u03c3 ({mean-2*std:.4f} \u2013 {mean+2*std:.4f})")
        ax.axvspan(mean - std, mean + std,
                   alpha=0.15, color=GH["purple"],
                   label=f"1\u03c3 ({mean-std:.4f} \u2013 {mean+std:.4f})")

        # Stats lines
        ax.axvline(mean, color=GH["green"], linestyle="-",
                   linewidth=2, label=f"Mean = {mean:.4f}")
        ax.axvline(mc["p95"], color=GH["orange"], linestyle="--",
                   linewidth=1.5, label=f"95th = {mc['p95']:.4f}")
        ax.axvline(mc["p99"], color=GH["red"], linestyle="--",
                   linewidth=1.5, label=f"99th = {mc['p99']:.4f}")
        ax.axvline(mc["max_dh"], color="#ff7b72", linestyle=":",
                   linewidth=1.5, label=f"Max = {mc['max_dh']:.4f}")

        # Threshold
        thr = self._get_threshold()
        if thr is not None:
            ax.axvline(thr, color=GH["red"], linestyle="-",
                       linewidth=2.5, alpha=0.9,
                       label=f"Threshold = {thr}")

        self._style_ax(ax, f"Monte Carlo Sensitivity  (N={mc['n']})",
                       "Max dh/h\u2080 per run", "Count")
        self._legend(ax)
        self.fig.set_facecolor(PLOT_BG)
        self.fig.tight_layout()
        self.canvas.draw_idle()
        self._store_default_limits()

    # ─────────────────────────────────────────────────────────────────────
    #  PDF REPORT GENERATION
    # ─────────────────────────────────────────────────────────────────────
    def _generate_report(self):
        if not self.sweep_results:
            messagebox.showinfo("Report",
                                "No results to report. Run a simulation first.")
            return

        path = filedialog.asksaveasfilename(
            title="Save PDF Report", defaultextension=".pdf",
            filetypes=[("PDF", "*.pdf"), ("All", "*.*")])
        if not path:
            return

        try:
            with PdfPages(path) as pdf:
                # Title page
                fig_t = Figure(figsize=(11, 8.5), dpi=100, facecolor=PLOT_BG)
                ax = fig_t.add_subplot(111)
                ax.set_facecolor(PLOT_BG)
                ax.text(0.5, 0.6, "LIQLEV SIMULATION REPORT",
                        transform=ax.transAxes, ha="center", va="center",
                        fontsize=28, fontweight="bold", color=GH["blue"],
                        family="Consolas")
                ax.text(0.5, 0.45,
                        "CRYOGENIC LIQUID LEVEL RISE SIMULATOR",
                        transform=ax.transAxes, ha="center", va="center",
                        fontsize=14, color=GH["text2"],
                        family="Consolas")
                ax.text(0.5, 0.32,
                        f"GENERATED: {datetime.datetime.now():%Y-%m-%d %H:%M}",
                        transform=ax.transAxes, ha="center", va="center",
                        fontsize=10, color=GH["muted"],
                        family="Consolas")
                ax.set_xticks([])
                ax.set_yticks([])
                for sp in ax.spines.values():
                    sp.set_visible(False)
                pdf.savefig(fig_t, facecolor=PLOT_BG)
                _plt.close(fig_t)

                # Input parameters page
                fig_ip = Figure(figsize=(11, 8.5), dpi=100, facecolor=PLOT_BG)
                ax = fig_ip.add_subplot(111)
                ax.set_facecolor(PLOT_BG)
                params_lines = [
                    "INPUT PARAMETERS",
                    "=" * 50,
                    f"Fluid:              {self.fluid_var.get()}",
                    f"Initial Pressure:   {self.pinit_var.get()} psia",
                    f"Final Pressure:     {self.pfinal_var.get()} psia",
                    f"Tank Diameter:      {self.dtank_var.get()} ft",
                    f"Tank Height:        {self.htank_var.get()} ft",
                    f"Duration:           {self.duration_var.get()} s",
                    f"Time Step (dt):     {self.dt_var.get()} s",
                    f"Vent Rate(s):       {self.vent_var.get()} lbm/s",
                    f"Fill Fraction(s):   {self.fill_var.get()}",
                    f"Epsilon Mode:       {self.eps_mode_var.get()}",
                    f"Ramp Duration:      {self.ramp_dur_var.get()} s",
                    f"Ramp Factor:        {self.ramp_factor_var.get()}",
                    f"Gravity Mode:       {self.grav_mode_var.get()}",
                    "",
                    f"Units Display:      {'SI' if self._si_mode else 'Imperial'}",
                    f"Threshold dh/h0:    {self._threshold_var.get() or 'None'}",
                    "",
                    f"Xmlzro Override:    {self.xmlzro_override_var.get() or 'computed'} lbm",
                    f"Tinit Override:     {self.tinit_override_var.get() or 'computed'} \u00b0R",
                ]
                if self.grav_mode_var.get() == "Constant":
                    params_lines.append(
                        f"Gravity Level:      {self.grav_const_val_var.get()} g's")
                elif self.grav_mode_var.get() == "Function of Time":
                    params_lines.append(
                        f"g(t) =              {self.grav_func_var.get()}")
                elif self.grav_mode_var.get() == "CSV Profile":
                    params_lines.append(
                        f"CSV File:           {self.grav_csv_var.get()}")
                    params_lines.append(
                        f"Hold G:             {self.hold_g_var.get()} g's")
                ax.text(0.05, 0.95, "\n".join(params_lines),
                        transform=ax.transAxes, va="top", ha="left",
                        fontsize=9, color=GH["text"],
                        family="monospace", linespacing=1.5)
                ax.set_xticks([])
                ax.set_yticks([])
                for sp in ax.spines.values():
                    sp.set_visible(False)
                pdf.savefig(fig_ip, facecolor=PLOT_BG)
                _plt.close(fig_ip)

                # Summary page (text)
                fig_s = Figure(figsize=(11, 8.5), dpi=100, facecolor=PLOT_BG)
                ax = fig_s.add_subplot(111)
                ax.set_facecolor(PLOT_BG)
                summary = self.summary_text.get("1.0", "end").strip()
                ax.text(0.02, 0.98, summary,
                        transform=ax.transAxes, va="top", ha="left",
                        fontsize=7, color=GH["text"],
                        family="monospace", linespacing=1.4)
                ax.set_xticks([])
                ax.set_yticks([])
                for sp in ax.spines.values():
                    sp.set_visible(False)
                pdf.savefig(fig_s, facecolor=PLOT_BG)
                _plt.close(fig_s)

                # One page per scenario x plot type
                all_scenarios = list(self.sweep_results.items())
                plot_types_for_report = [
                    ("Liquid Level Rise (dh/h0)", self._plot_level_rise),
                    ("Liquid Level Rise (inches)", self._plot_level_inches),
                    ("Tank Pressure", self._plot_pressure),
                    ("Vent Rate Profile", self._plot_vent_rate),
                    ("Epsilon", self._plot_epsilon),
                    ("Diagnostics", self._plot_diagnostics),
                    ("Vapor Generation", self._plot_vapor_gen),
                ]

                for pname, pfunc in plot_types_for_report:
                    fig_p = Figure(figsize=(11, 8.5), dpi=100,
                                  facecolor=PLOT_BG)
                    # Temporarily swap self.fig
                    orig_fig = self.fig
                    self.fig = fig_p
                    try:
                        pfunc(all_scenarios)
                        fig_p.set_facecolor(PLOT_BG)
                        pdf.savefig(fig_p, facecolor=PLOT_BG)
                    finally:
                        self.fig = orig_fig
                    _plt.close(fig_p)

                # MC histogram page
                if hasattr(self, '_mc_results') and self._mc_results:
                    fig_mc = Figure(figsize=(11, 8.5), dpi=100,
                                   facecolor=PLOT_BG)
                    orig_fig = self.fig
                    self.fig = fig_mc
                    try:
                        self._plot_mc_histogram()
                        fig_mc.set_facecolor(PLOT_BG)
                        pdf.savefig(fig_mc, facecolor=PLOT_BG)
                    finally:
                        self.fig = orig_fig
                    _plt.close(fig_mc)

            messagebox.showinfo("Report Saved",
                                f"PDF report saved to:\n{path}")
        except Exception as e:
            messagebox.showerror("Report Error", str(e))

    # ─────────────────────────────────────────────────────────────────────
    #  EXPORT / CONFIG SAVE-LOAD
    # ─────────────────────────────────────────────────────────────────────
    def _export_results(self):
        if not self.sweep_results:
            messagebox.showinfo("Export",
                                "No results to export. Run a simulation first.")
            return

        folder = filedialog.askdirectory(title="Select Export Folder")
        if not folder:
            return

        count = 0
        for key, s in self.sweep_results.items():
            for df, rate in zip(s["dfs"], s["vent_rates"]):
                safe = key.replace(" ", "_").replace(",", "").replace("%", "pct")
                fname = f"liqlev_{safe}_vent{rate}.csv"
                df.to_csv(os.path.join(folder, fname), index=False)
                count += 1

        messagebox.showinfo("Export Complete",
                            f"Exported {count} CSV file(s) to:\n{folder}")

    def _export_summary_csv(self):
        """Export a one-row-per-scenario summary CSV."""
        if not self.sweep_results:
            messagebox.showinfo("Export",
                                "No results to export. Run a simulation first.")
            return
        path = filedialog.asksaveasfilename(
            title="Save Summary CSV", defaultextension=".csv",
            filetypes=[("CSV", "*.csv"), ("All", "*.*")])
        if not path:
            return

        si = self._si_mode
        thr = self._get_threshold()
        rows = []
        for key, s in self.sweep_results.items():
            htank = s["htank"]
            for df, rate in zip(s["dfs"], s["vent_rates"]):
                if df.empty:
                    continue
                mr = df.loc[df["Hratio"].idxmax()]
                h0 = df["Height"].iloc[0]
                dh_in = (mr["Height"] - h0) * 12
                exceeds = (df["Height"] >= htank).any()
                conv_fails = int(df["Conv Failed"].sum()) if "Conv Failed" in df.columns else 0
                row = {
                    "Scenario": key,
                    "Vent Rate (lbm/s)": rate,
                    "Max dh/h0": round(mr["Hratio"], 6),
                    "Max dh (in)": round(dh_in, 4),
                    "Time to Peak (s)": round(mr["Time"], 2),
                    "Final Pressure (psia)": round(df["Press"].iloc[-1], 4),
                    "Exceeds Tank": "YES" if exceeds else "No",
                    "Convergence Failures": conv_fails,
                }
                if thr is not None:
                    row["Threshold PASS/FAIL"] = (
                        "FAIL" if mr["Hratio"] > thr else "PASS")
                rows.append(row)

        pd.DataFrame(rows).to_csv(path, index=False)
        messagebox.showinfo("Summary Exported", f"Summary CSV saved to:\n{path}")

    def _clear_results(self):
        if self.sweep_results or self._mc_results:
            if not messagebox.askyesno("Confirm Clear",
                                       "Clear all results? This cannot be undone."):
                return
        self.sweep_results.clear()
        self.scenario_keys.clear()
        self._mc_results = None
        self.scenario_menu.configure(values=["(run a simulation)"])
        self.scenario_var.set("(run a simulation)")
        self._show_welcome_plot()
        self.tree.delete(*self.tree.get_children())
        self.tree["columns"] = []
        self.summary_text.configure(state="normal")
        self.summary_text.delete("1.0", "end")
        self.summary_text.configure(state="disabled")
        self.log_text.configure(state="normal")
        self.log_text.delete("1.0", "end")
        self.log_text.configure(state="disabled")
        self.progress.set(0)

    def _save_config(self):
        path = filedialog.asksaveasfilename(
            title="Save Configuration", defaultextension=".json",
            filetypes=[("JSON", "*.json"), ("All", "*.*")])
        if not path:
            return

        config = {
            "dark_mode": self._dark_mode,
            "si_mode": self._si_mode,
            "fluid": self.fluid_var.get(),
            "pinit": self.pinit_var.get(),
            "pfinal": self.pfinal_var.get(),
            "dtank": self.dtank_var.get(),
            "htank": self.htank_var.get(),
            "duration": self.duration_var.get(),
            "delta_t": self.dt_var.get(),
            "vent_rates": self.vent_var.get(),
            "fill_fractions": self.fill_var.get(),
            "eps_mode": self.eps_mode_var.get(),
            "eps_custom": self.eps_custom_var.get(),
            "ramp_duration": self.ramp_dur_var.get(),
            "ramp_factor": self.ramp_factor_var.get(),
            "grav_mode": self.grav_mode_var.get(),
            "const_grav": self.grav_const_val_var.get(),
            "grav_func_expr": self.grav_func_var.get(),
            "grav_csv": self.grav_csv_var.get(),
            "hold_g": self.hold_g_var.get(),
            "threshold": self._threshold_var.get(),
            "mc_n": self.mc_n_var.get(),
            "mc_vent_min": self.mc_vent_min_var.get(),
            "mc_vent_max": self.mc_vent_max_var.get(),
            "mc_fill_min": self.mc_fill_min_var.get(),
            "mc_fill_max": self.mc_fill_max_var.get(),
            "mc_grav_min": self.mc_grav_min_var.get(),
            "mc_grav_max": self.mc_grav_max_var.get(),
            "xmlzro_override": self.xmlzro_override_var.get(),
            "tinit_override": self.tinit_override_var.get(),
        }

        with open(path, "w") as f:
            json.dump(config, f, indent=2)
        messagebox.showinfo("Saved", f"Configuration saved to:\n{path}")

    def _load_config(self):
        path = filedialog.askopenfilename(
            title="Load Configuration",
            filetypes=[("JSON", "*.json"), ("All", "*.*")])
        if not path:
            return

        with open(path, "r") as f:
            config = json.load(f)

        field_map = {
            "fluid": self.fluid_var,
            "pinit": self.pinit_var,
            "pfinal": self.pfinal_var,
            "dtank": self.dtank_var,
            "htank": self.htank_var,
            "duration": self.duration_var,
            "delta_t": self.dt_var,
            "vent_rates": self.vent_var,
            "fill_fractions": self.fill_var,
            "eps_mode": self.eps_mode_var,
            "eps_custom": self.eps_custom_var,
            "ramp_duration": self.ramp_dur_var,
            "ramp_factor": self.ramp_factor_var,
            "grav_mode": self.grav_mode_var,
            "const_grav": self.grav_const_val_var,
            "grav_func_expr": self.grav_func_var,
            "grav_csv": self.grav_csv_var,
            "hold_g": self.hold_g_var,
            "threshold": self._threshold_var,
            "mc_n": self.mc_n_var,
            "mc_vent_min": self.mc_vent_min_var,
            "mc_vent_max": self.mc_vent_max_var,
            "mc_fill_min": self.mc_fill_min_var,
            "mc_fill_max": self.mc_fill_max_var,
            "mc_grav_min": self.mc_grav_min_var,
            "mc_grav_max": self.mc_grav_max_var,
            "xmlzro_override": self.xmlzro_override_var,
            "tinit_override": self.tinit_override_var,
        }

        # Config values are stored in whatever unit mode was active when saved.
        # First, match our current mode to the config's mode so values load
        # correctly.  Then switch to the saved mode.
        # Restore theme mode
        config_dark = config.get("dark_mode", True)
        if self._dark_mode != config_dark:
            self._toggle_theme()

        config_si = config.get("si_mode", False)
        # Temporarily switch to the config's unit mode so values load directly
        if self._si_mode != config_si:
            self._toggle_units()  # match saved mode

        for key, var in field_map.items():
            if key in config:
                var.set(str(config[key]))

        # Handle legacy configs that used the old boolean toggle
        if "use_const_grav" in config and "grav_mode" not in config:
            if config["use_const_grav"]:
                self.grav_mode_var.set("Constant")
            else:
                self.grav_mode_var.set("CSV Profile")

        self._on_gravity_mode_change(self.grav_mode_var.get())

        self._on_eps_mode_change(self.eps_mode_var.get())
        self._update_computed()
        self._update_sweep_count()
        messagebox.showinfo("Loaded", f"Configuration loaded from:\n{path}")


# ── Stdout redirector ──
class _LogRedirector:
    def __init__(self, callback):
        self._cb = callback

    def write(self, text):
        if text.strip():
            self._cb(text.rstrip())

    def flush(self):
        pass


# ═════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    app = LIQLEVApp()
    app.mainloop()
