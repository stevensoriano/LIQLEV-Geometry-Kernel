# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec for the modern (PySide6) LIQLEV Cryovent Analysis Console.

Build with:
    pyinstaller LIQLEV-Modern.spec --clean --noconfirm

Outputs:
    dist/LIQLEV/LIQLEV.exe       — Windows executable
    dist/LIQLEV/<deps>           — bundled DLLs, Qt plugins, data files

Distribute to teammates by zipping the entire `dist/LIQLEV/` folder.
Recipients unzip and double-click LIQLEV.exe — no Python install required.

Note: The legacy CustomTkinter app spec is in `LIQLEV.spec`; this file is
the equivalent for the modern PySide6 entry point at `liqlev/ui_qt/app.py`.
"""
import sys
from pathlib import Path

from PyInstaller.utils.hooks import collect_all

# If the build env happens to be (or inherit from) a conda Python, its
# stdlib .pyd modules and PySide6 DLLs dynamically link against shared
# libraries under <base_prefix>/Library/bin/ instead of being self-contained.
# Resolve via base_prefix so that pip venvs created from conda Pythons
# still pick up the shared libs they actually need at runtime.
ENV_ROOT = Path(sys.base_prefix)
CONDA_BIN = ENV_ROOT / "Library" / "bin"
CONDA_PLUGINS = ENV_ROOT / "Library" / "plugins"

# --------------------------------------------------------------------------
# Data files and source modules to bundle
# --------------------------------------------------------------------------
datas = [
    # Example gravity / vent CSV profiles teammates may need to load.
    ('data', 'data'),
    # Legacy physics modules at repo root — the modern package imports these
    # directly via `from core import ...` and `from thermo_utils import ...`.
    ('core.py', '.'),
    ('thermo_utils.py', '.'),
]

binaries = []
hiddenimports = []

# --------------------------------------------------------------------------
# Conda PySide6 native libs — Qt6 DLLs go *next to* the PySide6 .pyd files
# so Windows' default DLL search (same directory as the loading module)
# finds them. The shiboken6 DLL similarly sits next to its .pyd in
# shiboken6/. ICU and other shared deps go to the bundle root which is
# always on the DLL search path for OneDir bundles.
# --------------------------------------------------------------------------
if CONDA_BIN.is_dir():
    for dll in CONDA_BIN.glob("Qt6*.dll"):
        binaries.append((str(dll), "PySide6"))
    for dll in CONDA_BIN.glob("shiboken6*.dll"):
        binaries.append((str(dll), "shiboken6"))
        # Also drop one copy into PySide6/ so QtCore.pyd's load of
        # shiboken6 succeeds via the same-directory search.
        binaries.append((str(dll), "PySide6"))
    # Shared deps Qt itself pulls in. Bundle to root so all .pyds find them.
    for pattern in ("icudt*.dll", "icuin*.dll", "icuuc*.dll",
                    "libcrypto*.dll", "libssl*.dll",
                    "libEGL.dll", "libGLESv2.dll", "opengl32sw.dll",
                    "d3dcompiler_*.dll"):
        for dll in CONDA_BIN.glob(pattern):
            binaries.append((str(dll), "."))
    # Conda Python stdlib .pyd files dynamically link to a handful of
    # shared libs (libexpat, libffi, sqlite3, zlib, bz2, lzma, etc.) that
    # also live under Library/bin. PyInstaller's stdlib hooks expect these
    # statically linked into the .pyd, so on a conda Python we need to
    # carry them along explicitly.
    for pattern in (
        "libexpat*.dll", "libffi*.dll", "ffi*.dll",
        "sqlite3*.dll", "libsqlite*.dll",
        "zlib*.dll", "libzlib*.dll", "libz*.dll",
        "bz2*.dll", "libbz2*.dll",
        "lzma*.dll", "liblzma*.dll",
        "libpng*.dll", "libtiff*.dll", "libjpeg*.dll",
        "libxml*.dll", "libxslt*.dll",
        "libsharpyuv*.dll", "libwebp*.dll",
        "tcl*.dll", "tk*.dll",  # harmless if unused
        "libssh2*.dll", "iconv*.dll",
        "vcruntime*.dll", "msvcp*.dll", "concrt*.dll",
        "libintl*.dll", "libuv*.dll",
        "harfbuzz*.dll", "freetype*.dll", "graphite2*.dll",
        "libglib*.dll", "pcre*.dll",
    ):
        for dll in CONDA_BIN.glob(pattern):
            binaries.append((str(dll), "."))

# --------------------------------------------------------------------------
# Conda PySide6 plugins — Qt needs the platforms plugin (qwindows.dll) at
# minimum or the application can't open a window. Stash everything under
# `PySide6/plugins/` and let the entry script point QT_PLUGIN_PATH at it.
# --------------------------------------------------------------------------
if CONDA_PLUGINS.is_dir():
    for sub in ("platforms", "imageformats", "iconengines", "styles"):
        src_dir = CONDA_PLUGINS / sub
        if src_dir.is_dir():
            for dll in src_dir.glob("*.dll"):
                binaries.append((str(dll), f"PySide6/plugins/{sub}"))

# --------------------------------------------------------------------------
# Hidden imports — modules PyInstaller's static analysis tends to miss
# --------------------------------------------------------------------------
hiddenimports += [
    # The modern package itself — picked up via dependency walk, but listed
    # explicitly so submodules used only via string-named load paths still
    # make it in.
    'liqlev',
    'liqlev.io',
    'liqlev.io.config_json',
    'liqlev.io.profiles',
    'liqlev.model',
    'liqlev.model.builder',
    'liqlev.model.config',
    'liqlev.model.parsing',
    'liqlev.model.units',
    'liqlev.model.validation',
    'liqlev.runner',
    'liqlev.runner.monte_carlo',
    'liqlev.runner.progress',
    'liqlev.runner.single',
    'liqlev.runner.sweep',
    'liqlev.ui_qt',
    'liqlev.ui_qt.app',
    'liqlev.viz',
    'liqlev.viz.datasets',
    'liqlev.viz.summaries',

    # Numba / llvmlite — JIT-compiled solver loop.
    'numba',
    'numba.core',
    'numba.core.types',
    'numba.core.typing',
    'numba.cpython',
    'numba.np',
    'llvmlite',
    'llvmlite.binding',

    # CoolProp — real-fluid thermodynamic properties.
    'CoolProp',
    'CoolProp.CoolProp',

    # matplotlib backends used by the PDF report export.
    'matplotlib.backends.backend_pdf',
    'matplotlib.backends.backend_agg',

    # pyqtgraph — plotting widget.
    'pyqtgraph',

    # PySide6 Qt modules referenced.
    'PySide6.QtCore',
    'PySide6.QtGui',
    'PySide6.QtWidgets',
]

# --------------------------------------------------------------------------
# Native libs and data via collect_all — these packages ship DLLs or
# resource files PyInstaller needs to copy verbatim into the bundle.
# --------------------------------------------------------------------------
for pkg in ('CoolProp', 'numba', 'llvmlite', 'pyqtgraph'):
    pkg_datas, pkg_binaries, pkg_hidden = collect_all(pkg)
    datas += pkg_datas
    binaries += pkg_binaries
    hiddenimports += pkg_hidden

# --------------------------------------------------------------------------
# Analysis / build
# --------------------------------------------------------------------------
a = Analysis(
    [str(Path('liqlev') / 'ui_qt' / 'app.py')],
    pathex=[str(Path.cwd())],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        # Legacy GUI stack — not used by the modern app. Excluding shrinks
        # the bundle and removes a chunk of unused Tcl/Tk DLLs.
        'customtkinter',
        'tkinter',
        '_tkinter',
        'matplotlib.backends.backend_tkagg',
        # We never run a Jupyter notebook from the bundle.
        'IPython',
        'jupyter',
        'notebook',
        # Alternative Qt bindings — PyInstaller refuses to bundle more than
        # one Qt family. Spyder ships PyQt5 in the same env; we only need
        # PySide6 for the modern app, so exclude the rest.
        'PyQt5',
        'PyQt6',
        'PySide2',
        'shiboken2',
        # Numba's test and CUDA sub-trees are huge and unused at runtime —
        # excluding shaves hundreds of MB and avoids needing a CUDA toolkit.
        'numba.tests',
        'numba.cuda',
        'numba.cuda.tests',
        'numba.roc',
    ],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='LIQLEV',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,                # UPX not assumed installed; skip compression
    console=False,            # GUI app — no console window
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    # icon='liqlev.ico',      # add when an icon file is available
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name='LIQLEV',
)
