"""Write a JSON baseline for canonical LIQLEV physics cases."""

from __future__ import annotations

import argparse
import datetime as dt
import importlib.metadata
import json
import platform
import sys
from pathlib import Path
from typing import Dict, Iterable, Optional

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from validation.physics_cases import dataframe_summary, get_case, iter_cases, run_case  # noqa: E402


DEFAULT_OUTPUT = ROOT / "validation" / "baselines" / "physics_baseline.json"


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def selected_cases(names: Optional[Iterable[str]]):
    if names:
        return [get_case(name) for name in names]
    return list(iter_cases())


def dataframe_payload(df) -> Dict[str, object]:
    data = []
    arr = df.to_numpy(dtype=float)
    for row in arr:
        clean_row = []
        for value in row:
            if np.isnan(value):
                clean_row.append(None)
            elif np.isposinf(value):
                clean_row.append("Infinity")
            elif np.isneginf(value):
                clean_row.append("-Infinity")
            else:
                clean_row.append(float(value))
        data.append(clean_row)
    return {
        "columns": list(df.columns),
        "data": data,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Baseline JSON path. Default: {DEFAULT_OUTPUT}",
    )
    parser.add_argument(
        "--case",
        action="append",
        dest="cases",
        help="Case name to include. May be provided more than once.",
    )
    args = parser.parse_args()

    output = args.output
    if not output.is_absolute():
        output = ROOT / output
    output.parent.mkdir(parents=True, exist_ok=True)

    payload: Dict[str, object] = {
        "metadata": {
            "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": package_version("numpy"),
            "pandas": package_version("pandas"),
            "numba": package_version("numba"),
            "CoolProp": package_version("CoolProp"),
        },
        "cases": {},
    }

    for case in selected_cases(args.cases):
        print(f"Running baseline case: {case.name}")
        df = run_case(case)
        payload["cases"][case.name] = {
            "description": case.description,
            "summary": dataframe_summary(df),
            "frame": dataframe_payload(df),
        }
        summary = payload["cases"][case.name]["summary"]
        print(
            "  rows={rows} max_dh_h0={max_dh_h0:.8g} final_pressure={final_pressure_psia:.8g}".format(
                **summary
            )
        )

    with output.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, allow_nan=False)

    print(f"Wrote baseline: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

