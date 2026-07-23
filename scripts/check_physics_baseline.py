"""Compare current LIQLEV physics output against a saved JSON baseline."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Iterable, Optional

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from validation.physics_cases import dataframe_summary, get_case, run_case  # noqa: E402


DEFAULT_BASELINE = ROOT / "validation" / "baselines" / "physics_baseline.json"


def selected_case_names(baseline_cases, requested: Optional[Iterable[str]]):
    if requested:
        return list(requested)
    return list(baseline_cases.keys())


def decode_frame(frame):
    data = []
    for row in frame["data"]:
        clean = []
        for value in row:
            if value is None:
                clean.append(np.nan)
            elif value == "Infinity":
                clean.append(np.inf)
            elif value == "-Infinity":
                clean.append(-np.inf)
            else:
                clean.append(float(value))
        data.append(clean)
    return frame["columns"], np.array(data, dtype=float)


def compare_arrays(expected, actual, rtol, atol):
    if expected.shape != actual.shape:
        return False, f"shape mismatch expected {expected.shape}, got {actual.shape}"

    close = np.isclose(actual, expected, rtol=rtol, atol=atol, equal_nan=True)
    if np.all(close):
        return True, "ok"

    bad = np.argwhere(~close)
    diffs = np.abs(actual - expected)
    finite = np.isfinite(diffs)
    max_abs = float(np.nanmax(diffs[finite])) if np.any(finite) else float("nan")
    row, col = bad[0]
    return (
        False,
        f"{len(bad)} values differ; max_abs={max_abs:.6g}; first difference at row={row}, col={col}, expected={expected[row, col]:.12g}, actual={actual[row, col]:.12g}",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--baseline",
        type=Path,
        default=DEFAULT_BASELINE,
        help=f"Baseline JSON path. Default: {DEFAULT_BASELINE}",
    )
    parser.add_argument("--case", action="append", dest="cases", help="Case name to check.")
    parser.add_argument("--rtol", type=float, default=1e-9, help="Relative tolerance.")
    parser.add_argument("--atol", type=float, default=1e-8, help="Absolute tolerance.")
    args = parser.parse_args()

    baseline_path = args.baseline
    if not baseline_path.is_absolute():
        baseline_path = ROOT / baseline_path

    if not baseline_path.exists():
        print(f"Baseline not found: {baseline_path}")
        print("Create it first with: python scripts/write_physics_baseline.py")
        return 2

    with baseline_path.open("r", encoding="utf-8") as f:
        baseline = json.load(f)

    failures = []
    for name in selected_case_names(baseline["cases"], args.cases):
        print(f"Checking baseline case: {name}")
        if name not in baseline["cases"]:
            failures.append(f"{name}: not present in baseline")
            print("  FAIL case is not present in the selected baseline")
            continue

        case = get_case(name)
        expected_columns, expected_data = decode_frame(baseline["cases"][name]["frame"])

        df = run_case(case)
        actual_columns = list(df.columns)
        actual_data = df.to_numpy(dtype=float)

        if actual_columns != expected_columns:
            failures.append(f"{name}: columns changed")
            print(f"  FAIL columns changed\n  expected={expected_columns}\n  actual={actual_columns}")
            continue

        ok, message = compare_arrays(expected_data, actual_data, args.rtol, args.atol)
        if not ok:
            failures.append(f"{name}: {message}")
            print(f"  FAIL {message}")
            print(f"  current summary: {dataframe_summary(df)}")
            print(f"  baseline summary: {baseline['cases'][name]['summary']}")
        else:
            summary = dataframe_summary(df)
            print(
                "  PASS rows={rows} max_dh_h0={max_dh_h0:.8g} final_pressure={final_pressure_psia:.8g}".format(
                    **summary
                )
            )

    if failures:
        print("\nPhysics baseline check failed:")
        for failure in failures:
            print(f" - {failure}")
        return 1

    print("\nPhysics baseline check passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
