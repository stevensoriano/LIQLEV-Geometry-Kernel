"""CSV profile loading and validation helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from liqlev.model.builder import G_TO_FT_S2


@dataclass(frozen=True)
class ProfileData:
    time_s: np.ndarray
    values: np.ndarray
    source: str

    @property
    def point_count(self) -> int:
        return int(len(self.time_s))


def load_gravity_profile_csv(
    path: str | Path, duration_s: float, hold_g: float
) -> ProfileData:
    """Load the legacy gravity CSV profile and convert g units to ft/s^2."""
    source = Path(path)
    data = pd.read_csv(source)
    missing = {"normalized_time", "az_positive"} - set(data.columns)
    if missing:
        raise ValueError(
            f"Gravity CSV missing required columns: {', '.join(sorted(missing))}"
        )

    time_s = data["normalized_time"].to_numpy(dtype=float)
    gravity_g = data["az_positive"].to_numpy(dtype=float)
    if len(time_s) == 0:
        raise ValueError("Gravity CSV contains no rows.")

    last_t = float(time_s[-1])
    if duration_s > last_t:
        time_s = np.append(time_s, [last_t + 1e-9, duration_s])
        gravity_g = np.append(gravity_g, [hold_g, hold_g])

    return ProfileData(time_s=time_s, values=gravity_g * G_TO_FT_S2, source=str(source))


def load_vent_rate_profile_csv(path: str | Path) -> ProfileData:
    """Load a vent-rate CSV using the first two columns, preserving legacy behavior."""
    source = Path(path)
    data = pd.read_csv(source)
    if data.shape[1] < 2:
        raise ValueError("Vent rate CSV must contain at least two columns.")
    if data.empty:
        raise ValueError("Vent rate CSV contains no rows.")

    return ProfileData(
        time_s=data.iloc[:, 0].to_numpy(dtype=float),
        values=data.iloc[:, 1].to_numpy(dtype=float),
        source=str(source),
    )
