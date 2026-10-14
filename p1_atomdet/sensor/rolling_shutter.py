"""Rolling-shutter timing relative to the illumination.

Row j starts integrating at t0 + j * t_line and integrates for t_exp. Illumination (imaging
light, hence atom signal and imaging-beam background) is on during [t_on, t_on + t_img].
The fraction of the illumination captured by row j is overlap/t_img. Dark current and
stray light that do not depend on the imaging beam scale with t_exp instead.
"""
from __future__ import annotations

import numpy as np


def row_capture_fraction(rows, t_line_s: float, t_exp_s: float, t_on_s: float, t_img_s: float,
                         t0_s: float = 0.0, shutter: str = "rolling") -> np.ndarray:
    rows = np.asarray(rows, float)
    start = t0_s + (rows * t_line_s if shutter == "rolling" else 0.0)
    end = start + t_exp_s
    ov = np.clip(np.minimum(end, t_on_s + t_img_s) - np.maximum(start, t_on_s), 0, None)
    return ov / t_img_s


def min_exposure_for_full_capture(n_rows_roi: int, t_line_s: float, t_img_s: float) -> float:
    """Shortest exposure such that every row of an ROI sees the whole pulse: t_img + (n-1) t_line."""
    return t_img_s + (n_rows_roi - 1) * t_line_s


def frame_time(n_rows_roi: int, t_line_s: float, t_exp_s: float) -> float:
    """Rolling-shutter frame time for an ROI of n rows (exposure + readout of all rows)."""
    return t_exp_s + n_rows_roi * t_line_s
