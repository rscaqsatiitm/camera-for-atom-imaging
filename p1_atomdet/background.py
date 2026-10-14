"""Background electrons per pixel per frame.

Three contributions, all Poisson and therefore additive in the mean:
  (i)  uniform stray light (and trap-light leakage through the filters), b_stray [e-/px/frame];
  (ii) imaging-beam light scattered near the site, optionally spatially shaped (b_local map);
  (iii) dark current D(T) * t_exp [e-/px].
The quantities are detected electrons (QE already applied), which is how the sweep range
0 to 1 e-/pixel/frame is defined.
"""
from __future__ import annotations

import numpy as np


def background_map(box: int, b_stray: float, t_exp_s: float, dark_e_per_s: float = 0.0,
                   b_local: np.ndarray | None = None) -> np.ndarray:
    b = np.full((box, box), float(b_stray) + dark_e_per_s * t_exp_s)
    if b_local is not None:
        b = b + np.asarray(b_local, float)
    return b
