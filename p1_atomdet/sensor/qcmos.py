"""qCMOS photon-number-resolving (PNR) output.

In PNR mode the camera reports round(n + noise) in electrons. This is the CMOS chain of
cmos.py with an ADC step of exactly 1 e-; with read noise sigma << 0.5 e- the output
reproduces the integer electron number with probability ~ 1 - 2 Phi(-0.5/sigma).
"""
from __future__ import annotations

import numpy as np
from scipy import special

from ..counts.lattice import Lattice
from .cmos import electrons_pmf, pixel_lattice


def pnr_pixel(lam_signal: float, lam_bg: float, sigma_e: float, prnu: float = 0.0) -> Lattice:
    n, pn = electrons_pmf(lam_signal, lam_bg, prnu)
    return pixel_lattice(n, pn, [(1.0, 0.0, sigma_e)], step=1.0)


def miscount_probability(sigma_e: float) -> float:
    """Probability that rounding assigns the wrong integer: P(|noise| > 1/2) = erfc(1/(2 sqrt(2) sigma))."""
    return float(special.erfc(0.5 / (np.sqrt(2) * sigma_e)))
