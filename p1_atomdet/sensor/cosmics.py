"""Cosmic-ray (muon) hits on a site's pixel box.

Hits on a box of area A during an exposure t_eff are Poisson with mean rate * A * t_eff.
A hit deposits dE/dx * path electrons; for a sensor of thickness d and zenith angle th the
path is d/cos(th) and the track crosses ~ d tan(th)/pitch pixels. For atom-detection
thresholds of tens of electrons any hit on the box is effectively a false-positive (or a
masked) event, so we model a hit as an additive deposit drawn from a broad distribution.
"""
from __future__ import annotations

import numpy as np

from ..counts.lattice import Lattice, mixture


def hit_probability(rate_per_cm2_s: float, box_area_cm2: float, t_eff_s: float) -> float:
    return float(1.0 - np.exp(-rate_per_cm2_s * box_area_cm2 * t_eff_s))


def add_cosmics(stat: Lattice, p_hit: float, deposit_e: float, weight: float = 1.0) -> Lattice:
    """Mixture: with probability p_hit the statistic is shifted by weight*deposit (single-pixel deposit)."""
    if p_hit <= 0:
        return stat
    shift = int(round(weight * deposit_e / stat.h))
    shifted = Lattice(stat.p.copy(), stat.x0 + shift * stat.h, stat.h)
    return mixture([stat, shifted], [1.0 - p_hit, p_hit])
