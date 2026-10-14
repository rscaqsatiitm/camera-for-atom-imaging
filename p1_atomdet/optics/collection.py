"""Collection fraction of dipole emission into a cone of numerical aperture NA (in vacuum).

The objective axis is z. Let c = cos(theta_max), sin(theta_max) = NA. Normalised emission
patterns and their integrals over the cone (DERIVED, report Sec. 2.2):

  pi dipole along z:        (3/8pi) sin^2(theta)        ->  1/2 - (3/4) c + (1/4) c^3
  sigma+- about z, or a
  linear dipole in the
  focal plane (phi-avg):    (3/16pi)(1 + cos^2 theta)   ->  (3/8)(1 - c) + (1/8)(1 - c^3)
  isotropic:                1/4pi                        ->  (1 - c)/2

All three tend to 1/2 at NA = 1 and to zero as NA^2 (isotropic, sigma) or NA^4 (pi along z).
"""
from __future__ import annotations

import numpy as np


def _c(na):
    na = np.asarray(na, float)
    if np.any((na < 0) | (na > 1)):
        raise ValueError("NA must be in [0, 1] (vacuum)")
    return np.sqrt(1.0 - na**2)


def collection_fraction(na, pattern: str = "isotropic"):
    c = _c(na)
    if pattern == "isotropic":
        return 0.5 * (1.0 - c)
    if pattern in ("sigma", "sigma_z", "pi_perp"):
        return 3.0 / 8.0 * (1.0 - c) + 1.0 / 8.0 * (1.0 - c**3)
    if pattern in ("pi", "pi_z"):
        return 0.5 - 0.75 * c + 0.25 * c**3
    raise ValueError(pattern)


def collection_fraction_numeric(na, pattern: str = "isotropic", n: int = 4000) -> float:
    """Direct quadrature of the emission pattern, used only to test the closed forms."""
    th = np.linspace(0, np.arcsin(na), n)
    if pattern == "isotropic":
        f = np.full_like(th, 1 / (4 * np.pi))
    elif pattern in ("sigma", "sigma_z", "pi_perp"):
        f = 3 / (16 * np.pi) * (1 + np.cos(th) ** 2)
    else:
        f = 3 / (8 * np.pi) * np.sin(th) ** 2
    return float(np.trapezoid(f * 2 * np.pi * np.sin(th), th))
