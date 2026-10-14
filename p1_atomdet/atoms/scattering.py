"""Two-level photon scattering rate."""
from __future__ import annotations

import numpy as np


def scattering_rate(gamma_s: float, s: float | np.ndarray, delta_s: float | np.ndarray = 0.0):
    """R = (Gamma/2) s / (1 + s + 4 Delta^2/Gamma^2)   [photons/s].

    Here, gamma_s is the decay rate Gamma (s^-1), s = I/I_sat the saturation parameter and
    delta_s the detuning Delta in angular units (s^-1), same units as Gamma.
    """
    s = np.asarray(s, float)
    d = np.asarray(delta_s, float)
    return 0.5 * gamma_s * s / (1.0 + s + 4.0 * d**2 / gamma_s**2)


def saturation_for_rate(gamma_s: float, rate_s: float, delta_s: float = 0.0) -> float:
    """Invert scattering_rate for s at fixed detuning (used when only R is known from a paper)."""
    x = 2.0 * rate_s / gamma_s
    if x >= 1.0:
        raise ValueError("rate exceeds Gamma/2")
    return x * (1.0 + 4.0 * delta_s**2 / gamma_s**2) / (1.0 - x)
