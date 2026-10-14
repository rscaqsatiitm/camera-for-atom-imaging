"""Scalar PSFs: Airy (uniform pupil) and Gaussian (empirical, e.g. a measured waist)."""
from __future__ import annotations

import numpy as np
from scipy import optimize, special


def psf_airy(r_m, wavelength_m: float, na: float):
    """Normalised Airy PSF (1/m^2): (k NA)^2/(2 pi) [2 J1(v)/v]^2, v = k NA r."""
    k = 2 * np.pi / wavelength_m
    v = k * na * np.asarray(r_m, float)
    with np.errstate(invalid="ignore", divide="ignore"):
        a = np.where(v == 0, 1.0, 2 * special.j1(v) / np.where(v == 0, 1.0, v))
    return (k * na) ** 2 / (2 * np.pi) * a**2


def psf_gaussian(r_m, waist_m: float):
    """Normalised Gaussian PSF with 1/e^2 intensity radius `waist_m`: 2/(pi w^2) exp(-2 r^2/w^2)."""
    r = np.asarray(r_m, float)
    return 2 / (np.pi * waist_m**2) * np.exp(-2 * r**2 / waist_m**2)


def airy_fwhm(wavelength_m: float, na: float) -> float:
    """FWHM of the Airy disc: 2 v_half/(k NA) with [2J1(v)/v]^2 = 1/2 at v_half = 1.6163..."""
    vh = optimize.brentq(lambda v: (2 * special.j1(v) / v) ** 2 - 0.5, 0.5, 3.0)
    return 2 * vh / (2 * np.pi / wavelength_m * na)


def radial_fwhm(psf_fn, r_max: float, n: int = 2001) -> float:
    """FWHM of a monotonically decreasing (near the core) radial PSF by root finding."""
    r = np.linspace(0, r_max, n)
    p = psf_fn(r)
    half = 0.5 * p[0]
    i = int(np.argmax(p < half))
    return 2 * optimize.brentq(lambda x: psf_fn(np.array([x]))[0] - half, r[i - 1], r[i])


def encircled_energy(psf_fn, radius: float, n: int = 4001) -> float:
    """EE(R) = int_0^R psf 2 pi r dr (Gauss-Legendre in r)."""
    x, w = np.polynomial.legendre.leggauss(n)
    r = 0.5 * radius * (x + 1)
    return float(np.sum(0.5 * radius * w * psf_fn(r) * 2 * np.pi * r))
