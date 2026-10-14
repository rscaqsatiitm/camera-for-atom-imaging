"""Pixel-integrated PSF, sub-pixel site placement and neighbour crosstalk.

Coordinates are object-plane metres; a camera pixel of pitch p maps to p/M in the object
plane (M = magnification). Pixel fractions are computed on a fine grid with s x s sub-samples
per pixel (midpoint rule, error O(1/s^2), convergence tested), optionally convolved with a
Gaussian blur (atom thermal motion, residual aberration, charge spread).
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import partial

import numpy as np
from scipy import signal

from .psf_airy import airy_fwhm, psf_airy, psf_gaussian, radial_fwhm
from .psf_vectorial import psf_vectorial


@dataclass
class PSFSpec:
    kind: str = "vectorial"          # vectorial | airy | gaussian
    wavelength_m: float = 852.347e-9
    na: float = 0.5
    dipole: str = "isotropic"
    waist_m: float | None = None     # for kind == gaussian (1/e^2 intensity radius, object plane)
    blur_sigma_m: float = 0.0        # extra Gaussian blur (object plane, rms per axis)
    n_quad: int = 64

    def radial(self):
        if self.kind == "vectorial":
            return partial(psf_vectorial, wavelength_m=self.wavelength_m, na=self.na, dipole=self.dipole,
                           n_quad=self.n_quad)
        if self.kind == "airy":
            return partial(psf_airy, wavelength_m=self.wavelength_m, na=self.na)
        if self.kind == "gaussian":
            return partial(psf_gaussian, waist_m=self.waist_m)
        raise ValueError(self.kind)

    def fwhm(self) -> float:
        """FWHM of the unblurred PSF (object plane)."""
        if self.kind == "airy":
            return airy_fwhm(self.wavelength_m, self.na)
        if self.kind == "gaussian":
            return self.waist_m * np.sqrt(2 * np.log(2))
        return radial_fwhm(self.radial(), r_max=2 * self.wavelength_m / self.na)


def pitch_for_pixels_per_fwhm(psf: PSFSpec, pixels_per_fwhm: float) -> float:
    """Object-plane pixel size giving the requested sampling of the PSF FWHM."""
    return psf.fwhm() / pixels_per_fwhm


def pixel_fractions(psf: PSFSpec, pitch_m: float, box: int = 7, offset_px=(0.0, 0.0), supersample: int = 8,
                    margin_px: int | None = None) -> np.ndarray:
    """Fraction of collected photons landing in each pixel of a box x box region.

    The box is centred on the pixel nearest the nominal site; `offset_px` = (dx, dy) is the
    site position relative to that pixel centre, in pixels (sub-pixel placement in [-0.5, 0.5)).
    Rows index y, columns index x.
    """
    s = int(supersample)
    blur_px = psf.blur_sigma_m / pitch_m
    if margin_px is None:
        margin_px = int(np.ceil(4 * blur_px)) + 1 if blur_px > 0 else 0
    n = box + 2 * margin_px
    h = 1.0 / s
    c = (np.arange(n * s) + 0.5) * h - n / 2.0  # sub-sample centres in pixel units, box centred at 0
    X, Y = np.meshgrid(c - offset_px[0], c - offset_px[1])
    r = np.hypot(X, Y) * pitch_m
    img = psf.radial()(r) * (h * pitch_m) ** 2  # probability per sub-cell
    if blur_px > 0:
        k = int(np.ceil(5 * blur_px * s))
        g = np.exp(-0.5 * ((np.arange(-k, k + 1) * h) / blur_px) ** 2)
        g /= g.sum()
        img = signal.fftconvolve(signal.fftconvolve(img, g[None, :], mode="same"), g[:, None], mode="same")
        img = np.clip(img, 0, None)
    px = img.reshape(n, s, n, s).sum(axis=(1, 3))
    if margin_px:
        px = px[margin_px:-margin_px, margin_px:-margin_px]
    return px


def neighbour_fractions(psf: PSFSpec, pitch_m: float, spacing_m: float, box: int = 7, offset_px=(0.0, 0.0),
                        supersample: int = 8, n_shell: int = 1) -> dict:
    """Pixel fractions in the central box due to atoms on a square lattice of the given spacing.

    Returns {(i, j): fractions} for lattice neighbours with max(|i|,|j|) <= n_shell, (0,0) excluded.
    """
    out = {}
    a = spacing_m / pitch_m
    for i in range(-n_shell, n_shell + 1):
        for j in range(-n_shell, n_shell + 1):
            if i == 0 and j == 0:
                continue
            out[(i, j)] = pixel_fractions(psf, pitch_m, box, (offset_px[0] + i * a, offset_px[1] + j * a), supersample)
    return out
