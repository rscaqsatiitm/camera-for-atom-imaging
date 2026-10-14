"""Vectorial (Richards-Wolf type) image of a radiating dipole through an aplanatic objective.

Object-side half-angle theta_max = arcsin(NA); the image side is paraxial (magnification >> 1),
so the image field is the Fourier transform of the pupil field. With the sine condition, the
pupil amplitude carries the factor 1/sqrt(cos theta), which with the area element gives the
integrand weight sin(theta) sqrt(cos theta). In object-plane coordinates r (DERIVED, report Sec. 2.2):

  I0(r) = int_0^tm  sin t sqrt(cos t) (1 + cos t)/2  J0(k r sin t) dt
  I2(r) = int_0^tm  sin t sqrt(cos t) (1 - cos t)/2  J2(k r sin t) dt
  I1(r) = int_0^tm  sin t sqrt(cos t)  sin t         J1(k r sin t) dt

Dipole in the focal plane (or sigma+- about z), azimuthally averaged:  |E|^2 ~ I0^2 + I2^2
Dipole along the optical axis:                                       |E|^2 ~ I1^2
Isotropic emitter (x, y, z incoherent):                                |E|^2 ~ 2(I0^2 + I2^2) + I1^2

Normalisation is exact by Parseval, so the PSF integrates to 1 over the image plane
(the collected power, not the emitted power; collection efficiency is handled separately):

  psf_inplane = k^2 (I0^2 + I2^2) / (2 pi Jx),          Jx = int (1 - sin^2 t / 2) sin t dt
  psf_axial   = k^2 I1^2 / (2 pi Jz),                    Jz = int sin^3 t dt
  psf_iso     = k^2 [2(I0^2 + I2^2) + I1^2] / (2 pi (2 Jx + Jz))
"""
from __future__ import annotations

import numpy as np
from scipy import special


def _j2(x):
    """J2 via the recurrence 2 J1(x)/x - J0(x) (much faster than jv), series for small x."""
    x = np.asarray(x, float)
    small = np.abs(x) < 1e-2
    xs = np.where(small, 1.0, x)
    out = 2 * special.j1(xs) / xs - special.j0(xs)
    return np.where(small, x**2 / 8 - x**4 / 96, out)


def _integrals(r_m: np.ndarray, wavelength_m: float, na: float, n_quad: int):
    k = 2 * np.pi / wavelength_m
    tm = np.arcsin(na)
    x, w = np.polynomial.legendre.leggauss(n_quad)
    t = 0.5 * tm * (x + 1.0)
    w = 0.5 * tm * w
    st, ct = np.sin(t), np.cos(t)
    base = w * st * np.sqrt(ct)
    arg = k * np.asarray(r_m, float)[..., None] * st
    I0 = np.sum(base * 0.5 * (1 + ct) * special.j0(arg), axis=-1)
    I2 = np.sum(base * 0.5 * (1 - ct) * _j2(arg), axis=-1)
    I1 = np.sum(base * st * special.j1(arg), axis=-1)
    Jx = np.sum(w * (1 - 0.5 * st**2) * st)
    Jz = np.sum(w * st**3)
    return k, I0, I1, I2, Jx, Jz


def psf_vectorial(r_m, wavelength_m: float, na: float, dipole: str = "isotropic", n_quad: int = 64):
    """Normalised radial PSF (1/m^2) in object-plane coordinates."""
    k, I0, I1, I2, Jx, Jz = _integrals(r_m, wavelength_m, na, n_quad)
    if dipole == "isotropic":
        return k**2 * (2 * (I0**2 + I2**2) + I1**2) / (2 * np.pi * (2 * Jx + Jz))
    if dipole in ("inplane", "sigma", "sigma_z"):
        return k**2 * (I0**2 + I2**2) / (2 * np.pi * Jx)
    if dipole in ("axial", "pi_z"):
        return k**2 * I1**2 / (2 * np.pi * Jz)
    raise ValueError(dipole)
