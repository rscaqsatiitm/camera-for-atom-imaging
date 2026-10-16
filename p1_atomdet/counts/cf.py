"""Characteristic-function (CF) engine for linear statistics S = sum_i w_i x_i.

For linear classifiers (box sum, matched filter) the CF of S factorises over pixels and the
atom signal enters only through s = tau * m, so

  phi_S(k | s) = exp( A(k) + mu s B(k) ) N(k)
    A(k) = sum_i b_i [psi_i(k) - 1],    B(k) = sum_i f_i [psi_i(k) - 1],
    psi_i(k) = E exp(i k w_i X_1)  (one detected electron through the gain stage),
    N(k) = prod_i phi_noise,i(k w_i)  (read noise, independent of the signal),

and the atom-present law is the mixture sum_s W_s phi_S(k | s). One inverse FFT then gives
the distribution on a lattice of step h. psi = exp(i k w) for CMOS/sCMOS/qCMOS (no gain
stage); for an EMCCD in analog mode the electron is replaced by a gamma-distributed register
output, psi = (1 - i k w theta)^-a with a = 1/(F^2 - 1), theta = F^2 - 1 (Erlang for F^2 = 2).

Scope and accuracy: no ADC rounding, no PRNU, ensemble read noise (fixed pixel maps need the
lattice engine), serial CIC neglected. The lattice values are h * density at the lattice
points (midpoint rule); h is chosen as sigma_S / n_per_sigma, sigma_S the total read-noise
width of S, which bounds the discretisation error of tail sums by ~ (z/n_per_sigma)^2/24
relative at z standard deviations (0.5% at z = 5 for the default 30). Tests compare with the
lattice engine. Used for sweeps, where it is ~10^3 times faster.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import fft as sfft
from scipy import special, stats

from .lattice import Lattice
from .pmf import SiteModel


@dataclass
class CFNoise:
    comps: list                 # ensemble Gaussian mixture [(w, mean_e, sigma_e)], input-referred electrons
    gain_kind: str = "none"     # none | gamma (EMCCD analog)
    f2: float = 2.0             # excess noise factor squared for gain_kind == gamma

    @classmethod
    def from_classes(cls, classes, **kw):
        comps = [(cw * w, m, s) for cw, cl in classes for w, m, s in cl]
        return cls(comps=comps, **kw)

    @property
    def var(self) -> float:
        m1 = sum(w * m for w, m, _ in self.comps)
        return sum(w * (m * m + s * s) for w, m, s in self.comps) - m1 * m1

    @property
    def max_abs(self) -> float:
        return max(abs(m) + 9 * s for _, m, s in self.comps)

    def psi(self, kw):
        if self.gain_kind == "none":
            return np.exp(1j * kw)
        a = 1.0 / (self.f2 - 1.0)
        th = self.f2 - 1.0
        return (1.0 - 1j * kw * th) ** (-a)

    def phi_noise(self, kw):
        out = np.zeros_like(kw, dtype=complex)
        for w, m, s in self.comps:
            out += w * np.exp(1j * kw * m - 0.5 * (kw * s) ** 2)
        return out


def _support(site: SiteModel, w, noise: CFNoise, s_max: float, nsd: float = 14.0):
    f, b = site.frac.ravel(), site.bg.ravel()
    jump_var = noise.f2 if noise.gain_kind == "gamma" else 1.0  # E[X^2] per electron
    lam_hi = b + site.mu * s_max * f
    mean_hi = np.sum(w * lam_hi)
    var_hi = np.sum(w**2 * (lam_hi * jump_var + noise.var))
    mean_lo = np.sum(w * b)
    var_lo = np.sum(w**2 * (b * jump_var + noise.var))
    lo = min(mean_lo, 0.0) - nsd * np.sqrt(var_lo) - np.sum(np.abs(w)) * noise.max_abs / np.sqrt(max(np.count_nonzero(w), 1)) - 1.0
    hi = mean_hi + nsd * np.sqrt(var_hi) + 30.0 * np.max(np.abs(w)) + 5.0
    return lo, hi


def cf_statistic(site: SiteModel, noise: CFNoise, weights=None, present: bool = True, h: float | None = None,
                 n_per_sigma: float = 30.0, lo_hi=None) -> Lattice:
    """Distribution of S on a lattice (x0 + j h); weights None means a full-box sum."""
    w = np.ones(site.frac.size) if weights is None else np.asarray(weights, float).ravel()
    f, b = site.frac.ravel(), site.bg.ravel()
    if site.neighbours or site.class_map is not None:
        raise NotImplementedError("CF engine: single site, ensemble read noise")
    if present:
        svals, sw = site.signal_multiplier_nodes()
    else:
        svals, sw = np.array([0.0]), np.array([1.0])
    sig_S = np.sqrt(np.sum(w**2) * noise.var)
    if h is None:
        if sig_S > 0:
            h = sig_S / n_per_sigma
        elif noise.gain_kind == "none" and np.allclose(w, np.round(w)):
            h = 1.0
        else:
            raise ValueError("zero read noise with non-integer weights: use the lattice engine")
    lo, hi = lo_hi if lo_hi is not None else _support(site, w, noise, float(np.max(svals)))
    lo = np.floor(lo / h) * h
    N = sfft.next_fast_len(int(np.ceil((hi - lo) / h)) + 1, real=True)
    m = np.arange(N // 2 + 1)
    k = 2 * np.pi * m / (N * h)
    kw = k[None, :] * w[:, None]
    psim1 = noise.psi(kw) - 1.0
    A = b @ psim1
    B = f @ psim1
    Nk = np.prod(noise.phi_noise(kw), axis=0)
    mixA = np.zeros(k.size, complex)
    for s, a in zip(svals, sw):
        mixA += a * np.exp(A + site.mu * s * B)
    phi = mixA * Nk * np.exp(-1j * k * lo)
    p = sfft.irfft(np.conj(phi), n=N)
    p = np.clip(p, 0.0, None)
    p[p < 1e-16 * p.max()] = 0.0
    return Lattice(p, lo, h)


def cf_pair(site: SiteModel, noise: CFNoise, weights=None, n_per_sigma: float = 30.0):
    """(P0, P1) on a common lattice, ready for fidelity.optimal_threshold."""
    w = np.ones(site.frac.size) if weights is None else np.asarray(weights, float).ravel()
    svals, _ = site.signal_multiplier_nodes()
    lo_hi = _support(site, w, noise, float(np.max(svals)))
    sig_S = np.sqrt(np.sum(w**2) * noise.var)
    if noise.gain_kind == "gamma":  # the register output is itself continuous: resolve its width
        sig_S = np.sqrt(np.sum(w**2 * (noise.var + site.bg.ravel() * noise.f2)))
    h = sig_S / n_per_sigma if sig_S > 0 else 1.0
    p0 = cf_statistic(site, noise, w, False, h, lo_hi=lo_hi)
    p1 = cf_statistic(site, noise, w, True, h, lo_hi=lo_hi)
    return p0, p1


def matched_weights_analytic(site: SiteModel, noise: CFNoise) -> np.ndarray:
    """Fisher weights w_i = delta_i / ((v0_i + v1_i)/2) from first and second moments."""
    f, b = site.frac.ravel(), site.bg.ravel()
    jm = noise.f2 if noise.gain_kind == "gamma" else 1.0
    d = site.mu * f
    v0 = b * jm + noise.var
    v1 = (b + d) * jm + noise.var
    w = d / (0.5 * (v0 + v1))
    return (w / w.max()).reshape(site.frac.shape)


def box_mask_for_snr(site: SiteModel, noise: CFNoise) -> np.ndarray:
    """Box classifier mask: pixels sorted by signal fraction, keep the prefix maximising the
    Gaussian discriminability d' = mu sum f / sqrt(var0 + var1)/sqrt(2) (cheap, moment-based)."""
    f, b = site.frac.ravel(), site.bg.ravel()
    jm = noise.f2 if noise.gain_kind == "gamma" else 1.0
    o = np.argsort(-f)
    cs = np.cumsum(f[o]) * site.mu
    v0 = np.cumsum(b[o] * jm + noise.var)
    v1 = v0 + cs * jm
    d = cs / np.sqrt(0.5 * (v0 + v1))
    n = int(np.argmax(d)) + 1
    mask = np.zeros(f.size)
    mask[o[:n]] = 1.0
    return mask.reshape(site.frac.shape)
