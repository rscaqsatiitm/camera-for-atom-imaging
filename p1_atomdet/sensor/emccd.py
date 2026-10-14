"""EMCCD output model: stochastic multiplication register with clock-induced charge (CIC).

References (retrieved and read 2026-10-14):
  Basden, Haniff & Mackay, MNRAS 345, 985 (2003), arXiv:astro-ph/0307305v1:
    Eq. (1)/(A5): output for n input electrons ~ x^(n-1) exp(-x/g)/(g^n (n-1)!)  (Erlang, large g, large r)
    Eq. (5): E^2 = 1/g + 2/g^(1/r) - 2/g^(1+1/r) for a register of r elements (Matsuo et al. 1985)
    Appendix A: stage recursion with per-stage multiplication probability P, mean gain g = (1+P)^r.

All semi-analytic distributions here are in INPUT-referred electron units (output / mean gain).
`register_model`:
  'erlang'    Gamma(shape n, scale 1): Basden Eq. (A5), excess noise factor F^2 = 2 exactly.
  'stagewise' Gamma moment-matched to the exact finite-register mean and variance (Basden Eq. 5);
              validated against the brute-force stage-by-stage Monte Carlo `register_mc`.
CIC: a fraction (1 - f_reg) is generated before the register (behaves as an input electron);
a fraction f_reg is generated at a stage j uniform in 1..r and sees gain (1+P)^(r-j) only.
"""
from __future__ import annotations

from dataclasses import dataclass

import numba
import numpy as np
from scipy import special, stats

from ..counts.lattice import Lattice, convolve, mixture


def enf2_finite(g: float, r: int) -> float:
    """Excess noise factor squared, Basden et al. (2003) Eq. (5); -> 2 - 1/g as P -> 0."""
    return 1.0 / g + 2.0 / g ** (1.0 / r) - 2.0 / g ** (1.0 + 1.0 / r)


@dataclass
class EMCCDReadout:
    em_gain: float
    n_stages: int
    cic: float                 # events/pixel/frame
    cic_register_fraction: float
    sigma_ro_e: float          # read noise of the output amplifier, output-referred electrons rms
    register_model: str = "erlang"
    mode: str = "analog"       # analog | photon_counting
    pc_threshold_sigma: float = 5.0
    h_fine: float = 0.02       # lattice step, input-referred electrons

    @property
    def f2(self) -> float:
        return 2.0 if self.register_model == "erlang" else enf2_finite(self.em_gain, self.n_stages)


def _gamma_bins(shape: float, scale: float, edges: np.ndarray) -> np.ndarray:
    return np.diff(special.gammainc(shape, np.clip(edges, 0, None) / scale))


def register_output_lattice(pn: np.ndarray, n: np.ndarray, ro: EMCCDReadout, x_max: float) -> Lattice:
    """Input-referred register output for an input electron PMF pn(n), without CIC and read noise."""
    h = ro.h_fine
    nb = int(np.ceil(x_max / h)) + 2
    edges = (np.arange(nb + 1) - 0.5) * h
    out = np.zeros(nb)
    out[0] += pn[n == 0].sum()
    k = ro.f2 - 1.0
    for ni, pi in zip(n, pn):
        if ni == 0 or pi < 1e-300:
            continue
        out += pi * _gamma_bins(ni / k, k, edges)
    return Lattice(out, 0.0, h)


def serial_cic_lattice(ro: EMCCDReadout, x_max: float) -> Lattice:
    """Input-referred output of register-generated CIC: Bernoulli(c_reg) x Exp(mean m^-j), j ~ U{1..r}."""
    h = ro.h_fine
    c_reg = ro.cic * ro.cic_register_fraction
    nb = int(np.ceil(x_max / h)) + 2
    edges = (np.arange(nb + 1) - 0.5) * h
    m = ro.em_gain ** (1.0 / ro.n_stages)
    j = np.arange(1, ro.n_stages + 1)
    means = m ** (-j.astype(float))
    cdf = 1.0 - np.exp(-np.clip(edges, 0, None)[None, :] / means[:, None])
    dens = np.diff(cdf, axis=1).mean(axis=0)
    out = c_reg * dens
    out[0] += 1.0 - c_reg
    return Lattice(out, 0.0, h)


def gaussian_lattice(sigma: float, h: float, nsig: float = 9.0) -> Lattice:
    if sigma <= 0:
        return Lattice(np.array([1.0]), 0.0, h)
    k = int(np.ceil(nsig * sigma / h)) + 1
    edges = (np.arange(-k, k + 2) - 0.5) * h
    p = np.diff(stats.norm.cdf(edges / sigma))
    return Lattice(p / p.sum(), -k * h, h)


def emccd_pixel(lam_signal: float, lam_bg: float, ro: EMCCDReadout) -> Lattice:
    """Per-pixel output: analog (input-referred lattice) or photon-counting (lattice {0, 1}, step 1)."""
    lam = lam_signal + lam_bg + ro.cic * (1.0 - ro.cic_register_fraction)
    hi = int(stats.poisson.isf(1e-16, lam)) + 2 if lam > 0 else 1
    n = np.arange(0, hi + 1)
    pn = stats.poisson.pmf(n, lam)
    x_max = hi + 12 * np.sqrt(hi * (ro.f2 - 1) + 1) + 12
    lat = register_output_lattice(pn, n, ro, x_max)
    if ro.cic * ro.cic_register_fraction > 0:
        lat = convolve(lat, serial_cic_lattice(ro, min(x_max, 50.0)))
    lat = convolve(lat, gaussian_lattice(ro.sigma_ro_e / ro.em_gain, ro.h_fine))
    if ro.mode == "analog":
        return lat
    thr = ro.pc_threshold_sigma * ro.sigma_ro_e / ro.em_gain
    p1 = lat.sf(thr)
    return Lattice(np.array([1.0 - p1, p1]), 0.0, 1.0)


# ------------------------------------------------------------- brute force
@numba.njit(cache=True)
def _register_mc(n_in, P, r, c_reg_per_stage, seed):
    np.random.seed(seed)
    out = np.empty(n_in.size, np.float64)
    for i in range(n_in.size):
        e = n_in[i]
        for s in range(r):
            if e > 0:
                e += np.random.binomial(e, P)
            if c_reg_per_stage > 0.0 and np.random.random() < c_reg_per_stage:
                e += 1
        out[i] = e
    return out


def register_mc(n_in: np.ndarray, ro: EMCCDReadout, seed: int) -> np.ndarray:
    """Stage-by-stage Monte Carlo of the multiplication register (Basden App. A recursion).

    Returns output electrons (not input-referred). Register CIC is injected with probability
    c_reg/r per stage, so its total probability per pixel is c_reg to first order.
    """
    P = ro.em_gain ** (1.0 / ro.n_stages) - 1.0
    c = ro.cic * ro.cic_register_fraction / ro.n_stages
    return _register_mc(np.asarray(n_in, np.int64), P, int(ro.n_stages), c, int(seed))
