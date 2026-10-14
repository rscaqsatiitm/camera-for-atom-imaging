"""CMOS / sCMOS / qCMOS per-pixel output distributions (electron units, offset subtracted).

Chain for one pixel: Poisson(lambda) photo- and dark electrons -> read noise (mixture of
Gaussian components, one pixel class or the ensemble mixture of classes) -> ADC code with
gain g (e-/ADU), clipped at full well and at 2^bits - 1. The output is the distribution of
g * code, i.e. a lattice of step g in electron units. A qCMOS in photon-number-resolving
mode is the same chain with g = 1 e-/code (rounding to integer electrons).

Without an ADC (ideal analog readout, used for analytic limits) the continuous output is
bin-integrated onto a fine lattice of step h; the Poisson part is exact.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import special, stats

from ..counts.lattice import Lattice


def _phi(z):
    return 0.5 * special.erfc(-z / np.sqrt(2.0))


@dataclass
class PixelReadout:
    """Readout description shared by all pixels of a site (per-pixel maps enter via `classes`)."""

    classes: list                    # [(class_weight, [(w, mean_e, sigma_e), ...]), ...]
    gain_e_per_adu: float | None     # None: ideal analog output on lattice step h_fine
    h_fine: float = 0.05
    full_well_e: float = np.inf
    adc_max_code: float = np.inf
    prnu: float = 0.0
    n_prnu_nodes: int = 5
    extra: dict = field(default_factory=dict)

    @property
    def step(self) -> float:
        return self.gain_e_per_adu if self.gain_e_per_adu else self.h_fine

    def ensemble_components(self):
        """Flatten class mixture into one list of Gaussian components (ensemble over pixels)."""
        comps = []
        for cw, cl in self.classes:
            for w, m, s in cl:
                comps.append((cw * w, m, s))
        return comps


def poisson_support(lam: float, tail: float = 1e-16) -> np.ndarray:
    lam = max(lam, 0.0)
    q = stats.poisson.isf(tail, lam) if lam > 0 else 0.0
    if not np.isfinite(q):  # below double precision: generous Chernoff-type bound
        q = lam + 14 * np.sqrt(lam) + 40
    hi = int(q) + 2 if lam > 0 else 1
    return np.arange(0, max(hi, 1) + 1)


def electrons_pmf(lam_signal: float, lam_bg: float, prnu: float = 0.0, n_nodes: int = 5):
    """P(n) for n = Poisson(lam_signal*(1+d) + lam_bg), d ~ N(0, prnu) per pixel (Gauss-Hermite mixture)."""
    if prnu <= 0 or lam_signal <= 0:
        n = poisson_support(lam_signal + lam_bg)
        return n, stats.poisson.pmf(n, lam_signal + lam_bg)
    x, w = np.polynomial.hermite_e.hermegauss(n_nodes)
    w = w / w.sum()
    lams = np.clip(lam_signal * (1 + prnu * x) + lam_bg, 0, None)
    n = poisson_support(lams.max())
    p = np.sum(w[:, None] * stats.poisson.pmf(n[None, :], lams[:, None]), axis=0)
    return n, p


def pixel_lattice(n: np.ndarray, pn: np.ndarray, comps, step: float, full_well_e: float = np.inf,
                  max_code: float = np.inf, nsig: float = 9.0) -> Lattice:
    """Output distribution on lattice k*step given electron PMF pn(n) and Gaussian noise components.

    P(k) = sum_n pn(n) sum_c w_c [Phi(((k+1/2) step - n - mu_c)/s_c) - Phi(((k-1/2) step - n - mu_c)/s_c)]
    Electrons above full well are clipped before the noise is added; codes above max_code clipped after.
    """
    n = np.minimum(n, full_well_e)
    smax = max(s for _, _, s in comps)
    mmax = max(abs(m) for _, m, _ in comps)
    lo = int(np.floor((n.min() - mmax - nsig * smax) / step)) - 1
    hi = int(np.ceil((n.max() + mmax + nsig * smax) / step)) + 1
    k = np.arange(lo, hi + 1)
    edges = (np.arange(lo, hi + 2) - 0.5) * step  # bin edges
    out = np.zeros(k.size)
    for w, mu, s in comps:
        if w <= 0:
            continue
        if s <= 0:
            # noiseless component: deterministic rounding of n + mu
            idx = np.clip(np.floor((n + mu) / step + 0.5).astype(int) - lo, 0, k.size - 1)
            np.add.at(out, idx, w * pn)
            continue
        cdf = _phi((edges[None, :] - n[:, None] - mu) / s)  # (n, edges)
        out += w * (pn @ np.diff(cdf, axis=1))
    lat = Lattice(out, lo * step, step)
    if np.isfinite(max_code):
        kmax = int(max_code)
        x = lat.x / step
        top = x >= kmax
        if top.any():
            p = lat.p.copy()
            j = int(np.argmax(top))
            p[j] = p[top].sum()
            p = p[: j + 1]
            lat = Lattice(p, lat.x0, step)
    return lat


def cmos_pixel(lam_signal: float, lam_bg: float, ro: PixelReadout, pixel_class: int | None = None) -> Lattice:
    """Output lattice of one pixel. pixel_class=None uses the ensemble mixture of read-noise classes."""
    n, pn = electrons_pmf(lam_signal, lam_bg, ro.prnu, ro.n_prnu_nodes)
    comps = ro.ensemble_components() if pixel_class is None else ro.classes[pixel_class][1]
    return pixel_lattice(n, pn, comps, ro.step, ro.full_well_e, ro.adc_max_code)


def readout_from_sensor(sm, h_fine: float = 0.05, use_adc: bool = True) -> PixelReadout:
    """Build a PixelReadout from a camsim_common SensorModel."""
    from camsim_common.provenance import val
    g = float(val(sm.adc.gain_e_per_adu)) if use_adc else None
    if sm.qcmos_pnr and use_adc:
        g = 1.0
    bits = val(sm.adc.adc_bits)
    max_code = (2 ** int(bits) - 1 - float(val(sm.adc.offset_adu))) if bits else np.inf
    return PixelReadout(classes=sm.read_noise.pixel_classes(), gain_e_per_adu=g, h_fine=h_fine,
                        full_well_e=float(val(sm.adc.full_well_e)), adc_max_code=max_code if use_adc else np.inf,
                        prnu=float(val(sm.prnu) or 0.0))
