"""Semi-analytic distributions of the per-site detection statistic.

The statistic is S = sum_i g_i(x_i) over the pixels i of a box around the site, where x_i is
the pixel output (electron units, offset subtracted) and g_i the classifier's per-pixel score:
  box sum         g_i(x) = m_i x        (m_i in {0, 1}: box mask)
  matched filter  g_i(x) = w_i x
  Bayesian LR     g_i(x) = log P1_i(x) / P0_i(x)   (per-pixel likelihood ratio, noise-map aware)

Conditional on the common variables (loss time tau, rate multiplier m, neighbour occupancy),
pixel outputs are independent, so P(S | common) is the convolution of per-pixel score
distributions. The common variables are then mixed with their quadrature weights:

  P(S | atom)    = sum_{s, nb} w_s w_nb  conv_i P(g_i(x_i) | lam_i = b_i + f_i mu s + nb_i),   s = tau * m
  P(S | no atom) = sum_{nb}         w_nb            conv_i P(g_i(x_i) | lam_i = b_i + nb_i)

Loss time tau and rate multiplier m enter only through their product s, so the mixture is
taken over a uniform grid in s (SiteModel.signal_multiplier_nodes). All errors of this
construction are controlled by numerical parameters covered by the convergence study:
tau nodes, the s-grid step ds, the statistic lattice step h_stat and the Poisson/Gaussian
truncation (1e-16 per pixel).
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np

from ..atoms.loss_phenom import TauDistribution
from ..sensor.cmos import PixelReadout, cmos_pixel
from ..sensor.emccd import EMCCDReadout, emccd_pixel
from .lattice import Lattice, convolve_many, map_to_grid, mixture, scale_to_grid


@dataclass
class SiteModel:
    frac: np.ndarray                      # (box, box) fraction of the atom's detected light per pixel
    mu: float                             # detected electrons if present for the whole exposure (all pixels, before frac)
    bg: np.ndarray                        # (box, box) background electrons per pixel (stray + dark)
    tau: TauDistribution = field(default_factory=TauDistribution.no_loss)
    rate_rel_sd: float = 0.0              # relative rms of mu between shots/sites (Gaussian; m < 0 mapped to 0)
    n_rate_nodes: int = 9                 # only for quadrature="tensor" (Gauss-Hermite, test use)
    neighbours: list = field(default_factory=list)   # list of (box, box) fraction arrays of neighbour atoms
    neighbour_filling: float = 0.5
    class_map: np.ndarray | None = None   # per-pixel read-noise class (None: ensemble mixture)
    quadrature: str = "s_grid"            # s_grid (default) | tensor
    ds: float = 0.01                      # step of the signal-multiplier grid s = tau * m
    x_range: float = 8.5                  # rate distribution integrated over |x| <= x_range sigma
    rate_hist: tuple | None = None        # empirical rate law: (bin centres m_k, weights w_k), uniform within bins;
                                          # overrides rate_rel_sd when given (m normalised to mean 1)

    def rate_nodes(self):
        if self.rate_rel_sd <= 0:
            return np.array([1.0]), np.array([1.0])
        x, w = np.polynomial.hermite_e.hermegauss(self.n_rate_nodes)
        w = w / w.sum()
        m = 1.0 + self.rate_rel_sd * x
        m = np.clip(m, 0.0, None)
        return m, w

    def signal_multiplier_nodes(self):
        """(s, weight) for the atom signal multiplier s = tau * m, the only way (tau, m) enter.

        tensor : Gauss-Hermite in m x tau nodes (fast but inaccurate for tail-dominated FN).
        s_grid : the joint law of s deposited on a uniform grid of step ds by linear (CIC) splitting,
                 with m integrated on a fine x grid (x_range sigma, 4000 points) for every tau node.
                 Mean-preserving; error O(ds^2), covered by the convergence study.
        """
        taus, wt = self.tau.all_nodes()
        keep = wt > 0
        taus, wt = taus[keep], wt[keep]
        if self.quadrature == "tensor":
            ms, wm = self.rate_nodes()
            return (taus[:, None] * ms[None, :]).ravel(), (wt[:, None] * wm[None, :]).ravel()
        if self.rate_hist is not None:
            c, w = (np.asarray(a, float) for a in self.rate_hist)
            w = w / w.sum()
            c = c / np.sum(c * w)
            half = 0.5 * np.median(np.diff(c))
            sub = (np.arange(40) + 0.5) / 40 * 2 - 1
            m = (c[:, None] + half * sub[None, :]).ravel()
            pm = np.repeat(w / 40, 40)
            vals = (taus[:, None] * m[None, :]).ravel()
            wts = (wt[:, None] * pm[None, :]).ravel()
        elif self.rate_rel_sd <= 0:
            vals, wts = taus, wt
        else:
            x = np.linspace(-self.x_range, self.x_range, 4001)
            px = np.exp(-0.5 * x**2)
            px /= px.sum()
            m = np.clip(1.0 + self.rate_rel_sd * x, 0.0, None)
            vals = (taus[:, None] * m[None, :]).ravel()
            wts = (wt[:, None] * px[None, :]).ravel()
        h = self.ds
        n = int(np.floor(vals.max() / h)) + 2
        u = vals / h
        i = np.floor(u + 1e-12).astype(int)
        f = np.clip(u - i, 0.0, 1.0)
        f[f < 1e-9] = 0.0
        grid = np.zeros(n + 1)
        np.add.at(grid, i, wts * (1 - f))
        np.add.at(grid, i + 1, wts * f)
        nz = np.nonzero(grid > 1e-20)[0]  # drops < 1e-18 total mass
        return nz * h, grid[nz]

    def neighbour_configs(self):
        """(extra light per pixel in units of mu, weight) for every occupancy configuration."""
        if not self.neighbours:
            return [(np.zeros_like(self.frac), 1.0)]
        out = []
        F = self.neighbour_filling
        for occ in itertools.product((0, 1), repeat=len(self.neighbours)):
            extra = sum((o * nb for o, nb in zip(occ, self.neighbours)), np.zeros_like(self.frac))
            k = sum(occ)
            out.append((extra, F**k * (1 - F) ** (len(occ) - k)))
        return out


def pixel_output(lam_signal: float, lam_bg: float, readout, pixel_class=None) -> Lattice:
    if isinstance(readout, EMCCDReadout):
        return emccd_pixel(lam_signal, lam_bg, readout)
    if isinstance(readout, PixelReadout):
        return cmos_pixel(lam_signal, lam_bg, readout, pixel_class)
    raise TypeError(type(readout))


@dataclass
class Classifier:
    kind: str = "box"                     # box | matched | bayes
    weights: np.ndarray | None = None     # box mask (0/1) or matched-filter weights, shape (box, box)
    h_stat: float | None = None           # lattice step of the statistic (None: automatic)
    llr: list | None = None               # per-pixel (x values, llr values) for kind == bayes

    def is_identity(self, step: float) -> bool:
        if self.kind != "box":
            return False
        w = self.weights
        return w is None or np.all((w == 0) | (w == 1))


def _score_lattice(lat: Lattice, clf: Classifier, i: int, h_stat: float) -> Lattice | None:
    if clf.kind == "box":
        w = 1.0 if clf.weights is None else clf.weights.flat[i]
        if w == 0:
            return None
        return lat if w == 1 else scale_to_grid(lat, w, h_stat)
    if clf.kind == "matched":
        w = clf.weights.flat[i]
        if w == 0:
            return None
        return scale_to_grid(lat, w, h_stat)
    if clf.kind == "bayes":
        xs, ls = clf.llr[i]
        g = np.interp(lat.x, xs, ls)
        return map_to_grid(lat, g, h_stat)
    raise ValueError(clf.kind)


def _default_h(readout, clf: Classifier) -> float:
    step = readout.h_fine if isinstance(readout, EMCCDReadout) else readout.step
    if isinstance(readout, EMCCDReadout) and readout.mode == "photon_counting":
        step = 1.0
    if clf.kind == "bayes":
        return 0.01
    if clf.is_identity(step):
        return step
    # Rebinning adds at most h^2/4 = (step*max w)^2/64 variance per pixel while keeping the
    # lattice short; the effect on error rates is covered by the h_stat convergence study.
    return step * float(np.max(clf.weights)) / 4.0


def site_statistic(site: SiteModel, readout, clf: Classifier, present: bool) -> Lattice:
    """Distribution of the statistic S for a site with (present=True) or without an atom."""
    h = clf.h_stat or _default_h(readout, clf)
    npx = site.frac.size
    classes = site.class_map.ravel() if site.class_map is not None else [None] * npx

    def one(signal_px: np.ndarray) -> Lattice:
        parts = []
        for i in range(npx):
            lat = pixel_output(float(signal_px.flat[i]), float(site.bg.flat[i]), readout, classes[i])
            s = _score_lattice(lat, clf, i, h)
            if s is not None:
                parts.append(s)
        return convolve_many(parts)

    comps, weights = [], []
    for extra, w_nb in site.neighbour_configs():
        nb_signal = site.mu * extra
        if not present:
            comps.append(one(nb_signal))
            weights.append(w_nb)
            continue
        svals, sw = site.signal_multiplier_nodes()
        for sv, a in zip(svals, sw):
            comps.append(one(nb_signal + site.mu * sv * site.frac))
            weights.append(w_nb * a)
    return mixture(comps, weights)


def optimal_linear_weights(site: SiteModel, readout, normalise: bool = True) -> np.ndarray:
    """Fisher linear discriminant for independent pixels: w_i = delta_i / ((v0_i + v1_i)/2).

    Means/variances from the per-pixel output lattices at tau = 1, m = 1 (no neighbours).
    """
    w = np.zeros(site.frac.size)
    for i in range(site.frac.size):
        l0 = pixel_output(0.0, float(site.bg.flat[i]), readout)
        l1 = pixel_output(float(site.mu * site.frac.flat[i]), float(site.bg.flat[i]), readout)
        w[i] = (l1.mean() - l0.mean()) / (0.5 * (l0.var() + l1.var()))
    w = w.reshape(site.frac.shape)
    return w / w.max() if normalise else w


def bayes_llr_tables(site: SiteModel, readout, use_class_map: bool = True, floor: float = 1e-300) -> list:
    """Per-pixel log-likelihood ratio tables log P1(x)/P0(x) at tau = 1, m = 1.

    The classifier's model of H1 ignores loss and rate spread (they are not knowable per shot);
    S = sum LLR is therefore the naive-Bayes statistic and its threshold is optimised on the
    true mixture distributions like any other statistic.
    """
    tabs = []
    classes = site.class_map.ravel() if (use_class_map and site.class_map is not None) else [None] * site.frac.size
    for i in range(site.frac.size):
        l0 = pixel_output(0.0, float(site.bg.flat[i]), readout, classes[i])
        l1 = pixel_output(float(site.mu * site.frac.flat[i]), float(site.bg.flat[i]), readout, classes[i])
        lo = min(l0.x0, l1.x0)
        hi = max(l0.x[-1], l1.x[-1])
        xs = np.arange(lo, hi + 0.5 * l0.h, l0.h)
        p0 = np.interp(xs, l0.x, l0.p, left=0, right=0)
        p1 = np.interp(xs, l1.x, l1.p, left=0, right=0)
        llr = np.log(np.maximum(p1, floor)) - np.log(np.maximum(p0, floor))
        # Clip: outside both supports the ratio is meaningless; saturate at +/- 50 nats.
        tabs.append((xs, np.clip(llr, -50, 50)))
    return tabs
