"""Loss model A (phenomenological) and the shared 'time-before-loss' representation.

Both loss models (A here, B in loss_randomwalk.py) are reduced to the distribution of the
fraction tau in [0, 1] of the exposure during which the atom is present, given that it was
present at the start. The count models then mix over tau, so swapping A and B is a direct
substitution of one TauDistribution for another.

Model A: a constant hazard per unit exposure time. Two contributions:
  (i) photon-induced loss, probability p_ph per scattered photon: hazard p_ph * R * t_exp;
  (ii) vacuum (background-gas) loss with lifetime tau_vac: hazard t_exp / tau_vac.
The total pseudo-loss L = p_ph R t + t/tau_vac gives survival S = exp(-L) and a truncated
exponential density L exp(-L tau) for the loss time.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import special, stats


@dataclass
class TauDistribution:
    """Atom-present fraction tau: point mass p_survive at tau = 1, plus nodes/weights on [0, 1)."""

    p_survive: float
    nodes: np.ndarray
    weights: np.ndarray
    label: str = ""

    def __post_init__(self):
        self.nodes = np.asarray(self.nodes, float)
        self.weights = np.asarray(self.weights, float)
        tot = self.p_survive + self.weights.sum()
        if not np.isclose(tot, 1.0, rtol=0, atol=1e-9):
            raise ValueError(f"TauDistribution weights sum to {tot}")

    @property
    def p_loss(self) -> float:
        return 1.0 - self.p_survive

    def all_nodes(self):
        """(tau, weight) including the survivor point mass at tau = 1."""
        return np.concatenate([[1.0], self.nodes]), np.concatenate([[self.p_survive], self.weights])

    def mean_tau(self) -> float:
        t, w = self.all_nodes()
        return float(np.sum(t * w))

    @classmethod
    def no_loss(cls):
        return cls(1.0, np.zeros(0), np.zeros(0), "no loss")

    @classmethod
    def from_hazard(cls, L: float, n_nodes: int = 96):
        """Exact truncated-exponential loss time, integrated by Gauss-Legendre on [0, 1]."""
        if L <= 0:
            return cls.no_loss()
        x, w = np.polynomial.legendre.leggauss(n_nodes)
        tau = 0.5 * (x + 1.0)
        wt = 0.5 * w * L * np.exp(-L * tau)
        S = np.exp(-L)
        wt *= (1.0 - S) / wt.sum()  # remove the O(1e-16) quadrature error so weights sum exactly
        return cls(S, tau, wt, f"hazard L={L:.3e}")

    @classmethod
    def from_samples(cls, tau_samples: np.ndarray, n_bins: int = 64, label: str = "samples"):
        """Empirical distribution (model B): tau >= 1 counts as survival; lost times are binned."""
        t = np.asarray(tau_samples, float)
        n = t.size
        lost = t[t < 1.0]
        S = 1.0 - lost.size / n
        if lost.size == 0:
            return cls(1.0, np.zeros(0), np.zeros(0), label)
        h, e = np.histogram(lost, bins=n_bins, range=(0.0, 1.0))
        c = 0.5 * (e[1:] + e[:-1])
        m = h > 0
        return cls(S, c[m], h[m] / n, label)


def hazard(p_photon: float, rate_s: float, t_exp_s: float, tau_vac_s: float = np.inf) -> float:
    """Pseudo-loss L for model A."""
    return p_photon * rate_s * t_exp_s + (t_exp_s / tau_vac_s if np.isfinite(tau_vac_s) else 0.0)


def p_photon_from_survival(S: float, rate_s: float, t_exp_s: float, tau_vac_s: float = np.inf) -> float:
    """Fit model A to a published survival: loss per scattered photon after removing vacuum loss."""
    L = -np.log(S) - (t_exp_s / tau_vac_s if np.isfinite(tau_vac_s) else 0.0)
    return max(L, 0.0) / (rate_s * t_exp_s)


def lossy_poisson_pmf(n: np.ndarray, mu: float, bg: float, tau: TauDistribution) -> np.ndarray:
    """P(N = n | atom at start) for N = Poisson(bg + mu*tau), mixed over tau (independent-hazard convention)."""
    n = np.asarray(n)
    t, w = tau.all_nodes()
    lam = bg + mu * t[:, None]
    return np.sum(w[:, None] * stats.poisson.pmf(n[None, :], lam), axis=0)


def manetsch_lossy_poisson_pmf(n: np.ndarray, lam0: float, lam1: float, L: float) -> np.ndarray:
    """Lossy Poisson model of Manetsch et al., arXiv:2403.12021v4, Methods ("Imaging model and
    characterization"), derived independently here (report Sec. 2.1):

    Each atom photon triggers loss with probability q = L/lam1 (the triggering photon is counted).
      survivors:  exp(-L) * Poisson(n; lam0 + lam1 - L)
      lost:       q (1-q)^(n-1) exp(lam0 L/(lam1-L)) [Q(n, b lam0) - Q(n, b lam0 + lam1)],  b = 1/(1-q), n >= 1
    with Q the regularised upper incomplete gamma function.
    """
    n = np.asarray(n, float)
    surv = np.exp(-L) * stats.poisson.pmf(n, lam0 + lam1 - L)
    q = L / lam1
    b = 1.0 / (1.0 - q)
    out = surv.copy()
    m = n >= 1
    nn = n[m]
    lost = q * (1 - q) ** (nn - 1) * np.exp(lam0 * L / (lam1 - L)) * (
        special.gammaincc(nn, b * lam0) - special.gammaincc(nn, b * lam0 + lam1))
    out[m] += lost
    return out
