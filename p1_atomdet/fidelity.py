"""Classification error rates, optimal threshold, imaging time to a target, confidence intervals.

Definitions (as in Manetsch et al., arXiv:2403.12021v4, Methods):
  FP = P(S >= T | no atom),  FN = P(S < T | atom present at the start of the image),
  infidelity eps = (1 - F) FP + F FN,  with F the filling fraction (default 0.5, a stated assumption).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize, stats

from .counts.lattice import Lattice


@dataclass
class ErrorRates:
    threshold: float
    fp: float
    fn: float
    eps: float
    filling: float

    @property
    def fidelity(self) -> float:
        return 1.0 - self.eps

    def as_dict(self):
        return {"threshold": self.threshold, "fp": self.fp, "fn": self.fn, "eps": self.eps,
                "fidelity": self.fidelity, "filling": self.filling}


def error_curves(p0: Lattice, p1: Lattice):
    """Thresholds on the union lattice and FP(T), FN(T) for each."""
    h = p0.h
    lo = min(p0.x0, p1.x0)
    hi = max(p0.x[-1], p1.x[-1]) + h
    T = np.arange(lo, hi + 0.5 * h, h)
    # FP(T) = P0(S >= T): reverse cumulative sum, aligned
    def sf(lat):
        a = np.zeros(T.size)
        i0 = int(round((lat.x0 - lo) / h))
        a[i0:i0 + lat.p.size] = lat.p
        c = np.cumsum(a[::-1])[::-1]
        return c
    def cdf_lt(lat):
        a = np.zeros(T.size)
        i0 = int(round((lat.x0 - lo) / h))
        a[i0:i0 + lat.p.size] = lat.p
        c = np.concatenate([[0.0], np.cumsum(a)[:-1]])
        return c
    return T, sf(p0), cdf_lt(p1)


def optimal_threshold(p0: Lattice, p1: Lattice, filling: float = 0.5) -> ErrorRates:
    T, fp, fn = error_curves(p0, p1)
    eps = (1 - filling) * fp + filling * fn
    i = int(np.argmin(eps))
    return ErrorRates(float(T[i]), float(fp[i]), float(fn[i]), float(eps[i]), filling)


def rates_at(p0: Lattice, p1: Lattice, threshold: float, filling: float = 0.5) -> ErrorRates:
    fp = p0.sf(threshold)
    fn = p1.cdf_lt(threshold)
    return ErrorRates(threshold, fp, fn, (1 - filling) * fp + filling * fn, filling)


def time_to_target(eps_of_t, target: float, t_lo: float, t_hi: float, survival_of_t=None,
                   min_survival: float = 0.99, n_scan: int = 25):
    """Shortest imaging time with eps(t) <= target and survival(t) >= min_survival.

    eps_of_t need not be monotone (loss makes it rise at long t), so we scan a log grid,
    then refine the first crossing by Brent's method. Returns (t, info) or (None, info).
    """
    ts = np.geomspace(t_lo, t_hi, n_scan)
    eps = np.array([eps_of_t(t) for t in ts])
    ok = eps <= target
    if survival_of_t is not None:
        surv = np.array([survival_of_t(t) for t in ts])
        ok &= surv >= min_survival
    info = {"t_scan": ts, "eps_scan": eps}
    if not ok.any():
        return None, info
    j = int(np.argmax(ok))
    if j == 0:
        return float(ts[0]), info
    f = lambda t: np.log(eps_of_t(t)) - np.log(target)
    t = optimize.brentq(f, ts[j - 1], ts[j], xtol=1e-4 * ts[j], rtol=1e-4)
    return float(t), info


def clopper_pearson(k: int, n: int, cl: float = 0.95):
    a = 1 - cl
    lo = 0.0 if k == 0 else stats.beta.ppf(a / 2, k, n - k + 1)
    hi = 1.0 if k == n else stats.beta.ppf(1 - a / 2, k + 1, n - k)
    return float(lo), float(hi)
