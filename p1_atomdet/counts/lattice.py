"""Probability mass functions on uniform lattices x = x0 + i*h, and the operations on them.

All count distributions in P1 are built from these: per-pixel PMFs, their convolution over
pixels (sums of independent variables), mixtures (over loss time, rate inhomogeneity,
neighbour occupancy) and rebinning onto a common step for weighted statistics.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import signal

# Below this many output bins, direct convolution is used (exact to rounding); above, FFT.
_DIRECT_MAX = 4000
# FFT round-off floor relative to the peak; values below are zeroed after FFT convolution.
FFT_FLOOR = 1e-15


@dataclass
class Lattice:
    p: np.ndarray
    x0: float
    h: float

    @property
    def x(self) -> np.ndarray:
        return self.x0 + self.h * np.arange(self.p.size)

    def total(self) -> float:
        return float(self.p.sum())

    def mean(self) -> float:
        return float(np.sum(self.x * self.p) / self.p.sum())

    def var(self) -> float:
        m = self.mean()
        return float(np.sum((self.x - m) ** 2 * self.p) / self.p.sum())

    def trim(self, eps: float = 1e-300) -> "Lattice":
        nz = np.nonzero(self.p > eps)[0]
        if nz.size == 0:
            return Lattice(np.array([0.0]), self.x0, self.h)
        a, b = nz[0], nz[-1] + 1
        return Lattice(self.p[a:b].copy(), self.x0 + a * self.h, self.h)

    def sf(self, t: float) -> float:
        """P(X >= t)."""
        return float(self.p[self.x >= t - 1e-9 * self.h].sum())

    def cdf_lt(self, t: float) -> float:
        """P(X < t)."""
        return float(self.p[self.x < t - 1e-9 * self.h].sum())


def delta(x: float = 0.0, h: float = 1.0) -> Lattice:
    return Lattice(np.array([1.0]), x, h)


def _conv(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    n = a.size + b.size - 1
    if min(a.size, b.size) < 64 or n <= _DIRECT_MAX:
        return np.convolve(a, b)
    c = signal.fftconvolve(a, b)
    peak = c.max()
    c[c < FFT_FLOOR * peak] = 0.0
    return c


def convolve(a: Lattice, b: Lattice) -> Lattice:
    if not np.isclose(a.h, b.h, rtol=1e-12):
        raise ValueError(f"lattice steps differ: {a.h} vs {b.h}")
    return Lattice(_conv(a.p, b.p), a.x0 + b.x0, a.h)


def convolve_many(items: list[Lattice], trim_eps: float = 1e-300) -> Lattice:
    """Pairwise tree reduction (balanced sizes keep FFT error and cost low)."""
    items = [it.trim(trim_eps) for it in items]
    if not items:
        raise ValueError("nothing to convolve")
    while len(items) > 1:
        nxt = []
        for i in range(0, len(items) - 1, 2):
            nxt.append(convolve(items[i], items[i + 1]).trim(trim_eps))
        if len(items) % 2:
            nxt.append(items[-1])
        items = nxt
    return items[0]


def mixture(items: list[Lattice], weights) -> Lattice:
    """Weighted sum of lattices with a common step (aligned onto the union support)."""
    weights = np.asarray(weights, float)
    h = items[0].h
    for it in items:
        if not np.isclose(it.h, h, rtol=1e-12):
            raise ValueError("mixture requires a common step")
    lo = min(it.x0 for it in items)
    hi = max(it.x0 + h * (it.p.size - 1) for it in items)
    n = int(round((hi - lo) / h)) + 1
    out = np.zeros(n)
    for it, w in zip(items, weights):
        i0 = int(round((it.x0 - lo) / h))
        out[i0:i0 + it.p.size] += w * it.p
    return Lattice(out, lo, h)


def rebin(values: np.ndarray, probs: np.ndarray, x0: float, h: float, n: int | None = None) -> Lattice:
    """Deposit point masses (values, probs) onto lattice x0 + i h by linear (cloud-in-cell) splitting.

    Preserves total mass and the mean exactly; adds variance <= h^2/4 per deposit (convergence in h tested).
    """
    v = (np.asarray(values, float) - x0) / h
    if n is None:
        n = int(np.floor(v.max())) + 2
    i = np.floor(v).astype(int)
    f = v - i
    out = np.zeros(n)
    ok = (i >= 0) & (i < n - 1)
    np.add.at(out, i[ok], probs[ok] * (1 - f[ok]))
    np.add.at(out, i[ok] + 1, probs[ok] * f[ok])
    edge = i == n - 1
    np.add.at(out, i[edge], probs[edge])
    if np.any(i < 0) or np.any(i > n - 1):
        raise ValueError("rebin: values outside lattice")
    return Lattice(out, x0, h)


def scale_to_grid(lat: Lattice, weight: float, h_out: float) -> Lattice:
    """Distribution of weight * X on a lattice with step h_out (anchored at a multiple of h_out)."""
    vals = weight * lat.x
    x0 = np.floor(vals.min() / h_out) * h_out
    return rebin(vals, lat.p, x0, h_out)


def map_to_grid(lat: Lattice, g_values: np.ndarray, h_out: float) -> Lattice:
    """Distribution of g(X) for arbitrary per-point scores g_values (same length as lat.p)."""
    x0 = np.floor(np.min(g_values) / h_out) * h_out
    return rebin(g_values, lat.p, x0, h_out)
