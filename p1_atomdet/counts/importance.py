"""Exponentially tilted importance sampling for linear statistics S = sum_i w_i x_i (CMOS chain).

Latent variables per pixel: n_i ~ Poisson(lam_i), e_i ~ Gaussian mixture {a_c, mu_c, s_c}.
For a tilt theta, the tilted law is again Poisson x Gaussian mixture:
  n_i ~ Poisson(lam_i exp(theta w_i)),
  component c with probability ~ a_c exp(theta w_i mu_c + theta^2 w_i^2 s_c^2/2),
  e_i ~ N(mu_c + theta w_i s_c^2, s_c),
and the likelihood ratio is exp(sum_i [K_i(theta) - theta w_i (n_i + e_i)]) with K_i the
cumulant generating function. The event is evaluated on the ADC-rounded output, which is a
deterministic function of the latent variables, so the estimator stays unbiased.

For FN, the rare atom-loss branch is sampled defensively: the 'lost' branch is drawn with
probability q_lost (default 0.5) and its loss time from the exact truncated exponential.
theta is chosen by the saddle-point condition K'(theta) = T at tau = 1, m = 1.
This module shares only the SiteModel container with the lattice engine.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize

from ..sensor.cmos import PixelReadout
from .pmf import SiteModel


@dataclass
class ISResult:
    p: float
    rel_err: float
    ess: float
    n: int
    theta: float

    def ci95(self):
        return self.p * (1 - 1.96 * self.rel_err), self.p * (1 + 1.96 * self.rel_err)


def _comps(ro: PixelReadout):
    a, mu, s = [], [], []
    for cw, cl in ro.classes:
        for w, m, sd in cl:
            a.append(cw * w); mu.append(m); s.append(sd)
    a = np.asarray(a); return a / a.sum(), np.asarray(mu), np.asarray(s)


def cgf(theta, lam, w, comps):
    a, mu, s = comps
    tw = theta * w
    K_pois = lam * np.expm1(tw)
    K_noise = np.log(np.sum(a[None, :] * np.exp(tw[:, None] * mu[None, :] + 0.5 * tw[:, None] ** 2 * s[None, :] ** 2), axis=1))
    return K_pois + K_noise


def cgf_prime_total(theta, lam, w, comps):
    eps = 1e-6
    return (np.sum(cgf(theta + eps, lam, w, comps)) - np.sum(cgf(theta - eps, lam, w, comps))) / (2 * eps)


def saddle_theta(T, lam, w, comps):
    f = lambda th: cgf_prime_total(th, lam, w, comps) - T
    lo, hi = -5.0, 5.0
    return optimize.brentq(f, lo, hi)


def _sample_tilted(rng, theta, lam, w, comps, n):
    """lam: (n, npx). Returns latent sums (n,), rounded outputs (n, npx), log-weights (n,)."""
    a, mu, s = comps
    tw = theta * w[None, :]
    nn = rng.poisson(lam * np.exp(tw))
    logits = np.log(a)[None, None, :] + tw[..., None] * mu[None, None, :] + 0.5 * tw[..., None] ** 2 * s[None, None, :] ** 2
    pc = np.exp(logits - logits.max(-1, keepdims=True))
    pc /= pc.sum(-1, keepdims=True)
    u = rng.random(lam.shape)[..., None]
    ci = (u > np.cumsum(pc, -1)).sum(-1)
    ci = np.minimum(ci, len(a) - 1)
    e = mu[ci] + tw * s[ci] ** 2 + s[ci] * rng.standard_normal(lam.shape)
    K = lam * np.expm1(tw) + np.log(np.sum(a * np.exp(tw[..., None] * mu + 0.5 * tw[..., None] ** 2 * s**2), axis=-1))
    logw = np.sum(K - tw * (nn + e), axis=1)
    return nn + e, logw


def _round(x, ro: PixelReadout):
    x = np.minimum(x, ro.full_well_e) if np.isfinite(ro.full_well_e) else x
    g = ro.gain_e_per_adu
    return np.floor(x / g + 0.5) * g if g else x


def is_fp(site: SiteModel, ro: PixelReadout, weights, T: float, n: int, rng) -> ISResult:
    """P(S >= T | no atom) by tilted sampling (neighbours and PRNU must be absent)."""
    if site.neighbours or ro.prnu > 0:
        raise NotImplementedError("IS supports the single-site, PRNU-free case")
    w = np.ones(site.frac.size) if weights is None else np.asarray(weights, float).ravel()
    lam0 = site.bg.ravel()
    comps = _comps(ro)
    th = saddle_theta(T, lam0, w, comps)
    lam = np.broadcast_to(lam0, (n, lam0.size))
    x, logw = _sample_tilted(rng, th, lam, w, comps, n)
    S = _round(x, ro) @ w
    ind = S >= T - 1e-9 * max(1, abs(T))
    return _estimate(ind, logw, th)


def is_fn(site: SiteModel, ro: PixelReadout, weights, T: float, n: int, rng, q_lost: float = 0.5,
          rate_shift_sd: float = 3.0) -> ISResult:
    """P(S < T | atom).

    Proposal: the loss branch is drawn with probability q_lost (exact truncated-exponential loss
    time, no pixel tilt: lost atoms already fall below T without help). Survivors get the
    saddle-point pixel tilt and, when the rate spread is non-zero, a defensive two-component
    proposal for the rate multiplier m: half natural N(1, rs), half shifted N(1 - k rs, rs).
    """
    if site.neighbours or ro.prnu > 0 or site.rate_hist is not None:
        raise NotImplementedError("IS supports the single-site, PRNU-free, Gaussian-rate case")
    w = np.ones(site.frac.size) if weights is None else np.asarray(weights, float).ravel()
    comps = _comps(ro)
    f, b = site.frac.ravel(), site.bg.ravel()
    S_true = site.tau.p_survive
    L = site.tau.hazard
    if L is None and site.tau.p_loss > 0:
        raise NotImplementedError("IS needs a model-A hazard")
    ql = q_lost if site.tau.p_loss > 0 else 0.0
    lost = rng.random(n) < ql
    nl, ns = int(lost.sum()), int((~lost).sum())
    rs = site.rate_rel_sd
    logw_all = np.empty(n)
    S_all = np.empty(n)
    # survivors
    if ns:
        if rs > 0:
            k = rate_shift_sd
            pick = rng.random(ns) < 0.5
            x = rng.standard_normal(ns) - np.where(pick, k, 0.0)
            # weight phi(x) / (0.5 phi(x) + 0.5 phi(x + k))
            lw_m = -np.log(0.5 + 0.5 * np.exp(-k * x - 0.5 * k * k))
            m = np.clip(1.0 + rs * x, 0, None)
        else:
            lw_m = np.zeros(ns)
            m = np.ones(ns)
        th = min(saddle_theta(T, b + site.mu * f, w, comps), 0.0)
        lam = b[None, :] + site.mu * m[:, None] * f[None, :]
        xs, lw = _sample_tilted(rng, th, lam, w, comps, ns)
        S_all[~lost] = _round(xs, ro) @ w
        logw_all[~lost] = lw + lw_m + np.log(S_true) - np.log(1 - ql)
    else:
        th = 0.0
    # lost atoms
    if nl:
        v = rng.random(nl)
        tau = -np.log1p(-v * (1 - np.exp(-L))) / L
        m = np.clip(1.0 + rs * rng.standard_normal(nl), 0, None) if rs > 0 else np.ones(nl)
        lam = b[None, :] + site.mu * (m * tau)[:, None] * f[None, :]
        xs, lw = _sample_tilted(rng, 0.0, lam, w, comps, nl)
        S_all[lost] = _round(xs, ro) @ w
        logw_all[lost] = lw + np.log(1 - S_true) - np.log(ql)
    ind = S_all < T - 1e-9 * max(1, abs(T))
    return _estimate(ind, logw_all, th)


def _estimate(ind, logw, th) -> ISResult:
    wts = np.where(ind, np.exp(logw), 0.0)
    n = wts.size
    p = wts.mean()
    sd = wts.std(ddof=1) / np.sqrt(n)
    ess = wts.sum() ** 2 / max(np.sum(wts**2), 1e-300)
    return ISResult(float(p), float(sd / p) if p > 0 else np.inf, float(ess), n, float(th))
