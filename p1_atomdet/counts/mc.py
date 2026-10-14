"""Brute-force Monte Carlo of the per-site statistic, sharing no code with the lattice engine.

Each sample draws: atom loss time tau, rate multiplier m, neighbour occupancy, per-pixel
Poisson electrons (with PRNU), the pixel's read-noise class and component, Gaussian read
noise, full-well clipping and ADC rounding (CMOS); or the EMCCD register (gamma or
stage-by-stage) with CIC, read noise and optional photon-counting threshold. It is used to
cross-check the semi-analytic distributions at error rates ~1e-3 (a test, not a report line).
"""
from __future__ import annotations

import numpy as np

from ..sensor.cmos import PixelReadout
from ..sensor.emccd import EMCCDReadout, register_mc
from .pmf import SiteModel


def _read_noise(rng, ro: PixelReadout, shape):
    cw = np.array([c[0] for c in ro.classes])
    cls = rng.choice(len(ro.classes), size=shape, p=cw / cw.sum())
    out = np.empty(shape)
    for k, (_, comps) in enumerate(ro.classes):
        m = cls == k
        nk = int(m.sum())
        if nk == 0:
            continue
        w = np.array([c[0] for c in comps])
        ci = rng.choice(len(comps), size=nk, p=w / w.sum())
        mu = np.array([c[1] for c in comps])[ci]
        sd = np.array([c[2] for c in comps])[ci]
        out[m] = mu + sd * rng.standard_normal(nk)
    return out


def _lam(site: SiteModel, rng, n: int, present: bool, prnu: float):
    npx = site.frac.size
    f = site.frac.ravel()[None, :]
    b = site.bg.ravel()[None, :]
    sig = np.zeros((n, npx))
    if present:
        tau = site.tau.sample(rng, n)
        m = np.clip(1.0 + site.rate_rel_sd * rng.standard_normal(n), 0, None) if site.rate_rel_sd > 0 else np.ones(n)
        sig = site.mu * (m * tau)[:, None] * f
    for nb in site.neighbours:
        occ = rng.random(n) < site.neighbour_filling
        sig = sig + site.mu * occ[:, None] * nb.ravel()[None, :]
    if prnu > 0:
        sig = sig * np.clip(1.0 + prnu * rng.standard_normal(sig.shape), 0, None)
    return sig + b


def simulate_pixels(site: SiteModel, readout, n: int, present: bool, rng: np.random.Generator) -> np.ndarray:
    """(n, npx) pixel outputs in electron units (CMOS) or input-referred units / clicks (EMCCD)."""
    if isinstance(readout, PixelReadout):
        lam = _lam(site, rng, n, present, readout.prnu)
        e = np.minimum(rng.poisson(lam).astype(float), readout.full_well_e)
        x = e + _read_noise(rng, readout, e.shape)
        g = readout.gain_e_per_adu
        if g:
            code = np.floor(x / g + 0.5)
            if np.isfinite(readout.adc_max_code):
                code = np.minimum(code, readout.adc_max_code)
            x = code * g
        return x
    if isinstance(readout, EMCCDReadout):
        lam = _lam(site, rng, n, present, 0.0) + readout.cic * (1 - readout.cic_register_fraction)
        n_in = rng.poisson(lam)
        if readout.register_model == "stagewise":
            seed = int(rng.integers(2**62))
            out = register_mc(n_in.ravel(), readout, seed).reshape(n_in.shape) / readout.em_gain
        else:
            k = readout.f2 - 1.0
            out = np.where(n_in > 0, rng.gamma(np.maximum(n_in, 1) / k, k), 0.0)
            c_reg = readout.cic * readout.cic_register_fraction
            if c_reg > 0:
                hit = rng.random(n_in.shape) < c_reg
                j = rng.integers(1, readout.n_stages + 1, size=n_in.shape)
                mstage = readout.em_gain ** (1.0 / readout.n_stages)
                out = out + hit * rng.exponential(1.0, size=n_in.shape) * mstage ** (-j.astype(float))
        out = out + rng.standard_normal(out.shape) * readout.sigma_ro_e / readout.em_gain
        if readout.mode == "photon_counting":
            return (out > readout.pc_threshold_sigma * readout.sigma_ro_e / readout.em_gain).astype(float)
        return out
    raise TypeError(type(readout))


def statistic_samples(x: np.ndarray, kind: str = "box", weights=None, llr=None) -> np.ndarray:
    if kind in ("box", "matched"):
        w = np.ones(x.shape[1]) if weights is None else np.asarray(weights, float).ravel()
        return x @ w
    if kind == "bayes":
        return sum(np.interp(x[:, i], *llr[i]) for i in range(x.shape[1]))
    raise ValueError(kind)


def mc_error_rates(site, readout, threshold, n: int, rng, kind="box", weights=None, llr=None, chunk=200_000):
    """Return dict with FP, FN counts and samples for a fixed threshold."""
    out = {}
    for present, key in ((False, "fp"), (True, "fn")):
        k = 0
        done = 0
        while done < n:
            m = min(chunk, n - done)
            s = statistic_samples(simulate_pixels(site, readout, m, present, rng), kind, weights, llr)
            # same tie convention as the lattice: S >= T counts as 'atom'
            tol = 1e-9 * max(1.0, abs(threshold))
            k += int(np.sum(s >= threshold - tol)) if not present else int(np.sum(s < threshold - tol))
            done += m
        out[key] = k
        out[f"n_{key}"] = n
    return out
