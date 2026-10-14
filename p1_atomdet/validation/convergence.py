"""Convergence appendix: every numerical parameter of the P1 pipeline, varied about the V1 baseline.

  (a) statistic lattice step h_stat (matched filter)        (b) signal-multiplier grid ds
  (c) number of loss-time nodes n_tau                        (d) pixel super-sampling s
  (e) vectorial PSF quadrature order n_quad                  (f) importance-sampling sample count
  (g) Langevin time step dt (model B)                        (h) Poisson truncation tail
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from ..atoms.loss_phenom import TauDistribution
from ..atoms.loss_randomwalk import RandomWalkParams, simulate_loss_times
from ..constants import CS133
from ..counts.importance import is_fp
from ..counts.pmf import Classifier, SiteModel, optimal_linear_weights, site_statistic
from ..fidelity import optimal_threshold, rates_at
from ..optics.array import PSFSpec, pixel_fractions
from ..sensor import cmos
from ..sensor.cmos import PixelReadout
from .common import md_table


def _base(cfg, **over):
    b = {**cfg["base"], **over}
    pitch = b["pitch_obj_um"] * 1e-6
    fr = pixel_fractions(PSFSpec("gaussian", waist_m=b["waist_px"] * pitch), pitch, b["box"], supersample=b.get("supersample", 8))
    tau = TauDistribution.from_hazard(b["L"], b.get("n_tau", 48))
    site = SiteModel(frac=fr, mu=b["mu"], bg=np.full(fr.shape, b["lam0"] / fr.size), tau=tau, rate_rel_sd=b["rs"],
                     ds=b.get("ds", 0.01))
    ro = PixelReadout(classes=[(1.0, [(1.0, 0.0, b["sigma"])])], gain_e_per_adu=b["gain"])
    return site, ro


def _box_rates(site, ro, T):
    p0 = site_statistic(site, ro, Classifier("box"), False)
    p1 = site_statistic(site, ro, Classifier("box"), True)
    return rates_at(p0, p1, T)


def run(cfg: dict, outdir: Path):
    res = {}
    T = cfg["threshold_e"]
    # (a) h_stat for the matched filter, relative to the default step*max(w)/4
    site, ro = _base(cfg)
    w = optimal_linear_weights(site, ro)
    h0 = ro.step * w.max() / 4
    rows = []
    for k in cfg["h_factors"]:
        c = Classifier("matched", weights=w, h_stat=h0 * k)
        r = optimal_threshold(site_statistic(site, ro, c, False), site_statistic(site, ro, c, True), cfg["filling"])
        rows.append({"h_factor": k, "h_stat": h0 * k, "fp": r.fp, "fn": r.fn, "eps": r.eps})
    res["h_stat"] = rows
    # (b) ds
    res["ds"] = []
    for ds in cfg["ds_values"]:
        site, ro = _base(cfg, ds=ds)
        r = _box_rates(site, ro, T)
        res["ds"].append({"ds": ds, "fp": r.fp, "fn": r.fn})
    # (c) n_tau
    res["n_tau"] = []
    for n in cfg["n_tau_values"]:
        site, ro = _base(cfg, n_tau=n)
        r = _box_rates(site, ro, T)
        res["n_tau"].append({"n_tau": n, "fp": r.fp, "fn": r.fn})
    # (d) supersampling of the pixel integral (vectorial PSF, NA 0.55, 2 px per FWHM: the sharpest case in V1-like use)
    psf = PSFSpec("vectorial", na=cfg["psf_na"], n_quad=64)
    pitch = psf.fwhm() / cfg["px_per_fwhm"]
    ref = pixel_fractions(psf, pitch, 7, (0.25, 0.1), supersample=max(cfg["supersample_values"]) * 2)
    res["supersample"] = []
    for s in cfg["supersample_values"]:
        f = pixel_fractions(psf, pitch, 7, (0.25, 0.1), supersample=s)
        res["supersample"].append({"s": s, "max_abs_diff": float(np.abs(f - ref).max()), "centre_rel": float(abs(f.max() / ref.max() - 1))})
    # (e) PSF quadrature order (NA 0.8, pixel fractions within 7x7 at 2 px/FWHM)
    res["n_quad"] = []
    psf_ref = PSFSpec("vectorial", na=0.8, n_quad=512)
    pitch8 = psf_ref.fwhm() / cfg["px_per_fwhm"]
    ref = pixel_fractions(psf_ref, pitch8, 7, supersample=8)
    for nq in cfg["n_quad_values"]:
        f = pixel_fractions(PSFSpec("vectorial", na=0.8, n_quad=nq), pitch8, 7, supersample=8)
        res["n_quad"].append({"n_quad": nq, "max_abs_diff": float(np.abs(f - ref).max())})
    # (f) IS sample count at a deep FP (threshold chosen where FP ~ 1e-7)
    site, ro = _base(cfg, rs=0.0, L=0.0)
    p0 = site_statistic(site, ro, Classifier("box"), False)
    Tdeep = float(p0.x[np.argmax(np.cumsum(p0.p[::-1])[::-1] < cfg["is_target"])])
    lat = p0.sf(Tdeep)
    rng = np.random.default_rng(cfg["seed"])
    res["is"] = {"T": Tdeep, "lattice": lat, "rows": []}
    for n in cfg["is_samples"]:
        r = is_fp(site, ro, None, Tdeep, n, rng)
        res["is"]["rows"].append({"n": n, "p": r.p, "rel_err": r.rel_err, "dev_sigma": (r.p - lat) / (r.rel_err * r.p),
                                  "ess": r.ess})
    # (g) Langevin time step (thermostat regime of V1)
    m, k = float(CS133.mass_kg), CS133.k
    res["dt"] = []
    for dt in cfg["dt_values"]:
        p = RandomWalkParams(mass_kg=m, k_m=k, **{**cfg["model_b"], "dt_s": dt})
        r = simulate_loss_times(p, cfg["model_b_atoms"], cfg["seed"])
        n = cfg["model_b_atoms"]
        res["dt"].append({"dt_us": dt * 1e6, "survival": r["survival"],
                          "survival_se": float(np.sqrt(r["survival"] * (1 - r["survival"]) / n)), "T_eff_uK": r["T_eff_K"] * 1e6})
    # (h) Poisson truncation tail
    site, ro = _base(cfg)
    res["poisson_tail"] = []
    old = cmos.poisson_support.__defaults__
    for tail in (1e-10, 1e-13, 1e-16):
        cmos.poisson_support.__defaults__ = (tail,)
        r = _box_rates(site, ro, T)
        res["poisson_tail"].append({"tail": tail, "fp": r.fp, "fn": r.fn})
    cmos.poisson_support.__defaults__ = old
    return res, {}


def tables(res: dict) -> dict:
    return {
        "h_stat": md_table(["h / default", "h_stat (e-)", "FP", "FN", "eps"],
                           [[r["h_factor"], r["h_stat"], r["fp"], r["fn"], r["eps"]] for r in res["h_stat"]]),
        "ds": md_table(["ds", "FP", "FN"], [[r["ds"], r["fp"], r["fn"]] for r in res["ds"]]),
        "n_tau": md_table(["n_tau", "FP", "FN"], [[r["n_tau"], r["fp"], r["fn"]] for r in res["n_tau"]]),
        "supersample": md_table(["s", "max |frac - ref|", "centre pixel rel. dev."],
                                [[r["s"], r["max_abs_diff"], r["centre_rel"]] for r in res["supersample"]]),
        "n_quad": md_table(["n_quad", "max |frac - ref(512)|"], [[r["n_quad"], r["max_abs_diff"]] for r in res["n_quad"]]),
        "is": f"Deep FP threshold T = {res['is']['T']:.2f} e-, lattice FP = {res['is']['lattice']:.4e}\n\n" + md_table(
            ["samples", "IS estimate", "rel. error", "(IS - lattice)/sigma", "ESS"],
            [[r["n"], r["p"], r["rel_err"], r["dev_sigma"], r["ess"]] for r in res["is"]["rows"]]),
        "dt": md_table(["dt (us)", "survival", "s.e.", "T_eff (uK)"],
                       [[r["dt_us"], r["survival"], r["survival_se"], r["T_eff_uK"]] for r in res["dt"]]),
        "poisson_tail": md_table(["Poisson tail cut", "FP", "FN"], [[r["tail"], r["fp"], r["fn"]] for r in res["poisson_tail"]]),
    }
