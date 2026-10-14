"""V3: analytic limits.

(1) Zero read noise and (zero or Poisson) background: the engine reproduces pure Poisson
    discrimination, FP = P(Pois(B) >= T), FN = P(Pois(B + mu) < T), to rounding error.
(2) EMCCD excess noise: F^2 -> 2 at high gain (Erlang lattice), and the stage-by-stage Monte
    Carlo reproduces the finite-register formula of Basden et al. (2003) Eq. (5).
(3) Dipole collection: closed forms vs direct quadrature, NA -> 0 slopes and NA = 1 -> 1/2.
(4) Airy vs vectorial PSF at NA 0.3, 0.5, 0.8: FWHM and encircled-energy residuals.
(5) Model B: energy conservation without scattering; recoil heating rate 2 E_r R.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy import stats

from ..atoms.loss_randomwalk import RandomWalkParams, simulate_loss_times
from ..constants import CS133, KB
from ..counts.pmf import Classifier, SiteModel, site_statistic
from ..fidelity import rates_at
from ..optics.collection import collection_fraction, collection_fraction_numeric
from ..optics.psf_airy import airy_fwhm, encircled_energy, psf_airy, radial_fwhm
from ..optics.psf_vectorial import psf_vectorial
from ..sensor.cmos import PixelReadout
from ..sensor.emccd import EMCCDReadout, emccd_pixel, enf2_finite, register_mc
from .common import md_table


def poisson_limit(cases):
    rows = []
    for mu, B, T, box in cases:
        fr = np.zeros((box, box))
        fr[box // 2, box // 2] = 1.0
        bg = np.full((box, box), B / box**2)
        ro = PixelReadout(classes=[(1.0, [(1.0, 0.0, 0.0)])], gain_e_per_adu=1.0)
        p0 = site_statistic(SiteModel(frac=fr, mu=mu, bg=bg), ro, Classifier("box"), False)
        p1 = site_statistic(SiteModel(frac=fr, mu=mu, bg=bg), ro, Classifier("box"), True)
        r = rates_at(p0, p1, T)
        fp_a = stats.poisson.sf(T - 1, B) if B > 0 else 0.0
        fn_a = stats.poisson.cdf(T - 1, B + mu)
        rows.append({"mu": mu, "B": B, "T": T, "box": box, "fp": r.fp, "fp_exact": fp_a, "fn": r.fn, "fn_exact": fn_a,
                     "rel_fp": abs(r.fp - fp_a) / fp_a if fp_a > 0 else abs(r.fp),
                     "rel_fn": abs(r.fn - fn_a) / fn_a})
    return rows


def emccd_enf(cfg, seed):
    rows = []
    lam = cfg["lambda_in"]
    for g in cfg["gains"]:
        ro = EMCCDReadout(em_gain=g, n_stages=cfg["n_stages"], cic=0.0, cic_register_fraction=0.0, sigma_ro_e=0.0,
                          register_model="erlang", h_fine=0.01)
        lat = emccd_pixel(lam, 0.0, ro)
        f2_erlang = lat.var() / lam
        ro_s = EMCCDReadout(em_gain=g, n_stages=cfg["n_stages"], cic=0.0, cic_register_fraction=0.0, sigma_ro_e=0.0,
                            register_model="stagewise")
        rng = np.random.default_rng(seed + int(g))
        n_in = rng.poisson(lam, cfg["mc_samples"])
        out = register_mc(n_in, ro_s, seed + 7 * int(g)) / g
        f2_mc = out.var(ddof=1) / lam
        # standard error of a variance ratio estimated from N samples (delta method, ~ sqrt(2/N) relative for near-Gaussian)
        se = f2_mc * np.sqrt((stats.kurtosis(out, fisher=False) - 1) / out.size)
        rows.append({"gain": g, "n_stages": cfg["n_stages"], "f2_formula_eq5": enf2_finite(g, cfg["n_stages"]),
                     "f2_erlang_lattice": f2_erlang, "f2_stagewise_mc": f2_mc, "f2_mc_se": se,
                     "mean_ratio_mc": out.mean() / lam})
    return rows


def collection_limits():
    rows = []
    for pat, slope_exp, k in (("isotropic", 1 / 4, 2), ("sigma", 3 / 8, 2), ("pi", 3 / 16, 4)):
        na_small = 1e-3
        slope = collection_fraction(na_small, pat) / na_small**k
        rows.append({"pattern": pat, "f(NA=1)": float(collection_fraction(1.0, pat)),
                     "small-NA coefficient": float(slope), "expected": slope_exp, "exponent": k,
                     "max |closed - quadrature| (NA 0.1..0.95)": max(
                         abs(collection_fraction(na, pat) - collection_fraction_numeric(na, pat, 20001))
                         for na in np.linspace(0.1, 0.95, 18))})
    return rows


def psf_compare(nas, lam=852.347e-9, n_quad=200):
    rows = []
    for na in nas:
        fv = lambda r, na=na: psf_vectorial(r, lam, na, "isotropic", n_quad)
        fa = lambda r, na=na: psf_airy(r, lam, na)
        r0 = 0.61 * lam / na
        rows.append({"NA": na, "FWHM_airy_nm": airy_fwhm(lam, na) * 1e9, "FWHM_vec_nm": radial_fwhm(fv, 2 * lam / na) * 1e9,
                     "EE_airy(r0)": encircled_energy(fa, r0), "EE_vec(r0)": encircled_energy(fv, r0),
                     "EE_airy(2r0)": encircled_energy(fa, 2 * r0), "EE_vec(2r0)": encircled_energy(fv, 2 * r0)})
    return rows


def model_b_checks(seed):
    m, k = float(CS133.mass_kg), CS133.k
    Er = CS133.recoil_energy_J
    base = dict(mass_kg=m, k_m=k, U0_K=1e-3, w0_m=0.8e-6, zR_m=2.4e-6, rate_s=0.0, gamma_c_s=0.0, T_eq_K=20e-6,
                T_init_K=20e-6, t_exp_s=0.01, dt_s=0.2e-6)
    r0 = simulate_loss_times(RandomWalkParams(**base), 2000, seed)
    r_short = simulate_loss_times(RandomWalkParams(**{**base, "t_exp_s": 2e-6}), 2000, seed)  # ~initial energy, same atoms
    R, t = 2e4, 0.02
    r1 = simulate_loss_times(RandomWalkParams(**{**base, "rate_s": R, "t_exp_s": t}), 4000, seed + 1)
    r1s = simulate_loss_times(RandomWalkParams(**{**base, "rate_s": R, "t_exp_s": 2e-6}), 4000, seed + 1)
    return {"energy_drift_rel": (r0["E_mean_K"] - r_short["E_mean_K"]) / r_short["E_mean_K"],
            "heating_uK": (r1["E_mean_K"] - r1s["E_mean_K"]) * 1e6, "heating_expected_uK": 2 * R * t * Er / KB * 1e6,
            "survival_no_scatter": r0["survival"]}


def run(cfg: dict, outdir: Path):
    res = {"poisson": poisson_limit([tuple(c) for c in cfg["poisson_cases"]]),
           "emccd": emccd_enf(cfg["emccd"], cfg["seed"]),
           "collection": collection_limits(),
           "psf": psf_compare(cfg["psf_na"]),
           "model_b": model_b_checks(cfg["seed"])}
    return res, {"psf_na": np.asarray(cfg["psf_na"], float)}


def tables(res: dict) -> dict:
    p = md_table(["mu", "B", "T", "box", "FP engine", "FP exact", "FN engine", "FN exact", "max rel. dev."],
                 [[r["mu"], r["B"], r["T"], r["box"], r["fp"], r["fp_exact"], r["fn"], r["fn_exact"],
                   max(r["rel_fp"], r["rel_fn"])] for r in res["poisson"]])
    e = md_table(["gain", "stages r", "F^2 Eq. (5)", "F^2 Erlang lattice", "F^2 stage MC", "MC s.e.", "(MC - Eq.5)/s.e."],
                 [[r["gain"], r["n_stages"], r["f2_formula_eq5"], r["f2_erlang_lattice"], r["f2_stagewise_mc"],
                   r["f2_mc_se"], (r["f2_stagewise_mc"] - r["f2_formula_eq5"]) / r["f2_mc_se"]] for r in res["emccd"]])
    c = md_table(["pattern", "f(NA = 1)", "small-NA coefficient", "expected", "NA exponent", "max |closed - quadrature|"],
                 [[r["pattern"], r["f(NA=1)"], r["small-NA coefficient"], r["expected"], r["exponent"],
                   r["max |closed - quadrature| (NA 0.1..0.95)"]] for r in res["collection"]])
    ps = md_table(["NA", "FWHM Airy (nm)", "FWHM vectorial (nm)", "ratio", "EE(r0) Airy", "EE(r0) vect.", "EE(2r0) Airy", "EE(2r0) vect."],
                  [[r["NA"], r["FWHM_airy_nm"], r["FWHM_vec_nm"], r["FWHM_vec_nm"] / r["FWHM_airy_nm"], r["EE_airy(r0)"],
                    r["EE_vec(r0)"], r["EE_airy(2r0)"], r["EE_vec(2r0)"]] for r in res["psf"]])
    b = res["model_b"]
    mb = md_table(["check", "value", "expected"], [["relative energy drift, no scattering (10 ms)", b["energy_drift_rel"], "0"],
                                                   ["survival, no scattering", b["survival_no_scatter"], "1"],
                                                   ["heating in 20 ms at R = 2e4 s^-1 (uK)", b["heating_uK"], b["heating_expected_uK"]]])
    return {"poisson": p, "emccd": e, "collection": c, "psf": ps, "model_b": mb}
