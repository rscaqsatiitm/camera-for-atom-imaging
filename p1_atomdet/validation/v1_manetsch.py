"""V1: reproduce Manetsch et al., arXiv:2403.12021v4 (Cs, ORCA-Quest, 80 ms, F = 0.9999374(8)).

Fixed from the paper: QE 0.44, 2.7% of scattered photons to the camera (excludes QE), 80 ms,
filling 0.512, 7x7 box, measured PSF waist 1.7 px, magnification 5.1, model-free survival
S = 0.999864(2) (sets the model-A hazard).
Taken from the paper's figures (DATA, digitised): background and atom histogram moments
(Ext. Data Fig. 4b), site-to-site scattering-rate spread (Ext. Data Fig. 4d).
Fitted: background mean lambda0, atom signal mu (both from the histogram means) and the
read noise sigma (from the background width). Everything else is a prediction.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from camsim_common import SensorModel
from camsim_common.provenance import val

from ..atoms.loss_phenom import TauDistribution, lossy_poisson_pmf, manetsch_lossy_poisson_pmf
from ..atoms.loss_randomwalk import RandomWalkParams, loss_time_shape, simulate_loss_times
from ..constants import CS133
from ..counts.importance import is_fn, is_fp
from ..counts.lattice import mixture
from ..counts.pmf import Classifier, SiteModel, bayes_llr_tables, optimal_linear_weights, site_statistic
from ..fidelity import optimal_threshold
from ..optics.array import PSFSpec, pixel_fractions
from ..sensor.cmos import PixelReadout
from .common import ROOT, md_table, solve_monotone, v


def _readout(sigma, gain):
    return PixelReadout(classes=[(1.0, [(1.0, 0.0, sigma)])], gain_e_per_adu=gain)


def _site(fr, mu, lam0, tau, rs, ds, rate_hist=None):
    return SiteModel(frac=fr, mu=mu, bg=np.full(fr.shape, lam0 / fr.size), tau=tau, rate_rel_sd=rs, ds=ds,
                     rate_hist=rate_hist)


def _rates(site, ro, kind, filling):
    if kind == "box":
        clf = Classifier("box")
    elif kind == "matched":
        clf = Classifier("matched", weights=optimal_linear_weights(site, ro))
    else:
        clf = Classifier("bayes", llr=bayes_llr_tables(site, ro))
    p0 = site_statistic(site, ro, clf, False)
    p1 = site_statistic(site, ro, clf, True)
    return optimal_threshold(p0, p1, filling), p0, p1, clf


def model_b_shapes(cfg_b: dict, seed: int, n_atoms: int):
    """Loss-time shapes from model B in two regimes with measurable loss (shape transferred, report Sec. 2.1)."""
    out = {}
    m, k = float(CS133.mass_kg), CS133.k
    for name, kw in cfg_b["regimes"].items():
        prm = {kk: v(vv) for kk, vv in {**cfg_b["base"], **kw}.items()}
        prm["cooled_axes"] = tuple(prm["cooled_axes"])
        p = RandomWalkParams(mass_kg=m, k_m=k, **prm)
        r = simulate_loss_times(p, n_atoms, seed)
        c, h, nlost = loss_time_shape(r["tau"], cfg_b.get("n_bins", 20))
        out[name] = {"centres": c, "density": h, "n_lost": nlost, "survival": r["survival"],
                     "T_eff_K": r["T_eff_K"], "params": prm}
    return out


def run(cfg: dict, outdir: Path) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    sm = SensorModel.from_yaml(ROOT / cfg["sensor_preset"])
    F, d, num = cfg["fixed"], cfg["data"], cfg["numerics"]
    filling = v(F["filling"])
    t_exp = v(F["t_exp_s"])
    qe, eta = v(F["qe"]), v(F["eta_cam"])
    gain = float(val(sm.adc.gain_e_per_adu))

    # ---- geometry: measured Gaussian PSF (waist 1.7 px) on 4.6 um pixels, magnification 5.1
    pitch_obj = float(val(sm.pixel_pitch_um)) * 1e-6 / v(F["magnification"])
    psf = PSFSpec("gaussian", waist_m=v(F["psf_waist_px"]) * pitch_obj)
    box = int(v(F["box_px"]))
    fr = pixel_fractions(psf, pitch_obj, box, supersample=num["supersample"])

    # ---- loss: model A hazard from the model-free imaging survival
    S = v(F["survival_imaging"])
    L = -np.log(S)
    tauA = TauDistribution.from_hazard(L, num["n_tau"])

    # ---- fits
    lam0 = v(d["bg_mean"])
    sd_bg = v(d["bg_sd"])
    def bg_sd(sig):
        st = SiteModel(frac=fr, mu=0.0, bg=np.full(fr.shape, lam0 / fr.size))
        return site_statistic(st, _readout(sig, gain), Classifier("box"), False).var() ** 0.5
    sigma_fit = solve_monotone(bg_sd, sd_bg, 0.08, 0.8, xtol=1e-3)
    sigma_note = "fitted to background width"
    if sigma_fit is None:
        sigma_fit, sigma_note = 0.08, "background width at or below the Poisson limit; sigma set to lower bound 0.08 e-"
    mu = (v(d["atom_mean"]) - lam0) / (fr.sum() * tauA.mean_tau())
    rs = v(d["rate_sd_kcps"]) / v(d["rate_mean_kcps"])
    rh = np.loadtxt(ROOT / "data" / "manetsch_ed4d_rate_hist_digitised.csv", delimiter=",", skiprows=4)
    rate_hist = (rh[:, 0], rh[:, 2])
    R_implied = mu / (t_exp * eta * qe)
    fitted = {"lambda0_e_per_box": lam0, "mu_e": mu, "sigma_read_e": sigma_fit, "sigma_note": sigma_note}

    # ---- predictions
    variants = {
        "baseline": dict(sig=sigma_fit, g=gain, rs=rs, tau=tauA, rh=rate_hist),
        "gaussian_rate": dict(sig=sigma_fit, g=gain, rs=rs, tau=tauA, rh=None),
        "no_rate_spread": dict(sig=sigma_fit, g=gain, rs=0.0, tau=tauA, rh=None),
        "sigma_0.43_spec": dict(sig=0.43, g=gain, rs=rs, tau=tauA, rh=rate_hist),
        "pnr_rounding": dict(sig=sigma_fit, g=1.0, rs=rs, tau=tauA, rh=rate_hist),
        "no_loss": dict(sig=sigma_fit, g=gain, rs=rs, tau=TauDistribution.no_loss(), rh=rate_hist),
    }
    shapes = model_b_shapes(cfg["model_b"], cfg["seed"], cfg["model_b"]["n_atoms"])
    for name, sh in shapes.items():
        variants[f"model_B_{name}"] = dict(sig=sigma_fit, g=gain, rs=rs, rh=rate_hist,
                                          tau=TauDistribution.with_shape(sh["centres"], sh["density"], S, name))
    rows, store = [], {}
    for name, p in variants.items():
        ro = _readout(p["sig"], p["g"])
        site = _site(fr, mu, lam0, p["tau"], p["rs"], num["ds"], p["rh"])
        kinds = ("box", "matched", "bayes") if name in ("baseline", "no_loss") else ("box",)
        for kind in kinds:
            r, p0, p1, clf = _rates(site, ro, kind, filling)
            rows.append({"variant": name, "classifier": kind, **r.as_dict(), "atom_sd": p1.var() ** 0.5})
            if name == "baseline" and kind == "box":
                store["mix_p0"], store["mix_p1"] = p0, p1
            if name == "gaussian_rate" and kind == "box":
                store["p0"], store["p1"], store["site"], store["ro"], store["r"] = p0, p1, site, ro, r

    # ---- independent cross-check of the Gaussian-rate box rates by tilted importance sampling
    rng = np.random.default_rng(cfg["seed"])
    T = store["r"].threshold
    isfp = is_fp(store["site"], store["ro"], None, T, num["is_samples"], rng)
    isfn = is_fn(store["site"], store["ro"], None, T, num["is_samples"], rng)

    # ---- lossy-Poisson formula of the paper vs our tau mixture (pure Poisson, same lambda0, lambda1, L)
    n = np.arange(0, 250)
    lam1 = mu * fr.sum()
    pm = manetsch_lossy_poisson_pmf(n, lam0, lam1, L)
    po = lossy_poisson_pmf(n, lam1, lam0, tauA)
    formula_maxdiff = float(np.max(np.abs(pm - po)))
    cdf_m, cdf_o = np.cumsum(pm), np.cumsum(po)
    formula_fn_diff = float(abs(cdf_m[int(T)] - cdf_o[int(T)]))

    # ---- histogram shape vs digitised inset (log density per photon)
    mix = mixture([store["mix_p0"], store["mix_p1"]], [1 - filling, filling])
    dig = np.loadtxt(ROOT / "data" / "manetsch_ed4b_inset_digitised.csv", delimiter=",", skiprows=5)
    model_logd = []
    for x, _ in dig:
        sel = (mix.x >= x - 0.5) & (mix.x < x + 0.5)
        model_logd.append(np.log10(max(mix.p[sel].sum(), 1e-300)))
    model_logd = np.array(model_logd)

    res = {
        "fitted": fitted,
        "fixed": {"qe": qe, "eta_cam": eta, "t_exp_s": t_exp, "filling": filling, "box": box,
                  "pitch_obj_um": pitch_obj * 1e6, "frac_in_box": float(fr.sum()), "L": L, "S": S,
                  "rate_rel_sd_data": rs, "adc_gain_placeholder": gain},
        "implied_rate_kcps": R_implied / 1e3, "rate_mean_data_kcps": v(d["rate_mean_kcps"]),
        "rows": rows,
        "is_check": {"T": T, "lattice_fp": store["r"].fp, "is_fp": isfp.p, "is_fp_rel": isfp.rel_err,
                     "lattice_fn": store["r"].fn, "is_fn": isfn.p, "is_fn_rel": isfn.rel_err,
                     "ess_fp": isfp.ess, "ess_fn": isfn.ess},
        "formula_check": {"max_abs_pmf_diff": formula_maxdiff, "fn_cdf_diff_at_T": formula_fn_diff},
        "model_b": {k: {kk: vv for kk, vv in s.items() if kk in ("n_lost", "survival", "T_eff_K", "params")}
                    for k, s in shapes.items()},
        "targets": {k: v(x) for k, x in cfg["targets"].items()},
        "data": {k: v(x) for k, x in d.items()},
    }
    arrays = {"inset_x": dig[:, 0], "inset_logd_data": dig[:, 1], "inset_logd_model": model_logd,
              "mix_x": mix.x, "mix_p": mix.p, "p0_x": store["mix_p0"].x, "p0_p": store["mix_p0"].p,
              "p1_x": store["mix_p1"].x, "p1_p": store["mix_p1"].p, "frac": fr,
              "rate_hist_kcps": rh[:, 0], "rate_hist_w": rh[:, 2]}
    for k, s in shapes.items():
        arrays[f"modelB_{k}_centres"] = s["centres"]
        arrays[f"modelB_{k}_density"] = s["density"]
    return res, arrays


def tables(res: dict) -> dict:
    t = res["targets"]
    rows = []
    for r in res["rows"]:
        rows.append([r["variant"], r["classifier"], f"{r['threshold']:.2f}", r["fp"], r["fn"], f"{r['fidelity']:.7f}",
                     f"{r['atom_sd']:.2f}" if r["classifier"] != "bayes" else "-"])
    pred = md_table(["variant", "classifier", "threshold", "FP", "FN", "fidelity", "atom sd (e-)"], rows)
    base = next(r for r in res["rows"] if r["variant"] == "baseline" and r["classifier"] == "box")
    resid = md_table(["quantity", "measured (SOURCE)", "model (box, baseline)", "model/measured"], [
        ["FP", t["fp"], base["fp"], base["fp"] / t["fp"]],
        ["FN", t["fn"], base["fn"], base["fn"] / t["fn"]],
        ["infidelity", 1 - t["fidelity"], base["eps"], base["eps"] / (1 - t["fidelity"])],
        ["fidelity", f"{t['fidelity']:.7f}", f"{base['fidelity']:.7f}", "-"],
        ["atom-peak sd (e-, DATA)", res["data"]["atom_sd"], base["atom_sd"], base["atom_sd"] / res["data"]["atom_sd"]],
        ["scattering rate (kc/s, DATA mean vs implied)", res["rate_mean_data_kcps"], res["implied_rate_kcps"],
         res["implied_rate_kcps"] / res["rate_mean_data_kcps"]],
    ])
    nl = next(r for r in res["rows"] if r["variant"] == "no_loss" and r["classifier"] == "box")
    resid += "\n" + md_table(["no-loss fidelity", "paper (lossy-Poisson fit, loss set to 0)", "model (box)"],
                              [["", f"{t['fidelity_no_loss']:.5f}(1)", f"{nl['fidelity']:.6f}"]])
    ic = res["is_check"]
    isc = md_table(["Gaussian-rate box, T = %.2f e-" % ic["T"], "lattice (SIM)", "importance sampling (SIM)", "IS rel. error", "ESS"], [
        ["FP", ic["lattice_fp"], ic["is_fp"], ic["is_fp_rel"], ic["ess_fp"]],
        ["FN", ic["lattice_fn"], ic["is_fn"], ic["is_fn_rel"], ic["ess_fn"]]])
    fc = res["formula_check"]
    form = md_table(["check", "value"], [["max |P_paper(n) - P_ours(n)|, pure Poisson", fc["max_abs_pmf_diff"]],
                                         ["|CDF difference| at threshold", fc["fn_cdf_diff_at_T"]]])
    fx = res["fitted"]
    fit = md_table(["parameter", "value", "how obtained"], [
        ["lambda0 (e-, 7x7 box)", fx["lambda0_e_per_box"], "fitted: background mean (DATA)"],
        ["mu (e-, atom present all exposure)", fx["mu_e"], "fitted: atom minus background mean (DATA)"],
        ["sigma_read (e- rms)", fx["sigma_read_e"], fx["sigma_note"]],
        ["rate distribution", "Fig. 4d histogram", "DATA (not fitted); Gaussian with rel. sd %.4f as variant" % res["fixed"]["rate_rel_sd_data"]],
        ["loss hazard L", res["fixed"]["L"], "DERIVED: -ln S, S = 0.999864 (SOURCE)"],
    ])
    return {"fit": fit, "predictions": pred, "residuals": resid, "is_check": isc, "formula_check": form}
