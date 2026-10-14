"""V2: reproduce Picken, Legaie & Pritchard, arXiv:1709.00946v2 (Cs, Zyla 5.5 sCMOS, 40 ms).

Fixed: QE 0.22, read noise 2.2 e- rms, 40 ms, detection efficiency 3.5% (excludes QE),
effective scattering rate 450 photons/ms, 3x3 box of 1 um object-plane pixels.
Fitted per trap: ADC gain g (e-/ADU) and background lambda0 (both from the background
peak mean and width of Table II) and the relative spread of the atom signal (from the atom
peak width). Predicted: the atom signal in electrons (compared with R t eta QE), the
false-positive rate at the paper's threshold, and the Table I SNR values.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy import special

from ..atoms.loss_phenom import TauDistribution
from ..counts.pmf import Classifier, SiteModel, site_statistic
from ..fidelity import optimal_threshold, rates_at
from ..optics.array import PSFSpec, pixel_fractions
from ..sensor.cmos import PixelReadout, cmos_pixel
from ..sensor.emccd import EMCCDReadout, emccd_pixel
from .common import md_table, solve_monotone, v


def gaussian_eps(p0, mu0, s0, c):
    """Picken et al. Eq. (4): eps = p0/2 [1 - erf((c - mu0)/(sqrt(2) s0))]."""
    return 0.5 * p0 * special.erfc((c - mu0) / (np.sqrt(2) * s0))


def gaussian_acc(mu1, s1, c):
    """Picken et al. Eq. (5)."""
    return 0.5 * special.erfc((c - mu1) / (np.sqrt(2) * s1))


def fit_gain_bg(mu0, s0, sigma, npx):
    """Solve mean = lam0/g, var = (lam0 + npx sigma^2 + npx g^2/12)/g^2 for (g, lam0)."""
    a = s0**2 - npx / 12.0
    g = (mu0 + np.sqrt(mu0**2 + 4 * a * npx * sigma**2)) / (2 * a)
    return g, mu0 * g


def snr_eq1(n, qe, F2, read, M):
    """Picken et al. Eq. (1) with n_b = 0."""
    return n * qe / np.sqrt(F2 * qe * n + (read / M) ** 2)


def snr_from_chain(n, det):
    """SNR = mean/sd of the simulated single-pixel output for n incident photons (our sensor chain)."""
    qe, read, M = det["qe"], det["read_e"], det["M"]
    if M > 1:
        h = max(0.02, n * qe / 4000.0)
        ro = EMCCDReadout(em_gain=M, n_stages=536, cic=0.0, cic_register_fraction=0.0, sigma_ro_e=read,
                          register_model="erlang", h_fine=h)
        lat = emccd_pixel(n * qe, 0.0, ro)
    else:
        ro = PixelReadout(classes=[(1.0, [(1.0, 0.0, read)])], gain_e_per_adu=None, h_fine=0.02)
        lat = cmos_pixel(n * qe, 0.0, ro)
    return lat.mean() / np.sqrt(lat.var())


def run(cfg: dict, outdir: Path):
    F, num = cfg["fixed"], cfg["numerics"]
    sigma = v(F["sigma_read_e"])
    npx = int(v(F["box_px"])) ** 2
    mu_pred = v(F["rate_s"]) * v(F["t_exp_s"]) * v(F["eta_det"]) * v(F["qe"])
    eta_chain = v(F["eta_collection"]) * v(F["eta_ito"]) * v(F["eta_filters"])

    box = int(v(F["box_px"]))
    pitch = v(F["pixel_obj_um"]) * 1e-6
    w = v(F["psf_um"]) * 1e-6 / np.sqrt(2 * np.log(2))  # FWHM -> 1/e^2 radius
    fr = pixel_fractions(PSFSpec("gaussian", waist_m=w), pitch, box, supersample=num["supersample"])
    L_vac = v(F["t_exp_s"]) / v(F["vacuum_lifetime_s"])

    traps = {}
    for tname in ("trap1", "trap2"):
        d = cfg["data"][tname]
        g, lam0 = fit_gain_bg(d["mu0"], d["sigma0"], sigma, npx)
        ro = PixelReadout(classes=[(1.0, [(1.0, 0.0, sigma)])], gain_e_per_adu=g)
        bg = np.full((box, box), lam0 / npx)
        p0 = site_statistic(SiteModel(frac=fr, mu=0.0, bg=bg), ro, Classifier("box"), False)
        # model check of the fitted moments (exact ADC rounding in the lattice)
        m0_adu, s0_adu = p0.mean() / g, np.sqrt(p0.var()) / g
        mu_meas = (d["mu1"] - d["mu0"]) * g / (fr.sum() * np.exp(-0.5 * L_vac))
        tau = TauDistribution.from_hazard(L_vac, num["n_tau"])

        def atom_sd(rs):
            s = SiteModel(frac=fr, mu=mu_meas, bg=bg, tau=tau, rate_rel_sd=rs, ds=num["ds"])
            return np.sqrt(site_statistic(s, ro, Classifier("box"), True).var()) / g
        rs = solve_monotone(atom_sd, d["sigma1"], 0.0, 0.5, xtol=1e-4)
        rs = 0.0 if rs is None else rs
        site = SiteModel(frac=fr, mu=mu_meas, bg=bg, tau=tau, rate_rel_sd=rs, ds=num["ds"])
        p1 = site_statistic(site, ro, Classifier("box"), True)
        c_e = d["cmin"] * g
        at = rates_at(p0, p1, c_e, filling=d["p1"])
        opt = optimal_threshold(p0, p1, filling=d["p1"])
        # the paper's own Gaussian estimate (DERIVED from Table II)
        eps_g = gaussian_eps(1 - d["p1"], d["mu0"], d["sigma0"], d["cmin"])
        acc_g = gaussian_acc(d["mu1"], d["sigma1"], d["cmin"])
        traps[tname] = {
            "g_e_per_adu": g, "lambda0_e": lam0, "model_mu0_adu": m0_adu, "model_sigma0_adu": s0_adu,
            "signal_e_measured": mu_meas * fr.sum(), "signal_e_predicted": mu_pred * fr.sum(),
            "signal_ratio": mu_meas / mu_pred, "rate_rel_sd_fit": rs,
            "cmin_adu": d["cmin"], "cmin_e": c_e,
            "fp_model_at_cmin": at.fp, "eps_model_at_cmin": (1 - d["p1"]) * at.fp,
            "eps_gauss_paper_eq4": eps_g, "acceptance_model": 1 - at.fn, "acceptance_gauss_eq5": acc_g,
            "opt_threshold_adu": opt.threshold / g, "opt_fp": opt.fp, "opt_fn": opt.fn, "opt_eps": opt.eps,
            "p1": d["p1"],
        }
        traps[tname]["_p0"], traps[tname]["_p1"], traps[tname]["_g"] = p0, p1, g

    # ---- Table I SNR
    snr_rows = []
    for name, det in cfg["snr_table"]["detectors"].items():
        for col, qe_key in (("best", "qe_best"), ("at852", "qe_852")):
            for n, ref in zip(cfg["snr_table"]["levels"], det[col]):
                dd = {"qe": det[qe_key], "read_e": det["read_e"], "M": det["M"]}
                s_eq = snr_eq1(n, det[qe_key], det["F2"], det["read_e"], det["M"])
                s_ch = snr_from_chain(n, dd)
                snr_rows.append({"detector": name, "qe": det[qe_key], "n": n, "paper": ref, "eq1": s_eq, "chain": s_ch})

    arrays = {"frac": fr}
    for tname, t in traps.items():
        p0, p1, g = t.pop("_p0"), t.pop("_p1"), t.pop("_g")
        arrays[f"{tname}_p0_adu"] = p0.x / g
        arrays[f"{tname}_p0"] = p0.p
        arrays[f"{tname}_p1_adu"] = p1.x / g
        arrays[f"{tname}_p1"] = p1.p
    res = {"fixed": {"mu_pred_e": mu_pred, "eta_chain_product": eta_chain, "L_vac": L_vac,
                     "frac_in_box": float(fr.sum())},
           "traps": traps, "snr": snr_rows, "targets": {k: v(x) for k, x in cfg["targets"].items()}}
    return res, arrays


def tables(res: dict) -> dict:
    tr = res["traps"]
    fx = res["fixed"]
    fit = md_table(["quantity", "trap 1", "trap 2"], [
        ["detection efficiency 5.4% x 79% x 83% (DERIVED) vs 3.5% stated", f"{fx['eta_chain_product']:.4f}", "-"],
        ["fitted ADC gain g (e-/ADU)", tr["trap1"]["g_e_per_adu"], tr["trap2"]["g_e_per_adu"]],
        ["fitted background lambda0 (e-, 3x3)", tr["trap1"]["lambda0_e"], tr["trap2"]["lambda0_e"]],
        ["model mu0 / sigma0 (ADU) [paper 206/29, 191/29]",
         f"{tr['trap1']['model_mu0_adu']:.1f} / {tr['trap1']['model_sigma0_adu']:.2f}",
         f"{tr['trap2']['model_mu0_adu']:.1f} / {tr['trap2']['model_sigma0_adu']:.2f}"],
        ["atom signal, measured (e-)", tr["trap1"]["signal_e_measured"], tr["trap2"]["signal_e_measured"]],
        ["atom signal, predicted R t eta QE (e-)", tr["trap1"]["signal_e_predicted"], tr["trap2"]["signal_e_predicted"]],
        ["measured / predicted", tr["trap1"]["signal_ratio"], tr["trap2"]["signal_ratio"]],
        ["fitted relative signal spread", tr["trap1"]["rate_rel_sd_fit"], tr["trap2"]["rate_rel_sd_fit"]],
    ])
    err = md_table(["quantity", "trap 1", "trap 2"], [
        ["threshold c_min (ADU) [paper]", tr["trap1"]["cmin_adu"], tr["trap2"]["cmin_adu"]],
        ["eps, paper Eq. (4) Gaussian (DERIVED from Table II)", tr["trap1"]["eps_gauss_paper_eq4"], tr["trap2"]["eps_gauss_paper_eq4"]],
        ["eps = p0 FP, Poisson x Gaussian model (SIM)", tr["trap1"]["eps_model_at_cmin"], tr["trap2"]["eps_model_at_cmin"]],
        ["model / Gaussian", tr["trap1"]["eps_model_at_cmin"] / tr["trap1"]["eps_gauss_paper_eq4"],
         tr["trap2"]["eps_model_at_cmin"] / tr["trap2"]["eps_gauss_paper_eq4"]],
        ["acceptance, paper Eq. (5) (DERIVED from Table II)", f"{tr['trap1']['acceptance_gauss_eq5']:.7f}", f"{tr['trap2']['acceptance_gauss_eq5']:.7f}"],
        ["acceptance, model incl. vacuum loss at tau_vac = 6 s (SIM)", f"{tr['trap1']['acceptance_model']:.5f}", f"{tr['trap2']['acceptance_model']:.5f}"],
        ["optimal threshold (ADU), model", tr["trap1"]["opt_threshold_adu"], tr["trap2"]["opt_threshold_adu"]],
        ["optimal FP / FN, model", f"{tr['trap1']['opt_fp']:.2e} / {tr['trap1']['opt_fn']:.2e}",
         f"{tr['trap2']['opt_fp']:.2e} / {tr['trap2']['opt_fn']:.2e}"],
    ])
    rows = []
    for r in res["snr"]:
        rows.append([r["detector"], r["qe"], r["n"], r["paper"], r["eq1"], r["chain"], r["chain"] - r["paper"]])
    snr = md_table(["detector", "QE", "photons/px", "Table I", "Eq. (1)", "sensor chain (SIM)", "chain - Table I"], rows)
    return {"fit": fit, "errors": err, "snr": snr}
