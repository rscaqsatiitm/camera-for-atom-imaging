"""Sweep evaluators: map one parameter point to scalar results.

`generic_point` covers the generic (QE, read-noise median, tail fraction, pixels per FWHM,
background, dark, imaging time, NA, loss) space with a CMOS chain and box / matched-filter
classifiers. It is used by the reduced local sweep and the resume test; production sweeps
(after the validation check-in) extend it with EMCCD, rolling shutter and imaging-time search.
"""
from __future__ import annotations

import numpy as np

from ..atoms.loss_phenom import TauDistribution
from ..constants import SPECIES
from ..counts.pmf import Classifier, SiteModel, optimal_linear_weights, site_statistic
from ..fidelity import optimal_threshold
from ..optics.array import PSFSpec, pitch_for_pixels_per_fwhm, pixel_fractions
from ..optics.collection import collection_fraction
from ..sensor.cmos import PixelReadout


def generic_point(p: dict, seed: int) -> dict:
    sp = SPECIES[p.get("species", "Cs133")]
    na = p["na"]
    t = p["t_img_s"]
    eta_coll = float(collection_fraction(na, p.get("dipole_pattern", "isotropic")))
    mu = p["rate_s"] * t * eta_coll * p["eta_opt"] * p["qe"]
    psf = PSFSpec("vectorial" if p.get("psf", "vectorial") == "vectorial" else "airy",
                  wavelength_m=float(sp.wavelength_m), na=na)
    pitch = pitch_for_pixels_per_fwhm(psf, p["px_per_fwhm"])
    box = int(p.get("box", 5))
    fr = pixel_fractions(psf, pitch, box, supersample=int(p.get("supersample", 6)))
    bg = np.full(fr.shape, p["bg_e_per_px"] + p.get("dark_e_per_px_s", 0.0) * t)
    L = p.get("loss_per_photon", 0.0) * p["rate_s"] * t + t / p.get("tau_vac_s", np.inf)
    tau = TauDistribution.from_hazard(L, 32)
    s = p["read_median_e"]
    f = p.get("tail_fraction", 0.0)
    A = p.get("tail_amplitude_e", 3.0 * s + 1.0)
    classes = [(1 - f, [(1.0, 0.0, s)])]
    if f > 0:
        classes.append((f, [(0.5, 0.0, s), (0.25, A, s), (0.25, -A, s)]))
    ro = PixelReadout(classes=classes, gain_e_per_adu=p.get("gain_e_per_adu", None), h_fine=p.get("h_fine", 0.05))
    site = SiteModel(frac=fr, mu=mu, bg=bg, tau=tau, rate_rel_sd=p.get("rate_rel_sd", 0.0), ds=p.get("ds", 0.02))
    out = {"mu_e": mu, "eta_coll": eta_coll, "frac_in_box": float(fr.sum()), "survival": float(np.exp(-L))}
    for kind in p.get("classifiers", ["box", "matched"]):
        clf = Classifier("box") if kind == "box" else Classifier("matched", weights=optimal_linear_weights(site, ro))
        r = optimal_threshold(site_statistic(site, ro, clf, False), site_statistic(site, ro, clf, True),
                              p.get("filling", 0.5))
        out.update({f"{kind}_eps": r.eps, f"{kind}_fp": r.fp, f"{kind}_fn": r.fn, f"{kind}_T": r.threshold})
    return out


def toy_point(p: dict, seed: int) -> dict:
    """Cheap evaluator for tests of the sweep machinery."""
    rng = np.random.default_rng(seed)
    return {"y": p["a"] ** 2 + p.get("b", 0.0), "noise": float(rng.random())}


# ------------------------------------------------------------------ production stage
def _with_file(p: dict) -> dict:
    """Merge a referenced scenario file (configs/production/*.yaml) under the point's own keys."""
    if "scenario_file" not in p:
        return p
    import yaml
    from ..validation.common import ROOT
    f = yaml.safe_load((ROOT / p["scenario_file"]).read_text())
    return {**f, **p}


def _scenario(p: dict):
    from ..production import Scenario
    from ..validation.common import v
    sc = Scenario()
    for k, x in p.get("scenario", {}).items():
        setattr(sc, k, v(x))
    for k in ("rate_s", "na", "p_loss_per_photon", "bg_rate_e_per_px_s", "rate_rel_sd"):
        if k in p:
            setattr(sc, k, p[k])
    if "bg_scale" in p:
        sc.bg_rate_e_per_px_s *= p["bg_scale"]
    return sc


def _times(sc, cam, pc, p, classifier, geo):
    from ..production import eps_any, time_to_any
    out = {}
    for target in p.get("targets", [1e-3, 1e-4]):
        t, _ = time_to_any(sc, cam, target, pc, classifier, geo=geo, t_lo=p.get("t_range_s", [1e-3, 0.3])[0],
                           t_hi=p.get("t_range_s", [1e-3, 0.3])[1], min_survival=p.get("min_survival", 0.99))
        out[f"{classifier}_t_{target:.0e}"] = t if t is not None else float("nan")
    for tt in p.get("eps_at_ms", [20, 80]):
        r = eps_any(sc, cam, geo, tt / 1000, pc, classifier)
        out.update({f"{classifier}_eps_{tt}ms": r.eps, f"{classifier}_fp_{tt}ms": r.fp, f"{classifier}_fn_{tt}ms": r.fn})
    return out


def production_generic_point(p: dict, seed: int) -> dict:
    """Generic CMOS-class camera in the reference scenario: QE, read-noise median and tail, sampling, background."""
    from ..production import Camera, geometry, mu_per_s, read_classes
    p = _with_file(p)
    sc = _scenario(p)
    cam = Camera(qe=p["qe"], read_classes=read_classes(p["read_median_e"], p.get("tail_fraction", 0.0),
                                                         p.get("tail_factor", 5.0)),
                 pixels_per_fwhm=p["px_per_fwhm"], dark_e_per_px_s=p.get("dark_e_per_px_s", 0.0))
    geo = geometry(sc, cam)
    out = {"box_px": geo.box, "mu_e_per_ms": mu_per_s(sc, cam) / 1000}
    for c in p.get("classifiers", ["matched", "box"]):
        out.update(_times(sc, cam, None, p, c, geo))
    return out


def production_preset_point(p: dict, seed: int) -> dict:
    """A SensorModel preset at a given pixel sampling in the reference scenario."""
    from camsim_common.schema import SensorModel
    from ..production import camera_from_sensor, geometry
    from ..validation.common import ROOT
    p = _with_file(p)
    sc = _scenario(p)
    sm = SensorModel.from_yaml(ROOT / "camsim_common" / "presets" / f"{p['preset']}.yaml")
    kw = {}
    if "cic_register_fraction" in p:
        kw["cic_register_fraction"] = p["cic_register_fraction"]
    if "pc_frame_rate_hz" in p:
        kw["pc_frame_rate_hz"] = p["pc_frame_rate_hz"]
    cam, pc = camera_from_sensor(sm, p["px_per_fwhm"], **kw)
    if "read_tail_fraction" in p and p["read_tail_fraction"] > 0:
        from ..production import read_classes
        med = cam.read_classes[0][1][0][2]
        cam.read_classes = read_classes(med, p["read_tail_fraction"], p.get("tail_factor", 5.0))
    if "qe_override" in p:
        cam.qe = p["qe_override"]
    geo = geometry(sc, cam)
    out = {"box_px": geo.box, "qe": cam.qe, "dark_e_per_px_s": cam.dark_e_per_px_s,
           "pitch_obj_um": geo.pitch_m * 1e6, "magnification": float(sm.pixel_pitch_um.value) / (geo.pitch_m * 1e6)}
    classifiers = ["box"] if pc is not None else p.get("classifiers", ["matched", "box"])
    for c in classifiers:
        out.update(_times(sc, cam, pc, p, c, geo))
    return out


SENSITIVITY_CASES = {
    "base_bsi": ("prime_bsi_express_cms", {}),
    "bsi_tail2pct": ("prime_bsi_express_cms", {"read_tail_fraction": 0.02}),
    "bsi_qe_lo": ("prime_bsi_express_cms", {"qe_override": 0.47}),
    "bsi_qe_hi": ("prime_bsi_express_cms", {"qe_override": 0.50}),
    "bsi_loss_high": ("prime_bsi_express_cms", {"loss": "high"}),
    "base_q2pnr": ("orca_quest2_uq_pnr", {}),
    "q2pnr_loss_high": ("orca_quest2_uq_pnr", {"loss": "high"}),
    "ixon_pc_freg0": ("ixon888_exf_pc", {"cic_register_fraction": 0.0}),
    "ixon_pc_freg1": ("ixon888_exf_pc", {"cic_register_fraction": 1.0}),
    "ixon_an_freg0": ("ixon888_exf_analog", {"cic_register_fraction": 0.0}),
    "ixon_an_freg1": ("ixon888_exf_analog", {"cic_register_fraction": 1.0}),
    "base_bfs": ("blackfly_s_bfs_u3_17s7m_hcg", {}),
    "bfs_loss_high": ("blackfly_s_bfs_u3_17s7m_hcg", {"loss": "high"}),
}


def production_sensitivity_point(p: dict, seed: int) -> dict:
    from ..validation.common import v
    preset, change = SENSITIVITY_CASES[p["case"]]
    q = {**p, "preset": preset, **{k: x for k, x in change.items() if k != "loss"}}
    if change.get("loss") == "high":
        q["p_loss_per_photon"] = v(p["p_loss_high"])
    return production_preset_point(q, seed)


SOBOL_LOG10_CAP_S = 0.0   # unreached targets are censored at t = 1 s (log10 t = 0)


def sobol_point(p: dict, seed: int) -> dict:
    """Generic CMOS-class point for the Sobol analysis: matched filter, t(1e-3) and t(1e-4) from one scan."""
    from ..production import Camera, eps_at, geometry, read_classes, survival, times_to_targets
    p = _with_file(p)
    sc = _scenario(p)
    sc.n_offsets = int(p.get("n_offsets", 2))
    sc.offset_rule = p.get("offset_rule", "midpoint")
    sc.ds = float(p.get("ds", 0.02))
    cam = Camera(qe=p["qe"], read_classes=read_classes(p["read_median_e"], p.get("tail_fraction", 0.0),
                                                         p.get("tail_factor", 5.0)),
                 pixels_per_fwhm=p["px_per_fwhm"], dark_e_per_px_s=p.get("dark_e_per_px_s", 0.0))
    geo = geometry(sc, cam)
    targets = p.get("targets", [1e-3, 1e-4])
    res, (ts, eps) = times_to_targets(lambda t: eps_at(sc, cam, geo, t, "matched").eps, targets,
                                      t_lo=p.get("t_range_s", [1e-3, 1.0])[0], t_hi=p.get("t_range_s", [1e-3, 1.0])[1],
                                      survival_of_t=lambda t: survival(sc, t), min_survival=p.get("min_survival", 0.99))
    out = {"box_px": geo.box, "eps_min_scan": float(np.min(eps))}
    for tg in targets:
        t = res[tg]
        out[f"t_{tg:.0e}"] = float("nan") if t is None else t
        out[f"log10t_{tg:.0e}"] = SOBOL_LOG10_CAP_S if t is None else float(np.log10(t))
        out[f"reached_{tg:.0e}"] = 0.0 if t is None else 1.0
    return out
