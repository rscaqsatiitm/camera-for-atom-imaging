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
