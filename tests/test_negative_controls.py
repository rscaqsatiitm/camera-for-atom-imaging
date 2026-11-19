"""Negative controls: a deliberately wrong threshold must fail each classifier's fidelity check.

Design fixed on 2026-11-19 before the first run (CLAUDE.md rules 10 to 12):
  - The fidelity check is the existing comparison of the reported error rates with an independent
    route: brute-force MC for the lattice engine (tests/test_counts.py), the lattice engine for the
    CF engine used in production (tests/test_cf_production.py). Same sites, readouts, seeds, sample
    sizes and tolerances as those tests.
  - The rates are reported at the optimal threshold T_opt, but the classifier applies
    T_wrong = T_opt +/- DELTA_FRAC x (mean score with atom - mean score without atom),
    DELTA_FRAC = 0.10 (tests/mc_compare.py). Both signs must fail the check.
  - A clean twin at T_opt, with the same random samples, must pass, so that a failure is caused by
    the threshold and not by the comparator rejecting everything.
"""
import numpy as np
import pytest

from p1_atomdet.atoms.loss_phenom import TauDistribution
from p1_atomdet.counts.cf import CFNoise, cf_pair, matched_weights_analytic
from p1_atomdet.counts.mc import mc_error_rates
from p1_atomdet.counts.pmf import Classifier, SiteModel, bayes_llr_tables, optimal_linear_weights, site_statistic
from p1_atomdet.fidelity import optimal_threshold, rates_at
from p1_atomdet.optics.array import PSFSpec, pixel_fractions
from p1_atomdet.sensor.cmos import PixelReadout
from p1_atomdet.sensor.emccd import EMCCDReadout

from mc_compare import DELTA_FRAC, mc_agrees, wrong_thresholds

# CF-versus-lattice tolerance of tests/test_cf_production.py::test_cf_matches_lattice_cmos.
CF_REL_TOL = 0.02


def _site(mu=22.0, bg=0.4, L=2e-3, rs=0.1, box=5, ds=0.01):
    """Same constructor as tests/test_counts.py::_site."""
    fr = pixel_fractions(PSFSpec("gaussian", waist_m=1.7 * 0.9e-6), 0.9e-6, box)
    return SiteModel(frac=fr, mu=mu, bg=np.full((box, box), bg), tau=TauDistribution.from_hazard(L, 48),
                     rate_rel_sd=rs, ds=ds)


def _mc_verdicts(site, ro, kind, clf, n, seed):
    """Comparator verdicts at T_opt (clean twin) and at the two wrong thresholds, identical MC samples."""
    p0 = site_statistic(site, ro, clf, False)
    p1 = site_statistic(site, ro, clf, True)
    r = optimal_threshold(p0, p1)
    out = {}
    for label, T in (("T_opt", r.threshold), *zip(("T_low", "T_high"), wrong_thresholds(p0, p1, r.threshold))):
        m = mc_error_rates(site, ro, T, n, np.random.default_rng(seed), kind=kind, weights=clf.weights, llr=clf.llr)
        out[label] = (mc_agrees(kind, r, p0, p1, m["fp"], m["fn"], n), T, m["fp"] / n, m["fn"] / n)
    return r, out


def _assert_control(r, out):
    msg = f"reported at T_opt={r.threshold:.3f}: FP={r.fp:.3e} FN={r.fn:.3e}; MC: " + "; ".join(
        f"{k}: T={v[1]:.3f} FP={v[2]:.3e} FN={v[3]:.3e} agrees={v[0]}" for k, v in out.items())
    assert out["T_opt"][0], "clean twin failed: " + msg
    assert not out["T_low"][0], "wrong threshold (low) passed the check: " + msg
    assert not out["T_high"][0], "wrong threshold (high) passed the check: " + msg


@pytest.mark.slow
@pytest.mark.parametrize("kind", ["box", "matched", "bayes"])
def test_wrong_threshold_fails_mc_check_cmos(kind):
    """Site, readout, seed and n as in tests/test_counts.py::test_lattice_vs_bruteforce_mc_at_1em3."""
    ro = PixelReadout(classes=[(0.97, [(1.0, 0.0, 0.8)]), (0.03, [(0.8, 0.0, 0.8), (0.1, 3.0, 0.8), (0.1, -3.0, 0.8)])],
                      gain_e_per_adu=0.5)
    site = _site(mu=30.0)
    if kind == "box":
        clf = Classifier("box")
    elif kind == "matched":
        clf = Classifier("matched", weights=optimal_linear_weights(site, ro))
    else:
        clf = Classifier("bayes", llr=bayes_llr_tables(site, ro))
    r, out = _mc_verdicts(site, ro, kind, clf, 1_500_000, 123)
    _assert_control(r, out)


def test_wrong_threshold_fails_mc_check_emccd_box():
    """Site, readout, seed and n as in tests/test_counts.py::test_emccd_lattice_vs_mc."""
    ro = EMCCDReadout(em_gain=300, n_stages=512, cic=0.005, cic_register_fraction=0.5, sigma_ro_e=40.0,
                      register_model="erlang", mode="photon_counting", pc_threshold_sigma=5.0, h_fine=0.01)
    site = _site(mu=8.0, bg=0.05, L=0.0, rs=0.0, box=3)
    r, out = _mc_verdicts(site, ro, "box", Classifier("box"), 400_000, 7)
    _assert_control(r, out)


@pytest.mark.parametrize("kind", ["box", "matched"])
def test_wrong_threshold_fails_cf_check(kind):
    """CF engine (production) against the lattice; setup of tests/test_cf_production.py::test_cf_matches_lattice_cmos."""
    classes = [(1.0, [(1.0, 0.0, 0.5)])]
    fr = pixel_fractions(PSFSpec("gaussian", waist_m=1.5e-6), 0.9e-6, 5)
    site = SiteModel(frac=fr, mu=40.0, bg=np.full((5, 5), 0.4), tau=TauDistribution.from_hazard(1e-3, 32),
                     rate_rel_sd=0.1, ds=0.01)
    ro = PixelReadout(classes=classes, gain_e_per_adu=None, h_fine=0.01)
    nz = CFNoise.from_classes(classes)
    w = None if kind == "box" else matched_weights_analytic(site, nz)
    clf = Classifier("box") if w is None else Classifier("matched", weights=w)
    rl = optimal_threshold(site_statistic(site, ro, clf, False), site_statistic(site, ro, clf, True))
    q0, q1 = cf_pair(site, nz, w)
    rc = optimal_threshold(q0, q1)

    def agrees(eps):  # the predicate of test_cf_matches_lattice_cmos: eps equal to rel 2%
        return abs(eps - rl.eps) <= CF_REL_TOL * abs(rl.eps)

    lo, hi = wrong_thresholds(q0, q1, rc.threshold)
    e_lo, e_hi = rates_at(q0, q1, lo).eps, rates_at(q0, q1, hi).eps
    msg = f"lattice eps={rl.eps:.4e}; CF eps at T_opt={rc.eps:.4e}, at T_low={e_lo:.4e}, at T_high={e_hi:.4e}"
    assert agrees(rc.eps), "clean twin failed: " + msg
    assert not agrees(e_lo), "wrong threshold (low) passed the check: " + msg
    assert not agrees(e_hi), "wrong threshold (high) passed the check: " + msg


def test_delta_frac_is_the_registered_value():
    """Guard against silent loosening of the perturbation (CLAUDE.md rule 10)."""
    assert DELTA_FRAC == 0.10
