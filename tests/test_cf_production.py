"""CF engine and production paths against the lattice engine."""
import numpy as np
import pytest

from p1_atomdet.atoms.loss_phenom import TauDistribution
from p1_atomdet.counts.cf import CFNoise, cf_pair, matched_weights_analytic
from p1_atomdet.counts.pmf import Classifier, SiteModel, site_statistic
from p1_atomdet.fidelity import optimal_threshold
from p1_atomdet.optics.array import PSFSpec, pixel_fractions
from p1_atomdet.production import (Camera, Scenario, eps_at, geometry, pc_click_prob, pc_pair, read_classes,
                                   time_to)
from p1_atomdet.sensor.cmos import PixelReadout
from p1_atomdet.sensor.emccd import EMCCDReadout, emccd_pixel


def _site(mu=40.0, bg=0.4, box=5, L=1e-3, rs=0.1):
    fr = pixel_fractions(PSFSpec("gaussian", waist_m=1.5e-6), 0.9e-6, box)
    return SiteModel(frac=fr, mu=mu, bg=np.full((box, box), bg), tau=TauDistribution.from_hazard(L, 32),
                     rate_rel_sd=rs, ds=0.01)


@pytest.mark.parametrize("classes", [[(1.0, [(1.0, 0.0, 0.5)])],
                                     [(0.95, [(1.0, 0.0, 0.3)]), (0.05, [(1.0, 0.0, 1.5)])]])
def test_cf_matches_lattice_cmos(classes):
    site = _site()
    ro = PixelReadout(classes=classes, gain_e_per_adu=None, h_fine=0.01)
    nz = CFNoise.from_classes(classes)
    for w in (None, matched_weights_analytic(site, nz)):
        clf = Classifier("box") if w is None else Classifier("matched", weights=w)
        rl = optimal_threshold(site_statistic(site, ro, clf, False), site_statistic(site, ro, clf, True))
        q0, q1 = cf_pair(site, nz, w)
        rc = optimal_threshold(q0, q1)
        assert rc.eps == pytest.approx(rl.eps, rel=0.02)


def test_cf_matches_lattice_emccd_analog():
    site = _site(mu=20.0, bg=0.05, box=3, L=0.0, rs=0.0)
    ro = EMCCDReadout(em_gain=500, n_stages=512, cic=0.0, cic_register_fraction=0.0, sigma_ro_e=40.0,
                      register_model="stagewise", h_fine=0.01)
    rl = optimal_threshold(site_statistic(site, ro, Classifier("box"), False), site_statistic(site, ro, Classifier("box"), True))
    nz = CFNoise(comps=[(1.0, 0.0, 40.0 / 500)], gain_kind="gamma", f2=ro.f2)
    q0, q1 = cf_pair(site, nz)
    rc = optimal_threshold(q0, q1)
    assert rc.eps == pytest.approx(rl.eps, rel=0.03)


def test_pc_click_probability_matches_lattice():
    ro = EMCCDReadout(em_gain=1000, n_stages=604, cic=0.0, cic_register_fraction=0.0, sigma_ro_e=40.0,
                      register_model="stagewise", mode="photon_counting", pc_threshold_sigma=5.0, h_fine=0.001)
    # the lattice converges to the quadrature value as h -> 0 (h = 0.004, 0.002, 0.001 tested by hand)
    for lam in (0.0, 0.3, 2.0):
        p_lat = emccd_pixel(lam, 0.0, ro).p[1]
        p_cf = pc_click_prob(lam, 1000, 604, ro.f2, 40.0, 5.0)[0]
        assert p_cf == pytest.approx(p_lat, rel=0.08 if lam == 0 else 2e-3)


@pytest.mark.slow
def test_pc_moment_matched_gamma_vs_stage_mc():
    from scipy import special
    from p1_atomdet.sensor.emccd import register_mc
    ro = EMCCDReadout(em_gain=1000, n_stages=604, cic=0.0, cic_register_fraction=0.0, sigma_ro_e=40.0,
                      register_model="stagewise")
    rng = np.random.default_rng(3)
    out = register_mc(np.ones(400000, int), ro, 7) + 40.0 * rng.standard_normal(400000)
    p_mc = np.mean(out > 200.0)
    a, th = 1 / (ro.f2 - 1), ro.f2 - 1
    x, w = np.polynomial.hermite_e.hermegauss(40)
    p_formula = np.sum(w / w.sum() * special.gammaincc(a, np.clip(0.2 - 0.04 * x, 0, None) / th))
    assert p_formula == pytest.approx(p_mc, abs=3 * np.sqrt(p_mc * (1 - p_mc) / 400000))


def test_pc_pair_single_frame_matches_lattice():
    ro = EMCCDReadout(em_gain=1000, n_stages=604, cic=0.0, cic_register_fraction=0.0, sigma_ro_e=40.0,
                      register_model="stagewise", mode="photon_counting", pc_threshold_sigma=5.0, h_fine=0.001)
    site = _site(mu=3.0, bg=0.02, box=3, L=0.0, rs=0.0)
    mask = np.ones((3, 3))
    P0, P1 = pc_pair(site, lambda lam: pc_click_prob(lam, 1000, 604, ro.f2, 40.0, 5.0), mask)
    l0 = site_statistic(site, ro, Classifier("box"), False)
    l1 = site_statistic(site, ro, Classifier("box"), True)
    assert np.allclose(P0.p, l0.p[:P0.p.size], rtol=3e-3, atol=1e-6)
    assert np.allclose(P1.p, l1.p[:P1.p.size], rtol=3e-3, atol=1e-6)


def test_time_to_target_monotone_in_qe():
    sc = Scenario(n_offsets=2, ds=0.04, n_tau=12)
    ts = []
    for qe in (0.4, 0.8):
        cam = Camera(qe=qe, read_classes=read_classes(0.5), pixels_per_fwhm=1.0)
        t, _ = time_to(sc, cam, 1e-3)
        ts.append(t)
    assert ts[1] < ts[0]
    # photon budget scaling: time at fixed infidelity scales roughly as 1/QE when read noise is small
    assert 0.4 < ts[1] / ts[0] < 0.7
