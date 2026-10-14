"""Semi-analytic count distributions vs exact limits, brute-force MC (at ~1e-3) and IS (deep tails)."""
import numpy as np
import pytest
from scipy import stats

from p1_atomdet.atoms.loss_phenom import TauDistribution
from p1_atomdet.counts.importance import is_fn, is_fp
from p1_atomdet.counts.mc import mc_error_rates, simulate_pixels, statistic_samples
from p1_atomdet.counts.pmf import Classifier, SiteModel, bayes_llr_tables, optimal_linear_weights, site_statistic
from p1_atomdet.fidelity import clopper_pearson, optimal_threshold, rates_at, time_to_target
from p1_atomdet.optics.array import PSFSpec, pixel_fractions
from p1_atomdet.sensor.cmos import PixelReadout
from p1_atomdet.sensor.emccd import EMCCDReadout


def _site(mu=22.0, bg=0.4, L=2e-3, rs=0.1, box=5, ds=0.01):
    fr = pixel_fractions(PSFSpec("gaussian", waist_m=1.7 * 0.9e-6), 0.9e-6, box)
    return SiteModel(frac=fr, mu=mu, bg=np.full((box, box), bg), tau=TauDistribution.from_hazard(L, 48),
                     rate_rel_sd=rs, ds=ds)


def test_poisson_discrimination_exact():
    fr = np.zeros((3, 3)); fr[1, 1] = 1.0
    ro = PixelReadout(classes=[(1.0, [(1.0, 0.0, 0.0)])], gain_e_per_adu=1.0)
    s = SiteModel(frac=fr, mu=12.0, bg=np.full((3, 3), 0.1))
    p0 = site_statistic(s, ro, Classifier("box"), False)
    p1 = site_statistic(s, ro, Classifier("box"), True)
    for T in (2, 4, 6):
        r = rates_at(p0, p1, T)
        assert r.fp == pytest.approx(stats.poisson.sf(T - 1, 0.9), rel=1e-12)
        assert r.fn == pytest.approx(stats.poisson.cdf(T - 1, 12.9), rel=1e-12)


def test_s_grid_matches_tensor_quadrature():
    ro = PixelReadout(classes=[(1.0, [(1.0, 0.0, 0.8)])], gain_e_per_adu=0.5)
    a = _site(ds=0.005)
    b = _site(); b.quadrature, b.n_rate_nodes = "tensor", 41
    p0 = site_statistic(a, ro, Classifier("box"), False)
    ra = rates_at(p0, site_statistic(a, ro, Classifier("box"), True), 20.0)
    rb = rates_at(p0, site_statistic(b, ro, Classifier("box"), True), 20.0)
    assert ra.fn == pytest.approx(rb.fn, rel=5e-3)


@pytest.mark.slow
@pytest.mark.parametrize("kind", ["box", "matched", "bayes"])
def test_lattice_vs_bruteforce_mc_at_1em3(kind):
    """Agreement within the 95% Clopper-Pearson interval at error rates ~1e-3 (CMOS chain, RTS tail, ADC)."""
    ro = PixelReadout(classes=[(0.97, [(1.0, 0.0, 0.8)]), (0.03, [(0.8, 0.0, 0.8), (0.1, 3.0, 0.8), (0.1, -3.0, 0.8)])],
                      gain_e_per_adu=0.5)
    site = _site(mu=30.0)
    if kind == "box":
        clf = Classifier("box")
    elif kind == "matched":
        clf = Classifier("matched", weights=optimal_linear_weights(site, ro))
    else:
        clf = Classifier("bayes", llr=bayes_llr_tables(site, ro))
    p0 = site_statistic(site, ro, clf, False)
    p1 = site_statistic(site, ro, clf, True)
    r = optimal_threshold(p0, p1)
    assert 1e-4 < r.eps < 2e-2
    rng = np.random.default_rng(123)
    n = 1_500_000
    m = mc_error_rates(site, ro, r.threshold, n, rng, kind=kind, weights=clf.weights, llr=clf.llr)
    # matched/bayes statistics are rebinned (linear splitting) on the lattice; MC uses exact values.
    # Compare at the lattice threshold with a half-step tolerance band.
    lo, hi = clopper_pearson(m["fp"], n)
    if kind == "box":
        assert lo <= r.fp <= hi
        lo, hi = clopper_pearson(m["fn"], n)
        assert lo <= r.fn <= hi
    else:
        h = p0.h
        fp_band = (p0.sf(r.threshold + h), p0.sf(r.threshold - h))
        assert fp_band[0] <= hi and fp_band[1] >= lo
        lo, hi = clopper_pearson(m["fn"], n)
        fn_band = (p1.cdf_lt(r.threshold - h), p1.cdf_lt(r.threshold + h))
        assert fn_band[0] <= hi and fn_band[1] >= lo


def test_emccd_lattice_vs_mc():
    ro = EMCCDReadout(em_gain=300, n_stages=512, cic=0.005, cic_register_fraction=0.5, sigma_ro_e=40.0,
                      register_model="erlang", mode="photon_counting", pc_threshold_sigma=5.0, h_fine=0.01)
    site = _site(mu=8.0, bg=0.05, L=0.0, rs=0.0, box=3)
    p0 = site_statistic(site, ro, Classifier("box"), False)
    p1 = site_statistic(site, ro, Classifier("box"), True)
    r = optimal_threshold(p0, p1)
    rng = np.random.default_rng(7)
    n = 400_000
    m = mc_error_rates(site, ro, r.threshold, n, rng)
    lo, hi = clopper_pearson(m["fn"], n)
    assert lo <= r.fn <= hi
    lo, hi = clopper_pearson(m["fp"], n)
    assert lo <= r.fp <= hi


@pytest.mark.slow
def test_importance_sampling_matches_lattice_deep_tail():
    ro = PixelReadout(classes=[(1.0, [(1.0, 0.0, 0.5)])], gain_e_per_adu=0.25)
    site = _site(mu=40.0, bg=0.5, L=1e-4, rs=0.1, box=5)
    p0 = site_statistic(site, ro, Classifier("box"), False)
    p1 = site_statistic(site, ro, Classifier("box"), True)
    T = 36.0
    lat_fp, lat_fn = p0.sf(T), p1.cdf_lt(T)
    assert 1e-9 < lat_fp < 1e-5
    rng = np.random.default_rng(99)
    a = is_fp(site, ro, None, T, 200_000, rng)
    b = is_fn(site, ro, None, T, 200_000, rng)
    assert a.rel_err < 0.05 and b.rel_err < 0.05
    assert abs(a.p - lat_fp) < 3 * a.rel_err * a.p
    assert abs(b.p - lat_fn) < 3 * b.rel_err * b.p


def test_time_to_target_on_synthetic_curve():
    eps = lambda t: 1e-2 * np.exp(-t / 0.01) + 1e-6 * t / 0.01
    t, info = time_to_target(eps, 1e-4, 1e-3, 0.2)
    assert eps(t) == pytest.approx(1e-4, rel=1e-3)
    t2, _ = time_to_target(eps, 1e-7, 1e-3, 0.2)
    assert t2 is None


def test_simulated_pixels_mean():
    ro = PixelReadout(classes=[(1.0, [(1.0, 0.0, 1.0)])], gain_e_per_adu=None)
    site = _site(L=0.0, rs=0.0)
    x = simulate_pixels(site, ro, 200000, True, np.random.default_rng(1))
    assert statistic_samples(x).mean() == pytest.approx(site.mu * site.frac.sum() + site.bg.sum(), rel=2e-3)
