import numpy as np
import pytest
from scipy import stats

from p1_atomdet.atoms.loss_phenom import (TauDistribution, lossy_poisson_pmf, manetsch_lossy_poisson_pmf,
                                          p_photon_from_survival)
from p1_atomdet.atoms.scattering import saturation_for_rate, scattering_rate
from p1_atomdet.counts.lattice import Lattice, convolve, convolve_many, mixture, rebin


def test_convolution_of_poissons_is_poisson():
    n = np.arange(60)
    a = Lattice(stats.poisson.pmf(n, 3.0), 0.0, 1.0)
    b = Lattice(stats.poisson.pmf(n, 5.5), 0.0, 1.0)
    c = convolve(a, b)
    assert np.allclose(c.p[:60], stats.poisson.pmf(n, 8.5), atol=1e-15)
    many = convolve_many([Lattice(stats.poisson.pmf(np.arange(40), 0.4), 0.0, 1.0)] * 49)
    assert many.sf(40) == pytest.approx(stats.poisson.sf(39, 19.6), rel=1e-9)


def test_mixture_and_rebin_preserve_mass_and_mean():
    a = Lattice(np.array([0.5, 0.5]), 0.0, 0.5)
    b = Lattice(np.array([1.0]), 2.0, 0.5)
    m = mixture([a, b], [0.25, 0.75])
    assert m.total() == pytest.approx(1.0)
    assert m.mean() == pytest.approx(0.25 * 0.25 + 0.75 * 2.0)
    vals = np.array([0.13, 1.71, 2.02])
    pr = np.array([0.2, 0.5, 0.3])
    r = rebin(vals, pr, 0.0, 0.25)
    assert r.total() == pytest.approx(1.0)
    assert r.mean() == pytest.approx(np.sum(vals * pr))


def test_scattering_rate_limits():
    G = 32.889e6
    assert scattering_rate(G, 1e9) == pytest.approx(G / 2, rel=1e-6)
    assert scattering_rate(G, 1.0, 0.0) == pytest.approx(G / 4)
    s = saturation_for_rate(G, 7e4, -15.5 * G)
    assert scattering_rate(G, s, -15.5 * G) == pytest.approx(7e4)


@pytest.mark.parametrize("L", [1e-4, 0.05, 0.7])
def test_manetsch_formula_normalised_and_close_to_mixture(L):
    n = np.arange(0, 400)
    pm = manetsch_lossy_poisson_pmf(n, 16.0, 70.0, L)
    po = lossy_poisson_pmf(n, 70.0, 16.0, TauDistribution.from_hazard(L, 128))
    assert pm.sum() == pytest.approx(1.0, abs=1e-12)
    assert po.sum() == pytest.approx(1.0, abs=1e-12)
    # The two conventions differ by O(L) photons (survivor mean lambda1 - L; trigger photon counted).
    assert np.abs(np.cumsum(pm) - np.cumsum(po)).max() < 2 * L * 0.1 + 1e-12


def test_tau_distribution_hazard_and_sampling():
    L = 0.3
    t = TauDistribution.from_hazard(L, 64)
    assert t.p_survive == pytest.approx(np.exp(-L))
    rng = np.random.default_rng(0)
    s = t.sample(rng, 200000)
    assert np.mean(s >= 1.0) == pytest.approx(np.exp(-L), abs=0.005)
    lost = s[s < 1]
    # truncated exponential mean: 1/L - e^-L/(1-e^-L)
    assert lost.mean() == pytest.approx(1 / L - np.exp(-L) / (1 - np.exp(-L)), abs=0.005)
    # shape transfer
    sh = TauDistribution.with_shape(np.array([0.25, 0.75]), np.array([1.0, 3.0]), 0.99)
    assert sh.p_loss == pytest.approx(0.01) and sh.weights[1] == pytest.approx(0.0075)


def test_p_photon_from_survival_roundtrip():
    p = p_photon_from_survival(0.999864, 7e4, 0.08, tau_vac_s=1374.0)
    L = p * 7e4 * 0.08 + 0.08 / 1374.0
    assert np.exp(-L) == pytest.approx(0.999864)
