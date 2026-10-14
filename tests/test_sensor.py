import numpy as np
import pytest
from scipy import stats

from p1_atomdet.sensor.cmos import PixelReadout, cmos_pixel
from p1_atomdet.sensor.emccd import EMCCDReadout, emccd_pixel, enf2_finite, register_mc, serial_cic_lattice
from p1_atomdet.sensor.qcmos import miscount_probability, pnr_pixel
from p1_atomdet.sensor.rolling_shutter import min_exposure_for_full_capture, row_capture_fraction
from p1_atomdet.sensor.cosmics import add_cosmics, hit_probability
from p1_atomdet.counts.lattice import Lattice


def test_cmos_moments_analog():
    ro = PixelReadout(classes=[(1.0, [(1.0, 0.0, 1.3)])], gain_e_per_adu=None, h_fine=0.01)
    lat = cmos_pixel(3.0, 0.7, ro)
    assert lat.total() == pytest.approx(1.0, abs=1e-12)
    assert lat.mean() == pytest.approx(3.7, abs=1e-9)
    assert lat.var() == pytest.approx(3.7 + 1.69 + 0.01**2 / 12, rel=1e-4)


def test_cmos_adc_quantisation_variance():
    g = 0.8
    ro = PixelReadout(classes=[(1.0, [(1.0, 0.0, 2.0)])], gain_e_per_adu=g)
    lat = cmos_pixel(10.0, 0.0, ro)
    # Sheppard: rounding adds ~g^2/12 when sigma >> g
    assert lat.var() == pytest.approx(10.0 + 4.0 + g**2 / 12, rel=2e-3)


def test_cmos_full_well_clipping():
    ro = PixelReadout(classes=[(1.0, [(1.0, 0.0, 0.5)])], gain_e_per_adu=1.0, full_well_e=20)
    lat = cmos_pixel(40.0, 0.0, ro)
    assert lat.x[lat.p > 1e-6].max() <= 22


def test_rts_ensemble_mixture_is_mixture_of_classes():
    ro = PixelReadout(classes=[(0.9, [(1.0, 0.0, 0.5)]), (0.1, [(0.5, 0.0, 0.5), (0.25, 3.0, 0.5), (0.25, -3.0, 0.5)])],
                      gain_e_per_adu=None, h_fine=0.05)
    ens = cmos_pixel(1.0, 0.0, ro)
    c0, c1 = cmos_pixel(1.0, 0.0, ro, 0), cmos_pixel(1.0, 0.0, ro, 1)
    assert ens.var() == pytest.approx(0.9 * c0.var() + 0.1 * c1.var(), rel=1e-6)


def test_pnr_rounding():
    lat = pnr_pixel(2.0, 0.0, 0.2)
    assert np.allclose(lat.x[lat.p > 1e-3] % 1.0, 0.0)
    assert miscount_probability(0.2) == pytest.approx(2 * stats.norm.sf(2.5))


def test_emccd_erlang_enf_is_two():
    ro = EMCCDReadout(em_gain=1000, n_stages=512, cic=0.0, cic_register_fraction=0.0, sigma_ro_e=0.0, h_fine=0.01)
    lat = emccd_pixel(5.0, 0.0, ro)
    assert lat.mean() == pytest.approx(5.0, rel=1e-3)
    assert lat.var() / 5.0 == pytest.approx(2.0, rel=2e-3)


def test_enf_formula_limits():
    assert enf2_finite(1.0, 512) == pytest.approx(1.0)
    assert enf2_finite(1e6, 10**6) == pytest.approx(2.0, abs=1e-4)


@pytest.mark.slow
def test_stagewise_mc_matches_eq5():
    g, r = 300, 400
    ro = EMCCDReadout(em_gain=g, n_stages=r, cic=0.0, cic_register_fraction=0.0, sigma_ro_e=0.0, register_model="stagewise")
    out = register_mc(np.ones(100000, int), ro, seed=11) / g
    assert out.mean() == pytest.approx(1.0, abs=0.01)
    f2 = 1 + out.var(ddof=1)
    assert f2 == pytest.approx(enf2_finite(g, r), abs=0.03)


def test_photon_counting_dark_click_rate():
    ro = EMCCDReadout(em_gain=1000, n_stages=512, cic=0.0, cic_register_fraction=0.0, sigma_ro_e=50.0,
                      mode="photon_counting", pc_threshold_sigma=5.0, h_fine=0.002)
    lat = emccd_pixel(0.0, 0.0, ro)
    assert lat.p[1] == pytest.approx(stats.norm.sf(5.0), rel=0.2)
    lat1 = emccd_pixel(1.0, 0.0, ro)
    # Poisson(1) input, Erlang register: P(click) = sum_n Pois(n; 1) Q(n, T/g), T/g = 5*50/1000 = 0.25
    from scipy import special
    n = np.arange(1, 30)
    expected = np.sum(stats.poisson.pmf(n, 1.0) * special.gammaincc(n, 0.25))
    assert lat1.p[1] == pytest.approx(expected, rel=0.01)


def test_serial_cic_mass():
    ro = EMCCDReadout(em_gain=1000, n_stages=512, cic=0.01, cic_register_fraction=1.0, sigma_ro_e=0.0, h_fine=0.01)
    lat = serial_cic_lattice(ro, 20.0)
    assert lat.total() == pytest.approx(1.0, abs=1e-6)
    assert lat.p[0] > 0.99


def test_rolling_shutter_capture():
    rows = np.arange(100)
    f = row_capture_fraction(rows, 10e-6, 10e-3, 0.0, 10e-3)
    assert f[0] == pytest.approx(1.0) and f[-1] < 1.0
    t = min_exposure_for_full_capture(100, 10e-6, 10e-3)
    assert np.allclose(row_capture_fraction(rows, 10e-6, t, 99 * 10e-6, 10e-3), 1.0)
    assert np.allclose(row_capture_fraction(rows, 10e-6, 10e-3, 0.0, 10e-3, shutter="global"), 1.0)


def test_cosmics():
    assert hit_probability(0.0167, 1e-5, 0.08) == pytest.approx(1.336e-8, rel=1e-3)
    lat = add_cosmics(Lattice(np.array([1.0]), 0.0, 1.0), 0.1, 50.0)
    assert lat.sf(25) == pytest.approx(0.1)
