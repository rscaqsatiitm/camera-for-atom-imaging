import numpy as np
import pytest

from p1_atomdet.optics.array import PSFSpec, neighbour_fractions, pitch_for_pixels_per_fwhm, pixel_fractions
from p1_atomdet.optics.collection import collection_fraction, collection_fraction_numeric
from p1_atomdet.optics.psf_airy import airy_fwhm, encircled_energy, psf_airy, psf_gaussian, radial_fwhm
from p1_atomdet.optics.psf_vectorial import psf_vectorial

LAM = 852.347e-9


@pytest.mark.parametrize("pat", ["isotropic", "sigma", "pi"])
def test_collection_limits(pat):
    assert collection_fraction(1.0, pat) == pytest.approx(0.5)
    assert collection_fraction(0.0, pat) == pytest.approx(0.0)
    for na in (0.3, 0.55, 0.8):
        assert collection_fraction(na, pat) == pytest.approx(collection_fraction_numeric(na, pat, 20001), abs=1e-8)


def test_collection_small_na_slopes():
    na = 1e-3
    assert collection_fraction(na, "isotropic") / na**2 == pytest.approx(0.25, rel=1e-4)
    assert collection_fraction(na, "sigma") / na**2 == pytest.approx(3 / 8, rel=1e-4)
    assert collection_fraction(na, "pi") / na**4 == pytest.approx(3 / 16, rel=1e-3)


def test_airy_normalisation_and_first_ring():
    f = lambda r: psf_airy(r, LAM, 0.5)
    assert encircled_energy(f, 0.6098 * LAM / 0.5) == pytest.approx(0.8378, abs=2e-4)
    g = lambda r: psf_gaussian(r, 1e-6)
    assert encircled_energy(g, 5e-6) == pytest.approx(1.0, abs=1e-9)


@pytest.mark.parametrize("dipole", ["isotropic", "inplane", "axial"])
def test_vectorial_normalised(dipole):
    na = 0.5
    f = lambda r: psf_vectorial(r, LAM, na, dipole, n_quad=400)
    # EE(R) -> 1 - O(1/(k NA R)); at 40 lambda/NA the Airy deficit is 2/(pi v) ~ 2.5e-3
    ee = encircled_energy(f, 40 * LAM / na, 4001)
    assert 0.99 < ee < 1.0005


def test_vectorial_tends_to_airy_at_low_na():
    na = 0.1
    fv = lambda r: psf_vectorial(r, LAM, na)
    assert radial_fwhm(fv, 2 * LAM / na) / airy_fwhm(LAM, na) == pytest.approx(1.0, abs=2e-3)
    r = np.linspace(0, 1.5 * LAM / na, 50)
    assert np.allclose(fv(r), psf_airy(r, LAM, na), rtol=0, atol=5e-3 * psf_airy(0.0, LAM, na))


def test_pixel_fractions_converge_and_conserve():
    psf = PSFSpec("vectorial", na=0.55)
    pitch = pitch_for_pixels_per_fwhm(psf, 2.0)
    f8 = pixel_fractions(psf, pitch, 9, (0.2, -0.3), supersample=8)
    f32 = pixel_fractions(psf, pitch, 9, (0.2, -0.3), supersample=32)
    assert np.abs(f8 - f32).max() < 2e-3
    g = PSFSpec("gaussian", waist_m=1.0e-6)
    fg = pixel_fractions(g, 0.5e-6, 15, supersample=8)
    assert fg.sum() == pytest.approx(1.0, abs=1e-6)
    blurred = pixel_fractions(PSFSpec("gaussian", waist_m=1.0e-6, blur_sigma_m=0.3e-6), 0.5e-6, 15, supersample=8)
    assert blurred.sum() == pytest.approx(1.0, abs=1e-4)
    assert blurred.max() < fg.max()


def test_neighbour_crosstalk_decreases_with_spacing():
    psf = PSFSpec("gaussian", waist_m=1.5e-6)
    near = sum(v.sum() for v in neighbour_fractions(psf, 0.9e-6, 4e-6, 7).values())
    far = sum(v.sum() for v in neighbour_fractions(psf, 0.9e-6, 8e-6, 7).values())
    assert far < near
