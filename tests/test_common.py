import numpy as np
import pytest

from camsim_common import SensorModel, pv
from camsim_common.io import datasets_equal, load_h5, save_h5
from camsim_common.placeholders import open_placeholders, scan
from camsim_common.provenance import ProvValue, placeholder
from camsim_common.rng import stream
from camsim_common.schema import (DarkArrhenius, EmpiricalSigmaCDF, GaussianRTSMixture, QETable, RTSAmplitude,
                                  read_noise_median, read_noise_rms)
from camsim_common.presets.build_validation_presets import orca_quest_20up_manetsch, zyla55_picken


def test_pointer_required_for_sourced_values():
    with pytest.raises(ValueError):
        ProvValue(value=1.0, unit="", prov="SOURCE")
    assert placeholder(0.5, "e-").tag().startswith("[VALUE?: approx. 0.5")


def test_yaml_roundtrip(tmp_path):
    for sm in (orca_quest_20up_manetsch(), zyla55_picken()):
        p = tmp_path / "s.yaml"
        sm.to_yaml(p)
        assert SensorModel.from_yaml(p) == sm
        assert SensorModel.from_yaml(sm.to_yaml()) == sm


def test_rts_mixture_median_and_rms():
    rn = GaussianRTSMixture(sigma_core_e=pv(0.3, "e-", "DERIVED", "test"), f_tail=pv(0.02, "", "DERIVED", "test"),
                            p_jump=pv(0.5, "", "DERIVED", "test"),
                            amplitudes=[RTSAmplitude(amplitude_e=pv(2.0, "e-", "DERIVED", "test"))])
    # median pixel is a core pixel; rms^2 = s^2 + f p A^2
    assert read_noise_median(rn) == pytest.approx(0.3)
    assert read_noise_rms(rn) == pytest.approx(np.sqrt(0.09 + 0.02 * 0.5 * 4.0))


def test_empirical_cdf_classes():
    rn = EmpiricalSigmaCDF(sigma_e=[0.5, 1.0, 3.0], cdf=[0.6, 0.99, 1.0], prov=pv("t", "", "DERIVED", "test"))
    assert read_noise_median(rn) == pytest.approx(0.5)
    assert read_noise_rms(rn) == pytest.approx(np.sqrt(0.6 * 0.25 + 0.39 * 1.0 + 0.01 * 9.0))


def test_qe_interpolation_and_dark():
    q = QETable(wavelength_nm=[780, 852], temperature_K=[253, 293], values=[[0.6, 0.4], [0.7, 0.5]],
                prov=pv("t", "", "DERIVED", "test"))
    assert q(816, 293) == pytest.approx(0.6)
    assert q(852, 273) == pytest.approx(0.45)
    d = DarkArrhenius(D_ref=pv(1.0, "", "DERIVED", "t"), T_ref_K=pv(273.15, "K", "DERIVED", "t"),
                      Ea_eV=pv(0.6, "eV", "DERIVED", "t"), T_op_K=pv(253.15, "K", "DERIVED", "t"))
    assert d(273.15) == pytest.approx(1.0)
    assert d() < 1.0


def test_placeholder_scan_finds_preset_entries():
    rows = open_placeholders()
    assert any("read_noise.sigma_e" in r["path"] for r in rows)  # ORCA-Quest read noise is a placeholder
    assert all(r["prov"] == "PLACEHOLDER" for r in rows)
    assert len(scan()) > len(rows)


def test_h5_roundtrip_and_equality(tmp_path):
    cfg = {"a": 1, "b": [1, 2]}
    p1 = save_h5(tmp_path / "x.h5", {"g": {"y": np.arange(5)}, "z": np.ones(3)}, cfg)
    p2 = save_h5(tmp_path / "y.h5", {"g": {"y": np.arange(5)}, "z": np.ones(3)}, cfg)
    d, a = load_h5(p1)
    assert a["config"] == cfg and np.array_equal(d["g/y"], np.arange(5))
    assert datasets_equal(p1, p2)


def test_rng_streams_reproducible_and_distinct():
    a = stream(1, 3, "x").random(5)
    b = stream(1, 3, "x").random(5)
    c = stream(1, 4, "x").random(5)
    assert np.array_equal(a, b) and not np.array_equal(a, c)
