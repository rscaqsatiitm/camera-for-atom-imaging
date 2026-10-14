import numpy as np
import pytest

from p1_atomdet.atoms.loss_phenom import TauDistribution
from p1_atomdet.atoms.loss_randomwalk import RandomWalkParams, combine_vacuum, simulate_loss_times
from p1_atomdet.constants import CS133, KB

BASE = dict(mass_kg=float(CS133.mass_kg), k_m=CS133.k, U0_K=1e-3, w0_m=0.8e-6, zR_m=2.4e-6, rate_s=0.0,
            gamma_c_s=0.0, T_eq_K=20e-6, T_init_K=20e-6, t_exp_s=0.005, dt_s=0.2e-6)


def test_energy_conserved_without_scattering():
    a = simulate_loss_times(RandomWalkParams(**BASE), 500, 1)
    b = simulate_loss_times(RandomWalkParams(**{**BASE, "t_exp_s": 2e-6}), 500, 1)
    assert a["survival"] == 1.0
    assert abs(a["E_mean_K"] / b["E_mean_K"] - 1) < 1e-3


def test_recoil_heating_rate():
    R, t = 2e4, 0.01
    a = simulate_loss_times(RandomWalkParams(**{**BASE, "rate_s": R, "t_exp_s": t}), 3000, 2)
    b = simulate_loss_times(RandomWalkParams(**{**BASE, "rate_s": R, "t_exp_s": 2e-6}), 3000, 2)
    expected = 2 * R * t * CS133.recoil_energy_J / KB
    assert (a["E_mean_K"] - b["E_mean_K"]) == pytest.approx(expected, rel=0.08)


def test_hot_atoms_are_lost_and_tau_distribution_builds():
    r = simulate_loss_times(RandomWalkParams(**{**BASE, "U0_K": 5e-5, "T_init_K": 1e-4}), 400, 3)
    assert r["survival"] < 0.5
    td = TauDistribution.from_samples(r["tau"])
    assert td.p_survive == pytest.approx(r["survival"])


def test_vacuum_combination():
    rng = np.random.default_rng(0)
    tau = combine_vacuum(np.ones(200000), 0.08, 8.0, rng)
    assert np.mean(tau < 1) == pytest.approx(1 - np.exp(-0.01), rel=0.05)
