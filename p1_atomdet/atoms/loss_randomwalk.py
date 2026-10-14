"""Loss model B: semiclassical 3D recoil random walk in a Gaussian-beam tweezer.

Potential (focused Gaussian beam along y, the objective axis; waist w0, Rayleigh range zR):
    U(x, y, z) = -U0 / (1 + y^2/zR^2) * exp(-2 (x^2 + z^2) / (w0^2 (1 + y^2/zR^2)))
Photon scattering at rate R: each event gives a recoil kick hbar k along a randomly chosen
imaging-beam direction (absorption) and one along an isotropic random direction (emission).
Cooling: Ornstein-Uhlenbeck thermostat on the cooled axes, v -> v e^(-g dt) + sqrt(kB T_eq/m (1 - e^(-2 g dt))) xi.
Loss: total energy E = m v^2/2 + U > 0 (the atom has more than U0 above the trap bottom).
Integrator: BAOAB splitting (half kick, half drift, OU + recoil kicks, half drift, half kick).

The output is the loss time of each atom as a fraction tau of the exposure (1.0 = survived),
which TauDistribution.from_samples turns into the same object as model A. Vacuum loss is an
independent exponential and is combined outside the kernel.
"""
from __future__ import annotations

from dataclasses import dataclass

import numba
import numpy as np

from ..constants import HBAR, KB

HBAR_ = HBAR  # module globals frozen into the numba kernels at compile time
KB_ = KB


@dataclass
class RandomWalkParams:
    mass_kg: float
    k_m: float                 # 2 pi / lambda of the imaging light
    U0_K: float                # trap depth in kelvin (U0 / kB)
    w0_m: float
    zR_m: float
    rate_s: float              # scattering rate R (photons/s)
    gamma_c_s: float           # cooling damping rate (1/s); 0 = no cooling
    T_eq_K: float              # thermostat temperature
    T_init_K: float
    t_exp_s: float
    dt_s: float
    beam_dirs: tuple = ((1, 0, 0), (-1, 0, 0), (0, 0, 1), (0, 0, -1))   # 2D PGC in the x-z plane
    cooled_axes: tuple = (1, 0, 1)   # thermostat on x and z only (2D cooling)


@numba.njit(cache=True)
def _force(x, y, z, U0, w0, zR):
    q = 1.0 + (y / zR) ** 2
    rho2 = x * x + z * z
    e = np.exp(-2.0 * rho2 / (w0 * w0 * q))
    U = -U0 / q * e
    # dU/dx = U * (-4 x / (w0^2 q));  F = -dU/dx
    fx = U * 4.0 * x / (w0 * w0 * q)
    fz = U * 4.0 * z / (w0 * w0 * q)
    # dU/dy: U = -U0 q^-1 exp(-2 rho2 / (w0^2 q)), dq/dy = 2 y / zR^2
    dq = 2.0 * y / (zR * zR)
    dUdq = U * (-1.0 / q + 2.0 * rho2 / (w0 * w0 * q * q))
    fy = -dUdq * dq
    return fx, fy, fz, U


@numba.njit(parallel=True, cache=True)
def _walk(n_atoms, m, k, U0, w0, zR, R, gam, Teq, Tinit, t_exp, dt, dirs, cooled, seed):
    n_steps = int(np.ceil(t_exp / dt))
    out = np.ones(n_atoms)
    Efin = np.zeros(n_atoms)
    vrec = HBAR_ * k / m
    om_r = np.sqrt(4.0 * U0 / (m * w0 * w0))
    om_z = np.sqrt(2.0 * U0 / (m * zR * zR))
    c1 = np.exp(-gam * dt)
    c2 = np.sqrt(KB_ * Teq / m * (1.0 - c1 * c1))
    lam = R * dt
    nd = dirs.shape[0]
    for i in numba.prange(n_atoms):
        np.random.seed(seed + i)
        sx = np.sqrt(KB_ * Tinit / m)
        x = np.random.normal() * sx / om_r
        z = np.random.normal() * sx / om_r
        y = np.random.normal() * sx / om_z
        vx = np.random.normal() * sx
        vy = np.random.normal() * sx
        vz = np.random.normal() * sx
        fx, fy, fz, U = _force(x, y, z, U0, w0, zR)
        lost = False
        for s in range(n_steps):
            vx += 0.5 * dt * fx / m
            vy += 0.5 * dt * fy / m
            vz += 0.5 * dt * fz / m
            x += 0.5 * dt * vx
            y += 0.5 * dt * vy
            z += 0.5 * dt * vz
            if gam > 0.0:
                if cooled[0]:
                    vx = c1 * vx + c2 * np.random.normal()
                if cooled[1]:
                    vy = c1 * vy + c2 * np.random.normal()
                if cooled[2]:
                    vz = c1 * vz + c2 * np.random.normal()
            nsc = np.random.poisson(lam)
            for _ in range(nsc):
                d = np.random.randint(nd)
                vx += vrec * dirs[d, 0]
                vy += vrec * dirs[d, 1]
                vz += vrec * dirs[d, 2]
                ct = 2.0 * np.random.random() - 1.0
                ph = 2.0 * np.pi * np.random.random()
                st = np.sqrt(1.0 - ct * ct)
                vx += vrec * st * np.cos(ph)
                vy += vrec * st * np.sin(ph)
                vz += vrec * ct
            x += 0.5 * dt * vx
            y += 0.5 * dt * vy
            z += 0.5 * dt * vz
            fx, fy, fz, U = _force(x, y, z, U0, w0, zR)
            vx += 0.5 * dt * fx / m
            vy += 0.5 * dt * fy / m
            vz += 0.5 * dt * fz / m
            E = 0.5 * m * (vx * vx + vy * vy + vz * vz) + U
            if E > 0.0:
                out[i] = (s + 1) * dt / t_exp
                lost = True
                break
        Efin[i] = (0.5 * m * (vx * vx + vy * vy + vz * vz) + U + U0) if not lost else np.nan
    return out, Efin



def simulate_loss_times(p: RandomWalkParams, n_atoms: int, seed: int) -> dict:
    """Run the walk. Returns tau per atom, survival, and the mean final energy (in kelvin, from trap bottom)."""
    dirs = np.asarray(p.beam_dirs, float)
    dirs /= np.linalg.norm(dirs, axis=1)[:, None]
    tau, E = _walk(int(n_atoms), p.mass_kg, p.k_m, p.U0_K * KB, p.w0_m, p.zR_m, p.rate_s, p.gamma_c_s,
                   p.T_eq_K, p.T_init_K, p.t_exp_s, p.dt_s, dirs, np.asarray(p.cooled_axes, np.int64), int(seed))
    surv = float(np.mean(tau >= 1.0))
    # Mean energy above the bottom of survivors; for a 3D harmonic trap <E> = 3 kB T.
    e_mean = float(np.nanmean(E) / KB) if np.any(np.isfinite(E)) else float("nan")
    return {"tau": tau, "survival": surv, "E_mean_K": e_mean, "T_eff_K": e_mean / 3.0}


def combine_vacuum(tau: np.ndarray, t_exp_s: float, tau_vac_s: float, rng: np.random.Generator) -> np.ndarray:
    """Apply independent background-gas loss: tau_total = min(tau_heating, tau_vacuum)."""
    if not np.isfinite(tau_vac_s):
        return tau
    tv = rng.exponential(tau_vac_s / t_exp_s, size=tau.size)
    return np.minimum(tau, np.where(tv < 1.0, tv, 1.0))


def loss_time_shape(tau: np.ndarray, n_bins: int = 20):
    """Normalised histogram of loss times of lost atoms (the 'shape' compared between models A and B)."""
    lost = tau[tau < 1.0]
    h, e = np.histogram(lost, bins=n_bins, range=(0, 1), density=True)
    return 0.5 * (e[1:] + e[:-1]), h, lost.size
