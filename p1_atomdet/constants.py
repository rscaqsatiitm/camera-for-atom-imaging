"""Physical and atomic constants with provenance.

Atomic data: D. A. Steck, "Cesium D Line Data" and "Rubidium 87 D Line Data", revision 2.3.4
(8 August 2025), https://steck.us/alkalidata, retrieved and read 2026-10-14 (Tables 2, 3, 7).
"""
from __future__ import annotations

from dataclasses import dataclass

from camsim_common.provenance import ProvValue, pv

# CODATA 2018 exact SI constants (SOURCE: Steck Table 1 reproduces these; SI definitions are exact)
H = 6.62607015e-34          # J s
HBAR = H / (2 * 3.141592653589793)
KB = 1.380649e-23           # J/K
C = 299792458.0             # m/s
AMU = 1.66053906660e-27     # kg

STECK_CS = "Steck, Cesium D Line Data rev. 2.3.4 (2025), steck.us/alkalidata"
STECK_RB87 = "Steck, Rubidium 87 D Line Data rev. 2.3.4 (2025), steck.us/alkalidata"


@dataclass(frozen=True)
class Species:
    name: str
    mass_kg: ProvValue
    wavelength_m: ProvValue
    gamma_s: ProvValue          # decay rate Gamma (s^-1), natural linewidth = Gamma/2pi
    isat_iso_W_m2: ProvValue    # isotropic-polarisation effective saturation intensity (cycling hyperfine pair)
    isat_cycling_W_m2: ProvValue
    recoil_velocity_m_s: ProvValue

    @property
    def k(self) -> float:
        return 2 * 3.141592653589793 / float(self.wavelength_m)

    @property
    def recoil_energy_J(self) -> float:
        return 0.5 * float(self.mass_kg) * float(self.recoil_velocity_m_s) ** 2


CS133 = Species(
    name="Cs133",
    mass_kg=pv(2.20694695096e-25, "kg", "SOURCE", f"{STECK_CS}, Table 2"),
    wavelength_m=pv(852.34727582e-9, "m", "SOURCE", f"{STECK_CS}, Table 3 (vacuum)"),
    gamma_s=pv(32.889e6, "s^-1", "SOURCE", f"{STECK_CS}, Table 3 (32.889(84)e6 s^-1; 2pi x 5.234(13) MHz)"),
    isat_iso_W_m2=pv(27.119, "W/m^2", "SOURCE", f"{STECK_CS}, Table 7 (2.7119(69) mW/cm^2, F=4->F'=5 isotropic)"),
    isat_cycling_W_m2=pv(11.049, "W/m^2", "SOURCE", f"{STECK_CS}, Table 7 (1.1049(28) mW/cm^2, cycling)"),
    recoil_velocity_m_s=pv(3.5225e-3, "m/s", "SOURCE", f"{STECK_CS}, Table 3"),
)

RB87 = Species(
    name="Rb87",
    mass_kg=pv(1.44316089500e-25, "kg", "SOURCE", f"{STECK_RB87}, Table 2"),
    wavelength_m=pv(780.241209686e-9, "m", "SOURCE", f"{STECK_RB87}, Table 3 (vacuum)"),
    gamma_s=pv(38.117e6, "s^-1", "SOURCE", f"{STECK_RB87}, Table 3 (38.117(11)e6 s^-1; 2pi x 6.0666(18) MHz)"),
    isat_iso_W_m2=pv(35.771, "W/m^2", "SOURCE", f"{STECK_RB87}, Table 7 (3.5771(10) mW/cm^2, F=2->F'=3 isotropic)"),
    isat_cycling_W_m2=pv(16.6933, "W/m^2", "SOURCE", f"{STECK_RB87}, Table 7 (1.66933(49) mW/cm^2, cycling)"),
    recoil_velocity_m_s=pv(5.8845e-3, "m/s", "SOURCE", f"{STECK_RB87}, Table 3"),
)

SPECIES = {"Cs133": CS133, "Rb87": RB87}
