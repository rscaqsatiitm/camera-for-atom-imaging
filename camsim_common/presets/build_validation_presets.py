"""Generate the SensorModel presets used by the P1 validation stage.

Only the two cameras of the validation papers are built here (ORCA-Quest C15550-20UP as used by
Manetsch et al.; Andor Zyla 5.5 as used by Picken et al.). The eight production presets follow
after the validation check-in. Values not stated in a retrieved source are PLACEHOLDER and are
either irrelevant at the validation operating point (shown in the report) or fitted.

    python -m camsim_common.presets.build_validation_presets
"""
from __future__ import annotations

from pathlib import Path

from camsim_common.provenance import placeholder, pv
from camsim_common.schema import (ADCBlock, CosmicBlock, DarkArrhenius, GaussianCore, QETable, SensorModel,
                                  ShutterBlock)

HERE = Path(__file__).parent
MAN = "Manetsch et al., arXiv:2403.12021v4, Methods 'Single-atom imaging'"
PIC = "Picken et al., arXiv:1709.00946v2"
USER_TABLE = "USER verified-parameter table (P1 prompt, Oct 2026)"


def cosmic_placeholder(thickness_um: float):
    return CosmicBlock(
        rate_per_cm2_s=placeholder(0.0167, "cm^-2 s^-1", "sea-level muon flux ~1 cm^-2 min^-1, recalled; retrieve PDG value"),
        deposit_e_per_um=placeholder(80, "e-/um", "minimum-ionising deposit in Si, recalled"),
        thickness_um=placeholder(thickness_um, "um", "sensitive thickness not stated in a retrieved source"),
    )


def orca_quest_20up_manetsch() -> SensorModel:
    return SensorModel(
        name="orca_quest_c15550_20up_manetsch",
        camera="Hamamatsu ORCA-Quest C15550-20UP",
        mode="as used by Manetsch et al. (scan mode not stated)",
        sensor_type="qCMOS",
        notes="Validation preset V1. Output treated in photoelectron units, as in the paper's histograms.",
        qe=QETable(wavelength_nm=[852.35], values=[[0.44]],
                   prov=pv("single point", "", "SOURCE", f"{MAN}: '44% are detected on the sensor due to the quantum efficiency at 852 nm'")),
        qe_780=placeholder(None, "", "not used in V1"),
        qe_852=pv(0.44, "", "SOURCE", f"{MAN}"),
        prnu=placeholder(0.0, "", "not stated; effect on V1 tested in the sensitivity table"),
        read_noise=GaussianCore(sigma_e=placeholder(0.43, "e- rms", "ORCA-Quest 2 standard-scan value from USER table; C15550-20UP value not retrieved")),
        read_noise_median_e=placeholder(0.43, "e-", "as sigma"),
        read_noise_rms_e=placeholder(0.43, "e-", "as sigma"),
        adc=ADCBlock(gain_e_per_adu=placeholder(0.107, "e-/ADU", "recalled; V1 works in photoelectrons, quantisation step < 0.25 sigma"),
                     adc_bits=placeholder(16, "bit", "recalled"),
                     full_well_e=placeholder(7000, "e-", "recalled; V1 pixel signals < 20 e-"),
                     offset_adu=placeholder(200, "ADU", "recalled")),
        dark=DarkArrhenius(D_ref=placeholder(0.016, "e-/px/s", "recalled; dark over 49 px x 80 ms < 0.07 e-"),
                           T_ref_K=placeholder(238.15, "K", "recalled cooling set point"),
                           Ea_eV=placeholder(0.6, "eV", "generic Si; P4 to supply"),
                           T_op_K=placeholder(238.15, "K", "recalled")),
        pixel_pitch_um=pv(4.6, "um", "DERIVED", f"{MAN}: PSF waist 7.8 um on camera = 1.7 px -> 4.59 um/px; USER table lists 4.6 um for ORCA-Quest 2"),
        shutter=ShutterBlock(type="rolling", line_time_us=placeholder(7.2, "us", "recalled; irrelevant for 80 ms continuous imaging"),
                             frame_time_ms=placeholder(80.0, "ms", "set by the experiment"), n_rows=2304),
        cosmic=cosmic_placeholder(10.0),
    )


def zyla55_picken() -> SensorModel:
    return SensorModel(
        name="zyla55_picken_globalshutter_200MHz",
        camera="Andor Zyla 5.5 sCMOS",
        mode="global shutter, 200 MHz",
        sensor_type="sCMOS",
        notes="Validation preset V2. ADC gain is not stated by Picken et al.; it is fitted in V2.",
        qe=QETable(wavelength_nm=[852.0], values=[[0.22]],
                   prov=pv("single point", "", "SOURCE", f"{PIC}, Sec. II: 'At our imaging wavelength, 852 nm, QE = 22 %'")),
        qe_780=placeholder(None, "", "not used in V2"),
        qe_852=pv(0.22, "", "SOURCE", f"{PIC}, Sec. II"),
        prnu=placeholder(0.0, "", "not stated"),
        read_noise=GaussianCore(sigma_e=pv(2.2, "e- rms", "SOURCE", f"{PIC}, Sec. II: 'global shutter mode with 200 MHz readout speed, delta_ro = 2.2 e- rms'")),
        read_noise_median_e=placeholder(2.2, "e-", "paper gives rms only; Gaussian assumed"),
        read_noise_rms_e=pv(2.2, "e-", "SOURCE", f"{PIC}, Sec. II"),
        adc=ADCBlock(gain_e_per_adu=placeholder(0.38, "e-/ADU", "fitted in V2 from the background histogram"),
                     adc_bits=placeholder(16, "bit", "recalled"),
                     full_well_e=placeholder(30000, "e-", "recalled"),
                     offset_adu=placeholder(100, "ADU", "recalled")),
        dark=DarkArrhenius(D_ref=placeholder(0.1, "e-/px/s", "recalled; neglected by Picken et al. Sec. II"),
                           T_ref_K=placeholder(273.15, "K", ""), Ea_eV=placeholder(0.6, "eV", ""),
                           T_op_K=placeholder(273.15, "K", "")),
        pixel_pitch_um=placeholder(6.5, "um", "recalled; object-plane pixel 1 um is SOURCE"),
        shutter=ShutterBlock(type="global", line_time_us=placeholder(9.24, "us", "recalled"),
                             frame_time_ms=placeholder(40.0, "ms", ""), n_rows=2160),
        cosmic=cosmic_placeholder(10.0),
    )


def main():
    for sm in (orca_quest_20up_manetsch(), zyla55_picken()):
        sm.to_yaml(HERE / f"{sm.name}.yaml")
        print("wrote", sm.name)


if __name__ == "__main__":
    main()
