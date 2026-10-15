"""Generate the production SensorModel presets (P1 production stage).

    python -m camsim_common.presets.build_production_presets

Provenance conventions used here:
  USER     values from PDFs supplied by the PI on 2026-10-15 (SHA-256 prefixes in sources/README.md)
           or from the USER table of the P1 prompt.
  SOURCE   values quoted from vendor documents retrieved and read on 2026-10-15 (sources/README.md).
  DATA     QE at 780 / 852 nm digitised from vendor QE curves (method and uncertainty in the pointer);
           vendors print these as typical curves, not guaranteed values.
  DERIVED  computed from tagged values (formula in the pointer).
  PLACEHOLDER  anything else; swept or bounded in the production study.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from camsim_common.provenance import placeholder, pv
from camsim_common.schema import (ADCBlock, CosmicBlock, DarkArrhenius, EMCCDBlock, EmpiricalSigmaCDF, GaussianCore,
                                  QETable, SensorModel, ShutterBlock, K_B_EV)

HERE = Path(__file__).parent

# document keys (see sources/README.md)
BSI = "USER PDF Prime_BSI_Express_DS_web_datasheet.pdf (rev. 2024-08-20, sha256 1c867b6f)"
Q2 = "USER PDF SCAS0166E03 ORCA-Quest 2 C15550-22UP (Jun. 2026, sha256 8497f943)"
Q2IM = "SOURCE Hamamatsu C15550-22UP instruction manual B9200301-02 Ver.1.2 (Mar. 2025)"
FBT = "SOURCE Hamamatsu ORCA-Fusion BT instruction manual B0850301-05 Ver.1.5 (Oct. 2025)"
FBTBR = "SOURCE Hamamatsu SCAS0143E03 ORCA-Fusion BT brochure (Jun. 2026)"
IXON = "USER PDF iXon_Ultra_EMCCD_c.pdf (Andor, sha256 d0b88cd4); same values in iXonUltraSS 1026 R1"
MAR = "USER PDF andor-marana-scmos-specifications.pdf (MARANAFAMILY SS 0926 R1, sha256 7863aba3)"
ZYLA = "SOURCE Andor 'Zyla for Physical Sciences sCMOS Specifications' LZYLASPECTSS 1020 R1"
PIC = "SOURCE Picken et al., arXiv:1709.00946v2 Sec. II"
BFS = "SOURCE FLIR EMVA 1288 page BFS-U3-17S7 (11/15/2022), softwareservices.flir.com/BFS-U3-17S7/latest/EMVA/EMVA.html"
BFSS = "SOURCE FLIR spec page softwareservices.flir.com/BFS-U3-17S7/latest/Model/spec.html"
QHY = "SOURCE qhyccd.com/astronomical-camera-qhy268/ (modified 2026-06-16)"
QHYQE = "DATA digitised from QHY 'Absolute Response' curve, qhyccd.com/measuring-the-absolute-qe-of-the-qhy268m/ (+-0.5% reading; +-3% rel. scale)"

C2K = 273.15


def ea_two_point(D1, T1, D2, T2):
    """Arrhenius energy with T^1.5 prefactor from two points (DERIVED)."""
    return -K_B_EV * np.log((D2 / D1) * (T1 / T2) ** 1.5) / (1 / T2 - 1 / T1)


def dark(D_ref, T_ref_C, T_op_C, ptr, Ea=None, ea_ptr=""):
    return DarkArrhenius(
        D_ref=pv(D_ref, "e-/px/s", "USER" if ptr.startswith("USER") else "SOURCE", ptr),
        T_ref_K=pv(T_ref_C + C2K, "K", "USER" if ptr.startswith("USER") else "SOURCE", ptr),
        Ea_eV=(pv(round(float(Ea), 4), "eV", "DERIVED", ea_ptr) if Ea is not None
               else placeholder(0.6, "eV", "single vendor point; P4 to supply D(T)")),
        T_op_K=pv(T_op_C + C2K, "K", "USER" if ptr.startswith("USER") else "SOURCE", ptr))


def cosmic(th):
    return CosmicBlock(rate_per_cm2_s=placeholder(0.0167, "cm^-2 s^-1", "sea-level muon flux ~1 cm^-2 min^-1, recalled"),
                       deposit_e_per_um=placeholder(80, "e-/um", "MIP deposit in Si, recalled"),
                       thickness_um=placeholder(th, "um", "sensitive thickness not stated"))


def qe_tab(q780, q852, ptr):
    return QETable(wavelength_nm=[780.24, 852.35], values=[[q780, q852]], prov=pv("780/852 points", "", "DATA", ptr))


def gcore(med, rms, ptr, note=""):
    """Gaussian core if median == rms; else a two-class mixture with core = median and a tail class
    sized to reproduce the rms with tail pixels at 3x the median (shape assumed; P4 will replace)."""
    if rms <= med * 1.0001:
        return GaussianCore(sigma_e=pv(med, "e- rms", "SOURCE" if not ptr.startswith("USER") else "USER", ptr))
    k = 3.0
    f = (rms**2 - med**2) / ((k * med) ** 2 - med**2)
    return EmpiricalSigmaCDF(sigma_e=[med, k * med], cdf=[1 - f, 1.0],
                             prov=pv(f"core {med} e-, {f:.4f} of pixels at {k * med:.2f} e-", "", "DERIVED",
                                     f"matches median {med} and rms {rms} ({ptr}); tail shape assumed 3x median"))


def common(name, camera, mode, stype, **kw):
    return SensorModel(name=name, camera=camera, mode=mode, sensor_type=stype, **kw)


def presets():
    out = []
    # ---------------- Hamamatsu ORCA-Quest 2, standard and ultra-quiet/PNR
    q2qe = "digitised from QE curve p.6 [printed p.11] of SCAS0166E03 (300 dpi; anchors 460 nm 84.7/85, 900 nm 29.5/30); +-1.5% abs"
    ea_q2 = ea_two_point(0.006, -35 + C2K, 0.016, -20 + C2K)
    for mode, med, rms, lt, fr_fps, pnr in (("standard scan", 0.39, 0.43, 3.6, 120, False),
                                            ("ultra quiet scan, PNR", 0.25, 0.30, 16.97, 25.4, True)):
        out.append(common(
            f"orca_quest2_{'std' if not pnr else 'uq_pnr'}", "Hamamatsu ORCA-Quest 2 C15550-22UP", mode, "qCMOS",
            qe=qe_tab(0.480, 0.384, f"{Q2}, {q2qe}"),
            qe_780=pv(0.480, "", "DATA", f"{Q2}, {q2qe}"), qe_852=pv(0.384, "", "DATA", f"{Q2}, {q2qe}"),
            prnu=pv(0.001, "", "USER", f"{Q2}, spec table: PRNU <0.1% (upper bound used)"),
            read_noise=gcore(med, rms, f"{Q2}, spec table: {mode.split(',')[0]} {rms} e- rms, {med} e- median"),
            read_noise_median_e=pv(med, "e-", "USER", f"{Q2}, spec table"),
            read_noise_rms_e=pv(rms, "e-", "USER", f"{Q2}, spec table"),
            adc=ADCBlock(gain_e_per_adu=pv(1.0 if pnr else 0.107, "e-/ADU", "SOURCE",
                                           f"{Q2IM}: {'PNR one count per photoelectron' if pnr else 'conversion factor 0.107 e-/count'}"),
                         adc_bits=pv(16, "bit", "USER", f"{Q2}: digital output 16/12/8 bit"),
                         full_well_e=pv(200 if pnr else 7000, "e-", "SOURCE" if pnr else "USER",
                                        f"{Q2IM}: PNR max output 200 counts at 1x1" if pnr else f"{Q2}: full well 7000 e-"),
                         offset_adu=pv(200, "ADU", "SOURCE", f"{Q2IM}: dark offset 200 counts")),
            dark=dark(0.016, -20, -20, f"{Q2}: forced-air cooled -20 C, 0.016 e-/px/s", ea_q2,
                      "two-point fit 0.006 @ -35 C and 0.016 @ -20 C (SCAS0166E03)"),
            pixel_pitch_um=pv(4.6, "um", "USER", f"{Q2}: pixel size 4.6 um"),
            shutter=ShutterBlock(type="rolling", line_time_us=pv(lt, "us", "DERIVED",
                                 f"{Q2IM}: 1H = {7.2 if not pnr else 33.9} us per row pair -> {lt} us per row"),
                                 frame_time_ms=pv(round(1000 / fr_fps, 2), "ms", "USER", f"{Q2}: {fr_fps} fps full frame CoaXPress"),
                                 n_rows=2304),
            cosmic=cosmic(10.0), qcmos_pnr=pnr))
    # ---------------- Teledyne Photometrics Prime BSI Express, CMS and Speed
    bsiqe = "digitised from QE curve p.3 (300 dpi; anchors 95.4% peak vs 95%, 800 nm 58.5% vs flyer 60%); +-1.5% abs; etalon ripple 46.8-48.9% within +-5 nm of 852"
    ea_bsi = ea_two_point(0.5, -20 + C2K, 1.5, 0 + C2K)
    for mode, rn, lt, fw, fps, bits in (("CMS (100 MHz)", 1.0, 11.4, 1000, 43, 12), ("Speed (200 MHz)", 1.6, 5.14, 10000, 95, 11)):
        tag = "cms" if mode.startswith("CMS") else "speed"
        out.append(common(
            f"prime_bsi_express_{tag}", "Teledyne Photometrics Prime BSI Express (GSENSE2020BSI)", mode, "sCMOS",
            qe=qe_tab(0.638, 0.485, f"{BSI}, {bsiqe}"),
            qe_780=pv(0.638, "", "DATA", f"{BSI}, {bsiqe}"), qe_852=pv(0.485, "", "DATA", f"{BSI}, {bsiqe}"),
            prnu=placeholder(0.005, "", "not stated in the datasheet"),
            read_noise=GaussianCore(sigma_e=pv(rn, "e-", "USER", f"{BSI}: camera modes table, read noise {rn} e- (median/rms not stated)")),
            read_noise_median_e=pv(rn, "e-", "USER", f"{BSI}: '{rn} e-' (statistic not stated)"),
            read_noise_rms_e=placeholder(rn, "e-", "rms not stated; tail unknown; swept"),
            adc=ADCBlock(gain_e_per_adu=placeholder(0.25 if tag == "cms" else 5.0, "e-/ADU", "not stated; ~full well / 2^bits"),
                         adc_bits=pv(bits, "bit", "USER", f"{BSI}: bit depth {bits}-bit"),
                         full_well_e=pv(fw, "e-", "USER", f"{BSI}: full well {fw} e-"),
                         offset_adu=placeholder(100, "ADU", "not stated")),
            dark=dark(1.5, 0, 0, f"{BSI}: camera modes table, 1.5 e-/p/s at 0 C (air)", ea_bsi,
                      "two-point fit 0.5 @ -20 C (cooling options line) and 1.5 @ 0 C (modes table); different tables"),
            pixel_pitch_um=pv(6.5, "um", "USER", f"{BSI}: 6.5 x 6.5 um"),
            shutter=ShutterBlock(type="rolling", line_time_us=pv(lt, "us", "USER", f"{BSI}: line time {lt} us/line"),
                                 frame_time_ms=pv(round(1000 / fps, 2), "ms", "USER", f"{BSI}: {fps} fps full frame"), n_rows=2048),
            cosmic=cosmic(10.0)))
    # ---------------- Andor Zyla 5.5, global shutter 200 MHz (the V2 configuration)
    out.append(common(
        "zyla55_global_200MHz", "Andor Zyla 5.5 sCMOS", "global shutter, 200 MHz", "sCMOS",
        qe=qe_tab(0.404, 0.249, f"DATA digitised from {ZYLA} QE chart p.8 (+-1.5% abs)"),
        qe_780=pv(0.404, "", "DATA", f"{ZYLA}, QE chart p.8, digitised +-1.5%"),
        qe_852=pv(0.249, "", "DATA", f"{ZYLA}, QE chart p.8, digitised +-1.5% (Picken et al. state 0.22)"),
        prnu=placeholder(0.005, "", "not stated"),
        read_noise=gcore(2.3, 2.5, f"{ZYLA} p.6: global shutter @200 MHz median 2.3 [rms 2.5]"),
        read_noise_median_e=pv(2.3, "e-", "SOURCE", f"{ZYLA} p.6"), read_noise_rms_e=pv(2.5, "e-", "SOURCE", f"{ZYLA} p.6"),
        adc=ADCBlock(gain_e_per_adu=pv(0.375, "e-/ADU", "DERIVED", "V2 fit: mean of 0.381 and 0.368 (reports/p1_validation.md Table V2-a)"),
                     adc_bits=pv(16, "bit", "SOURCE", f"{ZYLA} p.6: 12-bit and 16-bit"),
                     full_well_e=pv(30000, "e-", "SOURCE", f"{ZYLA} p.6"), offset_adu=placeholder(100, "ADU", "")),
        dark=dark(0.10, 0, 0, f"{ZYLA} p.6: air cooled 0 C, 0.10 e-/px/s"),
        pixel_pitch_um=pv(6.5, "um", "SOURCE", f"{ZYLA} p.6"),
        shutter=ShutterBlock(type="global", line_time_us=placeholder(0.0, "us", "global shutter"),
                             frame_time_ms=pv(round(1000 / 49, 2), "ms", "SOURCE", f"{ZYLA} p.7: 49 fps global, Camera Link 10-tap"), n_rows=2160),
        cosmic=cosmic(10.0)))
    # ---------------- Andor iXon Ultra 888 EXF, analog and photon counting (10 MHz EM amplifier)
    ixqe = "digitised from QE chart (EXF curve), iXon Ultra spec sheet; +-1.5% abs; sensor-manufacturer curve at 25 C"
    for mode in ("analog", "photon_counting"):
        out.append(common(
            f"ixon888_exf_{'analog' if mode == 'analog' else 'pc'}", "Andor iXon Ultra 888 EXF (CCD201)",
            f"EM amplifier 10 MHz, gain 1000, {mode}", "EMCCD",
            qe=qe_tab(0.808, 0.608, f"{IXON}, {ixqe}"),
            qe_780=pv(0.808, "", "DATA", f"{IXON}, {ixqe}"), qe_852=pv(0.608, "", "DATA", f"{IXON}, {ixqe}"),
            prnu=placeholder(0.005, "", "not stated"),
            read_noise=GaussianCore(sigma_e=pv(40.0, "e- (output amplifier)", "USER", f"{IXON}: read noise without EM gain 40 e- at 10 MHz")),
            read_noise_median_e=pv(40.0, "e-", "USER", f"{IXON}: 10 MHz, without EM"),
            read_noise_rms_e=pv(40.0, "e-", "USER", f"{IXON}: 10 MHz, without EM"),
            adc=ADCBlock(gain_e_per_adu=placeholder(10.0, "e-/ADU (output)", "not stated"),
                         adc_bits=pv(16, "bit", "USER", f"{IXON}: digitization 16-bit"),
                         full_well_e=pv(80000, "e-", "USER", f"{IXON}: active area pixel well depth 80,000 e-"),
                         offset_adu=placeholder(100, "ADU", "")),
            emccd=EMCCDBlock(em_gain=pv(1000, "", "USER", f"{IXON}: EM gain 1-1000"),
                             cic=pv(0.005, "events/px/frame", "USER", f"{IXON}: spurious background 0.005 events/pix @1000x, -85 C (footnote: -95 C, 30 MHz, 10 ms); includes dark"),
                             n_stages=placeholder(604, "", "CCD201 register length not stated; recalled"),
                             cic_register_fraction=placeholder(0.5, "", "split of CIC between parallel and register; swept 0-1"),
                             register_model="stagewise", mode=mode, pc_threshold_sigma=5.0),
            dark=dark(0.00025, -80, -80, f"{IXON}: dark current 0.00025 e-/px/s @ -80 C"),
            pixel_pitch_um=pv(13.0, "um", "USER", f"{IXON}: 13 x 13 um"),
            shutter=ShutterBlock(type="global", line_time_us=placeholder(0.0, "us", "frame transfer"),
                                 frame_time_ms=pv(round(1000 / 26, 2), "ms", "USER", f"{IXON}: 26 fps full frame at 30 MHz"), n_rows=1024),
            cosmic=cosmic(10.0)))
    # ---------------- Hamamatsu ORCA-Fusion BT, ultra quiet scan
    fbqe = "digitised from QE curve p.5 of SCAS0143E03 (anchors 800 nm 58.9 vs printed 58); +-2% abs"
    ea_fb = ea_two_point(0.7, -15 + C2K, 1.0, -8 + C2K)
    out.append(common(
        "orca_fusion_bt_uq", "Hamamatsu ORCA-Fusion BT C15440-20UP", "ultra quiet scan", "sCMOS",
        qe=qe_tab(0.646, 0.464, f"{FBTBR}, {fbqe}"),
        qe_780=pv(0.646, "", "DATA", f"{FBTBR}, {fbqe}"), qe_852=pv(0.464, "", "DATA", f"{FBTBR}, {fbqe}"),
        prnu=pv(0.0006, "", "SOURCE", f"{FBT} p.80: PRNU 0.06% rms"),
        read_noise=gcore(0.6, 0.7, f"{FBT} p.80: ultra quiet 0.6 median, 0.7 rms"),
        read_noise_median_e=pv(0.6, "e-", "SOURCE", f"{FBT} p.80"), read_noise_rms_e=pv(0.7, "e-", "SOURCE", f"{FBT} p.80"),
        adc=ADCBlock(gain_e_per_adu=pv(0.24, "e-/ADU", "SOURCE", f"{FBT} p.80: conversion factor 0.24"),
                     adc_bits=pv(16, "bit", "SOURCE", f"{FBT} p.81: digital output 16 bit"),
                     full_well_e=pv(15000, "e-", "SOURCE", f"{FBT} p.80"), offset_adu=pv(100, "ADU", "SOURCE", f"{FBT} p.80")),
        dark=dark(0.7, -15, -15, f"{FBT} p.80: 0.7 e-/px/s at -15 C", ea_fb, "two-point fit 0.7 @ -15 C and 1.0 @ -8 C"),
        pixel_pitch_um=pv(6.5, "um", "SOURCE", f"{FBT} p.80"),
        shutter=ShutterBlock(type="rolling", line_time_us=pv(80.0, "us", "SOURCE", f"{FBT} p.36: 1H = 80.0 us (ultra quiet)"),
                             frame_time_ms=pv(184.4, "ms", "SOURCE", f"{FBT} p.81: readout 184.4 ms"), n_rows=2304),
        cosmic=cosmic(10.0)))
    # ---------------- Andor Marana 4.2B-6, low-noise 12-bit
    out.append(common(
        "marana_42b6_lownoise", "Andor Marana 4.2B-6", "low noise 12-bit, 180 MHz", "sCMOS",
        qe=qe_tab(0.645, 0.470, f"DATA digitised from {MAR} QE chart (4.2B-6 dashed curve), +-2% abs"),
        qe_780=pv(0.645, "", "DATA", f"{MAR}, QE chart, +-2%"), qe_852=pv(0.470, "", "DATA", f"{MAR}, QE chart, +-2%"),
        prnu=placeholder(0.005, "", "not stated"),
        read_noise=GaussianCore(sigma_e=pv(1.1, "e- median", "USER", f"{MAR}: read noise median 1.1 e- (Low Noise, 12-bit)")),
        read_noise_median_e=pv(1.1, "e-", "USER", f"{MAR}"), read_noise_rms_e=placeholder(1.1, "e-", "rms not stated"),
        adc=ADCBlock(gain_e_per_adu=placeholder(0.3, "e-/ADU", "not stated"), adc_bits=pv(12, "bit", "USER", f"{MAR}"),
                     full_well_e=pv(1250, "e-", "USER", f"{MAR}: 1250 e- low noise mode"), offset_adu=placeholder(100, "ADU", "")),
        dark=dark(0.15, -25, -25, f"{MAR}: 0.15 e-/px/s air cooled @ -25 C", ea_two_point(0.10, -45 + C2K, 0.15, -25 + C2K),
                  "two-point fit 0.10 @ -45 C (liquid) and 0.15 @ -25 C (air)"),
        pixel_pitch_um=pv(6.5, "um", "USER", f"{MAR}"),
        shutter=ShutterBlock(type="rolling", line_time_us=pv(11.4, "us", "DERIVED", f"{MAR}: 43 fps full frame (12-bit, USB) / 2046 rows"),
                             frame_time_ms=pv(23.26, "ms", "USER", f"{MAR}: 43 fps"), n_rows=2046),
        cosmic=cosmic(10.0)))
    # ---------------- Sony IMX571 cooled astro camera (QHY268M), high-gain mode with 2x CMS
    out.append(common(
        "qhy268m_imx571_hg2cms", "QHYCCD QHY268M (Sony IMX571 mono)", "High Gain-2CMS, gain 56", "CMOS",
        qe=qe_tab(0.314, 0.208, QHYQE), qe_780=pv(0.314, "", "DATA", QHYQE), qe_852=pv(0.208, "", "DATA", QHYQE),
        prnu=placeholder(0.005, "", "not stated"),
        read_noise=GaussianCore(sigma_e=pv(1.58, "e-", "SOURCE", f"{QHY}: read-noise chart label 1.58 e- (High Gain-2CMS, gain 56)")),
        read_noise_median_e=placeholder(1.58, "e-", "statistic not stated"), read_noise_rms_e=placeholder(1.58, "e-", ""),
        adc=ADCBlock(gain_e_per_adu=pv(0.44, "e-/ADU", "SOURCE", f"{QHY}: system-gain chart 0.44 e-/ADU (HG, gain 56)"),
                     adc_bits=pv(16, "bit", "SOURCE", f"{QHY}"), full_well_e=pv(28528, "e-", "SOURCE", f"{QHY}: full-well chart, HG gain 56"),
                     offset_adu=placeholder(100, "ADU", "")),
        dark=dark(0.0005, -20, -20, f"{QHY}: 0.0005 e-/px/s at -20 C", ea_two_point(0.0005, -20 + C2K, 0.001, -10 + C2K),
                  "two-point fit 0.0005 @ -20 C and 0.001 @ -10 C (QHY268 spec table)"),
        pixel_pitch_um=pv(3.76, "um", "SOURCE", f"{QHY}"),
        shutter=ShutterBlock(type="rolling", line_time_us=placeholder(40.0, "us", "~ 1/(6 fps x 4176 rows); not stated"),
                             frame_time_ms=pv(166.7, "ms", "SOURCE", f"{QHY}: 6 fps full frame 16-bit"), n_rows=4176),
        cosmic=cosmic(10.0)))
    # ---------------- FLIR Blackfly S BFS-U3-17S7M (IMX425), high conversion gain
    bfqe = f"digitised from EMVA spectral-response PNG, {BFS} (+-1.0% abs; 525 nm reads 71.8 vs 71.50 tabulated)"
    out.append(common(
        "blackfly_s_bfs_u3_17s7m_hcg", "Teledyne FLIR Blackfly S BFS-U3-17S7M (Sony IMX425)", "high conversion gain", "CMOS",
        qe=qe_tab(0.316, 0.186, bfqe), qe_780=pv(0.316, "", "DATA", bfqe), qe_852=pv(0.186, "", "DATA", bfqe),
        prnu=placeholder(0.005, "", "not stated"),
        read_noise=GaussianCore(sigma_e=pv(4.69, "e-", "SOURCE", f"{BFS}: temporal dark noise 4.69 e- (high conversion gain)")),
        read_noise_median_e=placeholder(4.69, "e-", "EMVA temporal dark noise (mean); median not given"),
        read_noise_rms_e=pv(4.69, "e-", "SOURCE", f"{BFS}"),
        adc=ADCBlock(gain_e_per_adu=pv(0.31, "e-/ADU", "SOURCE", f"{BFS}: gain K 0.31 e-/ADU (HCG)"),
                     adc_bits=pv(12, "bit", "SOURCE", f"{BFSS}: ADC 10/12-bit"),
                     full_well_e=pv(19557, "e-", "SOURCE", f"{BFS}: saturation capacity 19557.10 e- (HCG)"),
                     offset_adu=placeholder(0, "ADU", "")),
        dark=DarkArrhenius(D_ref=placeholder(5.0, "e-/px/s", "uncooled; not stated"), T_ref_K=placeholder(313.0, "K", "uncooled"),
                           Ea_eV=placeholder(0.6, "eV", ""), T_op_K=placeholder(313.0, "K", "uncooled board temperature")),
        pixel_pitch_um=pv(9.0, "um", "SOURCE", f"{BFSS}: pixel size 9 um"),
        shutter=ShutterBlock(type="global", line_time_us=placeholder(0.0, "us", "global shutter"),
                             frame_time_ms=pv(round(1000 / 196, 2), "ms", "SOURCE", f"{BFSS}: 196 fps"), n_rows=1080),
        cosmic=cosmic(5.0)))
    return out


def main():
    for sm in presets():
        sm.to_yaml(HERE / f"{sm.name}.yaml")
        chk = sm.read_noise_check()
        print(f"wrote {sm.name:38s} QE852 {float(sm.qe_852.value):.3f}  read median/rms model {chk['median_model']:.2f}/{chk['rms_model']:.2f}")


if __name__ == "__main__":
    main()
