# P1 interim technical report: single-atom detection versus camera parameters at 852 nm

CAQSAT lab, Department of Physics, IIT Madras. Status: interim, 2026-10-16. Covers `camsim_common`
and P1 through the validation stage, the camera presets and sensitivity study with the measured
V1 point-spread function (PSF), and two interim Sobol analyses. Supersedes the summary sections of
`reports/p1_production_note.md` (which uses a diffraction-limited PSF). Detailed results:
`reports/p1_validation.md`, `reports/p1_sobol_interim_n256.md`, `reports/p1_sobol_camera_n256.md`.
All results are SIM unless marked otherwise; every SIM value should be reviewed before submission.

## Summary

We built a semi-analytic model that maps a camera description (`SensorModel`) and the photon budget
and loss of a trapped Cs atom to the false-positive and false-negative rates of single-site
occupancy detection, and from them to the imaging time needed for a target infidelity. The model
reproduces the false-negative rate of Manetsch et al. (arXiv:2403.12021) to within 11% and the
readout statistics of Picken et al. (arXiv:1709.00946), but under-predicts the false-positive rate of
Manetsch et al. by a factor of 14, which we attribute to events outside a camera model.

In a reference scenario anchored on Manetsch et al. (NA 0.55, measured PSF of FWHM 1.70 um, 2.7% of
scattered photons reaching the camera), five low-noise CMOS configurations reach infidelity 1e-4 in
40 to 44 ms with a matched-filter classifier, and an EMCCD in multi-frame photon counting in 34 ms.
The Prime BSI Express in CMS mode needs 28.4 ms for 1e-3 and 41.3 ms for 1e-4, within 20% of the best
camera in the set. Cameras with QE at 852 nm below 0.25 need 2.5 to 4 times longer than the CMOS
group or do not reach 1e-4.

Over the camera-side parameter space, the imaging time is set by QE and read noise, which act
separably: at the reference background and 1 px per FWHM, t is proportional to 1/QE times a
read-noise factor that is 0.77 at 0.1 e-, 0.84 at 0.5 e-, 1 at 1.0 e- and 1.42 at 2.0 e- (for 1e-4).
Hence QE is the main lever for a Stage 1 sensor, and a read-noise median below 0.5 e- gains less
than 10%. Two inputs that carry little variance globally matter at the operating point: a 2% fraction
of pixels at five times the median read noise slows a 1 e- camera by 42%, and the loss per scattered
photon decides whether 1e-4 is reachable at all. NA and atom loss, which are properties of the
experiment, outweigh every camera parameter.

The main limitations are (i) the unexplained false-positive floor of the V1 experiment, (ii) an
assumed two-class read-noise tail in place of a measured distribution, and (iii) a single reference
scenario at 852 nm; Rb at 780 nm and array crosstalk are not yet treated.

## 1. Scope and what was run

The P1 brief asks three questions: (i) which sensor parameters set the infidelity and imaging time,
and with what sensitivity; (ii) what minimum specification reaches 1e-3 and 1e-4 within a given
imaging time and loss budget (the Stage 1 camera specification); and (iii) how the lab should operate
the Prime BSI Express. This report gives interim answers to (i) and (iii) and the elements of (ii).

We ran (a) validation stages V1 to V3 with a convergence study; (b) eleven camera presets at seven
pixel samplings, with the diffraction-limited and with the measured PSF; (c) a sensitivity study of
the preset results to their uncertain inputs; (d) a reduced generic CMOS-class sweep of 750 points
with the diffraction-limited PSF; (e) Sobol analyses at N = 256 of an 8-parameter design (camera
parameters, NA and loss) and of a 6-parameter camera-only design (NA and loss fixed). Table 7 lists
the cost of each.

## 2. Model and validation

Here, an atom present for the full exposure t yields mu = R t eta_cam QE detected photoelectrons,
with R the scattering rate, eta_cam = eta_coll(NA) eta_opt the fraction of scattered photons reaching
the sensor and QE the quantum efficiency at 852 nm. Loss during imaging truncates the exposure;
model A uses a constant loss probability per scattered photon, model B a 3D Langevin walk in the trap.
Per-pixel counts are Poisson in signal plus background, convolved with a Gaussian-mixture read-noise
model (core plus tail pixels) or, for an EMCCD, a gamma-distributed multiplication register with
clock-induced charge. Site classification uses a box sum, a matched filter with Fisher-optimal pixel
weights, or a per-pixel Bayesian likelihood ratio. Count distributions are computed exactly on a
lattice or, for linear statistics, from characteristic functions (within 2% of the lattice result in
optimal infidelity, about 10^3 times faster), and cross-checked against brute-force and importance-
sampled Monte Carlo at 1e-3. Thresholds are optimised per site and averaged over 3 x 3 sub-pixel site
positions (within 1.1% of a 9 x 9 average).

**Table 1. Validation results** (details and residuals in `reports/p1_validation.md`).

| stage | target (source) | model result | status |
|---|---|---|---|
| V1 | FN 5.50e-5, FP 7.01e-5, ORCA-Quest, 80 ms (Manetsch et al., arXiv:2403.12021v4) | FN 4.91e-5 (empirical rate spread) to 5.73e-5 (Gaussian); FP 4.85e-6 | FN reproduced; FP 14x low |
| V1 | dependence on loss-time profile at fixed survival 0.999864 | FN -42% (model B, late losses) to +10% (model B, early escape) relative to model A | quantified |
| V2 | Zyla 5.5, 40 ms, readout error < 8e-7 (Picken et al., arXiv:1709.00946v2) | fitted gains 0.381 and 0.368 e-/ADU (agree between traps); Poisson readout error 2.9 to 4.0 x their Gaussian estimate (9.55e-7 for trap 1) | reproduced; their bound marginal |
| V2 | sCMOS vs EMCCD SNR table (their Table I) | all cells to 0.05 except two EMCCD 852 nm cells, consistent with QE 0.55 instead of 0.60 | reproduced |
| V3 | Poisson limit; EMCCD excess noise F^2 -> 2 (Basden et al. 2003) | relative error 4.42e-14; F^2 = 1.972 at gain 1000, 512 stages | pass |

The V1 false-positive excess is visible directly in the digitised histogram: between 32 and 42
photons the measured background tail lies 2 to 3 times above any Poisson-plus-read-noise model
consistent with the background width. We read this as a floor set by events outside the camera, for
example atoms entering during imaging. This has a direct consequence for the present work: at the
1e-4 level, a real experiment may be limited by non-camera errors that no camera choice removes.

## 3. Reference scenario

`configs/production/reference_cs.yaml`, derived from the V1 fit: Cs D2 at 852 nm; scattering rate
72.3e3 s^-1; NA 0.55 with eta_opt = 0.327 (2.7% of scattered photons to the camera, as in V1); loss
2.35e-8 per scattered photon (V1 survival, model A); site-to-site rate spread 12.3% rms (V1 Fig. 4d);
background 4.76 e-/px/s at the V1 pixel (0.902 um in the object plane, QE 0.44), scaled with pixel
area and QE; filling 0.5. The PSF is a Gaussian of 1/e^2 radius 1.44 um at NA 0.55 (FWHM 1.70 um),
derived from the waist stated by Manetsch et al. with the pixel variance removed, and scaled as 1/NA
when NA is varied. Its FWHM is 2.1 times that of the diffraction-limited vectorial PSF (0.82 um).
Survival is not limiting in this scenario: the loss rate is 1.7e-3 s^-1, so survival stays above 0.99
up to 5.9 s (DERIVED: 0.01005 / (2.35e-8 x 72.3e3 s^-1)).

## 4. Camera presets

**Table 2. Presets at their best pixel sampling, measured PSF** (matched filter; photon counting uses
an SNR-optimal box). QE at 852 nm digitised from vendor typical curves (DATA, about +-1.5 to 2%).
Magnification for the V1 PSF; n/r: not reached within 300 ms with survival >= 0.99.

| camera | mode | QE 852 | read median/rms (e-) | best px/FWHM | magnification | t(1e-3) ms | t(1e-4) ms | t(1e-4) ms, diffraction-limited PSF |
|---|---|---|---|---|---|---|---|---|
| iXon Ultra 888 EXF | photon counting, 670 fps | 0.608 | EM x1000 | 1.4 | 10.7 | 21.9 | 34.3 | 29.6 |
| ORCA-Fusion BT | ultra quiet | 0.464 | 0.60/0.70 | 1.0 | 3.8 | 26.7 | 39.7 | 36.9 |
| ORCA-Quest 2 | ultra quiet (PNR) | 0.384 | 0.25/0.30 | 1.4 | 3.8 | 26.4 | 40.2 | 35.7 |
| Prime BSI Express | CMS (100 MHz) | 0.485 | 1.00/1.00 | 1.0 | 3.8 | 28.4 | 41.3 | 37.5 |
| ORCA-Quest 2 | standard | 0.384 | 0.39/0.43 | 1.0 | 2.7 | 27.9 | 42.6 | 37.9 |
| Marana 4.2B-6 | low noise 12-bit | 0.470 | 1.10/1.10 | 1.0 | 3.8 | 30.3 | 43.6 | 39.6 |
| iXon Ultra 888 EXF | analog, EM x1000 | 0.608 | EM x1000 | 3.0 | 23.0 | 30.2 | 48.1 | 42.5 |
| Prime BSI Express | Speed (200 MHz) | 0.485 | 1.60/1.60 | 0.7 | 2.7 | 34.0 | 49.3 | 46.3 |
| QHY268M (IMX571) | High Gain-2CMS, gain 100 | 0.208 | 1.19/1.19 | 0.7 | 1.6 | 68.8 | 104.2 | 94.4 |
| Zyla 5.5 | global shutter | 0.249 | 2.30/2.50 | 0.5 | 1.9 | 94.6 | 156.2 | 154.2 |
| Blackfly S BFS-U3-17S7M (IMX425) | high conversion gain | 0.186 | 4.69/4.69 | 0.35 | 1.9 | 179.1 | n/r | 299.6 |

As can be seen from Table 2, the measured PSF, although 2.1 times wider than the diffraction limit,
lengthens the imaging time by only 1 to 16%. This is because the best magnification changes little
(3.8 instead of 4.0 for the Prime BSI Express): the optimal sampling moves from about 0.5 to about
1 px per FWHM, and the wider PSF costs only the read noise and background of the additional pixels
in the classifier window. The ranking is unchanged except that the ORCA-Fusion BT and the ORCA-Quest 2
in PNR mode exchange places, a difference of 1.3%. Three remarks. (i) The iXon advantage requires
photon counting over many short frames (26 frames in 40 ms for a 128 x 128 crop); in analog mode the
excess noise factor (F^2 = 1.97) cancels its QE advantage. (ii) With the QE of 0.44 stated by
Manetsch et al. for their ORCA-Quest, instead of 0.384 from the Quest 2 typical curve, the Quest 2
times would fall by about 13% (t roughly proportional to 1/QE). (iii) Differences below about 5%
between presets are within the QE uncertainty and should not be used to rank them.

## 5. Sensitivity of the preset results

**Table 3. Changes to the preset results, measured PSF** (Prime BSI Express CMS at 1.0 px per FWHM
unless stated; the high loss, 5.6e-7 per photon, is the bound implied by Picken et al.).

| change | t(1e-3) ms | t(1e-4) ms | change in t(1e-4) |
|---|---|---|---|
| baseline | 28.4 | 41.3 | |
| QE 0.47 / 0.50 (digitisation range) | 29.3 / 27.5 | 42.6 / 40.0 | +3% / -3% |
| 2% of pixels at 5x median read noise | 37.9 | 58.7 | +42% |
| loss 5.6e-7 per photon (24 x reference) | 29.5 | n/r | 1e-4 not reachable |
| ORCA-Quest 2 PNR, loss 5.6e-7 per photon | 27.7 | n/r | 1e-4 not reachable |
| iXon PC, CIC all in image area / all in register | 23.7 / 23.1 | 37.7 / 36.0 | 5% between the limits |

Hence the read-noise tail, which no vendor datasheet in the set specifies, is the most consequential
camera input after QE and the median, and atom loss can make 1e-4 unreachable whatever the camera:
at the high loss, eps at 40 ms is 3.6e-4 for the Prime BSI Express and is set by atoms lost during
the exposure.

## 6. Global sensitivity

The brief asks for Sobol indices of the time to reach 1e-3 and 1e-4 at survival >= 0.99. We ran two
designs (SALib Saltelli sampling, no second-order terms, 95% bootstrap confidence intervals), each as
the first N = 256 base points of a nested N = 1024 design. The output is log10 t, with targets not
reached within the scan censored at the scan limit.

**8-parameter design** (`reports/p1_sobol_interim_n256.md`; 2560 points; scan to 1 s). Over QE 0.2 to
0.95, read-noise median 0.1 to 4 e-, tail fraction 0 to 5%, 0.35 to 2 px per FWHM, background 0.1 to
3 times the reference, dark current 1e-3 to 3 e-/px/s, NA 0.3 to 0.8 and loss 2e-9 to 6e-7 per
photon, NA dominates t(1e-3) (S_T = 0.56 +- 0.13), followed by read noise (0.22 +- 0.07) and QE
(0.14 +- 0.04). At 1e-4, 31% of points are censored, mostly at high loss and low NA, so the indices
mix the time to reach the target with whether it is reached; NA (0.51 +- 0.12), loss (0.34 +- 0.09)
and read noise (0.30 +- 0.08) are not separable at this N. Over 0.3 to 0.8, the collection fraction
for isotropic emission rises from 0.023 to 0.200 (DERIVED, `optics/collection.py`), a factor of 8.7,
against a factor of 4.75 across the QE range.

**Camera-only design** (`reports/p1_sobol_camera_n256.md`; 2048 points; NA 0.55 and loss fixed; scan
to 5 s, 3% censored at 1e-4). Table 4 gives the indices.

**Table 4. Camera-only Sobol indices of log10 t** (S1: first order; S_T: total; +- 95% CI).

| parameter | S1, t(1e-4) | S_T, t(1e-4) | S1, t(1e-3) | S_T, t(1e-3) |
|---|---|---|---|---|
| read-noise median | 0.60 +- 0.15 | 0.65 +- 0.12 | 0.54 +- 0.13 | 0.59 +- 0.09 |
| QE | 0.34 +- 0.10 | 0.38 +- 0.08 | 0.35 +- 0.09 | 0.34 +- 0.06 |
| background scale | 0.05 +- 0.05 | 0.06 +- 0.02 | 0.05 +- 0.05 | 0.06 +- 0.02 |
| pixels per FWHM | 0.00 +- 0.02 | 0.04 +- 0.02 | 0.00 +- 0.03 | 0.03 +- 0.01 |
| read-noise tail fraction | 0.02 +- 0.02 | 0.02 +- 0.02 | 0.03 +- 0.03 | 0.02 +- 0.01 |
| dark current | 0.00 | 0.00 | 0.00 | 0.00 |

The response is additive (sum of S1 = 0.97 at 1e-3 and 1.00 at 1e-4). However, the ranking of read
noise above QE reflects the sampled ranges: the read-noise range spans a factor of 40 and includes
cameras that are clearly unsuitable, and the log-uniform background range has its median (0.55 of the
reference) on the dark side, where read noise weighs more. For a specification, the local behaviour
at a realistic operating point is the relevant quantity (Table 5).

**Table 5. t(1e-3) / t(1e-4) in ms at the reference background**, 1.0 px per FWHM, no tail, dark
current 0.01 e-/px/s, same model as the Sobol designs.

| QE | 0.1 e- | 0.25 e- | 0.5 e- | 1.0 e- | 2.0 e- |
|---|---|---|---|---|---|
| 0.485 | 20.1 / 31.4 | 20.7 / 32.0 | 22.6 / 34.1 | 28.1 / 40.7 | 41.5 / 57.8 |
| 0.85 | 11.5 / 17.8 | 11.8 / 18.2 | 12.9 / 19.3 | 16.0 / 23.0 | 23.6 / 32.7 |

As can be seen from Table 5, the read-noise factor relative to 1.0 e- is the same at both QE values
to within 0.01 (for 1e-4: 0.77, 0.79, 0.84, 1, 1.42; for 1e-3: 0.72, 0.74, 0.80, 1, 1.48), and the QE
dependence is 1/QE to within 1% (DERIVED from Table 5). Hence, at this operating point,

t(1e-4) = 40.7 ms x (0.485 / QE) x f_4(sigma_r),

where sigma_r is the read-noise median and f_4 is the factor above. The same model reproduces the
preset results (Prime BSI Express 41.3 ms; ORCA-Quest 2 PNR 39.3 ms against 40.2 ms for the preset,
whose read-noise model includes a tail class).

## 7. Implications

**Stage 1 specification (provisional, reference scenario, matched filter).** Our choice of
requirements is a compromise between imaging time and what a purchased back-illuminated CMOS sensor
can deliver. First, QE at 852 nm sets the imaging time in proportion: each relative gain in QE
shortens the time by the same factor. Second, a read-noise median of 0.5 e- or below is within 10% of
an ideal sensor, whereas 1.0 e- costs about 30% and 2.0 e- about 80%. Third, pixels with read noise
far above the median must be rare: at a 1.0 e- median, 2% of pixels at 5 times the median cost 42%;
the acceptable fraction depends on the measured tail shape, which P4 must supply, and on whether such
pixels can be masked. Fourth, cooling need only keep dark current below about 1 e-/px/s, since the
dark charge at 40 ms is then below 0.04 e- per pixel (DERIVED); the Sobol designs show no effect of dark
current up to 3 e-/px/s. Fifth, the magnification should give 0.5 to 1.4 px per PSF FWHM (within 4%
of the best in the camera-only design). The quantitative specification (minimum QE and maximum read
noise for a given time budget) follows from the expression in Sec. 6 once the lab fixes the imaging
time budget; the generic sweep with the measured PSF (Sec. 9) will give the full maps.

**Prime BSI Express operating point (provisional).** CMS mode (1.0 e-), magnification giving about 0.7
to 1.0 px per PSF FWHM (2.7 to 3.8 for the V1 PSF of 1.70 um FWHM; to be rescaled with our own measured
PSF), matched filter with per-site weights and thresholds: t(1e-3) = 28.4 ms, t(1e-4) = 41.3 ms,
eps = 1.2e-4 at 40 ms in the reference scenario. Speed mode costs 19% (49.3 ms). A 2048 x 128 region
of interest reads out in 1.46 ms (128 rows at 11.4 us, USER datasheet line time); with continuous
illumination the exposure must exceed the imaging time by this amount. The read-noise tail of the
GSENSE2020BSI is not specified; a dark-frame stack from our camera is the most useful next input.

**Programme level.** At the operating point, raising QE from 0.485 to 0.85 shortens t(1e-4) by a
factor of 1.77 (Table 5). Raising NA from 0.55 to 0.70 increases the collection fraction by a factor
of 1.7 (0.083 to 0.143 for isotropic emission, DERIVED), a comparable gain from the optics. In our
opinion, the case for a custom thick-silicon sensor therefore rests on QE at 780 and 852 nm together
with cost and supply, and the gain in detection speed alone (about 1.8 times over a good commercial
back-illuminated sensor) should be weighed against an NA upgrade.

## 8. Limitations

(i) The V1 false-positive rate is under-predicted 14-fold; if real experiments carry such a floor,
camera improvements below it are invisible. (ii) Read-noise tails are a two-class model (tail pixels at
3 times the median for presets with a stated rms, 5 times in the generic and Sobol designs); P4 will
replace this with measured distributions. (iii) All production results are for one reference scenario:
Cs at 852 nm, NA 0.55, the measured V1 PSF, a uniform background, no neighbour crosstalk, filling 0.5.
(iv) The characteristic-function engine omits ADC rounding; for the ORCA-Quest 2 in PNR mode exact
rounding lowers eps by 7% at the operating point (1 to 2% in time), so the PNR results are slightly
conservative. (v) The Sobol results are at N = 256, with half-widths of 0.04 to 0.15 on the leading
indices. (vi) Model B (Langevin) is used only to bound the effect of the loss-time profile; the
production results use model A.

## 9. Remaining work and why

**Table 6. Outstanding computations and deliverables.**

| item | size | why it is needed |
|---|---|---|
| generic CMOS sweep with the measured PSF | 725 of 750 points left, at most 19 CPU-h | spec-to-fidelity maps in the (QE, read noise) plane consistent with Table 2; the existing maps use the diffraction-limited PSF |
| full 8-parameter Sobol, N = 1024 | 7680 points left, at most about 60 CPU-h | halves the confidence intervals; separates NA, loss and read noise at 1e-4 |
| full camera-only Sobol, N = 1024 | 6144 points left, at most about 30 CPU-h | final ranking of camera parameters for the specification |
| Rb at 780 nm, and the Cs/Rb compromise | not started | required by the brief; QE and PSF differ at 780 nm |
| neighbour crosstalk and array geometry | not started | required by the brief; matters when site spacing approaches the PSF width |
| Bayesian per-pixel classifier in production | implemented, validated, not yet in the sweeps | can exploit a measured read-noise map, the input that matters most after QE |
| EMCCD branch of the generic sweep | not started | EM gain, CIC and frame rate as free parameters |
| written Stage 1 specification and Prime BSI Express recipe | after the items above | the two deliverables of the brief |

The three runs are prepared for the IIT Madras HPC (`runs/production/slurm_submit_measpsf.sh`, 8 tasks
x 16 cores by default; partition and account [CONFIRM?]). Since this report, the time-to-target scan
stops after the last first crossing; re-evaluated points give identical times at 4 to 20 times lower
cost, so the sizes above are upper bounds.

## 10. Computational cost

**Table 7. CPU time used** (4-core cloud instance).

| study | points | CPU-h |
|---|---|---|
| presets, measured PSF | 77 | 0.66 |
| sensitivity, measured PSF | 26 | 0.17 |
| generic sweep, diffraction-limited PSF | 750 | 7.3 |
| 8-parameter Sobol, N = 256 | 2560 | 20.3 |
| camera-only Sobol, N = 256 | 2048 | 11.0 |

## Provenance ledger (principal values)

| value | location | type | pointer |
|---|---|---|---|
| FN 5.50e-5, FP 7.01e-5, fidelity 0.9999374(8), QE 0.44, 2.7%, 80 ms, NA 0.55 | Sec. 2, 3 | SOURCE | Manetsch et al., arXiv:2403.12021v4 (main text and Methods) |
| readout error < 8e-7, 40 ms, Zyla 5.5 | Sec. 2 | SOURCE | Picken et al., arXiv:1709.00946v2 |
| F^2 = 1/g + 2/g^{1/r} - 2/g^{1+1/r} | Sec. 2 | SOURCE | Basden, Haniff and Mackay 2003, Eq. 5 (`sources/README.md`) |
| FN 4.91e-5 to 5.73e-5, FP 4.85e-6, gains 0.381 / 0.368 e-/ADU, F^2 = 1.972 | Table 1 | SIM | `runs/validation/`, `reports/p1_validation.md` (review before submission) |
| rate 72.3e3 s^-1, eta_opt 0.327, loss 2.35e-8, background 4.76 e-/px/s | Sec. 3 | DERIVED / DATA | `configs/production/reference_cs.yaml` pointers |
| PSF 1/e^2 radius 1.44 um, FWHM 1.70 um; vectorial FWHM 0.82 um | Sec. 3 | DERIVED | `reference_cs.yaml`; `optics/array.py` PSFSpec.fwhm |
| survival >= 0.99 up to 5.9 s | Sec. 3 | DERIVED | 0.01005 / (2.35e-8 x 72.3e3 s^-1) |
| QE at 852 nm per preset | Table 2 | DATA / USER | `camsim_common/presets/*.yaml` (digitised vendor curves; USER-supplied PDFs) |
| preset times, magnifications | Table 2 | SIM | `runs/production/prod_presets/prod_presets.h5`, `prod_presets_difflim/` (review before submission) |
| sensitivity times | Table 3 | SIM | `runs/production/prod_sensitivity/prod_sensitivity.h5` (review before submission) |
| loss bound 5.6e-7 | Table 3 | DERIVED | Picken et al.: retention > 99% at 450 photons/ms, 40 ms (`configs/production/sensitivity.yaml`) |
| 8-parameter Sobol indices | Sec. 6 | SIM | `runs/production/prod_sobol/points/` rows 0 to 2559 (review before submission) |
| camera-only Sobol indices | Table 4 | SIM | `runs/production/prod_sobol_camera/points/` rows 0 to 2047 (review before submission) |
| operating-point times | Table 5 | SIM | `p1_atomdet.sweep.evaluate.sobol_point`, `configs/production/sobol_camera.yaml` fixed block (review before submission) |
| read-noise factors, 1/QE scaling | Sec. 6 | DERIVED | ratios of Table 5 entries |
| collection fraction 0.023, 0.083, 0.143, 0.200 at NA 0.3, 0.55, 0.7, 0.8 | Sec. 6, 7 | DERIVED | (1 - sqrt(1 - NA^2)) / 2, isotropic emission (`optics/collection.py`) |
| line time 11.4 us (CMS) | Sec. 7 | USER | Prime BSI Express datasheet rev. 2024-08-20 |
| CPU time | Table 7 | DATA | per-point cost records in the run directories |

## Open placeholders

The presets carry 116 open `[VALUE?]` entries (full list in `reports/p1_production_note.md`). Of these,
40 belong to the two validation presets and the model-B parameters and are discussed in
`reports/p1_validation.md`; 76 belong to the eleven production presets. The production engine does not
read 54 of the 76 (cosmic-ray parameters, PRNU, ADC offset, line time). The remaining 22 can affect the
results above: (i) read-noise median or rms where the vendor states a single value (6 entries; sets the
tail class of the preset); (ii) dark-current parameters (7; negligible by Sec. 6 unless dark current
exceeds about 1 e-/px/s); (iii) EMCCD register length and CIC split (4; the split is swept in Table 3);
(iv) ADC parameters (5; ADC rounding is omitted by the production engine, Sec. 8). In addition, the
SLURM partition and account are [CONFIRM?].
