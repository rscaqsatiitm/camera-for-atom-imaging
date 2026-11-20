# Contact log: validation datasets V1 (Manetsch) and V2 (Picken)

Date: 2026-11-19. Scope: every place in this repository where a model parameter, a model choice
or a tolerance was set by looking at V1 or V2 data, and what each contact removes from what that
dataset can validate. Method: reading of `p1_atomdet/validation/v1_manetsch.py`, `v2_picken.py`,
`report.py`, `configs/validation/*.yaml`, `configs/production/*.yaml`,
`camsim_common/presets/build_production_presets.py`, `p1_atomdet/production.py` and the reports,
at commit 0515d23. Rule (independence-bookkeeping skill): a value that influenced a parameter,
model choice, tolerance or stopping decision cannot certify the quantity it influenced, and the
log is append-only.

**Limit of this log.** The repository history starts at 152cc59 (2026-10-14); the V1 and V2
drivers were committed whole in 571ba3d and have not changed since. Choices made during
development before that commit are not recorded. Where the code does not show whether a choice
was made before or after the corresponding residual was seen, the entry says "order unknown" and
the contact is assumed (one-way rule: absence of a record is not independence).

The validation report (`reports/p1_validation.md`) is not modified by this log. Items where the
log disagrees with the report's wording are collected in Sec. 4 for the PI to decide.

## 1. V1: Manetsch et al., arXiv:2403.12021v4

Data used: Ext. Data Fig. 4b histogram moments and inset (DATA, digitised), Ext. Data Fig. 4d rate
histogram (DATA, digitised), and SOURCE numbers from the text (QE, 2.7%, 80 ms, filling, box, PSF
waist, survival S, measured FP, FN, F).

| ID | what was set | where | kind | consequence |
|---|---|---|---|---|
| V1-1 | background mean lambda0 = 18.66 e- per 7x7 box, from the Fig. 4b background mean | `v1_manetsch.py:90` | parameter (fit) | V1 cannot validate the background level or the background-peak position. |
| V1-2 | atom signal mu = 68.7 e-, from (atom mean - background mean) of Fig. 4b, divided by the mean present fraction | `v1_manetsch.py:99` | parameter (fit) | V1 cannot validate the photon budget (R t eta QE) or the atom-peak position. The implied rate check (72.3 against 69.7 kc/s, Table V1-c) is internal consistency of the paper's numbers, as the report states. |
| V1-3 | read noise sigma = 0.115 e- rms, solved so that the model background width equals the Fig. 4b width 4.40 | `v1_manetsch.py:95` | parameter (fit) | V1 cannot validate the background width or the camera read-noise model. Because 0.115 e- is negligible against Poisson noise, V1 does not test read-noise modelling at all. |
| V1-4 | Fig. 1 histogram comparison: peak positions and widths of the model come from V1-1 to V1-3, which were fitted to moments of the same Fig. 4b histogram whose inset is plotted | `v1_manetsch.py:90-99`, report Fig. 1 | shared input | The agreement of the first two moments of each peak in Fig. 1 is by construction. Only the shape beyond the second moment (tails, skew from loss and rate spread) is a prediction. The caption "The model is not fitted to this curve" is literally true but does not say this. |
| V1-5 | rate distribution: empirical Fig. 4d histogram as baseline, Gaussian of the same mean and sd as variant | `v1_manetsch.py:101-102, 108-109` | model choice, order unknown | The choice moved FN from 5.73e-5 (Gaussian) to 4.91e-5 (empirical) against 5.50e-5 measured, i.e. away from the measurement. This argues against tuning toward the FN target but does not prove the order. The FN comparison stands, but is quoted with both variants (as the report does). |
| V1-6 | loss hazard L = -ln S with S = 0.999864 from the same experiment's image sequence | `v1_manetsch.py:84-86` | parameter (SOURCE, same experiment) | Not fitted to FN, but S and the measured FN come from the same image sequences and the same model-free analysis. V1 FN therefore tests the camera-and-atom model given the measured loss, not the loss model. |
| V1-7 | interpretation of the FP residual (model 4.85e-6 against 7.01e-5): FP excess declared a non-camera floor; V1 then "regarded as validating the camera-and-atom part (FN, histogram shape, no-loss fidelity)" | report Sec. 3.3 | model choice / scope, made after the residual was seen | V1 cannot validate the claim that the FP floor is non-camera; that claim is a hypothesis with candidate causes (i) to (iii) untested. The restriction of V1's scope to FN was chosen after seeing FP, so V1 counts as a validation of FN only, not of the model's FP. |
| V1-8 | no-loss fidelity: our 0.999996 against the paper's 0.99999(1) | report Sec. 3.4, Table V1-c | comparison of two models on shared inputs | The paper's value is the output of its own lossy-Poisson fit to the same histogram, and our lambda0 and mu were fitted to that histogram (V1-1, V1-2). This is model-versus-model agreement on shared inputs, not validation against a measurement. |
| V1-9 | model B regimes: thermostat temperature "chosen to give measurable loss"; heating regime at 7 kc/s without cooling; base rate 72 kc/s from V1-2 | `configs/validation/v1_manetsch.yaml` `model_b` | model choice, order unknown | Model B is used only for loss-time shapes renormalised to S, and the report presents the -42% to +10% FN range as a sensitivity, not a validation. V1 cannot be used to select between model A and model B. |
| V1-10 | numerical defaults (ds = 0.01, n_tau = 48, h_stat, IS sample count) chosen from the convergence study run at the V1 fitted parameters (`convergence.yaml` copies them) | `configs/validation/convergence.yaml` | numerical tolerance | Contact with V1 fitted parameters, not with V1 targets. Counterfactual test: a different measured FN or FP would not have changed these defaults. Not a contamination of the FN comparison. |

**Comparisons in V1 that survive the counterfactual test** (no V1 target fed them): FN against
the measured 5.5(1)e-5 (given V1-5, V1-6); atom-peak sd 12.7 against 13.2 e- (Fig. 4b atom width
fed nothing); the tail shape of Fig. 1 beyond the second moments. Shared layers with the data:
the photon budget and background (V1-1, V1-2) and the loss (V1-6). No pass criterion for any of
these was written down before the comparison; the report gives residuals, not pass/fail.

## 2. V2: Picken et al., arXiv:1709.00946v2

Data used: Table II (background and atom peak mean and width per trap, filling, threshold c_min),
Table I (SNR), SOURCE numbers from the text.

| ID | what was set | where | kind | consequence |
|---|---|---|---|---|
| V2-1 | ADC gain g and background lambda0 per trap, closed-form from the Table II background mean and width | `v2_picken.py:78` | parameter (fit) | V2 cannot validate the background mean or width (the model reproduces 206.0/29.00 and 191.0/29.00 ADU by construction, as the report states). |
| V2-2 | counts assumed offset-subtracted, because mu0 = 206 counts "is below any plausible 9-pixel offset" | report Sec. 4.1, assumption (vi) | model choice made by looking at Table II | Table II cannot validate the offset assumption. |
| V2-3 | relative signal spread per trap (0.055, 0.085), fitted to the Table II atom-peak width | `v2_picken.py:90` | parameter (fit) | V2 cannot validate the atom-peak width, and the model acceptance (FN side) at c_min depends on this fit. |
| V2-4 | measured atom signal derived with the fitted gain g of V2-1 | `v2_picken.py:84` | shared input | The comparison with R t eta QE (1.05 and 1.32) is a genuine prediction of the signal from SOURCE numbers, but shares the gain fit. The trap-2 excess (32%) is explained only qualitatively, after the fact (standing wave); V2 does not validate that explanation. |
| V2-5 | readout error at c_min: model 9.55e-7 and 9.32e-8 against the paper's Gaussian estimate 3.32e-7 and 2.35e-8 | report Table V2-b | comparison of two models on shared inputs | Both use the Table II mean and width (by construction identical, V2-1). Picken et al. did not measure this error rate; their "< 8e-7" is their Gaussian estimate. V2 therefore contains no measured error rate against which the model's FP can be validated. |
| V2-6 | SNR table tolerance: cells with abs(chain - Table I) > 0.06 are listed as exceptions; the text states the rest agree "within rounding (0.05)" | `report.py:46-47`, report Sec. 4.3 and Summary item 2 | tolerance, order unknown | The cut 0.06 lies above the a priori rounding half-unit 0.05 of one-decimal printed values. Under 0.05 the Zyla QE 0.6, 100 photons/px cell (deviation 0.0513) is also an exception: three cells deviate, not two. The printed "0.05" is that 0.0513 rounded to two decimals. |
| V2-7 | "EMCCD 852 nm cells consistent with QE = 0.55", obtained by inverting Eq. (1) on the two deviating cells | report Sec. 4.3 | post-hoc explanation | Those cells cannot validate QE = 0.55; it is a hypothesis about a misprint. |

**Comparisons in V2 that survive the counterfactual test:** the agreement of the two
independently fitted gains (0.381 and 0.368 e-/ADU; each fit used only its own trap; shared layers:
the 2.2 e- read noise, the offset assumption V2-2 and the Poisson-plus-Gaussian model); the atom
signal against R t eta QE (shares g, V2-4); the chain against Eq. (1) (an implementation check of
our sensor chain against a formula, not against data). No pass criterion was written down before
these comparisons.

## 3. Downstream contacts (production stage)

V1 and V2 numbers were carried into the production configuration. For these quantities V1 or V2
cannot serve as validation of production results; re-running V1 through the production (CF)
engine at the reference scenario would be a self-consistency check of the engine, not a
validation.

| ID | production quantity | from | where |
|---|---|---|---|
| P-1 | scattering rate 72300 s^-1 | V1-2 (fit to Fig. 4b means) | `configs/production/reference_cs.yaml` `rate_s`; `production.py:31` |
| P-2 | optical transmission eta_opt = 0.327 | V1 SOURCE 2.7% (text, not a fit) | `reference_cs.yaml` `eta_opt` |
| P-3 | measured PSF waist 1.44 um | V1 SOURCE 1.7 px (text) | `reference_cs.yaml` `psf_waist_m` |
| P-4 | loss per photon 2.35e-8 | V1-6 and V1-2 | `reference_cs.yaml` `p_loss_per_photon` |
| P-5 | site-to-site rate spread 0.1226 | V1 Fig. 4d (DATA) | `reference_cs.yaml` `rate_rel_sd` |
| P-6 | background 4.76 e-/px/s | V1-1 (Fig. 4b background mean) | `reference_cs.yaml` `bg_rate_e_per_px_s` |
| P-7 | Zyla 5.5 production preset ADC gain 0.375 e-/ADU | V2-1 (mean of the two trap fits) | `build_production_presets.py:145` |
| P-8 | high-loss sensitivity case p_loss = 5.6e-7 | V2 SOURCE retention > 99% (text) | `configs/production/sensitivity.yaml` `p_loss_high` |

Hence V1 can no longer validate production predictions of the photon budget, background, rate
spread or loss in the reference scenario; P-2 and P-3 are SOURCE text values and do not burn the
V1 measured error rates. V2's Table II background cannot validate the Zyla production preset.

## 4. Result, for the PI

(i) Every fitted parameter in V1 and V2 is already identified as fitted in the validation report
(Table V1-a, Sec. 4.1); no fitted parameter is presented as a prediction.
(ii) V1 validates the FN rate (4.91e-5 to 5.73e-5 against 5.5(1)e-5) and the atom-peak width, given
the measured loss; it does not validate FP (V1-7), the read-noise model (V1-3), or, through the
production anchoring (P-1 to P-6), any production quantity derived from V1.
(iii) V2 contains no measured error rate (V2-5). It validates the consistency of the gain fits and
the photon budget of trap 1, not the error-rate model.
(iv) Points where this log disagrees with the report's wording, not changed: V1-4 (Fig. 1 caption
omits that the peak moments were fitted to the same histogram), V1-8 (no-loss fidelity is
model-versus-model), V2-6 (three SNR cells exceed the rounding tolerance, not two; the cut in
`report.py` is 0.06).
(v) No V1 or V2 comparison had a pass criterion written down in advance. Future validation
targets need one (CLAUDE.md rule 10).

Appending: add new entries below with a date; do not delete entries.
