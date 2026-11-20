# P1 status note: validation stage (check-in 1 of 2)

Date: 2026-10-14/15 (overnight session). Branch: `claude/epic-ritchie-g3ret8`.
Full report: `reports/p1_validation.md`.

## Done

`camsim_common` (provenance-tagged `SensorModel` schema, HDF5 I/O with config and git hash, run
logging, RNG streams, placeholder scanner) and P1 through V1 to V3 with a convergence appendix.
Tests: 49 fast and 5 slow, all passing; two V1 runs give identical outputs. Sweep runner
(resumable, CPU-hours per point), Sobol wrapper and SLURM templates written but not used.

## Results

- **V1 (Manetsch, ORCA-Quest, 80 ms).** FN reproduced: 4.9e-5 (empirical rate distribution) to
  5.7e-5 (Gaussian rates) against 5.5(1)e-5 measured. The log histogram is reproduced over four
  decades. FP is under-predicted 14-fold (4.9e-6 against 7.01(8)e-5), so predicted fidelity is
  0.999973 against 0.9999374(8). The digitised histogram shows a background tail 2 to 3 times
  heavier than any Poisson-plus-read-noise model at 32 to 42 photons. We read this as a non-camera
  floor (e.g. atoms entering during imaging), not a model error. No-loss fidelity 0.999996,
  consistent with the paper's 0.99999(1).
- **V2 (Picken, Zyla 5.5, 40 ms).** The ADC gains fitted independently for the two traps agree
  (0.381 and 0.368 e-/ADU). Trap 1's signal matches R t eta QE to 5%; trap 2 is 32% brighter.
  The Poisson background gives 2.9x and 4.0x the paper's Gaussian readout error, so trap 1 is at
  9.6e-7, just above their "< 8e-7". The SNR table is reproduced except the EMCCD 852 nm column,
  which implies QE 0.55, not the stated 0.60.
- **V3.** Poisson limit exact to 1e-13. EMCCD stage Monte Carlo matches Basden Eq. (5)
  (F^2 = 1.97 at 512 stages, gain 1000; tends to 2). Collection, PSF and Langevin limits pass.
  The Airy fallback is adequate only below NA ~ 0.3 (FWHM +6% at NA 0.8).

## Sensitivities that matter for the Stage 1 specification

At fixed survival, the loss-time profile moves FN by -42% to +10% (model B versus A). The
low-rate tail of the site-to-site scattering-rate distribution moves FN by ~17%. Both act at the
1e-5 level, comparable to the camera terms, so production results will be quoted as A/B brackets.

## Decisions needed from you

1. Approve the production stage: eight camera presets, the reduced local sweep, then the first
   production status note.
2. Accept the V1 FP excess as an experimental floor outside the camera model (or point us to
   data on background-atom events)?
3. Please upload the Prime BSI Express datasheet PDF and, if available, the ORCA-Quest
   C15550-20UP read-noise specification. photometrics.com and teledynevisionsolutions.com block
   automated access, so these values stay USER/PLACEHOLDER until then.
4. The plan named branch `claude/sharp-dirac-i7j2vp`; this session pushes to
   `claude/epic-ritchie-g3ret8` as assigned.

40 placeholders are open (listed at the end of the report). None affects a validation residual
except the ORCA-Quest read noise, which is fitted.

## Amendment (2026-11-20)

Three statements above are corrected; the text above is left as sent. Details:
`reports/contact_log_v1_v2.md` (entries V1-7, V1-8, V2-6 and the resolution of 2026-11-20).

1. The FP excess as a non-camera floor is our attribution, not a tested result. V1 validates the
   FN rate, not the model's FP rate; the scope was set after the FP residual was seen.
2. The no-loss comparison (0.999996 against 0.99999(1)) is between our model and the paper's
   lossy-Poisson model on the same histogram, not against a measurement.
3. The SNR table is reproduced within rounding (0.05) except three of 18 cells: the two EMCCD
   852 nm cells at 100 and 1000 photons/px (consistent with QE 0.55, a hypothesis) and the Zyla
   QE 0.6 cell at 100 photons/px (+0.051, unexplained). The earlier text named only the EMCCD
   column, and the report used a cut of 0.06 while stating 0.05.
