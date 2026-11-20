# CLAUDE.md

Simulators that map camera parameters to single-atom imaging fidelity (CAQSAT lab, IIT Madras).
Layout, install and reproduction commands: `README.md`. Validation results: `reports/p1_validation.md`.

## Rules this repository already follows

Each rule cites where it is stated or enforced. Do not weaken any of them.

1. **Provenance on every number.** Every value carries one of USER, SOURCE, SIM, DATA, DERIVED or
   PLACEHOLDER (`camsim_common/provenance.py`). Every tag except PLACEHOLDER requires a pointer
   (enforced by the `ProvValue` validator). A PLACEHOLDER is rendered `[VALUE?: approx. X, verify]`.
   Recalled values are PLACEHOLDER until a source is retrieved and read.
2. **Sources are recorded, not recalled.** Each source read has its version, URL and SHA-256 of the
   PDF in `sources/README.md`; PDFs are not committed. Figure digitisation is DATA, with method and
   resolution recorded there.
3. **The open `[VALUE?]` list is generated, never hand-written.** `python -m camsim_common.placeholders`
   scans `camsim_common/presets/` and `configs/`; the report generators call the same scan, so the
   lists cannot drift. As of 2026-11-19, `reports/p1_validation.md` lists the 40 entries of the
   validation-stage files (the two validation presets and `configs/validation/`), and
   `reports/p1_production_note.md` lists the 116 entries of the full default scan. Both are checked
   by `tests/test_placeholder_lists.py`; `p1_atomdet/validation/report.py` scans only the
   validation-stage files (`VALIDATION_PATHS`), so regenerating the report keeps the 40-entry list. Do not put `[VALUE?` text or `prov: PLACEHOLDER` into new
   files under `configs/` or `camsim_common/presets/` unless the open list is meant to grow.
4. **Report numbers come from run outputs.** `python -m p1_atomdet.validation.report` writes the
   validation report from `runs/validation/`; SIM entries in provenance ledgers are marked
   "review before submission".
5. **Reproducibility.** Seeds are fixed in the configs; per-task RNG streams come from
   `camsim_common.rng`; HDF5 outputs embed the resolved config and the git hash; `runlog.jsonl`
   records wall time, CPU time and peak memory; two runs give identical datasets
   (`camsim_common.io.datasets_equal`).
6. **Fitted, fixed and predicted quantities are kept apart.** The validation report states which
   parameters were fitted to which data (Table V1-a, Sec. 4.1) and treats everything else as a
   prediction. The contact log `reports/contact_log_v1_v2.md` records every place where V1 or V2
   data set a parameter, model choice or tolerance; a dataset does not count as validation for a
   quantity it helped set.
7. **Semi-analytic error rates are cross-checked by independent routes, as tests.** Brute-force
   Monte Carlo (agreement within the 95% Clopper-Pearson interval at error rates near 1e-3) and
   exponentially tilted importance sampling check the lattice engine; the lattice engine checks the
   CF engine used in production (`tests/test_counts.py`, `tests/test_cf_production.py`). Long checks
   carry `@pytest.mark.slow`.
8. **Convergence is shown, not assumed.** Numerical defaults (lattice step, s-grid step `ds`, loss
   quadrature nodes, super-sampling, PSF quadrature, IS sample count, Langevin step, Poisson
   truncation) are justified by the convergence study (report Appendix A).
9. **Production sweeps follow a check-in with the PI** (called rule 7 in the reports). Sweeps write
   one JSON checkpoint per point, are resumable, and record CPU time per point.

## Rules added for verification work (2026-11-19)

10. **Fix thresholds before running the check they apply to.** Write the tolerance or pass
    criterion down (in the test, config or a dated file, committed) before the first run. Never
    loosen a threshold after seeing a result; a check that needs its threshold moved has found
    something, so find the cause.
11. **Every new check is broken once on purpose and seen to fail.** Record the deliberate break,
    what failed, and the restore in `reports/check_breaks.md`.
12. **Negative controls.** Each classifier's fidelity check has a test in which a deliberately wrong
    threshold must fail it (`tests/test_negative_controls.py`). For integer scores (EMCCD photon
    counting) the shift is at least one score step, since a 10% shift is below one click (finding F1
    in `reports/check_breaks.md`).
13. **Reserved points before a sweep.** Before a production sweep, the reserved validation points
    and their pass criteria are written in a dated file and committed before the sweep starts
    (`configs/production/reserved_validation_*.yaml`).
14. **Time a small sweep and extrapolate before launching a production sweep.** Run a pilot with the
    same evaluator, config and machine, on points drawn from the same design (not its easy corner).
    Fit the per-point rate from the steady-state window (drop the first points, which carry numba
    compilation). Before launch, append to the sweep log one line with: pilot size, window, measured
    CPU-s per point, number of points, cores, projected wall time and the scaling assumption. If the
    projection exceeds the available slot (background-task limit or SLURM wall time), do not launch:
    split, shrink or re-plan first. Re-measure after any change of machine or parameter range
    (the per-point costs in `runs/production/slurm_submit_measpsf.sh` were measured in the cloud
    container, so they do not transfer to the cluster).

## Verification skills and what they mean here

`.claude/skills/` holds four skills vendored verbatim from BootLoops-ai/skills (record in
`.claude/skills/UPSTREAM.md`, licence in `NOTICE`): `planted-truth`, `acceptance-gate`,
`independence-bookkeeping`, `timing-discipline`. They were written for high-precision analytic
results. In this repository, read their terms as follows.

- **"Digits"** means stability of a quantity under doubled Monte Carlo samples or a refined count
  grid. For lattice or CF results: the relative change of FP, FN or infidelity when `h_stat`, `ds`
  or the loss quadrature is halved or doubled (report Appendix A). For MC and IS results: the
  change when the sample count is doubled, compared with the stated standard error. For Sobol
  indices: the change from the nested N to 2N Saltelli design, compared with the bootstrap
  interval. Quote the number, never an adjective. "Two-precision stable" means the agreement holds
  at both resolutions.
- **"Reserved points"** means held-out validation data: measurements or design points that entered
  no fit, tuning decision, model choice or tolerance. Check `reports/contact_log_v1_v2.md` before
  citing V1 or V2 as validation of a quantity. For sweeps, reserved points are fixed in
  `configs/production/reserved_validation_*.yaml` before the sweep.
- **"Independent route"** means a second method: brute-force Monte Carlo (`counts/mc.py`) or
  importance sampling (`counts/importance.py`) against the semi-analytic distributions (lattice:
  `counts/pmf.py`, `counts/lattice.py`; CF: `counts/cf.py`). The routes share these layers, and an
  error in them is invisible to the comparison: (i) the `SiteModel` inputs, including the pixel
  light fractions from `optics/array.py` and the PSF, the background and the signal `mu`;
  (ii) `TauDistribution` parameters (MC samples the exact truncated exponential when `hazard` is
  set, but samples the lattice's own nodes for model-B shapes); (iii) the readout parameter
  objects, including `EMCCDReadout.f2` from `enf2_finite` and the stage-wise `register_mc`;
  (iv) the classifier definition, since matched-filter weights (`optimal_linear_weights`) and Bayes
  LLR tables (`bayes_llr_tables`) are computed by the lattice code and passed to the MC, and the MC
  is evaluated at the lattice's optimal threshold; (v) numpy and scipy. The MC therefore checks
  the distribution of a given score at a given threshold, not the optimality of the weights or of
  the threshold, nor the optics. The CF-versus-lattice check shares (i), (ii) and the weights.
