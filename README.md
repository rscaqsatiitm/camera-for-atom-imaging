# camera-for-atom-imaging

Simulators that map camera parameters to single-atom imaging fidelity, for the CAQSAT lab
(Department of Physics, IIT Madras). Four linked projects share one camera description:

| package | content | status |
|---|---|---|
| `camsim_common` | `SensorModel` schema (YAML/JSON, every number a provenance-tagged `ProvValue`), HDF5 I/O with embedded config and git hash, run logging, reproducible RNG streams, placeholder scanner | done |
| `p1_atomdet` | end-to-end atom-detection model: fluorescence and loss (models A and B), collection and PSF, background, sensor chains (CMOS/sCMOS/qCMOS, EMCCD), semi-analytic count distributions, classifiers, fidelity | validation stage (V1 to V3) done; production sweeps await the check-in |
| `p2_bsi`, `p3_texture`, `p4_readout_thermal` | thick BSI pixel, light trapping, readout noise and thermal design | stubs |

## Install and test

```bash
pip install -e ".[test,sweep]"
pytest -q -m "not slow"     # fast suite (~1 min)
pytest -q -m slow           # brute-force MC and importance-sampling cross-checks (~10 min)
```

## Reproduce the validation report

```bash
python -m p1_atomdet.validation.run --config configs/validation/v1_manetsch.yaml \
    configs/validation/v2_picken.yaml configs/validation/v3_analytic.yaml configs/validation/convergence.yaml
python -m p1_atomdet.validation.report          # writes reports/p1_validation.md and figures
python -m camsim_common.placeholders            # open [VALUE?] list (must match the report)
```

Outputs go to `runs/validation/` (HDF5 with the resolved config and git hash, JSON tables,
`runlog.jsonl` with wall time, CPU time and peak memory). Seeds are fixed in the configs; two
runs give identical datasets (`camsim_common.io.datasets_equal`).

## How the P1 numbers are computed

Error rates of 1e-4 to 1e-6 are computed semi-analytically, not by brute-force sampling:
per-pixel output distributions (Poisson electrons, read-noise mixture, ADC rounding, or the
EMCCD register) are built on a lattice, mapped through the classifier score (box sum, matched
filter, or per-pixel log-likelihood ratio), convolved over the pixels of a site, and mixed over
the atom-signal multiplier s = (fraction of the exposure before loss) x (rate multiplier).
Brute-force Monte Carlo (agreement within the 95% Clopper-Pearson interval at ~1e-3) and
exponentially tilted importance sampling (deep tails) check the lattice engine; both are tests.

## Sweeps

`python -m p1_atomdet.sweep.grid --config <yaml> --workers N [--shard i/n]` evaluates a grid
with one JSON checkpoint per point (resumable; CPU-hours recorded per point).
`p1_atomdet/sweep/slurm/` holds SLURM array templates for the IIT Madras HPC; site-specific
fields are marked `[CONFIRM?]`. Sobol analysis: `p1_atomdet.sweep.sobol`.

## Provenance rules

Every number carries one of USER, SOURCE, SIM, DATA, DERIVED or PLACEHOLDER (rendered as
`[VALUE?: approx. X, verify]`). Sources read for the validation stage, with versions and file
hashes, are listed in `sources/README.md`.
