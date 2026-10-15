#!/bin/sh
# Production runs with the measured V1 PSF (resumable). Order: presets, sensitivity, Sobol.
set -e
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for c in presets_sweep sensitivity sobol; do
  n=$(python3 -c "import yaml;print(yaml.safe_load(open('configs/production/$c.yaml'))['name'])")
  python3 -m p1_atomdet.sweep.grid --config configs/production/$c.yaml --out runs/production/$n --workers 4
  python3 -m p1_atomdet.sweep.grid --config configs/production/$c.yaml --out runs/production/$n --collect
done
