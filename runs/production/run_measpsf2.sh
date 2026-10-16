#!/bin/sh
# Measured-PSF production runs, resumable. Waits for any running sweep to finish, then:
# presets (resume + collect), sensitivity, generic (reduced), Sobol.
set -e
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
while ps -eo args | grep -q "[p]1_atomdet.sweep.grid"; do sleep 20; done
for c in presets_sweep sensitivity generic_reduced sobol; do
  n=$(python3 -c "import yaml;print(yaml.safe_load(open('configs/production/$c.yaml'))['name'])")
  python3 -m p1_atomdet.sweep.grid --config configs/production/$c.yaml --out runs/production/$n --workers 4
  python3 -m p1_atomdet.sweep.grid --config configs/production/$c.yaml --out runs/production/$n --collect
done
