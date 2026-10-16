#!/bin/sh
# Resume measured-PSF runs after a container restart: Sobol, then generic (reduced).
set -e
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for c in sobol generic_reduced; do
  n=$(python3 -c "import yaml;print(yaml.safe_load(open('configs/production/$c.yaml'))['name'])")
  python3 -m p1_atomdet.sweep.grid --config configs/production/$c.yaml --out runs/production/$n --workers 4
  python3 -m p1_atomdet.sweep.grid --config configs/production/$c.yaml --out runs/production/$n --collect
done
