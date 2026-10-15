#!/bin/sh
# Evaluate the Sobol design only up to row NROWS (default 2560 = N=256 Saltelli prefix with D=8).
# Points are taken in index order and checkpointed, so rerunning after a stop resumes.
set -e
NROWS=${1:-2560}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
LEFT=$(python3 -c "
import sys; from pathlib import Path
d = Path('runs/production/prod_sobol/points')
print(sum(not (d / f'{i:07d}.json').exists() for i in range($NROWS)))")
echo "$(date -u +%FT%H:%M) sobol prefix $NROWS: $LEFT points left"
[ "$LEFT" -gt 0 ] && python3 -m p1_atomdet.sweep.grid --config configs/production/sobol.yaml \
  --out runs/production/prod_sobol --workers 4 --limit "$LEFT"
echo "$(date -u +%FT%H:%M) sobol prefix $NROWS done"
