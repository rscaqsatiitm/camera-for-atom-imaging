#!/bin/sh
# Evaluate a Sobol design only up to row NROWS (a nested Saltelli prefix: NROWS = N_base x (D + 2)).
# Points are taken in index order and checkpointed, so rerunning after a stop resumes.
#   sh runs/production/run_sobol_prefix.sh [NROWS] [CONFIG_STEM]
#   defaults: 2560 sobol  (N=256 of the 8-parameter design); camera-only N=256: 2048 sobol_camera
set -e
NROWS=${1:-2560}
STEM=${2:-sobol}
CFG=configs/production/$STEM.yaml
NAME=$(python3 -c "import yaml;print(yaml.safe_load(open('$CFG'))['name'])")
OUT=runs/production/$NAME
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
LEFT=$(python3 -c "
from pathlib import Path
d = Path('$OUT/points')
print(sum(not (d / f'{i:07d}.json').exists() for i in range($NROWS)))")
echo "$(date -u +%FT%H:%M) $NAME prefix $NROWS: $LEFT points left"
[ "$LEFT" -gt 0 ] && python3 -m p1_atomdet.sweep.grid --config "$CFG" --out "$OUT" --workers 4 --limit "$LEFT"
echo "$(date -u +%FT%H:%M) $NAME prefix $NROWS done"
