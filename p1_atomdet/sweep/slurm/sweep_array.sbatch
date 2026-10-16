#!/bin/bash
# SLURM array template for P1 sweeps on the IIT Madras HPC.
# Partition, account and module names are site-specific: [CONFIRM?] with the HPC documentation.
#
#   sbatch --array=0-63 p1_atomdet/sweep/slurm/sweep_array.sbatch configs/sweep/<name>.yaml [OUT_DIR]
#
# OUT_DIR defaults to runs/sweeps/<name> (production runs pass runs/production/<name>).
#
# Each array task evaluates every N-th grid point (N = array size) with $SLURM_CPUS_PER_TASK
# workers. Finished points are checkpointed as JSON, so resubmitting the same array resumes.
#SBATCH --job-name=p1sweep
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%A_%a.out
##SBATCH --partition=[CONFIRM?]
##SBATCH --account=[CONFIRM?]

set -euo pipefail
CONFIG=${1:?config yaml}
OUT=${2:-}
N=${SLURM_ARRAY_TASK_COUNT:-1}
I=${SLURM_ARRAY_TASK_ID:-0}
export NUMBA_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# module load python/3.11   # [CONFIRM?] site module name
# source ~/venvs/camsim/bin/activate
mkdir -p logs
python -m p1_atomdet.sweep.grid --config "$CONFIG" ${OUT:+--out "$OUT"} --workers "${SLURM_CPUS_PER_TASK:-1}" --shard "$I/$N"
