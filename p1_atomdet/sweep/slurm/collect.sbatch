#!/bin/bash
# Merge checkpointed sweep points into one HDF5 after the array job (submit with --dependency=afterok:<jobid>).
#SBATCH --job-name=p1collect
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
set -euo pipefail
python -m p1_atomdet.sweep.grid --config "${1:?config yaml}" --collect
