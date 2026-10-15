#!/bin/bash
# Submit the measured-PSF production runs (Sobol N=1024 and the reduced generic sweep) as SLURM arrays,
# each followed by a dependent collect job that writes runs/production/<name>/<name>.h5.
#
#   bash runs/production/slurm_submit_measpsf.sh            # from the repository root on the cluster
#   NTASKS=8 CPUS=16 bash runs/production/slurm_submit_measpsf.sh
#
# Before submitting: copy the local checkpoints so finished points are skipped (they are gitignored):
#   rsync -a runs/production/prod_sobol/points/ <cluster>:<repo>/runs/production/prod_sobol/points/
#   rsync -a runs/production/prod_generic_reduced/points/ <cluster>:<repo>/runs/production/prod_generic_reduced/points/
# and install the pinned environment: pip install -e '.[sweep]' -r runs/production/requirements-lock.txt
#
# Cost per point measured in the cloud container (4 workers): Sobol mean 28 CPU-s, max 221 CPU-s,
# peak RSS 1.8 GB per worker; generic (first 25 points) mean 92 CPU-s, peak RSS 0.43 GB.
# Totals: Sobol about 80 CPU-h, generic at most about 19 CPU-h. With the defaults below
# (8 tasks x 16 CPUs) each Sobol task needs about 40 min, so 4 h of wall time leaves a wide margin.
# Partition and account are site-specific: pass them via SBATCH_PARTITION / SBATCH_ACCOUNT.
set -euo pipefail
NTASKS=${NTASKS:-8}
CPUS=${CPUS:-16}
MEM=${MEM:-48G}          # 16 workers x 1.8 GB peak, with margin
TIME=${TIME:-04:00:00}
S=p1_atomdet/sweep/slurm
mkdir -p logs
for c in sobol generic_reduced; do
  cfg=configs/production/$c.yaml
  name=$(python3 -c "import yaml;print(yaml.safe_load(open('$cfg'))['name'])")
  out=runs/production/$name
  jid=$(sbatch --parsable --job-name="p1_$name" --array=0-$((NTASKS-1)) --cpus-per-task="$CPUS" \
        --mem="$MEM" --time="$TIME" "$S/sweep_array.sbatch" "$cfg" "$out")
  cid=$(sbatch --parsable --job-name="p1_${name}_collect" --dependency=afterok:"$jid" \
        "$S/collect.sbatch" "$cfg" "$out")
  echo "$name: array $jid, collect $cid"
done
