#!/bin/bash
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --job-name=mv_s2
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=/user/work/fh6520/meanvar/logs/s2_%A_%a.out
set -euo pipefail
SCALES=(raw log rint)
I=$((SLURM_ARRAY_TASK_ID - 1))
S=${SCALES[$((I / 40))]}
C=$(printf "%04d" $((I % 40)))
module load languages/R/4.5.1
Rscript --vanilla /user/work/fh6520/wt_pooled-guards/sim_meanvar/src/mv_stage2.R \
  /user/work/fh6520/meanvar/mc/raw1/"$S" "${S}_chunk${C}" \
  /user/work/fh6520/meanvar/mc/s2/"$S"
