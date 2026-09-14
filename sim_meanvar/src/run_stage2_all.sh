#!/bin/bash
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --job-name=mv_s2a
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=2 --mem=16G --time=2:00:00
#SBATCH --output=/user/work/fh6520/meanvar/logs/s2all_%j.out
set -euo pipefail
module load languages/R/4.5.1
SRC=/user/work/fh6520/wt_pooled-guards/sim_meanvar/src
for S in raw log rint; do
  for f in /user/work/fh6520/meanvar/mc/raw1/$S/${S}_chunk*.json; do
    TAG=$(basename "$f" .json)
    Rscript --vanilla $SRC/mv_stage2.R /user/work/fh6520/meanvar/mc/raw1/$S "$TAG" \
      /user/work/fh6520/meanvar/mc/s2/$S 2>&1 | grep -v "GLS tau_cov diagnostics" | tail -1
  done
done
Rscript --vanilla $SRC/mv_summarise.R /user/work/fh6520/meanvar/mc/s2 \
  /user/work/fh6520/meanvar/mc/summary
