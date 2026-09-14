#!/bin/bash
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --job-name=poolsalv
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH --output=/user/work/fh6520/wt_dispersion/sim/logs/salvage.%j.out
#SBATCH --error=/user/work/fh6520/wt_dispersion/sim/logs/salvage.%j.err
set -euo pipefail
module load languages/R/4.5.1
export POOLCAL_RLIB=/user/work/fh6520/poolcal/rlib
Rscript --vanilla /user/work/fh6520/wt_dispersion/sim/src/analyse.R \
  /user/work/fh6520/poolcal/rds \
  /user/work/fh6520/wt_dispersion/sim/out/mc
