#!/bin/bash
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --job-name=ns2
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH -o /user/work/fh6520/nscale_run/logs/s2_%A_%a.out
WS=/user/work/fh6520/nscale_run
SRC=/user/work/fh6520/wt_nscale/sim/nscale/src
module load languages/R/4.5.1
export OMP_NUM_THREADS=2
Rscript --vanilla $SRC/stage2_driver.R \
    $WS $WS/lib $SRC $SLURM_ARRAY_TASK_ID $WS/raw $WS/out 2500
