#!/bin/bash
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --job-name=cs2
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH -o /user/work/fh6520/covsource_run/logs/s2_%A_%a.out
WS=/user/work/fh6520/covsource_run
module load languages/R/4.5.1
export OMP_NUM_THREADS=2
Rscript --vanilla $WS/src/stage2_driver_cs.R $SLURM_ARRAY_TASK_ID $WS/raw $WS/out 2500
