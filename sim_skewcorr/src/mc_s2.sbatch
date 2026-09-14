#!/bin/bash
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --job-name=mcs2
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH -o /user/work/fh6520/skewcorr/logs/s2_%a.out
#SBATCH --array=0-15
WS=/user/work/fh6520/skewcorr
module load languages/R/4.5.1
Rscript --vanilla $WS/src/mc_stage2.R $SLURM_ARRAY_TASK_ID $WS/s1 $WS/s2 2000
