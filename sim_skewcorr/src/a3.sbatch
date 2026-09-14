#!/bin/bash
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH -o /user/work/fh6520/skewcorr/logs/a3_%a.out
#SBATCH --array=0-2
source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh
conda activate fungwas-stage1
S=(normal skew heavy)
python -u /user/work/fh6520/skewcorr/src/a3_power.py --shape ${S[$SLURM_ARRAY_TASK_ID]}
