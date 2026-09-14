#!/bin/bash
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH -o /user/work/fh6520/skewcorr/logs/a2_%a.out
#SBATCH --array=0-3
source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh
conda activate fungwas-stage1
S=(normal skew heavy mildskew)
SH=${S[$SLURM_ARRAY_TASK_ID]}
python /user/work/fh6520/skewcorr/src/a2_typeI.py --shape $SH --out /user/work/fh6520/skewcorr/out/a2_$SH.tsv
