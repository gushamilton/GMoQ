#!/usr/bin/env bash
#SBATCH --job-name=pc-mc
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=4:00:00
#SBATCH --output=/user/work/fh6520/gmoq_power_conserv/logs/mc_%A_%a.log
set -euo pipefail
source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh
conda activate fungwas-stage1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
cd /user/work/fh6520/gmoq_power_conserv/scripts
python mc_null.py "$MCMODE" "$SLURM_ARRAY_TASK_ID" "$NTASK" "$MPER"
