#!/bin/bash
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --job-name=mv_nbg
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=2 --mem=24G --time=1:00:00
#SBATCH --output=/user/work/fh6520/meanvar/logs/nbg_%j.out
set -euo pipefail
export OMP_NUM_THREADS=2
source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh
conda activate fungwas-stage1
python /user/work/fh6520/meanvar/tmp/nbg.py
