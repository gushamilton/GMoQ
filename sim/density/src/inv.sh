#!/bin/bash
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --time=0:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
set -e
source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh
conda activate fungwas-stage1
export PYTHONPATH=/user/work/fh6520/kde_lumpy/tree_bw/stage1 OMP_NUM_THREADS=1
cd /user/work/fh6520/kde_lumpy/src
python invariance.py
