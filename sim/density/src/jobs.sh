#!/bin/bash
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
set -e
source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh
conda activate fungwas-stage1
export PYTHONPATH=${TREE:-/user/work/fh6520/kde_lumpy/tree}/stage1
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd /user/work/fh6520/kde_lumpy/src
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")
echo "TASK $SLURM_ARRAY_TASK_ID : $LINE"
eval "python run_cell.py $LINE"
