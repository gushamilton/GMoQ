#!/bin/bash
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --job-name=mv_s1
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --time=2:00:00
#SBATCH --output=/user/work/fh6520/meanvar/logs/s1_%A_%a.out
set -euo pipefail
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
SCALES=(raw log rint)
I=$((SLURM_ARRAY_TASK_ID - 1))
S=${SCALES[$((I / 40))]}
C=$((I % 40))
source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh
conda activate fungwas-stage1
python /user/work/fh6520/wt_pooled-guards/sim_meanvar/src/mv_sim_stage1.py \
  --scale "$S" --N 50000 --chunk "$C" --reps 10 --nbg 200000 \
  --out /user/work/fh6520/meanvar/mc/raw1/"$S"
