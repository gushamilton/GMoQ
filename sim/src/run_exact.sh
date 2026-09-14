#!/bin/bash
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --job-name=fdrexact
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --array=1-36%36
#SBATCH --output=/user/work/fh6520/wt_dispersion/sim/logs/exact.%A_%a.out
#SBATCH --error=/user/work/fh6520/wt_dispersion/sim/logs/exact.%A_%a.err
set -euo pipefail
source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh
conda activate fungwas-stage1
cd /user/work/fh6520/wt_dispersion/sim/src
export OMP_NUM_THREADS=2

MAFS=(0.005 0.01 0.05 0.2 0.3 0.5)
SHAPES=(normal skew)
NS=(50000 100000 461320)

i=$((SLURM_ARRAY_TASK_ID - 1))
mi=$((i % 6)); i=$((i / 6))
si=$((i % 2)); i=$((i / 2))
ni=$((i % 3))

MAF=${MAFS[$mi]}; SH=${SHAPES[$si]}; N=${NS[$ni]}
OUT=/user/work/fh6520/wt_dispersion/sim/out/exact/z_${SH}_maf${MAF}_N${N}.tsv
mkdir -p $(dirname $OUT)
python fdr_exact.py --out $OUT --N $N --mafs $MAF --shapes $SH \
  --zladder 3,5,8,12,16,22,30,40,50,60
