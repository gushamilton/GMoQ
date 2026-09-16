#!/bin/bash
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --job-name=cmx1
#SBATCH --time=1:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH -o /user/work/fh6520/gmoq_covmixW/logs/s1_%A_%a.out
WS=/user/work/fh6520/gmoq_covmixW
source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh
conda activate fungwas-stage1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
NCH=10
I=$SLURM_ARRAY_TASK_ID
CELL=$((I / NCH))
CH=$((I % NCH))
REPS=$(python3 -c "import json;print(json.load(open(\"$WS/cells_mix.json\"))[$CELL][\"reps\"]//$NCH)")
python $WS/src/sim_stage1_mix.py --cell $CELL --chunk $CH --reps $REPS --out $WS/raw
