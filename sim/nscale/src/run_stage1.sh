#!/bin/bash
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --job-name=ns1
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH -o /user/work/fh6520/nscale_run/logs/s1_%A_%a.out
WS=/user/work/fh6520/nscale_run
SRC=/user/work/fh6520/wt_nscale/sim/nscale/src
TREE=/user/work/fh6520/wt_nscale/stage1
source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh
conda activate fungwas-stage1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
CELL=$SLURM_ARRAY_TASK_ID
NP=16
REPS=$(python3 -c "import json;print(json.load(open('$WS/cells.json'))[$CELL]['reps'])")
PER=$((REPS / NP))
mkdir -p $WS/raw
for i in $(seq 0 $((NP-1))); do
  python $SRC/sim_stage1.py --ws $WS --tree $TREE --cell $CELL \
      --chunk $i --reps $PER --out $WS/raw &
done
wait
echo "CELL $CELL done: $(ls $WS/raw/cell$(printf %02d $CELL)_chunk*.json | wc -l) chunks"
