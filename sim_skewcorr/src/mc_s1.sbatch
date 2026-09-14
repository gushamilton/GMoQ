#!/bin/bash
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --job-name=mcs1
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH -o /user/work/fh6520/skewcorr/logs/s1_%a.out
#SBATCH --array=0-159%40
WS=/user/work/fh6520/skewcorr
export STAGE1_PATH=/user/work/fh6520/gls_coverage/tree/GMoQ_v04/stage1
source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh
conda activate fungwas-stage1
export OMP_NUM_THREADS=2
CELL=$(( SLURM_ARRAY_TASK_ID / 10 ))
CHUNK=$(( SLURM_ARRAY_TASK_ID % 10 ))
python -u $WS/src/mc_stage1.py --cells $WS/cells.json --cell $CELL --chunk $CHUNK \
    --reps 2000 --out $WS/s1
