#!/bin/bash
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --job-name=mcprep
#SBATCH --time=1:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -o /user/work/fh6520/skewcorr/logs/mcprep.%j.out
WS=/user/work/fh6520/skewcorr
source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh
conda activate fungwas-stage1
python $WS/src/mc_cells.py $WS/cells.json
python -u $WS/src/mc_truth.py --cells $WS/cells.json --truthdir $WS/truth --out $WS/out/mc_exact.tsv
