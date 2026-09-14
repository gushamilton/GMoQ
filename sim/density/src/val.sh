#!/bin/bash
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --time=0:20:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
set -e
source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh
conda activate fungwas-stage1
PY=$CONDA_PREFIX/bin/python
export PYTHONPATH=/user/work/fh6520/kde_lumpy/tree/stage1
cd /user/work/fh6520/kde_lumpy/src
rm -f /user/work/fh6520/kde_lumpy/test_col/s2_*.bin
$PY keepset.py --dir /user/work/fh6520/kde_lumpy/test_col
module load languages/R/4.5.1
Rscript --vanilla stage2.R /user/work/fh6520/kde_lumpy/test_col full
Rscript --vanilla stage2.R /user/work/fh6520/kde_lumpy/test_col tiedrop /user/work/fh6520/kde_lumpy/test_col/keep.json
$PY analyse.py --dir /user/work/fh6520/kde_lumpy/test_col --out /user/work/fh6520/kde_lumpy/summ/TEST.json
