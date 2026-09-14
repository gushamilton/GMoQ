#!/bin/bash
#SBATCH --account=sscm013902
#SBATCH --partition=compute
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
set -e
source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh
conda activate fungwas-stage1
PY=$CONDA_PREFIX/bin/python
export PYTHONPATH=/user/work/fh6520/kde_lumpy/tree/stage1 OMP_NUM_THREADS=4
cd /user/work/fh6520/kde_lumpy/src
D=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")
NAME=$(basename "$D"); RUN=$(basename $(dirname "$D"))
OUT=/user/work/fh6520/kde_lumpy/col/${RUN}__${NAME}
echo "=== qdedup $D"
$PY collate.py --dir "$D" --out "$OUT"
$PY keepset_q.py --dir "$OUT"
module load languages/R/4.5.1
if [ -s "$OUT/keep_q.json" ]; then Rscript --vanilla stage2.R "$OUT" qdedup "$OUT/keep_q.json"; fi
$PY analyse.py --dir "$OUT" --out /user/work/fh6520/kde_lumpy/summ/${RUN}__${NAME}.json
rm -f "$OUT"/cov.bin "$OUT"/pooled.bin
echo "QDEDUP DONE $RUN/$NAME"
