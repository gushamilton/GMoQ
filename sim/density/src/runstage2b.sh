#!/bin/bash
#SBATCH --account=sscm013902
#SBATCH --partition=compute
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
set -e
source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh
conda activate fungwas-stage1
PY=$CONDA_PREFIX/bin/python
export PYTHONPATH=/user/work/fh6520/kde_lumpy/tree/stage1
export OMP_NUM_THREADS=4
cd /user/work/fh6520/kde_lumpy/src
D=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBLIST")
NAME=$(basename "$D")
RUN=$(basename $(dirname "$D"))
OUT=/user/work/fh6520/kde_lumpy/col/${RUN}__${NAME}
echo "=== $D -> $OUT"
mkdir -p /user/work/fh6520/kde_lumpy/summ
$PY collate.py --dir "$D" --out "$OUT"
$PY keepset.py --dir "$OUT"
module load languages/R/4.5.1
Rscript --vanilla stage2.R "$OUT" full
if [ -s "$OUT/keep.json" ]; then Rscript --vanilla stage2.R "$OUT" tiedrop "$OUT/keep.json" || true; fi
$PY analyse.py --dir "$OUT" --out /user/work/fh6520/kde_lumpy/summ/${RUN}__${NAME}.json
rm -f "$OUT"/cov.bin "$OUT"/pooled.bin
echo "S2 DONE $RUN/$NAME"
