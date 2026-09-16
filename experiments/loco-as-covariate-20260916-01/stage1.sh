#!/usr/bin/env bash
#SBATCH --job-name=locovar-s1
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/user/work/fh6520/gmoq_locovar/logs/s1_%j.log
set -euo pipefail

S=/user/work/fh6520/gmoq_locovar
DATA=/user/work/fh6520/fungwas_paper
REF=/user/work/fh6520/Gmoq_paper/.local/runs/paper-core-genome-20260915-02/testosterone/chr22

source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh
conda activate fungwas-stage1
T=${SLURM_CPUS_PER_TASK:-8}
export PYTHONPATH=$S/pkg:/user/work/fh6520/gmoq_finish_20260915/inputs/pydeps
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export OMP_NUM_THREADS=$T OPENBLAS_NUM_THREADS=$T MKL_NUM_THREADS=$T
TAUS=$(python -c "import numpy as np; print(','.join(f'{t:.17g}' for t in np.linspace(0.02,0.98,45)))")

BGEN=/bp1/mrcieu1/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr22.bgen
SAMPLE=/bp1/mrcieu1/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample
PHENO=$DATA/data/derived/paper/phenotypes_paper.txt
COVAR=$DATA/data/derived/paper/covars_no_sex_paper.txt
KEEP=$DATA/data/derived/paper/keep_files/keep_testosterone_series.txt
PRED=$DATA/results/paper/regenie_step1/full/testosterone/testosterone_step1_pred.list
SNPS=$DATA/inputs/data/w_hm3_chr22_only.txt

run_one () {
  MODE=$1
  OUT=$S/stage1/$MODE; rm -rf "$OUT"; mkdir -p "$OUT"
  echo "=== stage1 loco-mode=$MODE ==="
  start=$(date +%s)
  python -m fungwas_stage1.cli \
    --bgen "$BGEN" --sample "$SAMPLE" \
    --pheno "$PHENO" --pheno-col testosterone \
    --covar "$COVAR" --keep "$KEEP" --snps "$SNPS" \
    --loco-preds "$PRED" --loco-scale regression --loco-mode "$MODE" \
    --taus "$TAUS" --blocks 1000 --batch-size 100 --cov-dtype float32 \
    --seed 42 --threads $T \
    --out "$OUT/testosterone_chr22" >"$S/logs/stage1_$MODE.log" 2>&1
  echo "  wall: $(( $(date +%s) - start ))s"
  grep -E "LOCO|Matched|covariate" "$S/logs/stage1_$MODE.log" | head -8
}

run_one subtract
run_one covariate

echo
echo "=== subtract-mode reproduction check vs paper reference ==="
a=$REF/testosterone_chr22.stage1.tsv.gz
b=$S/stage1/subtract/testosterone_chr22.stage1.tsv.gz
if cmp -s <(gzip -dc "$a") <(gzip -dc "$b"); then echo "  IDENTICAL stage1.tsv.gz"; else echo "  DIFFERS stage1.tsv.gz"; fi
for f in sigma0.tsv.gz cov.gz; do
  if cmp -s <(gzip -dc "$REF/testosterone_chr22.$f") <(gzip -dc "$S/stage1/subtract/testosterone_chr22.$f"); then
    echo "  IDENTICAL $f"; else echo "  DIFFERS   $f"; fi
done
echo
echo "=== covariate vs subtract stage1 differ (they must) ==="
if cmp -s <(gzip -dc "$S/stage1/subtract/testosterone_chr22.stage1.tsv.gz") \
          <(gzip -dc "$S/stage1/covariate/testosterone_chr22.stage1.tsv.gz"); then
  echo "  IDENTICAL -- SOMETHING IS WRONG"; else echo "  DIFFER (expected)"; fi
python - <<'PY'
import json
for m in ("subtract","covariate"):
    d=json.load(open(f"/user/work/fh6520/gmoq_locovar/stage1/{m}/testosterone_chr22.sigma0.meta.json"))
    print(m, "n=",d["n_samples"], "rank=",d.get("covariate_basis_rank"),
          "loco_mode=",d.get("loco_mode"), "scale=",d.get("loco_scale_factor"),
          "chrs=",d.get("chromosomes"))
PY
echo DONE
