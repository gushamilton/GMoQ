#!/usr/bin/env bash
#SBATCH --job-name=pc-gmoq21
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH --output=/user/work/fh6520/gmoq_power_conserv/logs/gmoq21_%j.log
set -euo pipefail

S=/user/work/fh6520/gmoq_power_conserv
DATA=/user/work/fh6520/fungwas_paper
CHR=21

source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh
conda activate fungwas-stage1
T=${SLURM_CPUS_PER_TASK:-8}
export PYTHONPATH=/user/work/fh6520/gmoq_locovar/pkg:/user/work/fh6520/gmoq_finish_20260915/inputs/pydeps
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export OMP_NUM_THREADS=$T OPENBLAS_NUM_THREADS=$T MKL_NUM_THREADS=$T
TAUS=$(python -c 'import numpy as np; print(",".join("%.17g" % t for t in np.linspace(0.02,0.98,45)))')
echo "TAUS=$TAUS"

BGEN=/bp1/mrcieu1/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr$CHR.bgen
SAMPLE=/bp1/mrcieu1/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample
PHENO=$DATA/data/derived/paper/phenotypes_paper.txt
COVAR=$DATA/data/derived/paper/covars_no_sex_paper.txt
KEEP=$DATA/data/derived/paper/keep_files/keep_testosterone_series.txt
PRED=$DATA/results/paper/regenie_step1/full/testosterone/testosterone_step1_pred.list
SNPS=$DATA/inputs/data/w_hm3_chr${CHR}_only.txt

for MODE in subtract covariate; do
  OUT=$S/stage1/$MODE; rm -rf "$OUT"; mkdir -p "$OUT"
  echo "=== stage1 chr$CHR loco-mode=$MODE ==="
  python -m fungwas_stage1.cli \
    --bgen "$BGEN" --sample "$SAMPLE" \
    --pheno "$PHENO" --pheno-col testosterone \
    --covar "$COVAR" --keep "$KEEP" --snps "$SNPS" \
    --loco-preds "$PRED" --loco-scale regression --loco-mode "$MODE" \
    --taus "$TAUS" --blocks 1000 --batch-size 100 --cov-dtype float32 \
    --seed 42 --threads $T \
    --out "$OUT/testosterone_chr$CHR" >"$S/logs/stage1_${MODE}_chr$CHR.log" 2>&1
  echo "  done $MODE"
done
conda deactivate || true

source /etc/profile.d/modules.sh
module load languages/R/4.5.1
export R_LIBS=/user/work/fh6520/gmoq_verify_work/deplib
RLIB=/user/work/fh6520/gmoq_k4_scratch/Rlib
BUILDER=/user/work/fh6520/gmoq_noloco/prepare_matched_weights_noloco.R
WRAP=/user/work/fh6520/Gmoq_paper/scripts/stage2/run_stage2_chr.R

do_stage2 () {
  MODE=$1; shift
  stage1=$S/stage1/$MODE/testosterone_chr$CHR.stage1.tsv.gz
  W=$S/weights/$MODE; rm -rf "$W"; mkdir -p "$W" "$S/stage2/$MODE"
  Rscript "$BUILDER" --trait testosterone --chr $CHR --stage1 "$stage1" \
    --pheno "$DATA/data/derived/paper/phenotypes_paper.txt" \
    --covar "$DATA/data/derived/paper/covars_no_sex_paper.txt" \
    --keep "$DATA/data/derived/paper/keep_files/keep_testosterone_series.txt" \
    --bgen-sample "$SAMPLE" --family skew --out-dir "$W" --r-lib "$RLIB" --seed 20260915 "$@" 2>&1 | tail -4
  cat "$W/weights_manifest.tsv"
  Rscript "$WRAP" --trait testosterone --model skew --chr $CHR --stage1 "$stage1" \
    --weights "$W/weights_skew.rds" --out "$S/stage2/$MODE/chr$CHR" \
    --cov-source pooled --estimator ols \
    --sigma0-file "${stage1%.stage1.tsv.gz}.sigma0.tsv.gz" \
    --threads 8 --r-lib "$RLIB" 2>&1 | tail -3
}
do_stage2 covariate
do_stage2 subtract --loco "$PRED" --loco-scale regression
echo ALLDONE
