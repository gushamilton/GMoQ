#!/usr/bin/env bash
#SBATCH --job-name=locovar-s2
#SBATCH --account=sscm013902
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/user/work/fh6520/gmoq_locovar/logs/s2_%j.log
set -euo pipefail
source /etc/profile.d/modules.sh
module load languages/R/4.5.1
export R_LIBS=/user/work/fh6520/gmoq_verify_work/deplib

S=/user/work/fh6520/gmoq_locovar
DATA=/user/work/fh6520/fungwas_paper
RLIB=/user/work/fh6520/gmoq_k4_scratch/Rlib
BUILDER=/user/work/fh6520/gmoq_noloco/prepare_matched_weights_noloco.R
WRAP=/user/work/fh6520/Gmoq_paper/scripts/stage2/run_stage2_chr.R
SAMPLE=/bp1/mrcieu1/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample
PRED=$DATA/results/paper/regenie_step1/full/testosterone/testosterone_step1_pred.list

do_mode () {
  MODE=$1; shift
  stage1=$S/stage1/$MODE/testosterone_chr22.stage1.tsv.gz
  W=$S/weights/$MODE
  rm -rf "$W"; mkdir -p "$W" "$S/stage2/$MODE"
  echo "=== weights for $MODE ($*) ==="
  Rscript "$BUILDER" \
    --trait testosterone --chr 22 --stage1 "$stage1" \
    --pheno "$DATA/data/derived/paper/phenotypes_paper.txt" \
    --covar "$DATA/data/derived/paper/covars_no_sex_paper.txt" \
    --keep "$DATA/data/derived/paper/keep_files/keep_testosterone_series.txt" \
    --bgen-sample "$SAMPLE" \
    --family skew --out-dir "$W" --r-lib "$RLIB" --seed 20260915 "$@" 2>&1 | tail -4
  cat "$W/weights_manifest.tsv"
  echo "=== stage2 for $MODE ==="
  Rscript "$WRAP" \
    --trait testosterone --model skew --chr 22 --stage1 "$stage1" \
    --weights "$W/weights_skew.rds" --out "$S/stage2/$MODE/chr22" \
    --cov-source pooled --estimator ols \
    --sigma0-file "${stage1%.stage1.tsv.gz}.sigma0.tsv.gz" \
    --threads 8 --r-lib "$RLIB" 2>&1 | tail -3
}

# (d) LOCO as covariate: Y unadjusted, so the Stage 2 baseline is raw Y -> no --loco
do_mode covariate
# (a') same-code subtract control: baseline is Y - c*L(22)
do_mode subtract --loco "$PRED" --loco-scale regression
echo DONE
