#!/usr/bin/env bash
# Run the paper's empirical log_bmi input definition through the shipped Stage 1
# CLI on the same bounded 24-SNP chr22 list used by real_chr22_smoke.py.
#
# The paper scripts define log_bmi = log(bmi), use age/sex/PC1..PC10, keep IDs
# complete for those covariates and log_bmi, and pass the phenotype-specific
# REGENIE pred.list with --loco-scale regression (the CLI default, made explicit
# here).  This script records that wiring without recreating or modifying inputs.
set -euo pipefail

OUT_ROOT="${1:?usage: run_real_chr22_empirical.sh OUT_ROOT [SNPS]}"
BGEN="/bp1/mrcieu1/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr22.bgen"
SAMPLE="/bp1/mrcieu1/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample"
PAPER_ROOT="/user/work/fh6520/fungwas_paper"
PHENO="${PAPER_ROOT}/data/derived/paper/phenotypes_paper.txt"
COVAR="${PAPER_ROOT}/data/derived/paper/covars_with_sex_paper.txt"
KEEP="${PAPER_ROOT}/data/derived/paper/keep_files/keep_log_bmi.txt"
LOCO="${PAPER_ROOT}/results/paper/regenie_step1/full/log_bmi/log_bmi_step1_pred.list"
SNPS="${2:-$(dirname "${OUT_ROOT}")/selected_snps.txt}"
PREFIX="${OUT_ROOT}/stage1/log_bmi_chr22_smoke"

mkdir -p "${OUT_ROOT}/stage1"
export PYTHONPATH="/user/work/fh6520/gmoq_finish_20260915/inputs/pydeps:/user/work/fh6520/gmoq_finish_20260915/inputs/stage1"
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2

/usr/bin/time -v /user/home/fh6520/.conda/envs/fungwas-stage1/bin/python -m fungwas_stage1.cli \
  --bgen "${BGEN}" --sample "${SAMPLE}" \
  --pheno "${PHENO}" --pheno-col log_bmi \
  --covar "${COVAR}" --keep "${KEEP}" \
  --loco-preds "${LOCO}" --loco-scale regression \
  --snps "${SNPS}" --threads 2 --batch-size 500 --blocks 65 \
  --cov-dtype int8 --out "${PREFIX}"

printf 'stage1_prefix=%s\nphenotype=log_bmi=log(bmi)\ncovariates=age,sex,PC1..PC10\nkeep=%s\nloco=%s\nloco_scale=regression\n' \
  "${PREFIX}" "${KEEP}" "${LOCO}"
