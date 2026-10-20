#!/usr/bin/env bash
# Stage 1: RIF quantile GWAS, 45 quantiles, one call per phenotype.
set -euo pipefail
cd "$(dirname "$0")"; source config.sh
for ph in null shift variance mixture; do
  gmoq-stage1 --bgen "$WORK/geno.bgen" --pheno "$WORK/pheno.txt" --pheno-col "$ph" \
    --covar "$WORK/covar.txt" --out "$WORK/s1_$ph" --threads "$THREADS" --blocks 100 --seed 1 \
    > "$WORK/s1_$ph.log" 2>&1
  echo "stage1 $ph done"
done
