#!/usr/bin/env bash
# Whole example: genotypes -> phenotypes -> Stage 1 -> Stage 2 -> checks. About a minute.
set -euo pipefail
cd "$(dirname "$0")"; source config.sh
[ -s "$WORK/geno.bgen" ] && [ -s "$WORK/geno.raw" ] || bash 01_make_bgen.sh
Rscript --vanilla 02_simulate.R "$WORK" "${SEED:-2026}"
bash 03_stage1.sh
Rscript --vanilla 04_stage2.R "$WORK"
Rscript --vanilla 05_check.R "$WORK"
