#!/bin/bash
# ---------------------------------------------------------------------------
# tests/tier4/t4_stage1.sh -- Run Stage 1 on real chr22 data with simulated phenotypes.
#
# Usage:
#   bash tests/tier4/t4_stage1.sh <case> <work_dir> <out_prefix>
# ---------------------------------------------------------------------------
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "usage: bash tests/tier4/t4_stage1.sh <case> <work_dir> <out_prefix>" >&2
  echo "  case: null|location|scale|mixed" >&2
  exit 1
fi

case_name="$1"
work_dir="$2"
out_prefix="$3"
out_dir=$(dirname "$out_prefix")

# Resolve paths to real data
chr22_bgen="/bp1/mrcieu1/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr22.bgen"
sample_file="/bp1/mrcieu1/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample"
snp_file="/user/work/fh6520/fungwas_paper/inputs/data/w_hm3_chr22_only.txt"

# Use paper's real covariate and keep files
covar_file="/user/work/fh6520/fungwas_paper/data/derived/paper/covars_with_sex_paper.txt"
keep_file="/user/work/fh6520/fungwas_paper/data/derived/paper/keep_files/keep_bmi.txt"

mkdir -p "$out_dir"

# Generate phenotype using tier4 simulator
pheno_file="$work_dir/pheno_${case_name}.txt"
if [[ ! -f "$pheno_file" ]]; then
  cat >&2 << EOF
ERROR: Expected phenotype file not found: $pheno_file
  Run t4_make_pheno.R first to generate simulated phenotypes.
EOF
  exit 1
fi

cat >&2 << EOF
[t4_stage1] Running Stage 1 on chr22
[t4_stage1] case: $case_name
[t4_stage1] bgen: $chr22_bgen
[t4_stage1] pheno: $pheno_file
[t4_stage1] out: $out_prefix
EOF

# Activate Stage 1 environment
source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh 2>/dev/null || true
conda activate fungwas-stage1 2>/dev/null || true
export PYTHONPATH="/user/work/fh6520/wt_testsuite/stage1:${PYTHONPATH:-}"

# Run Stage 1
python -m fungwas_stage1.cli \
  --bgen "$chr22_bgen" \
  --sample "$sample_file" \
  --pheno "$pheno_file" \
  --pheno-col "phenotype" \
  --covar "$covar_file" \
  --keep "$keep_file" \
  --snps "$snp_file" \
  --threads "${SLURM_CPUS_PER_TASK:-8}" \
  --batch-size 500 \
  --blocks 65 \
  --cov-dtype int8 \
  --out "$out_prefix"

cat >&2 << EOF
[t4_stage1] Stage 1 complete
[t4_stage1] outputs:
EOF
ls -lh "${out_prefix}"* >&2

echo "$out_prefix"
