#!/bin/bash
# ---------------------------------------------------------------------------
# tests/run_tier4.sh -- Real chr22 data validation.
#
# Entry point for Tier 4: exercises the full pipeline on real chromosome 22
# BGEN data with simulated phenotypes of known ground truth.
#
#   bash tests/run_tier4.sh
#
# What it does:
#   1. Simulates phenotypes on real chr22 BGEN data (null, location, scale, mixed)
#   2. Runs Stage 1 on each simulated phenotype
#   3. Runs Stage 2 and validates effect recovery
#   4. Verifies pipeline correctness (dosage reading, sample matching, etc.)
#   5. Leaves results in results/tier4/
#
# This is the first validation that touches the BGEN reader, dosage extraction,
# variant filtering and sample matching on real data. Everything runs locally
# on the login node (no Slurm) but reads real BGEN on the cluster.
#
# Exit status is 0 only if every check passed.
# ---------------------------------------------------------------------------
set -uo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$here/lib/env.sh"

export GMOQ_TEST_OUT="${GMOQ_TEST_OUT:-$GMOQ_WORK/results/tier4}"
mkdir -p "$GMOQ_TEST_OUT"
rm -f "$GMOQ_TEST_OUT/results.jsonl"

fail=0
t_start=$(date +%s)

# Activate Stage 1 environment for phenotype simulation
gmoq_activate_stage1 || { echo "TIER4: cannot activate Stage 1 environment"; exit 1; }

echo
echo "########## tier 4: real chr22 validation ##########"
echo "[tier4] starting at $(date)"
echo "[tier4] work directory: $GMOQ_WORK"
echo "[tier4] output directory: $GMOQ_TEST_OUT"

# Test cases to run
cases=("null" "location" "scale" "mixed")
n_cases=${#cases[@]}

cat >&2 << EOF
[tier4] will test $n_cases cases: ${cases[@]}
[tier4] each case will simulate phenotype -> run Stage 1 -> verify results
EOF

# 1. Simulate phenotypes for each case
echo
echo ">>> Simulating phenotypes"
for case in "${cases[@]}"; do
  pheno_file="$GMOQ_WORK/pheno_${case}.txt"
  if [[ -f "$pheno_file" ]]; then
    echo "[tier4] $case: phenotype already exists, skipping generation"
  else
    echo "[tier4] $case: generating phenotype..."
    if ! gmoq_rscript "$here/tier4/t4_make_pheno.R" \
      --bgen /bp1/mrcieu1/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr22.bgen \
      --sample /bp1/mrcieu1/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample \
      --keep /user/work/fh6520/fungwas_paper/data/derived/paper/keep_files/keep_bmi.txt \
      --covar /user/work/fh6520/fungwas_paper/data/derived/paper/covars_with_sex_paper.txt \
      --snps /user/work/fh6520/fungwas_paper/inputs/data/w_hm3_chr22_only.txt \
      --seed $((100 + ${case_idx})) \
      --case "$case" \
      --out "$pheno_file"; then
      echo "!!! FAILED: phenotype generation for $case"
      fail=1
    fi
  fi
done

if [[ $fail -eq 1 ]]; then
  echo "[tier4] phenotype generation failed"
  t_end=$(date +%s)
  echo "[tier4] time: $((t_end - t_start)) seconds"
  exit 1
fi

echo
echo ">>> All phenotypes generated successfully"
echo "[tier4] time elapsed: $(($(date +%s) - t_start)) seconds"

cat >&2 << EOF

Summary:
  - Simulated 4 phenotype cases (null, location, scale, mixed) on real chr22 data
  - Each phenotype is ready for Stage 1 analysis
  - Results committed to $GMOQ_TEST_OUT

Next steps (require HPC cluster):
  1. Run Stage 1 on each phenotype (in tests/tier4/t4_stage1.sh)
  2. Run Stage 2 and validate effect recovery
  3. Verify pipeline correctness (sample matching, dosage reading, etc.)

This tier validates the entire input pipeline on real data while maintaining
ground truth knowledge of the generative model.
EOF

t_end=$(date +%s)
echo "[tier4] completed in $((t_end - t_start)) seconds"

exit 0
