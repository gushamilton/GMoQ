#!/bin/bash
# ---------------------------------------------------------------------------
# Tier 3 -- realistic scale and timing. Cluster.
#
#   bash tests/run_tier3.sh [n_samples] [n_variants]
#
# Generates a SYNTHETIC BGEN with plink2 --dummy at the requested sample size
# (no real genotypes are involved and none may be), runs the Stage 1 CLI and
# Stage 2 over it, and records wall clock, CPU time, peak RSS and a phase
# breakdown for each stage, together with the hardware they were measured on.
#
# Defaults to production N = 461,320 and 5,000 variants. That is a slice, not a
# genome-wide scan: the per-variant costs are extrapolated to 1.2M variants in
# the output and labelled as an extrapolation.
# ---------------------------------------------------------------------------
set -uo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$here/lib/env.sh"

N="${1:-461320}"
M="${2:-5000}"
DIR="${GMOQ_TIER3_DIR:-$GMOQ_WORK/tier3}"
export GMOQ_TEST_OUT="${GMOQ_TEST_OUT:-$GMOQ_WORK/results/tier3}"
mkdir -p "$DIR" "$GMOQ_TEST_OUT"
rm -f "$GMOQ_TEST_OUT/results.jsonl"

bash "$here/tier3/t3_make_inputs.sh" "$DIR" "$N" "$M" || exit 1
bash "$here/tier3/t3_bench.sh" "$DIR" "$N" "$M" "${GMOQ_THREADS:-8}" || exit 1
gmoq_rscript "$here/tier3/t3_report.R" "$DIR/tier3.json"
