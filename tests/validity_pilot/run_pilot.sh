#!/usr/bin/env bash
# Run the bounded known-truth pilot on the BluePebble compute node.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
out="${1:?usage: run_pilot.sh <output-dir> [reps] [N] [blocks] [threads]}"
reps="${2:-10}"
N="${3:-5000}"
blocks="${4:-1000}"
threads="${5:-8}"

mkdir -p "$out"
python="${GMOQ_PYTHON:-}"
if [[ -z "$python" ]]; then
  echo 'GMOQ_PYTHON must point to the Stage 1 environment interpreter' >&2
  exit 2
fi

"$python" "$here/pilot_stage1.py" --repo "$repo" --out "$out" \
  --reps "$reps" --N "$N" --blocks "$blocks" --threads "$threads"

Rscript --vanilla "$here/pilot_stage2.R" "$out" "$out/pilot_rows.csv"
Rscript --vanilla "$here/summarise_pilot.R" "$out"
Rscript --vanilla "$here/report_pilot.R" "$out"
