#!/bin/bash
# ---------------------------------------------------------------------------
# Tier 1 -- integration. About a minute. No Slurm.
#
#   bash tests/run_tier1.sh
#
# Needs the Stage 1 conda environment (stage1/environment.yml) because Stage 1
# is Python and C++; point GMOQ_CONDA_SH / GMOQ_CONDA_ENV at your installation.
# It does NOT need a cluster, a BGEN, or any real data.
#
#   1. builds the fixtures by running the SHIPPED Stage 1 on simulated data
#   2. Stage 1 -> Stage 2 end to end against a known theta
#   3. every refusal path
#   4. the eight model families and the tournament
# ---------------------------------------------------------------------------
set -uo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$here/lib/env.sh"

export GMOQ_TEST_OUT="${GMOQ_TEST_OUT:-$GMOQ_WORK/results/tier1}"
export GMOQ_TIER1_DIR="${GMOQ_TIER1_DIR:-$GMOQ_WORK/tier1}"
mkdir -p "$GMOQ_TEST_OUT" "$GMOQ_TIER1_DIR"
rm -f "$GMOQ_TEST_OUT/results.jsonl"

fail=0
t_start=$(date +%s)

if [ "${GMOQ_SKIP_INSTALL:-0}" != "1" ]; then
  gmoq_install_stage2 || { echo "TIER1: package failed to install"; exit 1; }
fi

echo
echo "########## tier 1: building fixtures ##########"
gmoq_activate_stage1 || exit 1
t_fix0=$(date +%s)
if [ "${GMOQ_SKIP_FIXTURES:-0}" != "1" ]; then
  "$GMOQ_PYTHON" "$here/tier1/t1_make_fixtures.py" --out "$GMOQ_TIER1_DIR" \
      --threads "${GMOQ_THREADS:-1}" || { echo "fixture build failed"; exit 1; }
  "$GMOQ_PYTHON" "$here/tier1/t1_families.py" --out "$GMOQ_TIER1_DIR" \
      --threads "${GMOQ_THREADS:-1}" || { echo "family fixture build failed"; exit 1; }
fi
t_fix=$(( $(date +%s) - t_fix0 ))
echo "fixtures built in ${t_fix}s"
printf '{"tier":"tier1","group":"meta","id":"T1.META.fixture_seconds","status":"INFO","severity":"report","establishes":"Wall clock to build every tier 1 fixture, including the shipped Stage 1 runs.","stat":%d,"criterion":"reported"}\n' \
  "$t_fix" >>"$GMOQ_TEST_OUT/results.jsonl"

echo
echo "########## tier 1: checks ##########"
for f in "$here"/tier1/t1_*.R; do
  echo
  echo ">>> $(basename "$f")"
  if ! gmoq_rscript "$f"; then
    echo "!!! FAILED: $(basename "$f")"
    fail=1
  fi
done

t_end=$(date +%s)
echo
echo "########## tier 1 finished in $((t_end - t_start))s ##########"
printf '{"tier":"tier1","group":"meta","id":"T1.META.runtime_seconds","status":"%s","severity":"report","establishes":"Wall clock for the whole of tier 1.","stat":%d,"criterion":"reported"}\n' \
  "$([ $fail -eq 0 ] && echo INFO || echo INFO)" "$((t_end - t_start))" >>"$GMOQ_TEST_OUT/results.jsonl"

exit $fail
