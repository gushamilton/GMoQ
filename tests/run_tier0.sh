#!/bin/bash
# ---------------------------------------------------------------------------
# Tier 0 -- unit tests. Seconds. No cluster, no Slurm, no Stage 1, no network.
#
#   bash tests/run_tier0.sh
#
# What it does:
#   1. installs stage2/ into a private R library (skip with GMOQ_SKIP_INSTALL=1)
#   2. runs every tests/tier0/*.R
#   3. runs the pre-existing stage2/tests unit scripts that are standalone
#   4. leaves a machine-readable results.jsonl for tests/report.R
#
# Exit status is 0 only if every check passed.
# ---------------------------------------------------------------------------
set -uo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$here/lib/env.sh"

export GMOQ_TEST_OUT="${GMOQ_TEST_OUT:-$GMOQ_WORK/results/tier0}"
mkdir -p "$GMOQ_TEST_OUT"
rm -f "$GMOQ_TEST_OUT/results.jsonl"

fail=0
t_start=$(date +%s)

if [ "${GMOQ_SKIP_INSTALL:-0}" != "1" ]; then
  gmoq_install_stage2 || { echo "TIER0: package failed to install"; exit 1; }
fi

echo
echo "########## tier 0 ##########"
for f in "$here"/tier0/t0_*.R; do
  echo
  echo ">>> $(basename "$f")"
  if ! gmoq_rscript "$f"; then
    echo "!!! FAILED: $(basename "$f")"
    fail=1
  fi
done

# The pre-existing Stage 2 unit scripts. They are written to run standalone
# (Rscript --vanilla <path>, pkg_root derived from --file=), which is the same
# convention the tier scripts use, so they are run unchanged.
echo
echo "########## pre-existing stage2 unit scripts ##########"
LEGACY="test-namespace-dynlib.R test-dispersion-guards.R test-q-negative-eigenvalues.R
        test-q-calibration.R test-cov-ids-mandatory.R test-pooled-covariance.R
        test-multi-taucov.R test-gls-param-gwas.R"
# These scripts call library(fungwasStage2) with no path handling of their own,
# so the library is supplied through R_LIBS, which --vanilla still honours
# (--no-environ skips .Renviron FILES, not the environment variable).
export R_LIBS="$GMOQ_TEST_LIB${GMOQ_DEP_LIB:+:$GMOQ_DEP_LIB}"
for name in $LEGACY; do
  f="$GMOQ_REPO/stage2/tests/$name"
  [ -f "$f" ] || continue
  echo
  echo ">>> stage2/tests/$name"
  if gmoq_rscript "$f" >"$GMOQ_TEST_OUT/legacy_${name%.R}.log" 2>&1; then
    echo "  ok"
    st=PASS
  else
    echo "  FAIL (see $GMOQ_TEST_OUT/legacy_${name%.R}.log)"
    tail -n 15 "$GMOQ_TEST_OUT/legacy_${name%.R}.log"
    st=FAIL
    fail=1
  fi
  printf '{"tier":"tier0","group":"legacy","id":"T0.LEGACY.%s","status":"%s","establishes":"Pre-existing Stage 2 unit script %s still passes.","criterion":"exit status 0"}\n' \
    "${name%.R}" "$st" "$name" >>"$GMOQ_TEST_OUT/results.jsonl"
done

t_end=$(date +%s)
echo
echo "########## tier 0 finished in $((t_end - t_start))s ##########"
printf '{"tier":"tier0","group":"meta","id":"T0.META.runtime_seconds","status":"%s","establishes":"Wall clock for the whole of tier 0.","stat":%d,"criterion":"reported, not asserted"}\n' \
  "$([ $fail -eq 0 ] && echo PASS || echo FAIL)" "$((t_end - t_start))" >>"$GMOQ_TEST_OUT/results.jsonl"

exit $fail
