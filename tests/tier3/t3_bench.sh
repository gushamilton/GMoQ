#!/bin/bash
# ---------------------------------------------------------------------------
# Tier 3: run Stage 1 and Stage 2 at realistic scale and measure them.
#
#   bash tests/tier3/t3_bench.sh <dir> <n_samples> <n_variants> <threads>
#
# Wall clock, CPU time and peak RSS come from /usr/bin/time -v, separately for
# each stage. The PHASE breakdown comes from Stage 1's own log lines, which
# already instrument RIF build, Sigma_0, batch compute and batch write; the
# BGEN read phase is what is left of the batch loop after compute and write,
# which is the honest way to get it without editing the shipped code.
#
# The hardware is recorded alongside, because a timing with no machine attached
# to it is not reproducible and should not be quoted.
# ---------------------------------------------------------------------------
set -uo pipefail
dir="${1:?work directory}"
N="${2:-461320}"
M="${3:-5000}"
TH="${4:-${GMOQ_THREADS:-8}}"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "$here/lib/env.sh"

mkdir -p "$dir"
out="$dir/out"
mkdir -p "$out"

# --- hardware provenance ----------------------------------------------------
{
  echo "host=$(hostname)"
  echo "date=$(date -Is)"
  echo "cpu_model=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ *//')"
  echo "sockets=$(lscpu 2>/dev/null | awk -F: '/Socket\(s\)/{gsub(/ /,"",$2);print $2}')"
  echo "cores_per_socket=$(lscpu 2>/dev/null | awk -F: '/Core\(s\) per socket/{gsub(/ /,"",$2);print $2}')"
  echo "threads_used=$TH"
  echo "mem_total_kB=$(awk '/MemTotal/{print $2}' /proc/meminfo)"
  echo "slurm_job=${SLURM_JOB_ID:-none}"
  echo "slurm_node=${SLURM_NODELIST:-none}"
  echo "n_samples=$N"
  echo "n_variants=$M"
} >"$dir/hardware.txt"
cat "$dir/hardware.txt"

gmoq_activate_stage1

# --- Stage 1 ----------------------------------------------------------------
echo
echo "== Stage 1: $M variants x 45 taus x $N samples, $TH threads"
/usr/bin/time -v -o "$dir/stage1.time" \
  python -m fungwas_stage1.cli \
    --bgen "$dir/geno.bgen" \
    --pheno "$dir/pheno.txt" --pheno-col simpheno \
    --covar "$dir/covar.txt" \
    --out "$out/bench" \
    --blocks 1000 --threads "$TH" --batch-size 500 \
    --allow-multiallelic-rsids --no-seek-variants \
    --cov-dtype int8 \
  >"$dir/stage1.log" 2>&1
s1_rc=$?
tail -n 25 "$dir/stage1.log"
[ $s1_rc -eq 0 ] || { echo "Stage 1 failed, see $dir/stage1.log"; exit 1; }

# --- Stage 2 ----------------------------------------------------------------
echo
echo "== Stage 2"
/usr/bin/time -v -o "$dir/stage2.time" \
  "$GMOQ_RSCRIPT" --vanilla "$here/tier3/t3_stage2.R" "$out/bench" "$dir/stage2.json" \
  >"$dir/stage2.log" 2>&1
s2_rc=$?
tail -n 20 "$dir/stage2.log"
[ $s2_rc -eq 0 ] || { echo "Stage 2 failed, see $dir/stage2.log"; exit 1; }

# --- assemble ---------------------------------------------------------------
python "$here/tier3/t3_parse.py" --dir "$dir" --out "$dir/tier3.json"
echo
echo "tier 3 measurements in $dir/tier3.json"
