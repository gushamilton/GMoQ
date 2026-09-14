#!/bin/bash
# ---------------------------------------------------------------------------
# Tier 2 -- statistical validity. Cluster. Counted coverage and Type-I error
# against known truth, on a preregistered cell grid.
#
#   bash tests/run_tier2.sh submit     # build the grid and submit the array
#   bash tests/run_tier2.sh status     # how far along
#   bash tests/run_tier2.sh collate    # counts -> rates -> verdicts
#   bash tests/run_tier2.sh quick      # 2 cells x 2 chunks, locally, ~3 min
#
# The grid lives in tests/tier2/cells.py and the pass criteria in
# tests/criteria.json. Both are committed; neither is derived from a run.
#
# Cost at the shipped grid: 33 cells, 306 array tasks. Measured on the Bristol
# BluePebble `short` partition, see docs/TEST_SUITE.md.
# ---------------------------------------------------------------------------
set -uo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$here/lib/env.sh"

export GMOQ_TIER2_DIR="${GMOQ_TIER2_DIR:-$GMOQ_WORK/tier2}"
export GMOQ_TEST_OUT="${GMOQ_TEST_OUT:-$GMOQ_WORK/results/tier2}"
mkdir -p "$GMOQ_TIER2_DIR" "$GMOQ_TEST_OUT"
CELLS="$GMOQ_TIER2_DIR/cells.json"

ACCOUNT="${GMOQ_SLURM_ACCOUNT:-sscm013902}"
PARTITION="${GMOQ_SLURM_PARTITION:-short}"
CPUS="${GMOQ_SLURM_CPUS:-4}"
MEM="${GMOQ_SLURM_MEM:-24G}"
TIME="${GMOQ_SLURM_TIME:-3:00:00}"
THROTTLE="${GMOQ_SLURM_THROTTLE:-40}"

build_cells() {
  gmoq_activate_stage1
  python "$here/tier2/cells.py" "$CELLS"
}

one_task() {
  # Run one array task end to end: Stage 1 for a chunk, Stage 2 over it, then
  # drop the bulky intermediate. Keeping the chunk binaries would be ~5 GB.
  local idx="$1"
  gmoq_activate_stage1
  python "$here/tier2/t2_stage1.py" --cells "$CELLS" --task "$idx" \
      --out "$GMOQ_TIER2_DIR" --threads "${GMOQ_THREADS:-$CPUS}" || return 1
  local tag
  tag=$(python - "$CELLS" "$idx" <<'PY'
import json, sys
s = json.load(open(sys.argv[1])); t = s["tasks"][int(sys.argv[2])]
print(f"{s['cells'][t['cell']]['name']}__chunk{t['chunk']:03d}")
PY
)
  gmoq_rscript "$here/tier2/t2_stage2.R" "$tag" "$GMOQ_TIER2_DIR" || return 1
  rm -f "$GMOQ_TIER2_DIR/$tag".{beta,se,cov,dresid,sigma0}.bin
}

case "${1:-submit}" in

  submit)
    build_cells
    n=$(python -c "import json,sys;print(json.load(open('$CELLS'))['n_tasks'])")
    echo "submitting $n tasks (throttle $THROTTLE)"
    cat >"$GMOQ_TIER2_DIR/array.sh" <<EOF
#!/bin/bash
#SBATCH --account=$ACCOUNT
#SBATCH --partition=$PARTITION
#SBATCH --time=$TIME
#SBATCH --cpus-per-task=$CPUS
#SBATCH --mem=$MEM
#SBATCH --job-name=gmoq_t2
#SBATCH --output=$GMOQ_TIER2_DIR/logs/t2_%A_%a.out
#SBATCH --array=0-$((n - 1))%$THROTTLE
export GMOQ_WORK="$GMOQ_WORK"
export GMOQ_TEST_LIB="$GMOQ_TEST_LIB"
export GMOQ_DEP_LIB="$GMOQ_DEP_LIB"
export GMOQ_CONDA_SH="$GMOQ_CONDA_SH"
export GMOQ_CONDA_ENV="$GMOQ_CONDA_ENV"
export GMOQ_TIER2_DIR="$GMOQ_TIER2_DIR"
export GMOQ_THREADS=$CPUS
${GMOQ_MODULE_CMD:-module load languages/R/4.5.1}
bash "$here/run_tier2.sh" task "\$SLURM_ARRAY_TASK_ID"
EOF
    mkdir -p "$GMOQ_TIER2_DIR/logs"
    sbatch "$GMOQ_TIER2_DIR/array.sh"
    ;;

  task)
    one_task "$2"
    ;;

  status)
    n=$(python -c "import json;print(json.load(open('$CELLS'))['n_tasks'])" 2>/dev/null || echo "?")
    done_n=$(ls "$GMOQ_TIER2_DIR"/*.s2.rds 2>/dev/null | wc -l)
    echo "tier 2: $done_n / $n chunks complete"
    ;;

  collate)
    gmoq_rscript "$here/tier2/t2_collate.R" "$GMOQ_TIER2_DIR"
    ;;

  quick)
    # A smoke test of the whole tier on two cells, runnable without Slurm.
    build_cells
    export GMOQ_TIER2_DIR="$GMOQ_WORK/tier2_quick"
    export GMOQ_TEST_OUT="$GMOQ_WORK/results/tier2_quick"
    mkdir -p "$GMOQ_TIER2_DIR" "$GMOQ_TEST_OUT"
    rm -f "$GMOQ_TEST_OUT/results.jsonl"
    python - "$CELLS" "$GMOQ_TIER2_DIR/cells.json" <<'PY'
import json, sys
s = json.load(open(sys.argv[1]))
keep = [i for i, c in enumerate(s["cells"])
        if c["name"] in ("null_maf0.3_normal", "locscale_maf0.3_normal",
                         "pin_vqtl_maf0.01")]
cells = [dict(s["cells"][i]) for i in keep]
for c in cells:
    c["reps"] = 4; c["chunk_reps"] = 2; c["N"] = 10000; c["n_taus"] = 15
    c["blocks"] = 200
tasks = [{"cell": i, "chunk": k} for i, c in enumerate(cells)
         for k in range((c["reps"] + c["chunk_reps"] - 1) // c["chunk_reps"])]
json.dump(dict(cells=cells, tasks=tasks, n_cells=len(cells), n_tasks=len(tasks)),
          open(sys.argv[2], "w"), indent=2)
print(f"quick grid: {len(cells)} cells, {len(tasks)} tasks")
PY
    CELLS="$GMOQ_TIER2_DIR/cells.json"
    n=$(python -c "import json;print(json.load(open('$CELLS'))['n_tasks'])")
    rc=0
    for ((i = 0; i < n; i++)); do one_task "$i" || rc=1; done
    gmoq_rscript "$here/tier2/t2_collate.R" "$GMOQ_TIER2_DIR" || rc=1
    exit $rc
    ;;

  *)
    echo "usage: run_tier2.sh {submit|status|collate|quick}" >&2
    exit 2
    ;;
esac
