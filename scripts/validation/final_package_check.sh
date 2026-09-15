#!/usr/bin/env bash
# Build the final candidate on BluePebble, preserving every output in a new run.
set -euo pipefail
src="${1:?usage: final_package_check.sh <repo> <new-output-dir> <source-commit> [test.R ...]}"
out="${2:?new output directory required}"
source_commit="${3:?source commit required}"
shift 3
case "$(hostname)" in bp1-compute*) ;; *) echo 'Requires a BluePebble compute node' >&2; exit 2;; esac
mkdir "$out"
module load languages/R/4.5.1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 MAKEFLAGS=-j2
export R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null
export R_LIBS_USER="${GMOQ_DEP_LIB:-/user/work/fh6520/r_packages}"
printf '%s\n' "$source_commit" > "$out/source_commit.txt"
hostname > "$out/host.txt"
printf '%s\n' "${SLURM_JOB_ID:-unknown}" > "$out/slurm_job.txt"
Rscript --vanilla -e 'sessionInfo()' > "$out/environment.txt"
cd "$out"
R CMD build "$src/stage2" > build.log 2>&1
archives=(fungwasStage2_*.tar.gz)
[[ ${#archives[@]} -eq 1 && -f "${archives[0]}" ]]
R CMD check --no-manual --no-tests "${archives[0]}" > check.log 2>&1
mkdir lib
R CMD INSTALL --library="$out/lib" "${archives[0]}" > install.log 2>&1
export R_LIBS_USER="$out/lib:$R_LIBS_USER" GMOQ_TEST_LIB="$out/lib"
for test_file in "$@"; do
  Rscript --vanilla "$src/stage2/tests/$test_file" > "$out/$test_file.log" 2>&1
done
tail -10 fungwasStage2.Rcheck/00check.log
