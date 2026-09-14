#!/bin/bash
# ---------------------------------------------------------------------------
# tests/lib/env.sh -- one place where every tier resolves its environment.
#
# Source it, do not run it.  It sets:
#
#   GMOQ_REPO       repo root (derived from this file's location)
#   GMOQ_WORK       scratch root for installs, results and simulated data.
#                   NEVER inside the repo: nothing the suite produces is
#                   committed, and no simulated phenotype file should ever end
#                   up staged by accident.
#   GMOQ_TEST_LIB   private R library the package under test is installed into
#   GMOQ_DEP_LIB    library holding the CRAN dependencies
#   GMOQ_TEST_OUT   per-tier results directory (set by the caller)
#
# Everything is overridable from the environment so a developer on a laptop and
# a Slurm job can use the same scripts.
# ---------------------------------------------------------------------------

_env_self="${BASH_SOURCE[0]}"
GMOQ_REPO="$(cd "$(dirname "${_env_self}")/../.." && pwd)"
export GMOQ_REPO

export GMOQ_WORK="${GMOQ_WORK:-${TMPDIR:-/tmp}/gmoq_testsuite}"
export GMOQ_TEST_LIB="${GMOQ_TEST_LIB:-$GMOQ_WORK/lib}"
export GMOQ_DEP_LIB="${GMOQ_DEP_LIB:-}"
mkdir -p "$GMOQ_WORK" "$GMOQ_TEST_LIB"

# Stage 1 lives in a conda environment (stage1/environment.yml). Tiers 1-3 need
# it; tier 0 does not. GMOQ_CONDA_SH and GMOQ_CONDA_ENV let a site point at its
# own installation without editing any script.
export GMOQ_CONDA_SH="${GMOQ_CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
export GMOQ_CONDA_ENV="${GMOQ_CONDA_ENV:-fungwas-stage1}"

# `Rscript --vanilla` is mandatory: a bare .libPaths() in ~/.Rprofile REPLACES
# the search path, so R_LIBS_USER is silently ignored and a stale package is
# loaded. tests/lib/harness.R sets the path from GMOQ_TEST_LIB and then asserts
# the package really came from there.
export GMOQ_RSCRIPT="${GMOQ_RSCRIPT:-Rscript}"

gmoq_rscript() {
  "$GMOQ_RSCRIPT" --vanilla "$@"
}

gmoq_activate_stage1() {
  if [ -f "$GMOQ_CONDA_SH" ]; then
    # shellcheck disable=SC1090
    source "$GMOQ_CONDA_SH"
    conda activate "$GMOQ_CONDA_ENV"
  fi
  export PYTHONPATH="$GMOQ_REPO/stage1:${PYTHONPATH:-}"
}

# Install the package under test into the private library, and prove the
# compiled kernel loads from there afterwards.
gmoq_install_stage2() {
  echo "== installing fungwasStage2 into $GMOQ_TEST_LIB"
  "$GMOQ_RSCRIPT" --vanilla -e "
    lib <- Sys.getenv('GMOQ_TEST_LIB'); dep <- Sys.getenv('GMOQ_DEP_LIB')
    .libPaths(unique(c(lib, dep[nzchar(dep)], .libPaths())))
    install.packages(file.path(Sys.getenv('GMOQ_REPO'), 'stage2'),
                     repos = NULL, type = 'source', lib = lib)
    .libPaths(unique(c(lib, dep[nzchar(dep)], .libPaths())))
    suppressPackageStartupMessages(library(fungwasStage2))
    p <- find.package('fungwasStage2')
    if (!startsWith(normalizePath(p), normalizePath(lib)))
      stop('loaded from ', p, ' not ', lib)
    if (!('fungwasStage2' %in% names(getLoadedDLLs())))
      stop('NAMESPACE lost useDynLib: the compiled kernel did not load')
    cat('install ok:', p, '\n')
  " || return 1
}
