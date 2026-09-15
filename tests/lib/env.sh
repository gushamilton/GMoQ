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

# Resolve the Stage 1 interpreter EXPLICITLY rather than trusting PATH.
#
# `module load languages/R/4.5.1` puts an R conda environment first on PATH, and
# its python has numpy but not scipy. If `conda activate` then fails quietly --
# which it does whenever conda.sh has not defined the shell function in this
# particular subshell -- every Stage 1 call runs under the wrong interpreter and
# dies at `import scipy`, several minutes into a Slurm array. Naming the binary
# removes the whole class of problem, and the import check below turns a silent
# wrong-interpreter into an immediate, explicit failure.
gmoq_activate_stage1() {
  if [ -z "${GMOQ_PYTHON:-}" ] && [ -f "$GMOQ_CONDA_SH" ]; then
    local croot
    croot="$(cd "$(dirname "$GMOQ_CONDA_SH")/../.." 2>/dev/null && pwd)"
    if [ -n "$croot" ] && [ -x "$croot/envs/$GMOQ_CONDA_ENV/bin/python" ]; then
      GMOQ_PYTHON="$croot/envs/$GMOQ_CONDA_ENV/bin/python"
    fi
  fi
  if [ -f "$GMOQ_CONDA_SH" ]; then
    # shellcheck disable=SC1090
    source "$GMOQ_CONDA_SH" 2>/dev/null || true
    conda activate "$GMOQ_CONDA_ENV" 2>/dev/null || true
  fi
  export GMOQ_PYTHON="${GMOQ_PYTHON:-python}"
  export PYTHONPATH="$GMOQ_REPO/stage1:${PYTHONPATH:-}"
  if ! "$GMOQ_PYTHON" -c "import numpy, scipy" 2>/dev/null; then
    echo "ERROR: $GMOQ_PYTHON cannot import numpy and scipy." >&2
    echo "  Point GMOQ_PYTHON at the Stage 1 environment's interpreter, or set" >&2
    echo "  GMOQ_CONDA_SH / GMOQ_CONDA_ENV so it can be found." >&2
    return 1
  fi
}

gmoq_python() { "$GMOQ_PYTHON" "$@"; }

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
