# fungwasStage2 0.6.2 (and fungwas-stage1 0.6.2)

Stage 1 only; no change to any estimate.

* Installs without conda. `setup.py` finds Armadillo via `ARMADILLO_DIR`, conda,
  Homebrew or the system prefixes, embeds its library path so the extension loads,
  and stops with install instructions if Armadillo is missing. If the extension
  cannot be loaded, the error now gives the real import error and the right
  reinstall command.
* Warnings instead of silence: Stage 1 reports how many phenotype rows it dropped
  and why (missing or non-numeric phenotype or covariate, no covariate row for
  the FID), warns when a phenotype or covariate column contains `-9` (analysed as
  a number), and warns when the extension was built without OpenMP so `--threads`
  has no effect. Which rows are analysed is unchanged.

# fungwasStage2 0.6.1 (and fungwas-stage1 0.6.1)

* MIT `LICENSE` at the repository root and in `stage1/`.
* The repository is now the two packages, their tests, CI, user docs and the
  1000 Genomes example. Simulation campaigns (`sim/`), run records
  (`experiments/`), the tier 1-4 suites, BluePebble/UK Biobank runners, one-off
  diagnostics, demonstration scripts, the built pkgdown site and dated audit
  reports moved to the tag `archive/research-records-20261016`. Links to them
  from the remaining docs point there. No BluePebble paths remain on `main`.
* Removed dead code: two C++ files that were never compiled
  (`stage1/src/stage1_core.cpp`, `stage1_cpp.cpp`), an unreachable branch in
  the Stage 1 CLI and the unreachable body of `quantile_gwas(covar_mode =
  "include")`, which still stops (now saying it is not implemented).

# fungwasStage2 0.6.0 (and fungwas-stage1 0.6.0)

Both packages now carry the same version. Stage 1 had been left at 0.1.0.

## Behaviour change

* **The mixture weight builders default to `density = "hybrid"` whenever the
  weight matrix has a membership (`gamma`) column**, and to `"model"`
  otherwise. This applies to `make_weights_normal_mixture(include_membership =
  TRUE)`, `make_weights_normal_mixture_vqtl()`, and `make_weights_generic()`
  when a gradient is named in `stage1_columns`. Hybrid needs the Stage 1 density,
  so these calls now require `stage1_density` (for example
  `stage1_kde_density(y, q_tau)`) and stop with a message saying so if it is
  missing. Pass `density = "model"` to reproduce 0.5.0 exactly. On 1000 Genomes
  with a simulated, well-separated mixture, `"model"` recovered about 70% of the
  true membership effect and `"hybrid"` about 90%, within its standard errors
  (`docs/examples/1kg_chr22/`).
* `param_gwas(transform = "two_normal")` passes `density` and `stage1_density`
  from `transform_args` through to the weight builder.

## Fixes

* Stage 1 installs on macOS, including Apple Silicon. `setup.py` no longer
  passes `-fopenmp` unconditionally (Apple clang rejects it). It probes for
  working OpenMP flags (`-fopenmp`, or `-Xpreprocessor -fopenmp` with a conda or
  Homebrew libomp), and if none work it builds single-threaded with a warning.
  `environment.yml` now includes `cxx-compiler` and `llvm-openmp`. CI runs
  Stage 1 on macOS as well as Linux.
* Stage 1 no longer needs a `bgenix` `.bgi` index for a whole-file scan. Without
  one, the multiallelic-rsID filter counts rsIDs from the BGEN variant headers;
  the output is byte-identical.
* Stage 1's memory log reported macOS figures 1024 times too large
  (`ru_maxrss` is in bytes there and in kilobytes on Linux).
* The Stage 1 wheel no longer installs the repository's `tests/` folder as a
  top-level `tests` package. Stage 1 packaging metadata lives only in
  `pyproject.toml`; `setup.py` had disagreed with it on dependencies and on the
  minimum Python.
* `param_gwas_from_file(rtau_file = "*.gz")` now reads gzip without `R.utils`,
  like every other Stage 2 reader.
* `src/Makevars` no longer uses GNU make extensions, and `R CMD check` passes
  with no errors, warnings or notes. The standalone `tests/*.R` scripts are
  excluded from the tarball and run directly by CI and tier 0.

## Documentation

* `beta_sigma2` is documented as the per-allele effect on the variance (it
  always was); the Stage 2 README had called it the log-variance.
* The README states that the `fourth_order_gaussian` kernel was not taken
  forward (issue #1).
* Corrected: the W orientation (T x K), the Stage 1 output table (`.cov.ids.tsv.gz`,
  `d_resid`, `.cov.gz` absent under pooled modes), the SLURM template, and the
  non-existent CLI flags in `docs/CHR22_VALIDATION.md`.
* New: `docs/examples/1kg_chr22/`, a one-minute end-to-end check on 1000
  Genomes chr22 genotypes with simulated null, shift, variance and mixture
  phenotypes, ending in PASS/FAIL against the truth.

## Tests

* Standalone Stage 2 test scripts now test the installed package
  (`tests/load-installed.R`) instead of recompiling `src/` with
  `Rcpp::sourceCpp()`, which failed where `omp.h` was not on the default path.
* Optional dependencies (`bigsnpr`, `fitdistrplus`) are skipped cleanly when
  absent.
* Cross-build floating-point checks use a 1e-12 tolerance instead of bitwise
  equality, because clang on arm64 contracts multiply-adds and differs by about 1 ulp.
* Stage 2 CI had been failing on `main`: dependencies failed to build from
  source (libuv), and a test read `.gz` with `fread()`, which needs `R.utils`.
  CI now installs binary packages, the test uses the package's own gzip reader,
  and CI also runs `R CMD check`.

# fungwasStage2 0.5.0

* The mixture and generic weight builders (`make_weights_normal_mixture()`,
  `make_weights_normal_mixture_vqtl()`, `make_weights_generic()`) take
  `density = c("model", "hybrid", "stage1")`. With `"hybrid"` the
  component-share columns stay on the model density and still sum to 1. The
  membership (`gamma`) column is divided by the Stage 1 density, because the
  Stage 1 coefficient carries that density's scaling. The default `"model"` is
  unchanged: its output is identical to 0.4.0. `"stage1"` reproduces the old
  all-KDE weights and warns. See `docs/MIXTURE_WEIGHT_DENSITY.md`.
* New: `stage1_kde_density()`, which is the Stage 1 Silverman KDE at `q_tau`;
  `read_stage1_density()`, which reads an `fhat`/`taus`/`q_tau` JSON sidecar;
  and `weight_density_mode()`, which reads back the mode a weight matrix was
  built with.
* `param_gwas_from_file()` and `param_gwas_multi_from_file()` gain
  `require_density_mode`, which refuses weights built in another mode.
* `make_weights_generic()` gains `unit_sum` for checking identities such as
  `W_left + W_right = 1` or an additive column of 1.
* The package now has testthat tests (`tests/testthat/`). Tier 0 of the suite
  runs them.
