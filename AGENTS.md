# Notes for developers

This repository holds the two G-MoQ packages, their tests, continuous integration, the user
documentation and the 1000 Genomes example. Keep it that way: simulation campaigns, run records
and one-off diagnostics belong elsewhere.

- **Layout.** `stage1/` is the Python package `gmoq` (command `gmoq-stage1`) with a pybind11 C++
  kernel in `stage1/src/`. `stage2/` is the R package `gmoq`, with C++ via Rcpp. Both carry the
  same version number: bump `stage1/gmoq/__init__.py` and `stage2/DESCRIPTION` together and add an
  entry to `stage2/NEWS.md`.
- **Generated files.** `stage2/NAMESPACE`, `stage2/man/` and the Rcpp glue
  (`stage2/R/RcppExports.R`, `stage2/src/RcppExports.cpp`) are generated: edit the roxygen
  comments and run `Rcpp::compileAttributes("stage2")` and `roxygen2::roxygenise("stage2")`. The
  `useDynLib` directive is declared in `stage2/R/gmoq-package.R`.
- **Tests.** `python -m pytest stage1/tests` (from the repository root, against the installed
  package) and `testthat::test_dir("stage2/tests/testthat", package = "gmoq", load_package =
  "installed")`. Tests build their own simulated data. Add a test for each user-visible behaviour
  you change; do not add tests that pin past bugs or timings.
- **Numerical output.** Changes to Stage 1 or Stage 2 that are meant to be refactors should leave
  the numerical outputs unchanged. Check by running the 1000 Genomes example before and after and
  comparing the decompressed output files.
- **Defaults.** OLS with the pooled covariance is the primary analysis. The jackknife covariance,
  GLS and the alternative density options are kept as checks. Passing tests does not establish
  statistical calibration; the known limitations are in `docs/METHOD.md`.
- **Data.** Never commit individual-level data (genotypes, phenotypes, sample identifiers) or
  results computed on access-controlled cohorts. Examples and tests use simulated data or public
  1000 Genomes genotypes.
