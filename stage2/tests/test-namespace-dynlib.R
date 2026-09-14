# NAMESPACE must load the compiled kernel.
#
# Regression test for a defect that reached three commits: running the R suite
# invokes devtools::document(), which rewrites NAMESPACE from roxygen and drops
# anything hand-written. The C++ block was hand-written, so it vanished, and a
# package without it fails at run time on every path that reaches C++ -- with no
# build-time error to catch it.
#
# The directives now live in R/fungwasStage2-package.R so roxygen regenerates
# them. This test fails if they are ever lost again.

this_file <- tryCatch({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(f) else normalizePath(sys.frames()[[1]]$ofile)
}, error = function(e) normalizePath("tests/test-namespace-dynlib.R"))
pkg_root <- normalizePath(file.path(dirname(this_file), ".."))

ns <- readLines(file.path(pkg_root, "NAMESPACE"), warn = FALSE)

stopifnot(
  "NAMESPACE lost useDynLib(fungwasStage2, .registration = TRUE); the compiled kernel will not load. It is declared in R/fungwasStage2-package.R -- re-run devtools::document() rather than editing NAMESPACE." =
    any(grepl("useDynLib\\(\\s*fungwasStage2", ns)),
  "NAMESPACE lost importFrom(Rcpp, sourceCpp). Declared in R/fungwasStage2-package.R." =
    any(grepl("importFrom\\(\\s*Rcpp\\s*,\\s*sourceCpp\\s*\\)", ns))
)

# and the directives must be where roxygen can regenerate them
pkgdoc <- file.path(pkg_root, "R", "fungwasStage2-package.R")
stopifnot(
  "R/fungwasStage2-package.R is missing; without it document() will drop the C++ block again." =
    file.exists(pkgdoc)
)
src <- readLines(pkgdoc, warn = FALSE)
stopifnot(
  "R/fungwasStage2-package.R no longer carries @useDynLib." =
    any(grepl("@useDynLib\\s+fungwasStage2", src)),
  "R/fungwasStage2-package.R no longer carries @importFrom Rcpp sourceCpp." =
    any(grepl("@importFrom\\s+Rcpp\\s+sourceCpp", src))
)

cat("ok: NAMESPACE carries the C++ directives, and roxygen can regenerate them\n")
