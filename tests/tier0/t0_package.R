#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Tier 0 -- the package itself: it must install, it must load its kernel, and
# its documentation must not break installation.
#
# Three defects this exists to stop recurring:
#
#   1. devtools::document() rewrites NAMESPACE from roxygen and silently drops
#      anything hand-written. useDynLib lived there by hand, so it vanished and
#      the package installed cleanly while every C++ path failed at run time.
#      It now lives in R/fungwasStage2-package.R.
#   2. A NAMESPACE.backup was once committed alongside it.
#   3. An apostrophe inside \code{} in a generated .Rd file broke installation.
#      tools::checkRd is the tool that sees that class of problem without a
#      full R CMD check.
#
# The install itself is verified by tests/run_tier0.sh, which installs into a
# private library before running anything; this script then asserts that what
# loaded is what was installed.
#
# Run:  Rscript --vanilla tests/tier0/t0_package.R
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  .self <- local({
    a <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (!length(f)) stop("Run as: Rscript --vanilla tests/tier0/t0_package.R")
    normalizePath(f[1L])
  })
  source(file.path(dirname(.self), "..", "lib", "harness.R"))
}))
gmoq_init_roots(.self)
where <- gmoq_load_stage2()
gmoq_begin("tier0", "package")

src <- file.path(GMOQ_REPO_ROOT, "stage2")

# --- 1. the source NAMESPACE, and where its directives come from -----------
ns_src <- readLines(file.path(src, "NAMESPACE"), warn = FALSE)
gmoq_check("T0.PKG.namespace_usedynlib",
           "The source NAMESPACE carries useDynLib(fungwasStage2, .registration = TRUE).",
           any(grepl("useDynLib\\(\\s*fungwasStage2", ns_src)),
           criterion = "useDynLib present in stage2/NAMESPACE")
gmoq_check("T0.PKG.namespace_importfrom_rcpp",
           "The source NAMESPACE carries importFrom(Rcpp, sourceCpp).",
           any(grepl("importFrom\\(\\s*Rcpp\\s*,\\s*sourceCpp\\s*\\)", ns_src)),
           criterion = "importFrom(Rcpp, sourceCpp) present")

pkgdoc <- file.path(src, "R", "fungwasStage2-package.R")
doc <- if (file.exists(pkgdoc)) readLines(pkgdoc, warn = FALSE) else character(0)
gmoq_check("T0.PKG.directives_are_roxygen_regenerable",
           "The C++ directives live in R/fungwasStage2-package.R, so document() regenerates them.",
           file.exists(pkgdoc) &&
             any(grepl("@useDynLib\\s+fungwasStage2", doc)) &&
             any(grepl("@importFrom\\s+Rcpp\\s+sourceCpp", doc)),
           criterion = "both directives present as roxygen tags")

gmoq_check("T0.PKG.no_namespace_backup",
           "No NAMESPACE.backup is present in the tree -- one was committed once.",
           !any(file.exists(file.path(src, c("NAMESPACE.backup", "NAMESPACE.bak")))),
           criterion = "stage2/NAMESPACE.backup absent")

# --- 2. the INSTALLED package -----------------------------------------------
ns_inst <- readLines(file.path(where, "NAMESPACE"), warn = FALSE)
gmoq_check("T0.PKG.installed_namespace_usedynlib",
           "The installed NAMESPACE also carries useDynLib: the shipped artefact, not just the source.",
           any(grepl("useDynLib\\(\\s*fungwasStage2", ns_inst)),
           criterion = "useDynLib present in the installed NAMESPACE")

so <- list.files(file.path(where, "libs"), pattern = "\\.(so|dylib|dll)$", full.names = TRUE)
gmoq_check("T0.PKG.shared_object_present",
           "The compiled shared object was built and installed.",
           length(so) >= 1L, criterion = "at least one .so/.dylib/.dll in libs/")
gmoq_check("T0.PKG.dll_loaded",
           "The shared object is actually loaded after library(fungwasStage2).",
           "fungwasStage2" %in% names(getLoadedDLLs()),
           criterion = "getLoadedDLLs() names fungwasStage2")

# Every registered C++ entry point must be callable. A NAMESPACE that lost
# useDynLib passes every R-only test and fails here.
kernels <- c("compute_calibrated_se", "compute_param_cov", "compute_gls_param_gwas",
             "compute_q_calibration", "fit_multi_models", "fit_multi_models_tau_cov")
callable <- vapply(kernels, function(f) is.function(get0(f, envir = asNamespace("fungwasStage2"))),
                   logical(1))
gmoq_check("T0.PKG.kernels_present",
           "All six compiled kernels are reachable from the package namespace.",
           all(callable), criterion = paste(kernels, collapse = ", "),
           detail = if (!all(callable)) paste("missing:", paste(kernels[!callable], collapse = ", ")) else NULL)

# and one of them actually executes, so registration is real and not just a symbol
probe <- tryCatch({
  S <- diag(3); u <- c(S[1,1],S[1,2],S[1,3],S[2,2],S[2,3],S[3,3])
  A <- matrix(c(1,1,1)/3, nrow = 1)
  as.numeric(gmoq_internal("compute_calibrated_se")(u, 0, A, 3L))
}, error = function(e) NA_real_)
gmoq_check("T0.PKG.kernel_executes",
           "A compiled kernel runs and returns the right answer, proving .registration worked.",
           isTRUE(abs(probe - sqrt(3 * (1/3)^2)) < 1e-12),
           stat = probe, criterion = "compute_calibrated_se on I_3 returns 1/sqrt(3)")

# --- 3. the exported API surface --------------------------------------------
exports <- getNamespaceExports("fungwasStage2")
must_export <- c("quantile_gwas", "param_gwas", "param_gwas_from_file", "param_gwas_multi",
                 "param_gwas_multi_from_file", "make_weights_vqtl",
                 "make_weights_normal_mixture", "make_weights_normal_mixture_vqtl",
                 "make_weights_generic", "make_fd_grad", "read_sigma0",
                 "dispersion_skew_guard", "estimate_rav", "residualise_dispersion",
                 "dispersion_report", "mashr_inputs_from_cov", "prepare_regenie")
missing <- setdiff(must_export, exports)
gmoq_check("T0.PKG.api_surface",
           "The public API a release is expected to provide is exported.",
           length(missing) == 0L, criterion = sprintf("%d expected exports", length(must_export)),
           detail = if (length(missing)) paste("missing:", paste(missing, collapse = ", ")) else NULL)

# --- 4. the documentation cannot break installation ------------------------
rd_dir <- file.path(src, "man")
rds <- list.files(rd_dir, pattern = "\\.Rd$", full.names = TRUE)
bad <- character(0)
for (f in rds) {
  msgs <- tryCatch(
    utils::capture.output(tools::checkRd(f), type = "message"),
    error = function(e) paste0("ERROR: ", conditionMessage(e)))
  hard <- grep("^(Error|.*\\bError\\b)", msgs, value = TRUE)
  if (length(hard)) bad <- c(bad, paste0(basename(f), ": ", hard[1]))
}
gmoq_check("T0.PKG.rd_parses",
           "Every man page parses: an Rd that does not is what once broke installation.",
           length(bad) == 0L, stat = length(rds),
           criterion = sprintf("tools::checkRd raises no error on any of %d .Rd files", length(rds)),
           detail = if (length(bad)) paste(head(bad, 5), collapse = " | ") else NULL)

# Every export needs a man page, or R CMD check fails at release time.
alias <- unlist(lapply(rds, function(f) {
  x <- readLines(f, warn = FALSE)
  sub("^\\\\alias\\{(.*)\\}\\s*$", "\\1", grep("^\\\\alias\\{", x, value = TRUE))
}))
undoc <- setdiff(exports, alias)
gmoq_check("T0.PKG.exports_documented",
           "Every exported function has an \\alias in man/, so R CMD check will not fail at release.",
           length(undoc) == 0L, severity = "advisory",
           criterion = "no undocumented exports",
           detail = if (length(undoc)) paste("undocumented:", paste(undoc, collapse = ", ")) else NULL)

# Rd warnings are recorded but do not block: they are a documentation quality
# signal, not a correctness one.
warn_n <- 0L
for (f in rds) {
  msgs <- tryCatch(utils::capture.output(tools::checkRd(f), type = "message"),
                   error = function(e) character(0))
  warn_n <- warn_n + length(grep("checkRd", msgs))
}
gmoq_check("T0.PKG.rd_warnings",
           "Number of tools::checkRd warnings across man/ -- reported, not asserted.",
           TRUE, stat = warn_n, severity = "report",
           criterion = "reported")

# --- 5. version and provenance ---------------------------------------------
desc <- read.dcf(file.path(where, "DESCRIPTION"))
gmoq_check("T0.PKG.version_recorded",
           "The installed package reports a version, which the release report records.",
           nzchar(desc[1, "Version"]),
           criterion = "DESCRIPTION Version non-empty",
           detail = paste("version:", desc[1, "Version"]))

gmoq_finish()
