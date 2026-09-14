# ---------------------------------------------------------------------------
# tests/lib/harness.R -- shared machinery for the G-MoQ validation suite
#
# Every tier script is standalone:  Rscript --vanilla tests/tierN/xxx.R
# and derives its own location from --file=, the convention the existing
# stage2/tests scripts already use. Sourcing one via -e source(...) is NOT
# supported and never has been.
#
# Two things this file exists to enforce:
#
#   1. THE LIBRARY PATH. ~/.Rprofile on the analysis cluster contains a bare
#      .libPaths("/user/work/fh6520/r_packages"), which REPLACES the search
#      path rather than prepending to it, so R_LIBS_USER is silently ignored
#      and a stale fungwasStage2 gets loaded instead of the one under test.
#      --vanilla skips the profile; this file then sets the path explicitly
#      from GMOQ_TEST_LIB and ASSERTS afterwards that the package really came
#      from there. See docs/KNOWN_ISSUES.md section 6 trap 2.
#
#   2. MACHINE-READABLE RESULTS. Every check appends one JSON object to
#      $GMOQ_TEST_OUT/results.jsonl. tests/report.R turns those into the single
#      release table. A check that is never recorded cannot pass.
# ---------------------------------------------------------------------------

# --- where am I -------------------------------------------------------------

gmoq_this_file <- function(fallback = NULL) {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) return(normalizePath(f[1L]))
  if (!is.null(fallback)) return(normalizePath(fallback))
  stop("Cannot locate this script. Run it as `Rscript --vanilla <path>`, not via -e source().")
}

# tests/lib/harness.R -> repo root is two levels up from lib/
GMOQ_SUITE_ROOT <- NULL
GMOQ_REPO_ROOT  <- NULL

gmoq_init_roots <- function(script_path) {
  d <- dirname(script_path)                 # e.g. <repo>/tests/tier0
  GMOQ_SUITE_ROOT <<- normalizePath(file.path(d, ".."))   # <repo>/tests
  GMOQ_REPO_ROOT  <<- normalizePath(file.path(d, "..", ".."))
  invisible(TRUE)
}

# --- library path -----------------------------------------------------------

#' Put the package under test first on the search path and prove it took.
#'
#' GMOQ_TEST_LIB is the private library the package under test was installed
#' into. GMOQ_DEP_LIB (optional) holds the CRAN dependencies. Neither is read
#' from R_LIBS_USER, on purpose.
gmoq_use_lib <- function() {
  test_lib <- Sys.getenv("GMOQ_TEST_LIB", "")
  dep_lib  <- Sys.getenv("GMOQ_DEP_LIB", "")
  paths <- c(test_lib, dep_lib, .libPaths())
  paths <- paths[nzchar(paths) & dir.exists(paths)]
  if (length(paths)) .libPaths(unique(paths))
  invisible(.libPaths())
}

#' Load fungwasStage2 and assert it came from the library we asked for.
gmoq_load_stage2 <- function() {
  gmoq_use_lib()
  suppressPackageStartupMessages(library(fungwasStage2))
  where <- normalizePath(find.package("fungwasStage2"))
  test_lib <- Sys.getenv("GMOQ_TEST_LIB", "")
  if (nzchar(test_lib)) {
    want <- normalizePath(test_lib, mustWork = FALSE)
    if (substr(where, 1, nchar(want)) != want) {
      stop("fungwasStage2 loaded from ", where, " but GMOQ_TEST_LIB is ", want,
           ". A bare .libPaths() in ~/.Rprofile overrides R_LIBS_USER silently; ",
           "run with --vanilla and set GMOQ_TEST_LIB.")
    }
  }
  # The compiled kernel must actually be loaded: a NAMESPACE that lost
  # useDynLib installs cleanly and fails at the first C++ call.
  dlls <- names(getLoadedDLLs())
  if (!("fungwasStage2" %in% dlls)) {
    stop("fungwasStage2 is attached but its shared object is not loaded: ",
         "NAMESPACE has lost useDynLib(fungwasStage2, .registration = TRUE). ",
         "It is declared in R/fungwasStage2-package.R -- re-run devtools::document().")
  }
  message("fungwasStage2 loaded from: ", where)
  invisible(where)
}

#' Reach an internal (non-exported) function by name.
gmoq_internal <- function(name) {
  get(name, envir = asNamespace("fungwasStage2"))
}

# --- result recording -------------------------------------------------------

.gmoq_state <- new.env(parent = emptyenv())
.gmoq_state$n_pass <- 0L
.gmoq_state$n_fail <- 0L
.gmoq_state$n_skip <- 0L
.gmoq_state$n_advisory <- 0L
.gmoq_state$n_info <- 0L
.gmoq_state$group  <- NA_character_
.gmoq_state$tier   <- NA_character_
.gmoq_state$out    <- NULL

gmoq_begin <- function(tier, group) {
  .gmoq_state$tier  <- tier
  .gmoq_state$group <- group
  out <- Sys.getenv("GMOQ_TEST_OUT", "")
  if (nzchar(out)) {
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    .gmoq_state$out <- file.path(out, "results.jsonl")
  }
  cat("\n=== ", tier, " / ", group, " ===\n", sep = "")
  invisible(TRUE)
}

.gmoq_json_scalar <- function(x) {
  if (is.null(x) || length(x) == 0L) return("null")
  if (length(x) > 1L) {
    return(paste0("[", paste(vapply(x, .gmoq_json_scalar, ""), collapse = ","), "]"))
  }
  if (is.logical(x)) return(if (is.na(x)) "null" else if (x) "true" else "false")
  if (is.numeric(x)) {
    if (!is.finite(x)) return("null")
    return(format(x, digits = 10, scientific = TRUE, trim = TRUE))
  }
  s <- as.character(x)
  s <- gsub("\\\\", "\\\\\\\\", s)
  s <- gsub("\"", "\\\\\"", s)
  s <- gsub("[\r\n\t]", " ", s)
  paste0("\"", s, "\"")
}

.gmoq_json_obj <- function(lst) {
  paste0("{", paste(sprintf("%s:%s",
                            vapply(names(lst), .gmoq_json_scalar, ""),
                            vapply(lst, .gmoq_json_scalar, "")),
                    collapse = ","), "}")
}

#' Record one check.
#'
#' @param id Stable identifier, e.g. "T0.RIF.equivariance". Used as the key in
#'   the release table, so it must not change casually.
#' @param establishes One sentence: what a PASS lets you conclude.
#' @param pass TRUE/FALSE/NA (NA = skipped).
#' @param stat Optional measured statistic.
#' @param criterion Optional human-readable pass criterion.
#' @param mcse Optional Monte Carlo standard error on `stat`.
#' @param detail Optional free text.
#' @param severity One of:
#'   \describe{
#'     \item{blocking}{a failure means the release is not sound. The default.}
#'     \item{advisory}{a real defect that does not by itself invalidate results
#'       -- a documentation or hygiene problem. Listed in the report; does not
#'       change the verdict.}
#'     \item{report}{a measurement with no pass/fail at all. Recorded so the
#'       number is on the record, never counted.}
#'   }
gmoq_check <- function(id, establishes, pass, stat = NULL, criterion = NULL,
                       mcse = NULL, detail = NULL, severity = "blocking") {
  severity <- match.arg(severity, c("blocking", "advisory", "report"))
  status <- if (severity == "report") "INFO"
            else if (is.na(pass)) "SKIP"
            else if (isTRUE(pass)) "PASS" else "FAIL"
  if (status == "PASS") .gmoq_state$n_pass <- .gmoq_state$n_pass + 1L
  if (status == "FAIL") {
    if (severity == "blocking") .gmoq_state$n_fail <- .gmoq_state$n_fail + 1L
    else .gmoq_state$n_advisory <- .gmoq_state$n_advisory + 1L
  }
  if (status == "SKIP") .gmoq_state$n_skip <- .gmoq_state$n_skip + 1L
  if (status == "INFO") .gmoq_state$n_info <- .gmoq_state$n_info + 1L

  line <- .gmoq_json_obj(list(
    tier = .gmoq_state$tier, group = .gmoq_state$group, id = id,
    status = status, severity = severity, establishes = establishes,
    stat = if (is.null(stat)) NULL else as.numeric(stat),
    mcse = if (is.null(mcse)) NULL else as.numeric(mcse),
    criterion = criterion, detail = detail,
    time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  ))
  if (!is.null(.gmoq_state$out)) cat(line, "\n", sep = "", file = .gmoq_state$out, append = TRUE)

  marker <- switch(status, PASS = "  ok  ", FAIL = if (severity == "blocking") "  FAIL" else "  warn",
                   SKIP = "  skip", INFO = "  note")
  extra <- character(0)
  if (!is.null(stat)) extra <- c(extra, sprintf("stat=%.6g", as.numeric(stat)[1]))
  if (!is.null(mcse)) extra <- c(extra, sprintf("mcse=%.3g", as.numeric(mcse)[1]))
  if (!is.null(criterion)) extra <- c(extra, criterion)
  cat(marker, " ", id, if (length(extra)) paste0("  [", paste(extra, collapse = "; "), "]") else "",
      "\n", sep = "")
  if (status %in% c("FAIL", "INFO") && !is.null(detail)) cat("        ", detail, "\n", sep = "")
  invisible(isTRUE(pass))
}

#' Check that two numeric objects agree.
gmoq_close <- function(id, establishes, a, b, tol, criterion = NULL,
                       severity = "blocking") {
  a <- as.numeric(a); b <- as.numeric(b)
  if (length(a) != length(b)) {
    return(gmoq_check(id, establishes, FALSE, severity = severity,
                      detail = sprintf("length %d vs %d", length(a), length(b))))
  }
  d <- max(abs(a - b))
  gmoq_check(id, establishes, is.finite(d) && d <= tol, stat = d, severity = severity,
             criterion = if (is.null(criterion)) sprintf("max|diff| <= %.3g", tol) else criterion,
             detail = if (!is.finite(d) || d > tol) sprintf("max|diff| = %.6g > %.3g", d, tol) else NULL)
}

#' Check that an expression raises an error whose message matches a pattern.
#' This is how every refusal path is tested: a refusal that becomes a warning,
#' or that stops with a different reason, is a regression.
gmoq_refuses <- function(id, establishes, expr, pattern, severity = "blocking") {
  res <- tryCatch({ force(expr); list(ok = FALSE, msg = "no error raised") },
                  error = function(e) list(ok = TRUE, msg = conditionMessage(e)))
  hit <- res$ok && grepl(pattern, res$msg, perl = TRUE)
  gmoq_check(id, establishes, hit, severity = severity,
             criterion = sprintf("stops, message matches /%s/", pattern),
             detail = if (!hit) paste0("got: ", substr(res$msg, 1, 300)) else NULL)
}

gmoq_finish <- function() {
  cat(sprintf("\n-- %s / %s: %d passed, %d FAILED, %d advisory, %d skipped, %d reported\n",
              .gmoq_state$tier, .gmoq_state$group,
              .gmoq_state$n_pass, .gmoq_state$n_fail, .gmoq_state$n_advisory,
              .gmoq_state$n_skip, .gmoq_state$n_info))
  if (.gmoq_state$n_fail > 0L) quit(status = 1L, save = "no")
  invisible(TRUE)
}

# --- Monte Carlo helpers ----------------------------------------------------

#' Monte Carlo standard error of an observed rate.
gmoq_mcse <- function(rate, n) sqrt(pmax(rate * (1 - rate), 0) / n)

#' Wilson score interval, which behaves at rates near 0 and 1 where the
#' Wald interval does not. Used for every counted rate in tier 2.
gmoq_wilson <- function(k, n, conf = 0.99) {
  z <- stats::qnorm(1 - (1 - conf) / 2)
  p <- k / n
  d <- 1 + z^2 / n
  c(lo = (p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d,
    hi = (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d)
}
