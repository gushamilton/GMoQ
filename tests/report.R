#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# The single artefact a human reads for the historical regression checks.
#
#   Rscript --vanilla tests/report.R [results_root] [out_dir]
#
# Reads every tier's results.jsonl and writes
#
#   RELEASE_REPORT.json   machine readable: every check, its statistic, its
#                         Monte Carlo error, its criterion and its verdict
#   RELEASE_REPORT.md     the same thing as a table, with the verdict first
#
# THE VERDICT RULE, stated before any numbers are seen:
#
#   REGRESSION CHECKS PASSED   every blocking check passed, and tiers 0, 1 and
#                              2 all ran.
#   REGRESSION CHECKS PARTIAL  every blocking check that ran passed, but a tier
#                              did not run.
#   REGRESSION CHECKS FAILED   at least one blocking check failed.
#
# These are historical regression labels. They do not certify scientific
# validity, statistical calibration, or release readiness; use the paired
# known-truth validity pilot for that question.
#
# A tier that was not run is never counted as a pass. The report says what was
# not established as prominently as what was.
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  .self <- local({
    a <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (!length(f)) stop("Run as: Rscript --vanilla tests/report.R")
    normalizePath(f[1L])
  })
}))
suite <- dirname(.self)
repo <- normalizePath(file.path(suite, ".."))

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1L] else
  file.path(Sys.getenv("GMOQ_WORK", tempdir()), "results")
outdir <- if (length(args) >= 2) args[2L] else root
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  lib <- c(Sys.getenv("GMOQ_TEST_LIB"), Sys.getenv("GMOQ_DEP_LIB"))
  .libPaths(unique(c(lib[nzchar(lib)], .libPaths())))
}

TIERS <- "tier0"   # tiers 1-4 are at tag archive/research-records-20261016
TIER_MEANS <- c(
  tier0 = "Unit: closed-form identities. Establishes that the arithmetic is right. Establishes nothing about statistical validity.")

rows <- list()
present <- character(0)
for (tr in TIERS) {
  f <- file.path(root, tr, "results.jsonl")
  if (!file.exists(f)) next
  present <- c(present, tr)
  lines <- readLines(f, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  for (ln in lines) {
    o <- tryCatch(jsonlite::fromJSON(ln), error = function(e) NULL)
    if (is.null(o)) next
    gv <- function(k, d = NA) { v <- o[[k]]; if (is.null(v) || !length(v)) d else v[[1]] }
    rows[[length(rows) + 1L]] <- data.frame(
      tier = gv("tier", tr), group = gv("group", ""), id = gv("id", ""),
      status = gv("status", ""), severity = gv("severity", "blocking"),
      stat = suppressWarnings(as.numeric(gv("stat", NA))),
      mcse = suppressWarnings(as.numeric(gv("mcse", NA))),
      criterion = gv("criterion", ""), establishes = gv("establishes", ""),
      detail = gv("detail", ""), stringsAsFactors = FALSE)
  }
}
if (!length(rows)) stop("No results found under ", root,
                        ". Run tests/run_tier0.sh first.")
tab <- do.call(rbind, rows)
tab <- tab[!duplicated(paste(tab$tier, tab$id), fromLast = TRUE), , drop = FALSE]

blocking_fail <- tab[tab$status == "FAIL" & tab$severity == "blocking", , drop = FALSE]
advisory_fail <- tab[tab$status == "FAIL" & tab$severity == "advisory", , drop = FALSE]
skipped <- tab[tab$status == "SKIP", , drop = FALSE]
info <- tab[tab$status == "INFO", , drop = FALSE]
passed <- tab[tab$status == "PASS", , drop = FALSE]

missing_tiers <- setdiff(TIERS, present)
verdict <- if (nrow(blocking_fail) > 0) {
  "REGRESSION CHECKS FAILED"
} else if (length(missing_tiers)) {
  "REGRESSION CHECKS PARTIAL"
} else {
  "REGRESSION CHECKS PASSED"
}

pkg_version <- tryCatch({
  d <- read.dcf(file.path(repo, "stage2", "DESCRIPTION")); d[1, "Version"]
}, error = function(e) NA_character_)
commit <- tryCatch(trimws(system2("git", c("-C", repo, "rev-parse", "HEAD"),
                                  stdout = TRUE, stderr = FALSE)),
                   error = function(e) NA_character_)
if (!length(commit) || !nzchar(commit[1])) commit <- NA_character_

payload <- list(
  verdict = verdict,
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  package_version = pkg_version,
  commit = commit[1],
  tiers_run = present,
  tiers_not_run = setdiff(TIERS, present),
  tiers_not_run_affecting_verdict = missing_tiers,
  counts = list(pass = nrow(passed), blocking_fail = nrow(blocking_fail),
                advisory_fail = nrow(advisory_fail), skipped = nrow(skipped),
                reported = nrow(info), total = nrow(tab)),
  verdict_rule = paste(
    "REGRESSION CHECKS PASSED: every blocking check passed and tier 0 ran.",
    "REGRESSION CHECKS PARTIAL: every blocking check that ran passed, but a tier did not run.",
    "REGRESSION CHECKS FAILED: at least one blocking check failed.",
    "Advisory failures are real defects that do not by themselves invalidate results.",
    "A tier that did not run is never counted as a pass.",
    "These regression labels do not certify scientific validity; see docs/examples/1kg_chr22 for recovery of known effects."),
  checks = tab)
writeLines(jsonlite::toJSON(payload, auto_unbox = TRUE, na = "null", digits = 8,
                            pretty = TRUE),
           file.path(outdir, "RELEASE_REPORT.json"))

# ---- markdown --------------------------------------------------------------
fmt <- function(x) ifelse(is.na(x), "", formatC(x, digits = 5, format = "g"))
mdrow <- function(r) sprintf("| %s | %s | %s | %s | %s | %s |",
                             r$id, r$status,
                             if (identical(r$severity, "blocking")) "" else r$severity,
                             fmt(r$stat),
                             if (is.na(r$mcse)) "" else fmt(r$mcse),
                             gsub("|", "/", r$criterion, fixed = TRUE))

md <- c(
  sprintf("# G-MoQ historical regression checks: **%s**", verdict), "",
  sprintf("Generated %s. Package version `%s`. Commit `%s`.",
          payload$generated, ifelse(is.na(pkg_version), "?", pkg_version),
          ifelse(is.na(commit[1]), "?", substr(commit[1], 1, 12))), "",
  sprintf("%d checks: **%d passed**, **%d blocking failures**, %d advisory, %d skipped, %d reported-only.",
          nrow(tab), nrow(passed), nrow(blocking_fail), nrow(advisory_fail),
          nrow(skipped), nrow(info)), "",
  "**Verdict rule.** REGRESSION CHECKS PASSED means every blocking check passed and tiers 0-2 all ran.",
  "These labels are historical regression results and do not establish scientific validity.",
  "Advisory failures are real defects that do not by themselves invalidate a result.",
  "A tier that did not run is never counted as a pass; it is listed below as not established.",
  "")

absent <- setdiff(TIERS, present)
if (length(absent)) {
  md <- c(md, "## Not established (tier did not run)", "",
          paste0("- **", absent, "**",
                 ifelse(absent %in% missing_tiers, " (changes the verdict)",
                        " (does not change the verdict: it establishes cost, not correctness)"),
                 " -- ", TIER_MEANS[absent]), "")
}
if (nrow(blocking_fail)) {
  md <- c(md, "## Blocking failures", "",
          "| check | statistic | criterion | what a pass would have established |",
          "|---|---|---|---|",
          apply(blocking_fail, 1, function(r) sprintf("| `%s` | %s | %s | %s |",
                r[["id"]], fmt(as.numeric(r[["stat"]])),
                gsub("|", "/", r[["criterion"]], fixed = TRUE),
                gsub("|", "/", r[["establishes"]], fixed = TRUE))), "",
          "Details:", "",
          paste0("- `", blocking_fail$id, "`: ", blocking_fail$detail), "")
}
if (nrow(advisory_fail)) {
  md <- c(md, "## Advisory failures (do not change the verdict)", "",
          paste0("- `", advisory_fail$id, "` -- ", advisory_fail$establishes,
                 " (", advisory_fail$detail, ")"), "")
}

md <- c(md, "## Every check", "")
for (tr in present) {
  sub <- tab[tab$tier == tr, , drop = FALSE]
  md <- c(md, sprintf("### %s", tr), "", TIER_MEANS[[tr]], "",
          "| check | status | severity | statistic | MC error | criterion |",
          "|---|---|---|---|---|---|")
  for (i in seq_len(nrow(sub))) md <- c(md, mdrow(sub[i, ]))
  md <- c(md, "")
}

md <- c(md, "## Measurements reported without a pass/fail", "")
if (nrow(info)) {
  md <- c(md, "| measurement | value | what it is |", "|---|---|---|",
          apply(info, 1, function(r) sprintf("| `%s` | %s | %s |", r[["id"]],
                fmt(as.numeric(r[["stat"]])),
                gsub("|", "/", r[["establishes"]], fixed = TRUE))), "")
} else md <- c(md, "_none_", "")

md <- c(md, "## What this report does NOT establish", "",
        "Tier 0 checks arithmetic only: no real genotypes, no LD, no relatedness,",
        "no measurement at p = 5e-8 and no comparison against any other method.",
        "`docs/examples/1kg_chr22` checks recovery of known effects on real genotypes.",
        "")
writeLines(md, file.path(outdir, "RELEASE_REPORT.md"))

cat("\n=============================================================\n")
cat("  G-MoQ historical regression checks:", verdict, "\n")
cat(sprintf("  %d passed, %d blocking failures, %d advisory, %d skipped, %d reported\n",
            nrow(passed), nrow(blocking_fail), nrow(advisory_fail),
            nrow(skipped), nrow(info)))
if (length(missing_tiers)) cat("  NOT RUN:", paste(missing_tiers, collapse = ", "), "\n")
if (nrow(blocking_fail)) {
  cat("  blocking failures:\n")
  cat(paste0("    - ", blocking_fail$id, "\n"), sep = "")
}
cat("  ->", file.path(outdir, "RELEASE_REPORT.md"), "\n")
cat("=============================================================\n")
if (nrow(blocking_fail)) quit(status = 1L, save = "no")
