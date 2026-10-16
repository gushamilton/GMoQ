#!/usr/bin/env Rscript

## Build Stage 2 weights on exactly the baseline used by a fresh LOCO Stage 1
## chromosome scan.  This script deliberately writes only aggregate baseline
## parameters and the tau-by-parameter W matrix; it never writes phenotype,
## covariate, LOCO, or sample-ID arrays.

args <- commandArgs(trailingOnly = TRUE)
if (!grepl("^bp1-compute", Sys.info()[["nodename"]])) {
  stop("Run matched-weight preparation on a BP compute node (bp1-compute*).",
       call. = FALSE)
}
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) args[[i + 1L]] else default
}
r_lib <- get_arg("--r-lib")
if (!is.null(r_lib)) .libPaths(c(r_lib, .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
})

stop_usage <- function() {
  stop(paste(
    "Usage: Rscript prepare_matched_weights.R",
    "--trait <trait> --chr <1..22> --stage1 <file>",
    "--pheno <file> --covar <file> --keep <file> [--loco <file>]",
    "[--loco-mode <subtract|covariate>]  (default: matches Stage 1)",
    "--bgen-sample <file> --family <shifted_lognormal|skew|mixture>",
    "--out-dir <directory> [--loco-scale regression] [--seed 20260915]",
    "[--expected-n N --expected-mean M --expected-sd S --expected-scale A]",
    sep = "\n"
  ), call. = FALSE)
}

trait <- get_arg("--trait")
chr_raw <- get_arg("--chr")
stage1_file <- get_arg("--stage1")
pheno_file <- get_arg("--pheno")
covar_file <- get_arg("--covar")
keep_file <- get_arg("--keep")
loco_file <- get_arg("--loco")
## Mirrors Stage 1's --loco-mode. Under 'covariate' Stage 1 leaves Y unadjusted
## and puts L(chr) in the residualising basis, so the matched W baseline MUST be
## raw Y. Under 'subtract' it is Y - scale*L(chr), the historical behaviour.
loco_mode <- get_arg("--loco-mode", if (is.null(loco_file)) "covariate" else "subtract")
if (!loco_mode %in% c("subtract", "covariate")) {
  stop("--loco-mode must be subtract or covariate.", call. = FALSE)
}
bgen_sample_file <- get_arg("--bgen-sample")
family <- get_arg("--family")
out_dir <- get_arg("--out-dir")
loco_scale_mode <- get_arg("--loco-scale", "regression")
seed <- as.integer(get_arg("--seed", "20260915"))
optional_num <- function(flag, integer = FALSE) {
  raw <- get_arg(flag)
  if (is.null(raw)) return(if (integer) NA_integer_ else NA_real_)
  suppressWarnings(if (integer) as.integer(raw) else as.numeric(raw))
}
expected_n <- optional_num("--expected-n", integer = TRUE)
expected_mean <- optional_num("--expected-mean")
expected_sd <- optional_num("--expected-sd")
expected_scale <- optional_num("--expected-scale")

if (any(vapply(list(trait, chr_raw, stage1_file, pheno_file, covar_file,
                    keep_file, bgen_sample_file, family, out_dir),
               is.null, logical(1)))) stop_usage()
if (!trait %in% c(
  "bmi", "log_bmi", "weight",
  "testosterone", "testosterone_shrink_0.2", "testosterone_shrink_0.5",
  "testosterone_shrink_0.7", "jcv_log_all", "mcv_log_all"
)) stop("Unsupported paper trait: ", trait, call. = FALSE)

chr_num <- suppressWarnings(as.integer(chr_raw))
if (!is.finite(chr_num) || chr_num < 1L || chr_num > 22L) {
  stop("--chr must be an autosome from 1 through 22.", call. = FALSE)
}
chr <- as.character(chr_num)

if (is.null(family) || !family %in% c("shifted_lognormal", "skew", "mixture")) {
  stop("--family must be shifted_lognormal, skew, or mixture.", call. = FALSE)
}
expected_family <- if (trait %in% c("bmi", "log_bmi", "weight")) {
  "shifted_lognormal"
} else if (grepl("^testosterone", trait)) {
  "skew"
} else {
  "mixture"
}
if (!identical(family, expected_family)) {
  stop("Trait ", trait, " requires family ", expected_family,
       "; got ", family, ".", call. = FALSE)
}
if (!is.finite(seed)) seed <- 20260915L

if (identical(loco_mode, "subtract") && is.null(loco_file)) {
  stop("--loco-mode subtract requires --loco.", call. = FALSE)
}
required <- c(stage1_file, pheno_file, covar_file, keep_file,
              if (!is.null(loco_file)) loco_file else NULL,
              bgen_sample_file)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing required input(s): ", paste(missing, collapse = ", "),
                         call. = FALSE)

normalize_chr <- function(x) {
  s <- trimws(as.character(x))
  s <- sub("^(chr|CHR)", "", s)
  is_num <- grepl("^[0-9]+$", s)
  s[is_num] <- as.character(as.integer(s[is_num]))
  toupper(s)
}

read_json <- function(path) {
  if (!file.exists(path)) return(NULL)
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to validate Stage 1 covariance metadata: ", path,
         call. = FALSE)
  }
  jsonlite::fromJSON(path, simplifyVector = TRUE)
}

## Stage 1 serialises the 45-level grid in its beta_tau_* column names.  The
## values recovered here are the only tau values used by any weight builder.
stage1_header <- names(fread(stage1_file, nrows = 0L))
beta_cols <- grep("^beta_tau_", stage1_header, value = TRUE)
if (length(beta_cols) != 45L) {
  stop("Expected 45 serialized beta_tau_* columns in ", stage1_file,
       "; found ", length(beta_cols), ".", call. = FALSE)
}
tau_labels <- as.numeric(sub("^beta_tau_", "", beta_cols))
if (any(!is.finite(tau_labels)) || any(diff(tau_labels) <= 0)) {
  stop("Stage 1 tau columns are not finite and strictly increasing.", call. = FALSE)
}
## Stage 1 labels the canonical grid to three decimals, but the production
## runner may pass the underlying full-precision linspace.  Recover that
## canonical grid when the labels identify it; preserve custom grids at the
## precision available in their labels.  The labels remain the sole column
## indexing key, while the recovered values keep Stage 2 thresholds/densities
## on the same probabilities used by Stage 1.
expected_taus <- round(seq(0.02, 0.98, length.out = 45L), 3L)
if (!isTRUE(all.equal(tau_labels, expected_taus, tolerance = 5e-7))) {
  stop("Stage 1 tau grid is not the pinned 45-level grid; refusing to build W.",
       call. = FALSE)
}
taus_full <- seq(0.02, 0.98, length.out = 45L)
taus <- if (isTRUE(all.equal(round(taus_full, 3L), tau_labels, tolerance = 5e-7))) {
  taus_full
} else {
  tau_labels
}

## Keep files are two-column FID/IID lists without a header.  Stage 1 keys this
## filter on FID, so the preparation follows that rule exactly.
keep_dt <- fread(keep_file, header = FALSE, colClasses = "character", showProgress = FALSE)
if (!ncol(keep_dt)) stop("Keep file has no columns: ", keep_file, call. = FALSE)
keep_ids <- unique(as.character(keep_dt[[1L]]))

pheno_hdr <- names(fread(pheno_file, nrows = 0L))
if (!all(c("FID", "IID", trait) %in% pheno_hdr)) {
  stop("Phenotype file must contain FID, IID, and ", trait, ".", call. = FALSE)
}
ph <- fread(pheno_file, select = c("FID", "IID", trait), showProgress = FALSE)
setnames(ph, c("FID", "IID", trait), c("FID", "IID", "y"))
ph[, FID := as.character(FID)]
ph[, IID := as.character(IID)]

cov_hdr <- names(fread(covar_file, nrows = 0L))
if (!all(c("FID", "IID") %in% cov_hdr)) {
  stop("Covariate file must contain FID and IID: ", covar_file, call. = FALSE)
}
covar_cols <- setdiff(cov_hdr, c("FID", "IID"))
if (!length(covar_cols)) stop("Covariate file has no covariate columns: ", covar_file,
                              call. = FALSE)
cov <- fread(covar_file, select = c("FID", "IID", covar_cols), showProgress = FALSE)
cov[, FID := as.character(FID)]
if (anyDuplicated(cov$FID)) stop("Covariate FID values are not unique: ", covar_file,
                                call. = FALSE)

## Match covariates by FID in phenotype-file order, as Stage 1 does.  Keeping
## this order also makes LOCO prediction alignment explicit.
cov_idx <- match(ph$FID, cov$FID)
if (anyNA(cov_idx)) {
  ph <- ph[!is.na(cov_idx)]
  cov_idx <- cov_idx[!is.na(cov_idx)]
}
for (nm in covar_cols) ph[[nm]] <- cov[[nm]][cov_idx]
ph <- ph[FID %in% keep_ids]
for (nm in c("y", covar_cols)) ph <- ph[is.finite(as.numeric(get(nm)))]
if (!nrow(ph)) stop("No complete phenotype/covariate rows remain after keep filtering.",
                   call. = FALSE)
if (anyDuplicated(ph$FID)) {
  stop("Phenotype FID values are not unique; refusing ambiguous sample matching.",
       call. = FALSE)
}

## Read the BGEN sample IDs.  The UK Biobank sample file has two metadata
## lines; the first data column is ID_1, which is the identifier Stage 1 uses
## for BGEN matching.
bgen <- fread(bgen_sample_file, skip = 2L, header = FALSE,
              colClasses = "character", showProgress = FALSE)
if (!nrow(bgen) || !ncol(bgen)) stop("BGEN sample file has no sample rows: ",
                                    bgen_sample_file, call. = FALSE)
bgen_ids <- as.character(bgen[[1L]])
bgen_ids <- bgen_ids[nzchar(bgen_ids) & bgen_ids != "0"]
if (anyDuplicated(bgen_ids)) stop("BGEN sample IDs are not unique: ", bgen_sample_file,
                                 call. = FALSE)
bgen_match <- ph$FID %in% bgen_ids
if (!any(bgen_match)) stop("No phenotype rows match the BGEN sample file.", call. = FALSE)
if (mean(bgen_match) < 0.01) {
  stop("BGEN sample match is below 1%; refusing a likely identifier mismatch.",
       call. = FALSE)
}

## Resolve either a direct .loco path or the REGENIE *_pred.list master file.
## This mirrors Stage 1's exact label/canonical-label rule; substring matching
## is deliberately disallowed because bmi/log_bmi and testosterone prefixes
## overlap.
sanitize_pheno <- function(x) {
  x <- as.character(x)
  gsub("[^[:alnum:]_-]", "_", x, perl = TRUE)
}
resolve_loco_file <- function(path, phenotype) {
  if (grepl("\\.loco$", path, ignore.case = TRUE)) return(path)
  lines <- readLines(path, warn = FALSE)
  entries <- strsplit(trimws(lines[nzchar(trimws(lines))]), "[[:space:]]+")
  entries <- entries[vapply(entries, function(x) length(x) >= 1L, logical(1))]
  loco_entries <- entries[vapply(entries, function(x) {
    grepl("\\.loco$", x[[length(x)]], ignore.case = TRUE)
  }, logical(1))]
  if (!length(loco_entries)) stop("No .loco files found in LOCO pred.list: ", path,
                                 call. = FALSE)
  loco_paths <- vapply(loco_entries, function(x) x[[length(x)]], character(1))
  if (length(unique(loco_paths)) == 1L) return(loco_paths[[1L]])
  labels <- vapply(loco_entries, `[[`, character(1), 1L)
  candidates <- loco_paths[labels == phenotype]
  if (!length(candidates)) {
    candidates <- loco_paths[sanitize_pheno(labels) == sanitize_pheno(phenotype)]
  }
  candidates <- unique(candidates)
  if (length(candidates) != 1L) {
    stop("Could not resolve exactly one LOCO file for ", phenotype,
         " in ", path, "; matched ", length(candidates), ".", call. = FALSE)
  }
  candidates[[1L]]
}

## Read the LOCO file one chromosome line at a time.  A line contains the
## chromosome label followed by one prediction per LOCO header sample.  This
## avoids materialising a 22 x 486k data.table and retains only the predictions
## needed for the scale regression and requested chromosome.
read_loco <- function(path, sample_ids, sample_iid, target_chr) {
  con <- file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  header_line <- readLines(con, n = 1L, warn = FALSE)
  if (!length(header_line)) stop("Empty LOCO file: ", path, call. = FALSE)
  header <- scan(text = header_line, what = character(), quiet = TRUE)
  if (length(header) < 2L) stop("Malformed LOCO header: ", path, call. = FALSE)
  loco_ids <- header[-1L]

  ## Stage 1 accepts several aliases.  Resolve in the same order, while
  ## retaining one integer position per phenotype row.
  ## Exactly the order used by Stage 1's build_sample_aliases()/
  ## resolve_sample_column(): FID, IID, FID_IID, "FID IID", FID:IID.
  pos <- match(sample_ids, loco_ids)
  aliases <- list(
    as.character(sample_iid),
    paste(sample_ids, sample_iid, sep = "_"),
    paste(sample_ids, sample_iid, sep = " "),
    paste(sample_ids, sample_iid, sep = ":")
  )
  for (alias in aliases) {
    miss <- is.na(pos)
    if (any(miss)) pos[miss] <- match(alias[miss], loco_ids)
  }

  pred_by_chr <- list()
  target <- NULL
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) break
    line <- trimws(line)
    if (!nzchar(line)) next
    first_space <- regexpr("[[:space:]]", line, perl = TRUE)[1L]
    if (first_space < 1L) next
    key <- normalize_chr(substr(line, 1L, first_space - 1L))
    values_text <- sub("^[^[:space:]]+[[:space:]]+", "", line, perl = TRUE)
    pred_all <- scan(text = values_text, what = numeric(), quiet = TRUE)
    if (length(pred_all) != length(loco_ids)) {
      stop("LOCO chromosome ", key, " has ", length(pred_all),
           " predictions; expected ", length(loco_ids), ".", call. = FALSE)
    }
    pred <- rep(NA_real_, length(sample_ids))
    ok <- !is.na(pos)
    pred[ok] <- pred_all[pos[ok]]
    pred_by_chr[[key]] <- pred
    if (identical(key, target_chr)) target <- pred
  }
  if (!length(pred_by_chr)) stop("No chromosome predictions in LOCO file: ", path,
                                 call. = FALSE)
  if (is.null(target)) stop("LOCO file has no chromosome ", target_chr, ".", call. = FALSE)
  list(pred_by_chr = pred_by_chr, target = target, loco_valid = !is.na(pos))
}

use_subtract <- identical(loco_mode, "subtract")
if (!use_subtract) {
  message("--loco-mode covariate: W baseline is the UNADJUSTED phenotype, ",
          "matching Stage 1, which puts L(chr) in the basis rather than in Y.")
  loco <- list(loco_valid = rep(TRUE, nrow(ph)),
               target = rep(0, nrow(ph)),
               pred_by_chr = list(rep(0, nrow(ph))))
} else {
loco_path <- resolve_loco_file(loco_file, trait)
loco <- read_loco(loco_path, ph$FID, ph$IID, chr)
if (!any(loco$loco_valid)) stop("No phenotype rows match LOCO samples.", call. = FALSE)
ph <- ph[loco$loco_valid]
pred_target <- loco$target[loco$loco_valid]
pred_by_chr <- lapply(loco$pred_by_chr, function(x) x[loco$loco_valid])
}

## Match Stage 1's regression LOCO scale: pooled OLS of the original phenotype
## on all chromosome predictions, before the BGEN mask is applied.
resolve_scale <- function(mode, y, pred_by_chr) {
  if (identical(mode, "none")) return(1)
  if (identical(mode, "sd")) return(stats::sd(y))
  if (!identical(mode, "regression")) {
    value <- suppressWarnings(as.numeric(mode))
    if (!is.finite(value)) stop("--loco-scale must be regression, sd, none, or numeric.",
                               call. = FALSE)
    return(value)
  }
  num <- 0
  den <- 0
  for (pred in pred_by_chr) {
    ok <- is.finite(y) & is.finite(pred)
    if (!any(ok)) next
    yc <- y[ok] - mean(y[ok])
    pc <- pred[ok] - mean(pred[ok])
    num <- num + sum(pc * yc)
    den <- den + sum(pc * pc)
  }
  if (!is.finite(den) || den <= 0) stop("LOCO predictions have zero variance.", call. = FALSE)
  num / den
}

scale_factor <- if (use_subtract) resolve_scale(loco_scale_mode, as.numeric(ph$y), pred_by_chr) else 0
if (!is.finite(scale_factor)) stop("Non-finite LOCO scale factor.", call. = FALSE)

## The BGEN mask is applied after the LOCO scale, exactly as Stage 1 does.
ph <- ph[bgen_match[loco$loco_valid]]
pred_target <- pred_target[bgen_match[loco$loco_valid]]
if (!nrow(ph)) stop("No rows remain after LOCO and BGEN sample matching.", call. = FALSE)
y_adj <- as.numeric(ph$y) - scale_factor * pred_target

## Stage 1 constructs the RIF on the LOCO-adjusted phenotype and only then
## residualizes each influence-function column on X.  Therefore the baseline
## distribution for W is y_adj itself on the same complete, LOCO/BGEN-matched
## rows; residualizing y here would define a different estimand.  Covariates
## remain part of the complete-case mask and are recorded in provenance.
y_base <- as.numeric(y_adj)
if (any(!is.finite(y_base))) stop("Non-finite matched baseline values.", call. = FALSE)
if (is.finite(expected_n) && length(y_base) != expected_n) {
  stop("Matched N=", length(y_base), " does not match --expected-n=", expected_n,
       call. = FALSE)
}
if (is.finite(expected_mean) && !isTRUE(all.equal(mean(y_base), expected_mean,
                                                  tolerance = 1e-8))) {
  stop("Matched baseline mean=", format(mean(y_base), digits = 16),
       " does not match --expected-mean=", format(expected_mean, digits = 16),
       call. = FALSE)
}
if (is.finite(expected_sd) && !isTRUE(all.equal(stats::sd(y_base), expected_sd,
                                                tolerance = 1e-8))) {
  stop("Matched baseline sd=", format(stats::sd(y_base), digits = 16),
       " does not match --expected-sd=", format(expected_sd, digits = 16),
       call. = FALSE)
}
if (is.finite(expected_scale) && !isTRUE(all.equal(scale_factor, expected_scale,
                                                   tolerance = 1e-8))) {
  stop("LOCO scale=", format(scale_factor, digits = 16),
       " does not match --expected-scale=", format(expected_scale, digits = 16),
       call. = FALSE)
}
q_tau <- as.numeric(stats::quantile(y_base, probs = taus, type = 8, names = FALSE))
if (any(!is.finite(q_tau))) stop("Non-finite baseline quantiles.", call. = FALSE)

## Stage 1 records the density kernel in its covariance metadata. Read it
## before rebuilding f(q), so an optional fourth-order Stage 1 run cannot
## silently receive Gaussian weights. Legacy sidecars predate this field and
## are interpreted as the historical Gaussian estimator.
cov_meta_path <- sub("\\.stage1\\.tsv\\.gz$", ".cov.meta.json", stage1_file)
cov_meta <- read_json(cov_meta_path)
analysis <- if (is.list(cov_meta)) cov_meta$analysis else NULL
density_method <- if (!is.null(analysis) && !is.null(analysis$density_method) &&
                      length(analysis$density_method) &&
                      !is.na(analysis$density_method[[1L]])) {
  as.character(analysis$density_method[[1L]])
} else {
  "silverman_gaussian"
}
if (!density_method %in% c("silverman_gaussian", "fourth_order_gaussian")) {
  stop("Unsupported Stage 1 density method in covariance metadata: ", density_method,
       call. = FALSE)
}

## Stage 1's RIF divides by this exact kernel estimate evaluated at empirical
## type-8 quantiles. Keep this calculation in R byte-for-byte at the level of
## the Python Stage 1 formula: population SD, Silverman's 1.06 rule, and the
## same denominator floor. The fourth-order method uses
## K(u)=(1.5-0.5*u^2)*phi(u) with the same bandwidth.
stage1_kde_at <- function(y, x, density_floor = 1e-8,
                          density_method = "silverman_gaussian") {
  n <- length(y)
  if (n < 2L) stop("Too few baseline values for Stage 1 KDE.", call. = FALSE)
  if (!density_method %in% c("silverman_gaussian", "fourth_order_gaussian")) {
    stop("Unsupported Stage 1 density method: ", density_method, call. = FALSE)
  }
  sd_population <- sqrt(sum((y - mean(y))^2) / n)
  h <- 1.06 * sd_population * n^(-0.2)
  if (!is.finite(h) || h <= 0) stop("Invalid Stage 1 KDE bandwidth.", call. = FALSE)
  f <- vapply(as.numeric(x), function(q) {
    u <- (y - q) / h
    ker <- exp(-0.5 * u^2)
    if (identical(density_method, "fourth_order_gaussian")) {
      ker <- ker * (1.5 - 0.5 * u^2)
    }
    mean(ker) / (h * sqrt(2 * pi))
  }, numeric(1))
  pmax(f, density_floor)
}
f_stage1 <- stage1_kde_at(y_base, q_tau, density_method = density_method)

## Build shifted-lognormal W using the package's generic derivative API.  The
## parameter names are retained from the paper outputs: additive is the shift
## parameter and multiplicative is the log-scale location parameter.
build_lognormal <- function(y, q_tau, taus) {
  shift0 <- as.numeric(stats::quantile(y, probs = 0.01, type = 8, names = FALSE)) - 1e-3
  y_shift <- y - shift0
  y_shift <- y_shift[y_shift > 0]
  if (length(y_shift) < 10L) stop("Too few positive baseline values after shift.",
                                 call. = FALSE)
  params0 <- list(
    shift = shift0,
    meanlog = mean(log(y_shift)),
    sdlog = stats::sd(log(y_shift))
  )
  if (!is.finite(params0$sdlog) || params0$sdlog <= 0) {
    stop("Invalid shifted-lognormal baseline sdlog.", call. = FALSE)
  }
  dist_pdf <- function(x, p) {
    out <- numeric(length(x)); ok <- x > p$shift
    out[ok] <- dlnorm(x[ok] - p$shift, meanlog = p$meanlog, sdlog = p$sdlog)
    out
  }
  grad_shift <- function(x, p) {
    out <- numeric(length(x)); ok <- x > p$shift
    out[ok] <- -dlnorm(x[ok] - p$shift, meanlog = p$meanlog, sdlog = p$sdlog)
    out
  }
  grad_meanlog <- function(x, p) {
    out <- numeric(length(x)); ok <- x > p$shift
    yy <- x[ok] - p$shift
    z <- (log(yy) - p$meanlog) / p$sdlog
    out[ok] <- -dnorm(z) / p$sdlog
    out
  }
  W <- fungwasStage2::make_weights_generic(
    taus = taus, q_tau = q_tau,
    dist_cdf = function(x, p) {
      out <- numeric(length(x)); ok <- x > p$shift
      out[ok] <- plnorm(x[ok] - p$shift, meanlog = p$meanlog, sdlog = p$sdlog)
      out
    },
    dist_pdf = dist_pdf,
    params = params0,
    grad_funcs = list(additive = grad_shift, multiplicative = grad_meanlog)
  )
  W_add <- W[, "additive", drop = FALSE]
  W_shift <- W[, "multiplicative", drop = FALSE]
  W_logscale <- cbind(
    meanlog = as.numeric(W_shift[, 1L]),
    sigma2 = as.numeric(W_shift[, 1L]) * qnorm(taus) / (2 * params0$sdlog)
  )
  out <- list(
    additive = list(W = W_add, taus = taus, q_tau = q_tau, params0 = params0,
                    type = "matched_shifted_lognormal"),
    shifted = list(W = W_shift, taus = taus, q_tau = q_tau, params0 = params0,
                   type = "matched_shifted_lognormal"),
    both = list(W = W, taus = taus, q_tau = q_tau, params0 = params0,
                type = "matched_shifted_lognormal")
  )
  if (identical(trait, "bmi")) {
    out$logscale_vqtl <- list(W = W_logscale, taus = taus, q_tau = q_tau,
                              params0 = params0, type = "matched_logscale_vqtl")
  }
  out
}

## Fit the same two-component skew-normal mixture used by the historical paper
## builder, but on the matched LOCO-adjusted baseline.  Components are sorted by
## location so the parameter columns have stable meaning across chromosomes.
build_skew_mixture <- function(y, q_tau, taus, seed, f_stage1, density_method) {
  if (!requireNamespace("mixsmsn", quietly = TRUE) || !requireNamespace("sn", quietly = TRUE)) {
    stop("The skew/mixture weight builder requires mixsmsn and sn in freshRlib.",
         call. = FALSE)
  }
  set.seed(seed)
  fit <- mixsmsn::smsn.mix(
    y = y, g = 2, family = "Skew.normal", nu = 1,
    get.init = TRUE, iter.max = 500, error = 1e-4, calc.im = FALSE
  )
  ord <- order(as.numeric(fit$mu))
  pii <- as.numeric(fit$pii)[ord]
  mu <- as.numeric(fit$mu)[ord]
  sigma2 <- as.numeric(fit$sigma2)[ord]
  shape <- as.numeric(fit$shape)[ord]
  if (length(pii) != 2L || any(!is.finite(c(pii, mu, sigma2, shape))) ||
      any(pii <= 0) || any(sigma2 <= 0)) {
    stop("Skew-normal mixture fit returned invalid component parameters.", call. = FALSE)
  }
  p <- list(
    p1 = pii[1L], xi1 = mu[1L], omega1 = sqrt(sigma2[1L]), alpha1 = shape[1L],
    xi2 = mu[2L], omega2 = sqrt(sigma2[2L]), alpha2 = shape[2L]
  )
  cdf <- function(x, z) {
    z$p1 * sn::psn(x, xi = z$xi1, omega = z$omega1, alpha = z$alpha1) +
      (1 - z$p1) * sn::psn(x, xi = z$xi2, omega = z$omega2, alpha = z$alpha2)
  }
  pdf <- function(x, z) {
    z$p1 * sn::dsn(x, xi = z$xi1, omega = z$omega1, alpha = z$alpha1) +
      (1 - z$p1) * sn::dsn(x, xi = z$xi2, omega = z$omega2, alpha = z$alpha2)
  }
  ## Use the density that scaled the Stage 1 RIF for testosterone.  The
  ## historical mixture path for the other mixture traits remains unchanged
  ## until it gets a separately registered validation.  The fallback is only
  ## for defensive callers evaluating the closure away from q_tau; the normal
  ## path receives q_tau exactly and reuses the already computed vector.
  use_stage1_density <- grepl("^testosterone", trait)
  density_for_weights <- function(x, z) {
    if (use_stage1_density && length(x) == length(q_tau) &&
        isTRUE(all.equal(as.numeric(x), q_tau))) {
      return(f_stage1)
    }
    if (use_stage1_density) {
      stage1_kde_at(y, x, density_method = density_method)
    } else {
      pdf(x, z)
    }
  }
  grad_p <- fungwasStage2::make_fd_grad("p1", cdf)
  grad_gamma <- function(x, z) grad_p(x, z) * z$p1 * (1 - z$p1)
  W <- fungwasStage2::make_weights_generic(
    taus = taus, q_tau = q_tau, dist_cdf = cdf,
    dist_pdf = density_for_weights, params = p,
    grad_funcs = list(
      gamma = grad_gamma,
      beta_left = fungwasStage2::make_fd_grad("xi1", cdf),
      beta_right = fungwasStage2::make_fd_grad("xi2", cdf)
    ),
    # dist_pdf above already selects the density; keep the builder on it.
    density = "model"
  )
  f_model <- as.numeric(pdf(q_tau, p))
  list(W = W, taus = taus, q_tau = q_tau,
       density_stage1 = f_stage1,
       density_model = f_model,
       density_model_over_stage1 = f_model / f_stage1,
       density_method = density_method,
       components = data.table(component = c("lower", "upper"),
                               weight = pii, xi = mu,
                               omega = sqrt(sigma2), alpha = shape),
       type = "matched_skew_normal_mixture")
}

make_outputs <- function() {
  if (identical(family, "shifted_lognormal")) {
    return(build_lognormal(y_base, q_tau, taus))
  }
  mix <- build_skew_mixture(y_base, q_tau, taus, seed, f_stage1, density_method)
  if (identical(family, "skew")) return(list(skew = mix))
  list(mixture = mix)
}

## Validate against the fresh Stage 1 provenance when available. This catches
## accidental use of the raw-BMI baseline for log-BMI before any Stage 2 call.
if (!is.null(analysis)) {
  if (!is.null(analysis$phenotype) && !identical(as.character(analysis$phenotype), trait)) {
    stop("Stage 1 metadata phenotype is ", analysis$phenotype, "; expected ", trait,
         call. = FALSE)
  }
  if (!is.null(analysis$n_samples) && as.integer(analysis$n_samples) != nrow(ph)) {
    stop("Matched baseline N=", nrow(ph), " but Stage 1 metadata records N=",
         analysis$n_samples, "; refusing mismatched W.", call. = FALSE)
  }
  if (!is.null(analysis$covariates) &&
      !identical(as.character(analysis$covariates), as.character(covar_cols))) {
    stop("Covariate columns differ from Stage 1 metadata; refusing mismatched W.",
         call. = FALSE)
  }
  if (!is.null(analysis$loco_scale_mode) &&
      !identical(as.character(analysis$loco_scale_mode), as.character(loco_scale_mode))) {
    stop("LOCO scale mode differs from Stage 1 metadata; refusing mismatched W.",
         call. = FALSE)
  }
  if (!is.null(analysis$loco_scale_factor) &&
      !isTRUE(all.equal(as.numeric(analysis$loco_scale_factor), scale_factor,
                        tolerance = 1e-8))) {
    stop("LOCO scale factor differs from Stage 1 metadata; refusing mismatched W.",
         call. = FALSE)
  }
  if (!is.null(analysis$density_method) &&
      !identical(as.character(analysis$density_method), density_method)) {
    stop("Density method differs from Stage 1 metadata; refusing mismatched W.",
         call. = FALSE)
  }
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
manifest_path <- file.path(out_dir, "weights_manifest.tsv")
if (file.exists(manifest_path)) {
  stop("Refusing to overwrite existing weight manifest: ", manifest_path,
       ". Use a new exclusive run root.", call. = FALSE)
}
outputs <- make_outputs()
provenance <- list(
  phenotype = trait,
  chromosome = chr,
  family = family,
  n_samples = nrow(ph),
  taus = taus,
  tau_labels = tau_labels,
  tau_grid_source = if (identical(taus, taus_full)) {
    "canonical_full_precision_linspace_recovered_from_three_decimal_labels"
  } else {
    "serialized_three_decimal_labels"
  },
  covariates = covar_cols,
  loco_mode = loco_mode,
  loco_scale_mode = loco_scale_mode,
  loco_scale_factor = unname(scale_factor),
  baseline = "LOCO-adjusted phenotype, BGEN-matched (Stage 1 RIF then residualizes X)",
  density_stage1 = paste0(density_method,
                          " at empirical type-8 quantiles, population SD, floor 1e-8"),
  density_method = density_method,
  density_weight_mode = if (grepl("^testosterone", trait)) {
    "matched_to_stage1_kde"
  } else {
    "family_specific"
  },
  stage1_file = normalizePath(stage1_file, mustWork = TRUE),
  seed = seed
)

for (nm in names(outputs)) {
  obj <- outputs[[nm]]
  obj$provenance <- provenance
  path <- file.path(out_dir, paste0("weights_", nm, ".rds"))
  if (file.exists(path)) stop("Refusing to overwrite existing weight file: ", path,
                              ". Use a new exclusive run root.", call. = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  saveRDS(obj, tmp)
  if (!file.rename(tmp, path)) stop("Could not move temporary weight file to ", path,
                                   call. = FALSE)
}

summary <- data.table(
  trait = trait, chromosome = chr, family = family, n_samples = nrow(ph),
  n_taus = length(taus), loco_mode = loco_mode, loco_scale_mode = loco_scale_mode,
  loco_scale_factor = scale_factor,
  baseline_mean = mean(y_base), baseline_sd = stats::sd(y_base),
  stage1_file = normalizePath(stage1_file, mustWork = TRUE),
  output_dir = normalizePath(out_dir, mustWork = TRUE)
)
fwrite(summary, manifest_path, sep = "\t")
message("Wrote matched aggregate weights for ", trait, " chr", chr,
        " (N=", nrow(ph), ", ", length(taus), " taus) to ", out_dir)
