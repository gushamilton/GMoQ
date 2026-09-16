#!/usr/bin/env Rscript

## Blind LOCO diagnostic for G-MoQ.
##
## Detects a REGENIE step-1 whole-genome predictor that has learned the
## phenotype's latent mixture structure instead of polygenic mean signal --
## WITHOUT reading any membership label.
##
## The only structural input is the two-component fit G-MoQ's own Stage 2
## weight builder already produced (`components` inside weights_<model>.rds).
## From it we derive posterior membership weights w_k(y) per sample and ask two
## questions:
##
##   1. How much of the LOCO predictor is explained by the posterior weights?
##      (`r2_loco_on_posterior`).  High means the predictor is a membership
##      proxy, so subtracting it erodes exactly the contrast Stage 2 estimates.
##
##   2. How much of the phenotype does the predictor explain WITHIN a
##      posterior-weighted component (`r2_within_*`)?  Near zero means it is
##      buying no polygenic adjustment to pay for that damage.
##
## A predictor that is doing its job has low (1) and non-trivial (2).  A
## membership proxy has the reverse.  We also report the within-component
## variance ratio var(y - c*loco)/var(y), the quantity that actually determines
## whether the subtraction helps or hurts Stage 1.
##
## Everything printed is aggregate. No sample identifiers or per-person values
## are emitted.
##
## Usage:
##   Rscript scripts/loco_blind_diagnostic.R \
##     --weights <weights_skew.rds> --pheno <phenotypes.txt> \
##     --pheno-col <column> --keep <keep.txt> \
##     --loco <file.loco|*_pred.list> --chr 22 \
##     [--loco-label <pred.list label>] [--r-lib <dir>] [--out <tsv>] \
##     [--validate-col <label column>]
##
## --validate-col is OPTIONAL and is for checking that the blind diagnostic
## agrees with a known label on a trait where one happens to exist. It plays no
## part in the diagnostic itself and must never be required.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) args[[i + 1L]] else default
}
r_lib <- get_arg("--r-lib")
if (!is.null(r_lib)) .libPaths(c(r_lib, .libPaths()))
suppressPackageStartupMessages({ library(data.table) })
if (!requireNamespace("sn", quietly = TRUE)) {
  stop("The blind diagnostic needs the `sn` package to evaluate the fitted ",
       "skew-normal mixture.", call. = FALSE)
}

weights_file <- get_arg("--weights")
pheno_file   <- get_arg("--pheno")
pheno_col    <- get_arg("--pheno-col")
keep_file    <- get_arg("--keep")
loco_file    <- get_arg("--loco")
chr_raw      <- get_arg("--chr", "22")
loco_label   <- get_arg("--loco-label")
out_file     <- get_arg("--out")
validate_col <- get_arg("--validate-col")
tag          <- get_arg("--tag", "loco")

if (any(vapply(list(weights_file, pheno_file, pheno_col, keep_file, loco_file),
               is.null, logical(1)))) {
  stop("Usage: --weights --pheno --pheno-col --keep --loco [--chr 22]", call. = FALSE)
}
chr <- as.character(as.integer(chr_raw))

## ---------------------------------------------------------------------------
## Fitted mixture -> posterior membership weights.  Nothing else is used.
## ---------------------------------------------------------------------------
wobj <- readRDS(weights_file)
if (is.null(wobj$components)) {
  stop("Weight file carries no fitted mixture components; the blind diagnostic ",
       "needs a Stage 2 family that estimates membership: ", weights_file,
       call. = FALSE)
}
comp <- as.data.table(wobj$components)
if (nrow(comp) != 2L) stop("Expected a two-component fit.", call. = FALSE)
pii <- as.numeric(comp$weight); xi <- as.numeric(comp$xi)
omega <- as.numeric(comp$omega); alpha <- as.numeric(comp$alpha)
delta <- alpha / sqrt(1 + alpha^2)
mu_k <- xi + omega * delta * sqrt(2 / pi)
sd_k <- omega * sqrt(1 - 2 * delta^2 / pi)

posterior_w1 <- function(y) {
  d1 <- pii[1] * sn::dsn(y, xi = xi[1], omega = omega[1], alpha = alpha[1])
  d2 <- pii[2] * sn::dsn(y, xi = xi[2], omega = omega[2], alpha = alpha[2])
  den <- d1 + d2
  bad <- !is.finite(den) | den <= 0
  if (any(bad)) {
    d1[bad] <- as.numeric(abs(y[bad] - mu_k[1]) / sd_k[1] <
                            abs(y[bad] - mu_k[2]) / sd_k[2])
    den[bad] <- 1
  }
  d1 / den
}

## ---------------------------------------------------------------------------
## Phenotype
## ---------------------------------------------------------------------------
sel <- c("FID", "IID", pheno_col, if (!is.null(validate_col)) validate_col)
ph <- fread(pheno_file, select = sel, showProgress = FALSE)
setnames(ph, pheno_col, "y")
ph[, FID := as.character(FID)][, IID := as.character(IID)]
keep_ids <- unique(as.character(fread(keep_file, header = FALSE,
                                      colClasses = "character",
                                      showProgress = FALSE)[[1L]]))
ph <- ph[FID %in% keep_ids][is.finite(as.numeric(y))]

## ---------------------------------------------------------------------------
## LOCO predictor for the requested chromosome
## ---------------------------------------------------------------------------
resolve_loco_file <- function(path, label) {
  if (grepl("\\.loco$", path, ignore.case = TRUE)) return(path)
  lines <- readLines(path, warn = FALSE)
  e <- strsplit(trimws(lines[nzchar(trimws(lines))]), "[[:space:]]+")
  paths <- vapply(e, function(x) x[[length(x)]], character(1))
  labels <- vapply(e, `[[`, character(1), 1L)
  if (length(unique(paths)) == 1L) return(paths[[1L]])
  if (is.null(label)) stop("pred.list has several entries; pass --loco-label.",
                           call. = FALSE)
  cand <- unique(paths[labels == label])
  if (length(cand) != 1L) stop("Could not resolve one .loco for label ", label,
                               call. = FALSE)
  cand[[1L]]
}
normalize_chr <- function(x) {
  s <- sub("^(chr|CHR)", "", trimws(as.character(x)))
  i <- grepl("^[0-9]+$", s); s[i] <- as.character(as.integer(s[i])); toupper(s)
}

read_loco_chr <- function(path, fid, iid, target_chr) {
  con <- file(path, open = "rt"); on.exit(close(con), add = TRUE)
  header <- scan(text = readLines(con, n = 1L, warn = FALSE),
                 what = character(), quiet = TRUE)
  loco_ids <- header[-1L]
  pos <- match(fid, loco_ids)
  for (alias in list(iid, paste(fid, iid, sep = "_"),
                     paste(fid, iid, sep = " "), paste(fid, iid, sep = ":"))) {
    m <- is.na(pos); if (any(m)) pos[m] <- match(alias[m], loco_ids)
  }
  target <- NULL; by_chr <- list()
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) break
    line <- trimws(line); if (!nzchar(line)) next
    sp <- regexpr("[[:space:]]", line, perl = TRUE)[1L]; if (sp < 1L) next
    key <- normalize_chr(substr(line, 1L, sp - 1L))
    vals <- scan(text = sub("^[^[:space:]]+[[:space:]]+", "", line, perl = TRUE),
                 what = numeric(), quiet = TRUE)
    pred <- rep(NA_real_, length(fid)); ok <- !is.na(pos)
    pred[ok] <- vals[pos[ok]]
    by_chr[[key]] <- pred
    if (identical(key, target_chr)) target <- pred
  }
  if (is.null(target)) stop("LOCO file has no chromosome ", target_chr, call. = FALSE)
  list(target = target, by_chr = by_chr, valid = !is.na(pos))
}

lc <- read_loco_chr(resolve_loco_file(loco_file, loco_label), ph$FID, ph$IID, chr)
ph <- ph[lc$valid]
loco <- lc$target[lc$valid]
by_chr <- lapply(lc$by_chr, function(x) x[lc$valid])
y <- as.numeric(ph$y)

## The pipeline's own LOCO scale: pooled OLS of y on all chromosome predictions.
num <- 0; den <- 0
for (p in by_chr) {
  ok <- is.finite(y) & is.finite(p)
  yc <- y[ok] - mean(y[ok]); pc <- p[ok] - mean(p[ok])
  num <- num + sum(pc * yc); den <- den + sum(pc * pc)
}
scale_factor <- num / den

w1 <- posterior_w1(y); w2 <- 1 - w1

## ---------------------------------------------------------------------------
## Diagnostic 1: how much of the LOCO predictor IS the posterior?
## ---------------------------------------------------------------------------
r2_loco_on_posterior <- summary(stats::lm(loco ~ w1))$r.squared

## ---------------------------------------------------------------------------
## Diagnostic 2: what does it explain WITHIN a posterior-weighted component?
## ---------------------------------------------------------------------------
wstat <- function(w) {
  sw <- sum(w)
  my <- sum(w * y) / sw; mp <- sum(w * loco) / sw
  vy <- sum(w * (y - my)^2) / sw; vp <- sum(w * (loco - mp)^2) / sw
  cyp <- sum(w * (y - my) * (loco - mp)) / sw
  r2 <- if (vy > 0 && vp > 0) (cyp^2) / (vy * vp) else NA_real_
  r <- y - scale_factor * loco
  mr <- sum(w * r) / sw
  vr <- sum(w * (r - mr)^2) / sw
  list(r2 = r2, var_ratio = vr / vy, eff_n = sw, beta = cyp / vp)
}
s1 <- wstat(w1); s2 <- wstat(w2)

overall_r <- y - scale_factor * loco
var_ratio_overall <- var(overall_r) / var(y)

res <- data.table(
  tag = tag,
  chromosome = chr,
  n = length(y),
  loco_scale_factor = scale_factor,
  cor_y_loco = stats::cor(y, loco),
  ## (1) membership-proxy score
  r2_loco_on_posterior = r2_loco_on_posterior,
  ## (2) within-component usefulness
  r2_within_lower = s1$r2,
  r2_within_upper = s2$r2,
  ## what the subtraction actually does to Stage 1's working variance
  var_ratio_within_lower = s1$var_ratio,
  var_ratio_within_upper = s2$var_ratio,
  var_ratio_overall = var_ratio_overall,
  eff_n_lower = s1$eff_n,
  eff_n_upper = s2$eff_n,
  posterior_pi_lower = mean(w1),
  frac_ambiguous_w_05_95 = mean(w1 > 0.05 & w1 < 0.95)
)

## Verdict.  A predictor is flagged when it explains far more of the posterior
## than of the phenotype within a component, or when subtracting it inflates
## within-component variance.
worst_within <- max(s1$r2, s2$r2, na.rm = TRUE)
flag_proxy <- r2_loco_on_posterior > 10 * worst_within
flag_inflate <- max(s1$var_ratio, s2$var_ratio, na.rm = TRUE) > 1
res[, membership_proxy_ratio := r2_loco_on_posterior / worst_within]
res[, flag := flag_proxy || flag_inflate]

if (!is.null(validate_col)) {
  lab <- as.numeric(ph[[validate_col]])
  res[, validate_cor_label_loco := stats::cor(lab, loco)]
  res[, validate_cor_label_posterior := stats::cor(lab, w1)]
}

print(t(res))
if (!is.null(out_file)) {
  fwrite(res, out_file, sep = "\t")
  message("wrote ", out_file)
}
message(if (res$flag) {
  "VERDICT: FLAGGED -- this LOCO predictor is a membership proxy and/or inflates within-component variance. Do not subtract it."
} else {
  "VERDICT: clean -- the LOCO predictor is not a membership proxy and does not inflate within-component variance."
})
