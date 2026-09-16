#!/usr/bin/env Rscript

## Two-pass blind LOCO -- Pass 1 residualisation.
##
## Reads the Pass-1 (no-LOCO) G-MoQ Stage 2 skew-normal mixture fit, computes
## posterior membership weights for every sample from that fit alone, and writes
## a residualised phenotype in which the ESTIMATED mixture structure has been
## removed.  No membership label is read at any point: the only inputs are the
## phenotype column and the aggregate mixture parameters G-MoQ itself estimated.
##
## Writes participant-level phenotype files that stay on BluePebble.  All
## console output is aggregate only.

args <- commandArgs(trailingOnly = TRUE)
if (!grepl("^bp1-compute", Sys.info()[["nodename"]])) {
  stop("Run on a BP compute node (bp1-compute*).", call. = FALSE)
}
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) args[[i + 1L]] else default
}
r_lib <- get_arg("--r-lib")
if (!is.null(r_lib)) .libPaths(c(r_lib, .libPaths()))
suppressPackageStartupMessages({ library(data.table); library(sn) })

weights_file <- get_arg("--weights")
pheno_file   <- get_arg("--pheno")
pheno_col    <- get_arg("--pheno-col", "testosterone")
keep_file    <- get_arg("--keep")
out_pheno    <- get_arg("--out-pheno")
out_summary  <- get_arg("--out-summary")

stopifnot(!is.null(weights_file), !is.null(pheno_file), !is.null(keep_file),
          !is.null(out_pheno), !is.null(out_summary))

wobj <- readRDS(weights_file)
if (is.null(wobj$components)) {
  stop("Weight file carries no fitted mixture components: ", weights_file, call. = FALSE)
}
comp <- as.data.table(wobj$components)
if (nrow(comp) != 2L) stop("Expected a two-component mixture fit.", call. = FALSE)
message("Pass-1 mixture fit (from ", weights_file, "):")
print(comp)
message("Pass-1 LOCO scale factor: ", wobj$provenance$loco_scale_factor,
        "  (0 means the pass-1 fit used no LOCO, as required)")
if (!isTRUE(all.equal(as.numeric(wobj$provenance$loco_scale_factor), 0))) {
  stop("Pass-1 weights must come from a NO-LOCO run (loco_scale_factor 0).",
       call. = FALSE)
}

pii   <- as.numeric(comp$weight)
xi    <- as.numeric(comp$xi)
omega <- as.numeric(comp$omega)
alpha <- as.numeric(comp$alpha)

## Skew-normal component moments.
delta  <- alpha / sqrt(1 + alpha^2)
mu_k   <- xi + omega * delta * sqrt(2 / pi)
sd_k   <- omega * sqrt(1 - 2 * delta^2 / pi)
if (any(!is.finite(mu_k)) || any(!is.finite(sd_k)) || any(sd_k <= 0)) {
  stop("Degenerate component moments.", call. = FALSE)
}
message(sprintf("component means  : %.6f %.6f", mu_k[1], mu_k[2]))
message(sprintf("component sds    : %.6f %.6f", sd_k[1], sd_k[2]))
message(sprintf("separation S = |m1-m2| / sqrt(mean var) : %.4f",
                abs(diff(mu_k)) / sqrt(mean(sd_k^2))))

ph <- fread(pheno_file, select = c("FID", "IID", pheno_col), showProgress = FALSE)
setnames(ph, pheno_col, "y")
ph[, FID := as.character(FID)][, IID := as.character(IID)]
keep_ids <- unique(as.character(fread(keep_file, header = FALSE,
                                      colClasses = "character",
                                      showProgress = FALSE)[[1L]]))
ph <- ph[FID %in% keep_ids]
ph <- ph[is.finite(as.numeric(y))]
y <- as.numeric(ph$y)
message("N with phenotype in keep set: ", length(y))

## Posterior membership weights from the estimated mixture alone.
d1 <- pii[1] * sn::dsn(y, xi = xi[1], omega = omega[1], alpha = alpha[1])
d2 <- pii[2] * sn::dsn(y, xi = xi[2], omega = omega[2], alpha = alpha[2])
den <- d1 + d2
bad <- !is.finite(den) | den <= 0
if (any(bad)) {
  ## Fall back to the nearer component in the far tails where both densities
  ## underflow; affects a handful of extreme values at most.
  message("posterior underflow rows (assigned to nearer component): ", sum(bad))
  d1[bad] <- as.numeric(abs(y[bad] - mu_k[1]) / sd_k[1] <
                          abs(y[bad] - mu_k[2]) / sd_k[2])
  d2[bad] <- 1 - d1[bad]
  den[bad] <- 1
}
w1 <- d1 / den
w2 <- 1 - w1

message(sprintf("posterior mean w1 = %.6f (fitted pi1 = %.6f)", mean(w1), pii[1]))
message(sprintf("fraction with w1 in [0.05,0.95] (ambiguous) : %.5f",
                mean(w1 > 0.05 & w1 < 0.95)))
message(sprintf("fraction with w1 in [0.20,0.80]            : %.5f",
                mean(w1 > 0.20 & w1 < 0.80)))

## --- Residualisation -------------------------------------------------------
## resid_mean : remove the posterior-weighted component LOCATION only.
##              y - sum_k w_k(y) m_k.  Kills the between-component mean
##              contrast, which is the only thing a ridge MEAN model can learn,
##              while preserving every within-component deviation intact.
## resid_z    : remove posterior-weighted component LOCATION AND SCALE.
##              sum_k w_k(y) (y - m_k)/s_k.  Additionally equalises the two
##              components' variances so neither dominates the ridge objective.
mu_hat <- w1 * mu_k[1] + w2 * mu_k[2]
resid_mean <- y - mu_hat
resid_z <- w1 * (y - mu_k[1]) / sd_k[1] + w2 * (y - mu_k[2]) / sd_k[2]

message(sprintf("var(y)=%.6f  var(resid_mean)=%.6f  var(resid_z)=%.6f",
                var(y), var(resid_mean), var(resid_z)))
message(sprintf("cor(y, resid_mean)=%.6f  cor(y, resid_z)=%.6f",
                cor(y, resid_mean), cor(y, resid_z)))
## Within hard-assigned components, the transform should be close to affine:
hard <- w1 > 0.5
for (lab in c(TRUE, FALSE)) {
  i <- hard == lab
  message(sprintf("  hard-component %s (n=%d): cor(y,resid_mean)=%.6f cor(y,resid_z)=%.6f",
                  if (lab) "lower" else "upper", sum(i),
                  cor(y[i], resid_mean[i]), cor(y[i], resid_z[i])))
}

out <- data.table(FID = ph$FID, IID = ph$IID,
                  resid_mean = resid_mean, resid_z = resid_z)
fwrite(out, out_pheno, sep = " ", quote = FALSE, na = "NA")
message("wrote ", out_pheno, " (", nrow(out), " rows)")

## Aggregate-only summary, safe to copy off the cluster.
summ <- data.table(
  component = c("lower", "upper"),
  pi = pii, xi = xi, omega = omega, alpha = alpha,
  mean = mu_k, sd = sd_k
)
fwrite(summ, out_summary, sep = "\t")
message("wrote ", out_summary)
