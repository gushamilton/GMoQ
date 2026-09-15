#!/usr/bin/env Rscript
# Stage 2 side of the compact validity pilot.
#
# Fits the installed package on the binary output from pilot_stage1.py.  The
# two covariance arms are pooled+OLS and jackknife+OLS, paired on every
# dataset.  We also fit an oracle baseline W and an estimated W, so parameter
# coverage is not silently attributed to a model basis fixed by the simulator.
suppressPackageStartupMessages({
  .self <- local({
    a <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (!length(f)) stop("run as Rscript --vanilla tests/validity_pilot/pilot_stage2.R")
    normalizePath(f[1L])
  })
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("usage: pilot_stage2.R <pilot_dir> <output_csv>")
pilot_dir <- normalizePath(args[1L])
out_csv <- args[2L]
test_lib <- Sys.getenv("GMOQ_TEST_LIB", "")
dep_lib <- Sys.getenv("GMOQ_DEP_LIB", "")
.libPaths(unique(c(test_lib, dep_lib[nzchar(dep_lib)], .libPaths())))
suppressPackageStartupMessages(library(fungwasStage2))
if (!("fungwasStage2" %in% names(getLoadedDLLs()))) stop("compiled Stage 2 DLL is not loaded")

suite <- normalizePath(file.path(dirname(.self), ".."))
internal <- function(n) get(n, envir = asNamespace("fungwasStage2"))
pooled_cov_data <- internal(".pooled_cov_data")

truth_dir <- file.path(pilot_dir, "truth")
stage_dir <- file.path(pilot_dir, "stage1")
truth_files <- list.files(truth_dir, pattern = "\\.json$", full.names = TRUE)
if (!length(truth_files)) stop("no truth files in ", truth_dir)

read_arr <- function(path, what, dims, size = 8L) {
  x <- readBin(path, what = what, size = size, n = prod(dims), endian = "little")
  if (length(x) != prod(dims)) stop("short binary file: ", path)
  array(x, dim = dims)
}

rows <- list(); datasets <- list(); ir <- 0L; id <- 0L
append_row <- function(x) { ir <<- ir + 1L; rows[[ir]] <<- x }
append_dataset <- function(x) { id <<- id + 1L; datasets[[id]] <<- x }

for (tf in truth_files) {
  tr <- jsonlite::fromJSON(tf, simplifyVector = FALSE)
  cell <- tr$cell; name <- cell$name
  for (path in if (identical(cell$family, "cov_both")) c("current", "y_first") else "current") {
    pfx <- file.path(stage_dir, paste0(name, "__", path))
    jm <- jsonlite::fromJSON(paste0(pfx, ".json"), simplifyVector = FALSE)
    R <- as.integer(jm$reps); M <- as.integer(jm$M); T <- as.integer(jm$T)
    NE <- as.integer(jm$ne); N <- as.integer(jm$N); B <- as.integer(jm$blocks)
    taus <- as.numeric(unlist(tr$taus))
    beta <- read_arr(paste0(pfx, ".beta.bin"), "double", c(T, M, R))
    se <- read_arr(paste0(pfx, ".se.bin"), "double", c(T, M, R))
    covA <- read_arr(paste0(pfx, ".cov.bin"), "double", c(NE, M, R))
    dres <- read_arr(paste0(pfx, ".dresid.bin"), "double", c(M, R))
    s0A <- read_arr(paste0(pfx, ".sigma0.bin"), "double", c(T, T, R))
    qA <- read_arr(paste0(pfx, ".q_tau.bin"), "double", c(T, R))
    ymean <- as.numeric(read_arr(paste0(pfx, ".y_mean.bin"), "double", c(R)))
    ysd <- as.numeric(read_arr(paste0(pfx, ".y_sd.bin"), "double", c(R)))
    bproj <- matrix(as.numeric(unlist(tr$projection[[path]])), nrow = M, ncol = T, byrow = TRUE)
    sci_a <- as.numeric(unlist(tr$scientific_a))
    sci_b <- as.numeric(unlist(tr$scientific_variance_coefficient))
    sci_vr <- as.numeric(unlist(tr$scientific_variance_ratio))
    taus <- as.numeric(unlist(tr$taus))
    # Both W arms use the shipped Normal mean/variance basis.  The oracle arm
    # knows the DGM baseline; the estimated arm uses each dataset's empirical
    # quantiles and SD.  This keeps beta_sigma2's scientific target (a local
    # variance coefficient) consistent across arms.
    W_oracle <- make_weights_vqtl(taus, q_tau = qnorm(taus), mu = 0, sd = 1)

    for (r in seq_len(R)) {
      Qs <- beta[, , r, drop = FALSE]; dim(Qs) <- c(T, M)
      Ss <- se[, , r, drop = FALSE]; dim(Ss) <- c(T, M)
      colnames(Qs) <- colnames(Ss) <- paste0("v", seq_len(M))
      st <- list(taus = taus, q_tau = qA[, r], Q_slope = Qs, SE_tau = Ss)
      iu <- which(upper.tri(matrix(0, T, T), diag = TRUE))
      # Python writes row-major upper triangles. The shipped reader and the
      # package's pooled helper use this same ordering.
      cv_jk <- list(cov_vec = as.numeric(covA[, , r]), offsets = seq(0, by = NE * 4, length.out = M))
      S0 <- s0A[, , r, drop = FALSE]; dim(S0) <- c(T, T)
      cv_pl <- pooled_cov_data(list(ALL = S0), rep("1", M), dres[, r], T)
      Ws <- list(oracle = W_oracle,
                 estimated = make_weights_vqtl(taus, q_tau = qA[, r], mu = ymean[r], sd = ysd[r]))
      for (src in c("jackknife", "pooled")) {
        cvd <- if (src == "jackknife") cv_jk else list(cov_vec = cv_pl$cov_vec, offsets = cv_pl$offsets)
        for (wn in names(Ws)) {
          W <- Ws[[wn]]
          fit <- tryCatch(
            suppressWarnings(suppressMessages(param_gwas(
              st, transform = "custom_W", transform_args = list(W = W),
              se_mode = "tau_cov", tau_cov_estimator = "ols", cov_data = cvd,
              q_method = "saddlepoint",
              n_blocks = if (src == "jackknife") B else NULL,
              n_eff = if (src == "jackknife") dres[, r] else NULL,
              nu = if (src == "pooled") N else NULL,
              nu_source = if (src == "pooled") "pooled_sigma0_n_samples" else NULL
            ))), error = function(e) list(error = conditionMessage(e)))
          if (!is.null(fit$error)) {
            append_dataset(data.frame(cell = name, path = path, rep = r, source = src,
                                      W = wn, fit_ok = FALSE, fit_error = fit$error,
                                      null_type1 = NA, null_any = NA, causal_cover_proj = NA,
                                      causal_cover_sci = NA, stringsAsFactors = FALSE))
            next
          }
          A <- solve(crossprod(W), t(W))
          th_proj <- A %*% t(bproj)  # K x M
          est <- fit$params; se0 <- fit$SE_params
          pp <- fit$P_params
          pn <- rownames(est)
          # This pilot has a location/scale W. Report each variant/parameter,
          # retaining both target types rather than silently choosing one.
          for (j in seq_len(M)) for (k in seq_along(pn)) {
            pnm <- pn[k]
            sci <- if (pnm == "beta_mu") sci_a[j] else sci_b[j]
            sci_var_ratio <- if (pnm == "beta_mu") NA_real_ else sci_vr[j]
            proj <- th_proj[k, j]
            zproj <- (est[k, j] - proj) / se0[k, j]
            zsci <- (est[k, j] - sci) / se0[k, j]
            append_row(data.frame(cell = name, path = path, rep = r, source = src,
                                  W = wn, variant = j, parameter = pnm,
                                  estimate = est[k, j], se = se0[k, j],
                                  p = pp[k, j], target_projection = proj,
                                  target_scientific = sci, z_projection = zproj,
                                  target_scientific_variance_ratio = sci_var_ratio,
                                  z_scientific = zsci, Q = fit$Q[j],
                                  Q_p = if (is.null(fit$Q_pval)) NA else fit$Q_pval[j],
                                  Q_ref_dist = if (is.null(fit$Q_ref_dist)) NA else fit$Q_ref_dist[j],
                                  param_ref_dist = fit$param_ref_dist,
                                  param_ref_df = fit$param_ref_df[k, j],
                                  stringsAsFactors = FALSE))
          }
          causal_idx <- if (cell$family %in% c("cov_both", "null")) integer(0) else 1L
          null_idx <- if (length(causal_idx)) setdiff(seq_len(M), causal_idx) else seq_len(M)
          murow <- which(pn == "beta_mu"); srow <- which(pn == "beta_sigma2")
          null_p <- c(pp[murow, null_idx], pp[srow, null_idx])
          null_type1 <- mean(is.finite(null_p) & null_p < 0.05)
          null_any <- as.numeric(any(is.finite(null_p) & null_p < 0.05))
          causal_cover_proj <- if (length(causal_idx)) {
            crit <- qt(0.975, df = fit$param_ref_df[, causal_idx])
            mean(abs(est[, causal_idx] - th_proj[, causal_idx]) < crit * se0[, causal_idx])
          } else NA_real_
          causal_cover_sci <- if (length(causal_idx)) {
            sci <- c(sci_a[causal_idx], sci_b[causal_idx])
            crit <- qt(0.975, df = fit$param_ref_df[, causal_idx])
            mean(abs(c(est[, causal_idx]) - sci) < crit * c(se0[, causal_idx]))
          } else NA_real_
          append_dataset(data.frame(cell = name, path = path, rep = r, source = src,
                                    W = wn, fit_ok = TRUE, fit_error = "",
                                    null_type1 = null_type1, null_any = null_any,
                                    causal_cover_proj = causal_cover_proj,
                                    causal_cover_sci = causal_cover_sci,
                                    stringsAsFactors = FALSE))
        }
      }
    }
  }
}

detail <- if (length(rows)) do.call(rbind, rows) else data.frame()
ds <- if (length(datasets)) do.call(rbind, datasets) else data.frame()
write.csv(detail, out_csv, row.names = FALSE, na = "")
write.csv(ds, sub("\\.csv$", "_datasets.csv", out_csv), row.names = FALSE, na = "")
cat("wrote", nrow(detail), "variant rows and", nrow(ds), "dataset-arm rows to", out_csv, "\n")
