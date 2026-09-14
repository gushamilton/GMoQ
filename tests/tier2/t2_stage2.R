#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Tier 2, Stage 2 side: fit one chunk of replicates and accumulate COUNTS.
#
#   Rscript --vanilla tests/tier2/t2_stage2.R <chunk_tag> <dir>
#
# Nothing here decides anything. It counts events -- intervals that covered,
# tests that rejected -- and writes them to an .rds. tests/tier2/t2_collate.R
# turns counts into rates, attaches Monte Carlo error, and compares them with
# the criteria that were committed before the run.
#
# Four estimator arms are evaluated per replicate at no extra Stage 1 cost:
#   {jackknife, pooled} x {ols, gls}
# and, per arm, the location-scale model a real scan would fit plus (where the
# family has one) the parsimonious TRUE model.
#
# WHAT COVERAGE IS OF. theta_pop = A b(tau), the projection of Stage 1's own
# exact population estimand onto the fitted W. That is the quantity the
# estimator is entitled to be judged on. Coverage of the generative parameter
# is a different question and is reported separately, because under a
# misspecified W they are not the same number and conflating them is how a
# suite talks itself into a false conclusion.
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  .self <- local({
    a <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (!length(f)) stop("Run as: Rscript --vanilla tests/tier2/t2_stage2.R <tag> <dir>")
    normalizePath(f[1L])
  })
  source(file.path(dirname(.self), "..", "lib", "harness.R"))
}))
gmoq_init_roots(.self)
gmoq_load_stage2()
source(file.path(GMOQ_SUITE_ROOT, "lib", "models.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("usage: t2_stage2.R <chunk_tag> <dir>")
tag <- args[1L]; dir <- normalizePath(args[2L])

meta  <- jsonlite::fromJSON(file.path(dir, paste0(tag, ".json")), simplifyVector = FALSE)
cellnm <- meta$cell$name
truth <- jsonlite::fromJSON(file.path(dir, paste0(cellnm, ".truth.json")), simplifyVector = FALSE)

T  <- as.integer(meta$T); M <- as.integer(meta$M)
R  <- as.integer(meta$reps); NE <- as.integer(meta$n_elem)
taus <- as.numeric(unlist(meta$taus))
Nn <- as.integer(meta$cell$N); B <- as.integer(meta$cell$blocks)
causal <- as.integer(unlist(truth$causal))
nulls  <- setdiff(seq_len(M), causal)
maf <- as.numeric(unlist(truth$maf))

rd <- function(suffix, what, size, n, dims) {
  v <- readBin(file.path(dir, paste0(tag, suffix)), what = what, size = size, n = n,
               endian = "little")
  array(v, dim = dims)          # C order (last index fastest) -> reversed dims
}
beta <- rd(".beta.bin",   "double", 8, R * M * T,  c(T, M, R))
se   <- rd(".se.bin",     "double", 8, R * M * T,  c(T, M, R))
covA <- rd(".cov.bin",    "double", 4, R * M * NE, c(NE, M, R))
dres <- rd(".dresid.bin", "double", 8, R * M,      c(M, R))
s0A  <- rd(".sigma0.bin", "double", 8, R * T * T,  c(T, T, R))

b_true <- array(readBin(file.path(dir, paste0(cellnm, ".btrue.bin")), "double",
                        size = 8, n = M * T, endian = "little"), dim = c(T, M))

# --- the models fitted ------------------------------------------------------
y_sd <- mean(as.numeric(unlist(meta$y_sd)))
mixarg <- NULL
if (!is.null(truth$mix)) {
  mixarg <- list(q_tau = as.numeric(unlist(truth$q_pop)),
                 p1 = truth$mix$p1, mu1 = truth$mix$mu1, sd1 = truth$mix$sd1,
                 mu2 = truth$mix$mu2, sd2 = truth$mix$sd2)
}
all_models <- gmoq_model_set(taus, sigma = y_sd, mix = mixarg)
fam <- meta$cell$family
tmap <- GMOQ_TRUTH_MAP[[fam]]
want <- unlist(meta$cell$fit_models)
models <- list()
if ("LS" %in% want) models[["LS"]] <- all_models[["LS"]]
if ("TRUE" %in% want && !is.null(tmap)) models[["TRUE"]] <- all_models[[tmap$true]]
# also fit the models that CANNOT represent this family's truth, so the other
# half of the tournament claim -- that they are rejected -- is measured at scale
wrong <- if (is.null(tmap)) character(0) else intersect(tmap$wrong, names(all_models))
for (w in wrong) models[[paste0("WRONG_", w)]] <- all_models[[w]]

pooled_cov_data <- gmoq_internal(".pooled_cov_data")
chrv <- rep("1", M)
offsets <- seq(0, by = NE * 4, length.out = M)

acc <- list()
key <- function(...) paste(..., sep = "|")

bump <- function(k, field, by = 1) {
  if (is.null(acc[[k]])) acc[[k]] <<- list()
  acc[[k]][[field]] <<- (if (is.null(acc[[k]][[field]])) 0 else acc[[k]][[field]]) + by
}

alphas <- c(0.05, 1e-3, 1e-4)

for (r in seq_len(R)) {
  Qs <- beta[, , r, drop = TRUE]; dim(Qs) <- c(T, M)
  Ss <- se[, , r, drop = TRUE];   dim(Ss) <- c(T, M)
  colnames(Qs) <- colnames(Ss) <- paste0("v", seq_len(M))
  st <- list(taus = taus, q_tau = rep(0, T), Q_slope = Qs, SE_tau = Ss)

  cv_jk <- list(cov_vec = as.numeric(covA[, , r]), offsets = offsets)
  S0 <- s0A[, , r]; dim(S0) <- c(T, T)
  cd_pl <- pooled_cov_data(list(ALL = S0), chrv, dres[, r], T)
  cv_pl <- list(cov_vec = cd_pl$cov_vec, offsets = cd_pl$offsets)

  for (src in c("jackknife", "pooled")) {
    cvd <- if (src == "jackknife") cv_jk else cv_pl
    for (est in c("ols", "gls")) {
      for (mn in names(models)) {
        W <- models[[mn]]
        fit <- tryCatch(suppressWarnings(suppressMessages(param_gwas(
          st, transform = "custom_W", transform_args = list(W = W),
          se_mode = "tau_cov", tau_cov_estimator = est, cov_data = cvd,
          q_method = "saddlepoint",
          n_blocks = if (src == "jackknife") B else NULL,
          n_eff = if (src == "jackknife") dres[, r] else NULL,
          nu = if (src == "pooled") Nn else NULL,
          nu_source = if (src == "pooled") "pooled_sigma0_n_samples" else NULL))),
          error = function(e) NULL)
        if (is.null(fit)) { bump(key(src, est, mn, "_fit"), "n_fail"); next }

        A <- solve(crossprod(W), t(W))
        th_pop <- A %*% b_true                       # K x M
        pn <- rownames(fit$params)
        for (pi in seq_along(pn)) {
          est_v <- fit$params[pi, ]; se_v <- fit$SE_params[pi, ]
          pv <- if (is.null(fit$P_params)) rep(NA_real_, M) else fit$P_params[pi, ]
          z_pop <- (est_v - th_pop[pi, ]) / se_v
          for (grp in c("null", "causal")) {
            idx <- if (grp == "null") nulls else causal
            if (!length(idx)) next
            k <- key(src, est, mn, pn[pi], grp)
            ok <- is.finite(z_pop[idx]) & is.finite(se_v[idx]) & se_v[idx] > 0
            bump(k, "n", sum(ok))
            bump(k, "cover95", sum(abs(z_pop[idx][ok]) < 1.959964))
            bump(k, "sum_z", sum(z_pop[idx][ok]))
            bump(k, "sum_z2", sum(z_pop[idx][ok]^2))
            bump(k, "sum_est", sum(est_v[idx][ok]))
            bump(k, "sum_true", sum(th_pop[pi, idx][ok]))
            # Bias and empirical sd are of the DEVIATION from each variant's
            # own truth. Taking the spread of the estimates themselves would
            # include the spread of the true values across variants, which in
            # a causal group is most of it.
            dv <- est_v[idx][ok] - th_pop[pi, idx][ok]
            bump(k, "sum_d", sum(dv))
            bump(k, "sum_d2", sum(dv^2))
            bump(k, "sum_se", sum(se_v[idx][ok]))
            bump(k, "sum_se2", sum(se_v[idx][ok]^2))
            for (al in alphas) {
              bump(k, paste0("rej_", format(al, scientific = TRUE)),
                   sum(is.finite(pv[idx]) & pv[idx] < al))
            }
            bump(k, "n_p", sum(is.finite(pv[idx])))
          }
        }
        # Q lack of fit
        for (grp in c("null", "causal")) {
          idx <- if (grp == "null") nulls else causal
          if (!length(idx)) next
          k <- key(src, est, mn, "Q", grp)
          qp <- fit$Q_pval
          if (is.null(qp)) qp <- rep(NA_real_, M)
          bump(k, "n", sum(is.finite(qp[idx])))
          bump(k, "rej05", sum(is.finite(qp[idx]) & qp[idx] < 0.05))
          bump(k, "rej1e-3", sum(is.finite(qp[idx]) & qp[idx] < 1e-3))
        }

        # --- the false-vQTL pin: bucket the dispersion rejection by G -------
        if (mn == "LS" && "beta_mu" %in% pn && "beta_sigma" %in% pn) {
          z_mu <- fit$params["beta_mu", ] / fit$SE_params["beta_mu", ]
          Gv <- dispersion_skew_guard(z_mu, N = Nn, maf = maf)
          p_sig <- if (is.null(fit$P_params)) rep(NA_real_, M) else fit$P_params["beta_sigma", ]
          true_sig <- (A %*% b_true)[which(pn == "beta_sigma"), ]
          # only variants whose TRUE dispersion effect is zero can produce a
          # false dispersion finding, which is the whole point of the pin
          zero <- abs(true_sig) < 1e-12
          for (bk in list(c(-Inf, 0.05), c(0.05, 0.2), c(0.2, Inf))) {
            sel <- zero & is.finite(Gv) & abs(Gv) >= bk[1] & abs(Gv) < bk[2] &
              is.finite(p_sig)
            if (!any(sel)) next
            k <- key(src, est, "GBUCKET", sprintf("[%g,%g)", bk[1], bk[2]), "-")
            bump(k, "n", sum(sel))
            bump(k, "rej05", sum(p_sig[sel] < 0.05))
            bump(k, "rej1e-3", sum(p_sig[sel] < 1e-3))
            bump(k, "sum_absG", sum(abs(Gv[sel])))
            bump(k, "sum_absZmu", sum(abs(z_mu[sel])))
          }
        }

        if (!is.null(fit$gls_diag)) {
          k <- key(src, est, mn, "GLSDIAG", "-")
          bump(k, "n_snps", fit$gls_diag$n_snps)
          bump(k, "n_ridge_tried", fit$gls_diag$n_ridge_tried)
          bump(k, "n_sigma_not_spd", fit$gls_diag$n_sigma_not_spd)
        }
      }
    }
  }
}

out <- list(cell = meta$cell, tag = tag, reps = R, M = M, T = T,
            taus = taus, causal = causal, maf = maf, y_sd = y_sd,
            models = names(models), acc = acc,
            stage1_seconds = meta$stage1_seconds)
saveRDS(out, file.path(dir, paste0(tag, ".s2.rds")))
cat("wrote", file.path(dir, paste0(tag, ".s2.rds")), "with",
    length(acc), "accumulators\n")
