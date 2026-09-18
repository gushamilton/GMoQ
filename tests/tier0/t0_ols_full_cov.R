#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Tier 0 -- the OLS-fit + full-covariance residual score.
#
# For an OLS fit with hat matrix H = W (W'W)^-1 W' and residual r = b - W theta,
# the coherent lack-of-fit score weights the residual by the per-SNP full
# covariance:
#
#     Q_full = r' Sigma^-1 r,    score = Q_full + 2 rank(W)
#
# because the optimism is exactly 2 tr(Sigma^-1 H Sigma) = 2 tr(H) = 2 rank(W).
# The diagonal-Q + 2k pairing is not that: with a diagonal loss the optimism is
# 2 tr(D^-1 H Sigma), which is not 2k when the taus are correlated. These checks
# pin the arithmetic, the rank accounting, and the refusal paths on tiny inputs
# where every expected number can be written down.
#
# Run:  Rscript --vanilla tests/tier0/t0_ols_full_cov.R
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  .self <- local({
    a <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (!length(f)) stop("Run as: Rscript --vanilla tests/tier0/t0_ols_full_cov.R")
    normalizePath(f[1L])
  })
  source(file.path(dirname(.self), "..", "lib", "harness.R"))
}))
gmoq_init_roots(.self)
gmoq_load_stage2()
gmoq_begin("tier0", "ols_full_cov")

fit_ols_fc <- gmoq_internal("fit_multi_models_ols_full_cov")

## ---- tiny, exactly-known inputs: T = 3 -------------------------------------
taus <- c(0.25, 0.5, 0.75)
b <- c(1, 2, 3)
s2 <- 0.25

upper_rowmajor <- function(S) {
  n <- nrow(S); out <- numeric(n * (n + 1) / 2); i <- 1L
  for (r in seq_len(n)) for (cc in r:n) { out[i] <- S[r, cc]; i <- i + 1L }
  out
}

## Sigma = s^2 I (spherical) and an equicorrelated Sigma with rho = 0.5.
Sigma_sp <- diag(s2, 3)
Sigma_co <- matrix(0.5, 3, 3); diag(Sigma_co) <- 1

stage1_sp <- list(taus = taus, q_tau = rep(0, 3),
                  Q_slope = matrix(b, ncol = 1),
                  SE_tau = matrix(sqrt(s2), ncol = 1))
colnames(stage1_sp$Q_slope) <- "snp1"
colnames(stage1_sp$SE_tau)  <- "snp1"
stage1_co <- stage1_sp
stage1_co$SE_tau <- matrix(1, ncol = 1)

cov_sp <- list(cov_vec = upper_rowmajor(Sigma_sp), offsets = 0)
cov_co <- list(cov_vec = upper_rowmajor(Sigma_co), offsets = 0)

W1 <- matrix(1, nrow = 3, ncol = 1)
colnames(W1) <- "mu"
W2 <- cbind(mu = rep(1, 3), slope = c(-1, 0, 1))

## =====================================================================
## 1. Spherical covariance: Q_full must equal the diagonal Q = r' D^-1 r
## =====================================================================
res_sp <- fit_ols_fc(matrix(b, nrow = 1), cov_sp$cov_vec, cov_sp$offsets,
                     list(W1 = W1), 100, 1L)
r_sp <- b - mean(b)
diag_sp <- sum(r_sp^2 / diag(Sigma_sp))
gmoq_close("T0.OLSFC.spherical_equals_diagonal",
           "With Sigma = s^2 I, Q_full = r' Sigma^-1 r equals the diagonal Q = r' D^-1 r.",
           as.numeric(res_sp$Q_full[1, 1]), diag_sp, 1e-12)

## =====================================================================
## 2. Correlated covariance: an exact, hand-computed number
## =====================================================================
# For an equicorrelated 3x3 with rho = 0.5, Sigma^-1 = 2 I - 0.5 J.
# With W = 1 (intercept), r = b - mean(b) = (-1, 0, 1), so
#   Q_full = r' Sigma^-1 r = 2 (r'r) - 0.5 (sum r)^2 = 4 - 0 = 4.
res_co <- fit_ols_fc(matrix(b, nrow = 1), cov_co$cov_vec, cov_co$offsets,
                     list(W1 = W1), 100, 1L)
gmoq_close("T0.OLSFC.correlated_exact",
           "With correlated Sigma, Q_full equals the hand-computed r' Sigma^-1 r = 4.",
           as.numeric(res_co$Q_full[1, 1]), 4, 1e-12)
gmoq_check("T0.OLSFC.correlated_differs_from_diagonal",
           "The correlated Q_full (4) differs from the diagonal Q (2), so the two
            losses really are different.",
           isTRUE(abs(res_co$Q_full[1, 1] - sum(r_sp^2 / diag(Sigma_co))) > 1e-6),
           stat = abs(res_co$Q_full[1, 1] - sum(r_sp^2 / diag(Sigma_co))),
           criterion = "|Q_full - Q_diagonal| > 1e-6")

## =====================================================================
## 3. Penalty identity on nested models
## =====================================================================
# b = (1, 2, 4) is not exactly in the span of W2, so both models have a
# non-zero Q_full and the identity is not vacuous.
b_pen <- c(1, 2, 4)
res_pen <- fit_ols_fc(matrix(b_pen, nrow = 1), cov_co$cov_vec, cov_co$offsets,
                      list(W1 = W1, W2 = W2), 100, 1L)
q1 <- res_pen$Q_full[1, 1]; q2 <- res_pen$Q_full[1, 2]
a1 <- res_pen$AIC[1, 1];    a2 <- res_pen$AIC[1, 2]
b1 <- res_pen$BIC[1, 1];    b2 <- res_pen$BIC[1, 2]
gmoq_close("T0.OLSFC.penalty_identity_aic",
           "For nested models, AIC difference equals the Q_full difference plus 2 * rank difference.",
           a2 - a1, (q2 - q1) + 2 * (2 - 1), 1e-12)
gmoq_close("T0.OLSFC.penalty_identity_bic",
           "For nested models, BIC difference equals the Q_full difference plus log(N) * rank difference.",
           b2 - b1, (q2 - q1) + log(100) * (2 - 1), 1e-12)
gmoq_close("T0.OLSFC.aic_uses_rank",
           "AIC is Q_full + 2 rank(W) for each model.",
           as.numeric(res_pen$AIC[1, ]), as.numeric(res_pen$Q_full[1, ]) + 2 * res_pen$rank, 1e-12)

## =====================================================================
## 4. OLS coefficients unchanged relative to the existing OLS path
## =====================================================================
out_pen <- suppressWarnings(param_gwas_multi(
  list(taus = taus, q_tau = rep(0, 3),
       Q_slope = matrix(b_pen, ncol = 1), SE_tau = matrix(1, ncol = 1)),
  list(W2 = W2), N = 100, se_mode = "tau_cov", cov_data = cov_co,
  residual_loss = "full_cov", return_params = TRUE))
colnames(out_pen$params$W2) <- "snp1"

fit_ols <- suppressWarnings(param_gwas(
  list(taus = taus, q_tau = rep(0, 3),
       Q_slope = structure(matrix(b_pen, ncol = 1), dimnames = list(NULL, "snp1")),
       SE_tau = structure(matrix(1, ncol = 1), dimnames = list(NULL, "snp1"))),
  transform = "custom_W", transform_args = list(W = W2),
  se_mode = "tau_cov", cov_data = cov_co, tau_cov_estimator = "ols",
  q_method = "chisq_naive"))

gmoq_close("T0.OLSFC.coef_matches_param_gwas_ols",
           "The full-covariance path returns exactly the same OLS point estimates as param_gwas(tau_cov_estimator = 'ols').",
           as.numeric(out_pen$params$W2), as.numeric(fit_ols$params), 1e-10)
A_ols <- solve(crossprod(W2), t(W2))
gmoq_close("T0.OLSFC.coef_matches_normal_equations",
           "The full-covariance path returns exactly (W'W)^-1 W'b.",
           as.numeric(out_pen$params$W2), as.numeric(A_ols %*% b_pen), 1e-10)

## =====================================================================
## 5. Rank deficiency: reported, and scored with the true rank
## =====================================================================
Wdup <- cbind(a = rep(1, 3), b = rep(1, 3))
res_dup <- fit_ols_fc(matrix(b_pen, nrow = 1), cov_co$cov_vec, cov_co$offsets,
                      list(Wdup = Wdup), 100, 1L)
gmoq_check("T0.OLSFC.rank_deficiency_detected",
           "A duplicated-column W is detected as rank deficient rather than treated as full rank.",
           isTRUE(res_dup$rank[1] < ncol(Wdup)) && isTRUE(res_dup$rank_deficient[1]),
           stat = res_dup$rank[1], criterion = "rank(W) < ncol(W) and rank_deficient = TRUE")
gmoq_close("T0.OLSFC.rank_deficiency_penalty",
           "The rank-deficient model is penalised by its true rank (1), not its column count (2).",
           as.numeric(res_dup$AIC[1, 1]),
           as.numeric(res_dup$Q_full[1, 1]) + 2 * res_dup$rank[1], 1e-12)

## =====================================================================
## 6. GLS contrast: OLS and GLS coincide under spherical Sigma
## =====================================================================
gls_sp_res <- gmoq_internal("fit_multi_models_tau_cov")(matrix(b_pen, nrow = 1), cov_sp$cov_vec,
                                       cov_sp$offsets, list(W2 = W2), 100, 1L)
ols_sp_res <- fit_ols_fc(matrix(b_pen, nrow = 1), cov_sp$cov_vec,
                         cov_sp$offsets, list(W2 = W2), 100, 1L)
gmoq_close("T0.OLSFC.gls_spherical_coincide",
           "With Sigma = s^2 I the OLS full-covariance score and the GLS score coincide.",
           as.numeric(ols_sp_res$AIC[1, 1]), as.numeric(gls_sp_res[1, 1]), 1e-10)

# On the correlated case they need not coincide; record rather than assert.
gls_co_res <- gmoq_internal("fit_multi_models_tau_cov")(matrix(b_pen, nrow = 1), cov_co$cov_vec,
                                       cov_co$offsets, list(W2 = W2), 100, 1L)
ols_co_res <- fit_ols_fc(matrix(b_pen, nrow = 1), cov_co$cov_vec,
                         cov_co$offsets, list(W2 = W2), 100, 1L)
gmoq_check("T0.OLSFC.gls_correlated_report",
           "On the correlated case the OLS and GLS scores are reported, not required to agree.",
           isTRUE(is.finite(ols_co_res$AIC[1, 1]) && is.finite(gls_co_res[1, 1])),
           stat = as.numeric(ols_co_res$AIC[1, 1] - gls_co_res[1, 1]),
           criterion = "both finite", severity = "report")

## =====================================================================
## 7. SNP/covariance offset alignment and per-SNP scale
## =====================================================================
b2 <- c(1, 2, 4)
beta_2 <- rbind(b, b2)
n_elem <- length(upper_rowmajor(Sigma_co))
cov_2 <- c(upper_rowmajor(Sigma_co), upper_rowmajor(Sigma_sp))
off_2 <- c(0, n_elem * 4)
res_2 <- fit_ols_fc(beta_2, cov_2, off_2, list(W1 = W1), 100, 1L)
q_a <- fit_ols_fc(matrix(b, nrow = 1), upper_rowmajor(Sigma_co), 0, list(W1 = W1), 100, 1L)$Q_full[1, 1]
q_b <- fit_ols_fc(matrix(b2, nrow = 1), upper_rowmajor(Sigma_sp), 0, list(W1 = W1), 100, 1L)$Q_full[1, 1]
gmoq_close("T0.OLSFC.offsets_aligned",
           "Each SNP is scored against the covariance record its byte offset names, not a neighbour's.",
           as.numeric(res_2$Q_full[, 1]), c(q_a, q_b), 1e-12)

res_scaled <- fit_ols_fc(matrix(b, nrow = 1), upper_rowmajor(Sigma_co), 0,
                         list(W1 = W1), 100, 1L, cov_scales = 2)
gmoq_close("T0.OLSFC.cov_scales_applied",
           "A per-SNP positive scale multiplies Sigma, so Q_full scales by its inverse.",
           as.numeric(res_scaled$Q_full[1, 1]), as.numeric(res_co$Q_full[1, 1]) / 2, 1e-12)

## =====================================================================
## 8. Refusal paths return NA with a recorded reason
## =====================================================================
Sig_neg <- matrix(c(1, 2, 0, 2, 1, 0, 0, 0, 1), nrow = 3)   # eigenvalues -1, 1, 3
res_npd <- fit_ols_fc(matrix(b_pen, nrow = 1), upper_rowmajor(Sig_neg), 0,
                      list(W2 = W2), 100, 1L)
gmoq_check("T0.OLSFC.non_pd_refused",
           "A non-positive-definite Sigma returns NA with reason sigma_not_pd, with no silent repair.",
           is.na(res_npd$Q_full[1, 1]) && res_npd$reason_code[1, 1] == 5L,
           stat = res_npd$reason_code[1, 1], criterion = "Q_full is NA and reason_code = 5")

cov_bad <- cov_co$cov_vec; cov_bad[1] <- NA_real_
res_bad <- fit_ols_fc(matrix(b_pen, nrow = 1), cov_bad, 0, list(W2 = W2), 100, 1L)
gmoq_check("T0.OLSFC.nonfinite_cov_refused",
           "A non-finite covariance entry returns NA with reason cov_nonfinite.",
           is.na(res_bad$Q_full[1, 1]) && res_bad$reason_code[1, 1] == 4L,
           stat = res_bad$reason_code[1, 1], criterion = "Q_full is NA and reason_code = 4")

res_off <- fit_ols_fc(matrix(b_pen, nrow = 1), cov_co$cov_vec, 1e9,
                      list(W2 = W2), 100, 1L)
gmoq_check("T0.OLSFC.bad_offset_refused",
           "An out-of-range covariance offset returns NA with reason offset_bad.",
           is.na(res_off$Q_full[1, 1]) && res_off$reason_code[1, 1] == 2L,
           stat = res_off$reason_code[1, 1], criterion = "Q_full is NA and reason_code = 2")

res_na_off <- fit_ols_fc(matrix(b_pen, nrow = 1), cov_co$cov_vec, NA_real_,
                         list(W2 = W2), 100, 1L)
gmoq_check("T0.OLSFC.na_offset_refused",
           "A non-finite covariance offset returns NA with reason offset_bad.",
           is.na(res_na_off$Q_full[1, 1]) && res_na_off$reason_code[1, 1] == 2L,
           stat = res_na_off$reason_code[1, 1], criterion = "Q_full is NA and reason_code = 2")

## =====================================================================
## 9. The default path is unchanged
## =====================================================================
auto <- param_gwas_multi(
  list(taus = taus, q_tau = rep(0, 3),
       Q_slope = structure(matrix(b_pen, ncol = 1), dimnames = list(NULL, "snp1")),
       SE_tau = matrix(1, ncol = 1)),
  list(W1 = W1, W2 = W2), N = 100, se_mode = "tau_cov", cov_data = cov_co)
ref <- gmoq_internal("fit_multi_models_tau_cov")(matrix(b_pen, nrow = 1), cov_co$cov_vec,
                                cov_co$offsets, list(W1 = W1, W2 = W2), 100, 1L)
gmoq_close("T0.OLSFC.default_auto_unchanged",
           "With residual_loss = 'auto' the returned matrix is bit-for-bit the existing GLS tournament.",
           as.numeric(auto), as.numeric(ref), 1e-14)
gmoq_check("T0.OLSFC.default_columns_unchanged",
           "The default tournament still returns the AIC_/BIC_ column names.",
           identical(colnames(auto), c("AIC_W1", "BIC_W1", "AIC_W2", "BIC_W2")),
           criterion = "colnames are AIC_/BIC_ per model")

gmoq_finish()
