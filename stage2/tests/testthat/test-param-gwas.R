# Mapping quantile effects to parameters: OLS estimates, the sandwich
# covariance under a pooled Stage 1 covariance, p-values, the lack-of-fit Q,
# the GLS path with its Hotelling reference, and refusals.

set.seed(42)
taus <- seq(0.05, 0.95, length.out = 15)
Tt <- length(taus)
W <- make_weights_vqtl(taus, qnorm(taus), mu = 0, sd = 1)
# A pooled cross-quantile covariance (Sigma_0) with the shape of a RIF covariance:
# min(tau_s, tau_t) - tau_s tau_t over the product of the densities.
fq <- dnorm(qnorm(taus))
Sigma0 <- (outer(taus, taus, pmin) - outer(taus, taus)) / outer(fq, fq)
n_snp <- 400
D <- runif(n_snp, 2000, 20000)            # sum of squared residual genotypes
theta <- rbind(beta_mu = rep(c(0, 0.05), length.out = n_snp), beta_sigma2 = 0)
R <- chol(Sigma0)
b <- vapply(seq_len(n_snp), function(j)
  drop(W %*% theta[, j]) + drop(crossprod(R, rnorm(Tt))) / sqrt(D[j]), numeric(Tt))
colnames(b) <- paste0("rs", seq_len(n_snp))
se_tau <- sqrt(outer(diag(Sigma0), D, "/"))
stage1 <- list(taus = taus, q_tau = qnorm(taus), Q_slope = b, SE_tau = se_tau)
cov_data <- gmoq:::.pooled_cov_data(list(ALL = Sigma0), rep("1", n_snp), D, Tt)

fit_ols <- function(...) param_gwas(stage1, transform_args = list(W = W),
                                    se_mode = "tau_cov", cov_data = cov_data, ...)

test_that("OLS estimates and sandwich covariance match a direct computation", {
  fit <- suppressMessages(fit_ols(nu = 1e5, return_cov = "full"))
  A <- solve(crossprod(W), t(W))
  expect_equal(unname(fit$params), unname(A %*% b), tolerance = 1e-12)
  V1 <- A %*% (Sigma0 / D[1]) %*% t(A)
  expect_equal(unname(fit$param_cov[, , 1]), unname(V1), tolerance = 1e-10)
  se <- sqrt(outer(diag(A %*% Sigma0 %*% t(A)), D, "/"))
  expect_equal(unname(fit$SE_params), unname(se), tolerance = 1e-10)
  # two-sided t p-values on the supplied degrees of freedom
  expect_equal(unname(fit$P_params),
               unname(2 * pt(-abs(fit$params / fit$SE_params), 1e5)), tolerance = 1e-12)
  # estimates recover the simulated effects, and null SNPs are calibrated
  z <- (fit$params - theta) / fit$SE_params
  expect_lt(abs(mean(z)), 0.1)
  expect_gt(sd(z), 0.85); expect_lt(sd(z), 1.15)
})

test_that("Q is the diagonal lack-of-fit statistic with a calibrated p-value", {
  fit <- suppressMessages(fit_ols(nu = 1e5))
  r <- b - W %*% fit$params
  expect_equal(fit$Q, unname(colSums(r^2 / se_tau^2)), tolerance = 1e-10)
  expect_equal(fit$df, Tt - 2)
  expect_true(all(fit$Q_pval >= 0 & fit$Q_pval <= 1))
  expect_gt(mean(fit$Q_pval < 0.05), 0.02)       # correctly specified: about 5%
  expect_lt(mean(fit$Q_pval < 0.05), 0.09)

  # a quantile profile the Normal model cannot produce: skewness
  bad <- stage1
  bad$Q_slope <- b + 0.06 * (qnorm(taus)^2 - 1)
  fit_bad <- suppressMessages(param_gwas(bad, transform_args = list(W = W),
                                         se_mode = "tau_cov", cov_data = cov_data,
                                         nu = 1e5))
  expect_gt(mean(fit_bad$Q_pval < 0.05), 0.5)
})

test_that("GLS estimates and the Hotelling reference for Q", {
  nu <- 500
  fit <- suppressMessages(fit_ols(tau_cov_estimator = "gls", nu = nu))
  Si <- solve(Sigma0)
  A_gls <- solve(t(W) %*% Si %*% W, t(W) %*% Si)
  expect_equal(unname(fit$params), unname(A_gls %*% b), tolerance = 1e-8)
  r <- b - W %*% fit$params
  Q <- unname(colSums(r * (Si %*% r)) * D)
  expect_equal(fit$Q, Q, tolerance = 1e-8)
  m <- Tt - 2
  expect_equal(fit$Q_pval,
               pf(Q * (nu - m + 1) / (nu * m), m, nu - m + 1, lower.tail = FALSE),
               tolerance = 1e-8)
  expect_identical(fit$Q_ref_dist, "hotelling_f")
  expect_error(fit_ols(q_method = "hotelling"), "applies only")
})

test_that("unknown degrees of freedom give no p-value rather than a normal one", {
  expect_warning(fit <- suppressMessages(fit_ols()), "degrees of freedom|nu")
  expect_null(fit$P_params)
  expect_error(suppressMessages(fit_ols(unknown_nu = "error")))
})

test_that("singular and mis-sized designs are refused", {
  W_sing <- cbind(a = W[, 1], b = 2 * W[, 1])
  expect_error(param_gwas(stage1, transform_args = list(W = W_sing), se_mode = "diagonal"))
  expect_error(param_gwas(stage1, transform_args = list(W = W[-1, ]), se_mode = "diagonal"),
               "nrow\\(W\\)")
  expect_error(param_gwas(stage1, se_mode = "tau_cov"), "transform_args\\$W|cov_data")
  expect_error(param_gwas(stage1, transform_args = list(W = W), se_mode = "tau_cov"),
               "requires cov_data")
})

test_that("cached pooled SE, param_cov, and Q calibration match uncached per-SNP path", {
  # Build an uncached copy where every SNP has its own distinct offset in cov_vec.
  n_cov <- Tt * (Tt + 1L) / 2L
  s0_upper <- Sigma0[upper.tri(Sigma0, diag = TRUE)]
  # row-major upper triangle matching stage2_se_helper.cpp
  s0_rm <- numeric(n_cov)
  k <- 1L
  for (r in seq_len(Tt)) {
    for (c in r:Tt) {
      s0_rm[k] <- Sigma0[r, c]
      k <- k + 1L
    }
  }
  cov_uncached <- list(
    cov_vec = as.numeric(outer(s0_rm, 1 / D)),
    offsets = seq(0, by = n_cov * 4, length.out = n_snp),
    cov_scales = NULL
  )
  fit_cached   <- suppressMessages(fit_ols(nu = 1e5, return_cov = "upper"))
  fit_uncached <- suppressMessages(param_gwas(
    stage1, transform_args = list(W = W), se_mode = "tau_cov",
    cov_data = cov_uncached, nu = 1e5, return_cov = "upper"
  ))
  expect_equal(fit_cached$SE_params, fit_uncached$SE_params, tolerance = 1e-12)
  expect_equal(fit_cached$param_cov_upper, fit_uncached$param_cov_upper, tolerance = 1e-12)
  expect_equal(fit_cached$Q, fit_uncached$Q, tolerance = 1e-12)
  expect_equal(fit_cached$Q_pval, fit_uncached$Q_pval, tolerance = 1e-11)
  expect_equal(fit_cached$Q_df_eff, fit_uncached$Q_df_eff, tolerance = 1e-12)
  expect_equal(fit_cached$Q_scale, fit_uncached$Q_scale, tolerance = 1e-12)
})
