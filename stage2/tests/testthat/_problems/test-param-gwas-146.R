# Extracted from test-param-gwas.R:146

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "gmoq", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
set.seed(42)
taus <- seq(0.05, 0.95, length.out = 15)
Tt <- length(taus)
W <- make_weights_vqtl(taus, qnorm(taus), mu = 0, sd = 1)
fq <- dnorm(qnorm(taus))
Sigma0 <- (outer(taus, taus, pmin) - outer(taus, taus)) / outer(fq, fq)
n_snp <- 400
D <- runif(n_snp, 2000, 20000)
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

# test -------------------------------------------------------------------------
taus <- c(0.25, 0.5, 0.75)
set.seed(3)
S0 <- crossprod(matrix(rnorm(9), 3)) + diag(3)
d <- c(100, 400)
se <- outer(1 / sqrt(d), sqrt(diag(S0)))
beta <- matrix(rnorm(6), 2) * se
lab <- sprintf("%.3f", taus)
base <- data.table::data.table(chr = "22", d_resid = d)
bs <- data.table::copy(base)
zt <- data.table::copy(base)
for (k in seq_along(taus)) {
    bs[[paste0("beta_tau_", lab[k])]] <- beta[, k]
    bs[[paste0("se_tau_", lab[k])]] <- se[, k]
    zt[[paste0("z_tau_", lab[k])]] <- beta[, k] / se[, k]
  }
local_mocked_bindings(read_sigma0 = function(sigma0_file) list(meta = list(), matrices = list(ALL = S0)))
from_bs <- gmoq:::.stage1_effects(bs, "x.stage1.tsv.gz")
from_z <- gmoq:::.stage1_effects(zt, "x.stage1.tsv.gz", sigma0_file = "x.sigma0.tsv.gz")
expect_equal(from_z$taus, from_bs$taus)
expect_equal(from_z$SE_tau, from_bs$SE_tau, tolerance = 1e-12)
expect_equal(from_z$Q_slope, from_bs$Q_slope, tolerance = 1e-12)
