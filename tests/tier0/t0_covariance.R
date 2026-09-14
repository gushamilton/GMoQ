#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Tier 0 -- the covariance machinery against exact expressions on tiny inputs.
#
# Everything the C++ kernels do here is a few lines of linear algebra, so on a
# 5-tau grid with a hand-built Sigma the exact answer can be written down in R
# and compared to machine precision. Nothing is tolerance-tuned.
#
# The one non-obvious contract being pinned is the WIRE FORMAT: `cov_vec` holds
# the row-major upper triangle including the diagonal, and `offsets` are BYTE
# offsets that the kernels divide by 4. Getting either wrong silently gives
# every SNP the wrong covariance, which is why a permuted-offset case is here.
#
# Run:  Rscript --vanilla tests/tier0/t0_covariance.R
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  .self <- local({
    a <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (!length(f)) stop("Run as: Rscript --vanilla tests/tier0/t0_covariance.R")
    normalizePath(f[1L])
  })
  source(file.path(dirname(.self), "..", "lib", "harness.R"))
}))
gmoq_init_roots(.self)
gmoq_load_stage2()
gmoq_begin("tier0", "covariance")

# The compiled kernels are registered in the namespace but not exported, so
# they are reached the same way the package reaches them.
compute_calibrated_se  <- gmoq_internal("compute_calibrated_se")
compute_param_cov      <- gmoq_internal("compute_param_cov")
compute_gls_param_gwas <- gmoq_internal("compute_gls_param_gwas")

set.seed(4242)
Tt <- 5L; K_par <- 2L
taus <- c(0.1, 0.3, 0.5, 0.7, 0.9)

# --- a tiny, exactly-known Sigma -------------------------------------------
# Built as the Brownian-bridge form so it is genuinely positive definite and
# genuinely has the cross-tau structure the method lives with.
fpd <- dnorm(qnorm(taus))
Sig1 <- (outer(taus, taus, pmin) - outer(taus, taus)) / outer(fpd, fpd) / 1000
Sig2 <- Sig1 * 2.5 + diag(1e-4, Tt)
Sig3 <- Sig1 * 0.4 + diag(3e-5, Tt)
sigmas <- list(Sig1, Sig2, Sig3)

# Written out as an explicit loop so there is no doubt about the ordering.
upper_rowmajor <- function(S) {
  n <- nrow(S); out <- numeric(n * (n + 1) / 2); i <- 1L
  for (r in seq_len(n)) for (cc in r:n) { out[i] <- S[r, cc]; i <- i + 1L }
  out
}

n_elem <- Tt * (Tt + 1L) / 2L
cov_vec <- c(upper_rowmajor(Sig1), upper_rowmajor(Sig2), upper_rowmajor(Sig3))
offsets <- seq(0, by = n_elem * 4, length.out = 3L)

W <- cbind(beta_mu = rep(1, Tt), beta_sigma = qnorm(taus))
A <- solve(crossprod(W), t(W))            # K x T, the fixed OLS contrast

# =====================================================================
# 1. compute_calibrated_se -- the OLS delta-method standard errors
# =====================================================================
se_cpp <- compute_calibrated_se(cov_vec, offsets, A, Tt)   # n_snps x K
se_ref <- t(vapply(sigmas, function(S) sqrt(diag(A %*% S %*% t(A))), numeric(K_par)))
gmoq_close("T0.COV.calibrated_se",
           "compute_calibrated_se returns sqrt(diag(A Sigma A')) exactly, per SNP.",
           as.numeric(se_cpp), as.numeric(se_ref), 1e-12)

# =====================================================================
# 2. compute_param_cov -- the full parameter covariance upper triangle
# =====================================================================
pc_cpp <- compute_param_cov(cov_vec, offsets, A, Tt)       # n_snps x K(K+1)/2
pc_ref <- t(vapply(sigmas, function(S) upper_rowmajor(A %*% S %*% t(A)), numeric(3L)))
gmoq_close("T0.COV.param_cov",
           "compute_param_cov returns the row-major upper triangle of A Sigma A'.",
           as.numeric(pc_cpp), as.numeric(pc_ref), 1e-14)

# its diagonal must reproduce the SEs computed by the other kernel
gmoq_close("T0.COV.param_cov_matches_se",
           "The two covariance kernels agree: diag of A Sigma A' is se^2.",
           sqrt(pc_cpp[, c(1, 3)]), se_cpp, 1e-12)

# =====================================================================
# 3. BYTE offsets, and the ordering they imply
# =====================================================================
perm <- c(3L, 1L, 2L)
se_perm <- compute_calibrated_se(cov_vec, offsets[perm], A, Tt)
gmoq_close("T0.COV.offsets_are_bytes",
           "offsets are byte offsets into a float32 payload; a permutation reorders the SNPs.",
           as.numeric(se_perm), as.numeric(se_ref[perm, ]), 1e-12)

# A transposed Sigma would be indistinguishable if the triangle were packed
# column-major, so pin the ordering with an asymmetric probe.
Asym <- matrix(0, Tt, Tt); Asym[1, 2] <- 1; Asym <- Asym + t(Asym) + diag(Tt)
probe <- upper_rowmajor(Asym)
gmoq_check("T0.COV.upper_triangle_rowmajor",
           "The stored triangle is row-major with the diagonal first in each row.",
           abs(probe[2] - 1) < 1e-15 && abs(probe[1] - 1) < 1e-15 && abs(probe[3]) < 1e-15,
           criterion = "element 2 of the packed triangle is Sigma[1,2]")

# =====================================================================
# 4. compute_gls_param_gwas -- the whole GLS fit, written out in R
# =====================================================================
b <- matrix(c(0.031, 0.028, 0.030, 0.034, 0.041,
              0.010, 0.012, 0.009, 0.011, 0.013,
              -0.02, -0.01, 0.00, 0.01, 0.02), nrow = Tt)
fit <- compute_gls_param_gwas(b, cov_vec, offsets, W, return_cov = TRUE,
                              gls_ridge = 1e-8, n_threads = 1L)
th_ref <- se_ref_g <- matrix(NA_real_, K_par, 3L); Q_ref <- numeric(3L)
for (j in 1:3) {
  Si <- solve(sigmas[[j]])
  V  <- solve(t(W) %*% Si %*% W)
  th <- V %*% t(W) %*% Si %*% b[, j]
  r  <- b[, j] - W %*% th
  th_ref[, j] <- th; se_ref_g[, j] <- sqrt(diag(V))
  Q_ref[j] <- as.numeric(t(r) %*% Si %*% r)
}
gmoq_close("T0.COV.gls_theta",
           "GLS returns (W'S^-1 W)^-1 W'S^-1 b exactly.",
           as.numeric(fit$params), as.numeric(th_ref), 1e-10)
gmoq_close("T0.COV.gls_se",
           "GLS standard errors are sqrt(diag((W'S^-1 W)^-1)) exactly.",
           as.numeric(fit$se), as.numeric(se_ref_g), 1e-12)
gmoq_close("T0.COV.gls_Q",
           "The GLS lack-of-fit statistic is r' Sigma^-1 r with r the GLS residual.",
           fit$Q, Q_ref, 1e-9)
gmoq_close("T0.COV.gls_param_cov",
           "The GLS parameter covariance is (W'S^-1 W)^-1, upper triangle row-major.",
           as.numeric(fit$param_cov_upper[1, ]),
           upper_rowmajor(solve(t(W) %*% solve(Sig1) %*% W)), 1e-12)

# =====================================================================
# 5. OLS via param_gwas() -- the same numbers, through the public entry point
# =====================================================================
colnames(b) <- paste0("SNP", 1:3)
stage1 <- list(taus = taus, q_tau = qnorm(taus), Q_slope = b,
               SE_tau = vapply(sigmas, function(S) sqrt(diag(S)), numeric(Tt)))
fit_ols <- param_gwas(stage1, transform = "custom_W", transform_args = list(W = W),
                      se_mode = "tau_cov", tau_cov_estimator = "ols",
                      cov_data = list(cov_vec = cov_vec, offsets = offsets),
                      q_method = "chisq_naive", nu = 1000, nu_source = "test")
gmoq_close("T0.COV.ols_params",
           "The OLS point estimate is the fixed projection A b, independent of Sigma.",
           as.numeric(fit_ols$params), as.numeric(A %*% b), 1e-14)
gmoq_close("T0.COV.ols_se",
           "The OLS standard error is the sandwich sqrt(diag(A Sigma A')).",
           as.numeric(fit_ols$SE_params), as.numeric(t(se_ref)), 1e-12)

# =====================================================================
# 6. The pooled covariance identity: Sigma_v = Sigma_0 / d_resid
# =====================================================================
pooled_cov_data <- gmoq_internal(".pooled_cov_data")
d_resid <- c(1000, 4000, 250)
cd <- pooled_cov_data(list(ALL = Sig1), chr = c("1", "1", "1"),
                      d_resid = d_resid, n_taus = Tt)
exp_vec <- as.numeric(vapply(d_resid, function(d) upper_rowmajor(Sig1 / d), numeric(n_elem)))
gmoq_close("T0.COV.pooled_identity",
           "The pooled path forms Sigma_0/d_resid per variant, in the jackknife wire format.",
           cd$cov_vec, exp_vec, 1e-16)
gmoq_close("T0.COV.pooled_offsets",
           "The pooled path produces the same byte offsets the jackknife reader does.",
           cd$offsets, seq(0, by = n_elem * 4, length.out = 3L), 0)

# SE must then scale exactly as 1/sqrt(d_resid)
se_pool <- compute_calibrated_se(cd$cov_vec, cd$offsets, A, Tt)
gmoq_check("T0.COV.pooled_se_scales",
           "Pooled standard errors scale exactly as 1/sqrt(d_resid).",
           max(abs(se_pool[2, ] / se_pool[1, ] - sqrt(d_resid[1] / d_resid[2]))) < 1e-12,
           stat = max(abs(se_pool[2, ] / se_pool[1, ] - sqrt(d_resid[1] / d_resid[2]))),
           criterion = "max deviation < 1e-12")

# a non-positive d_resid must produce NA, not a silently huge covariance
cd_bad <- pooled_cov_data(list(ALL = Sig1), chr = rep("1", 2), d_resid = c(0, -3), n_taus = Tt)
gmoq_check("T0.COV.pooled_bad_d_resid",
           "A non-positive d_resid yields NA covariance and is counted, not silently scaled.",
           cd_bad$n_bad == 2L && all(is.na(cd_bad$cov_vec)),
           criterion = "n_bad == 2 and cov_vec all NA")

# per-chromosome keying: a variant must get its own chromosome's matrix
cd_chr <- pooled_cov_data(list(`1` = Sig1, `2` = Sig2), chr = c("1", "2"),
                          d_resid = c(1, 1), n_taus = Tt)
gmoq_close("T0.COV.pooled_per_chromosome",
           "Each variant gets the Sigma_0 of its own chromosome, keyed on the normalised label.",
           cd_chr$cov_vec, c(upper_rowmajor(Sig1), upper_rowmajor(Sig2)), 1e-16)
gmoq_refuses("T0.COV.pooled_missing_chromosome",
             "A chromosome with no Sigma_0 is refused, not substituted with another chromosome's.",
             pooled_cov_data(list(`1` = Sig1), chr = c("1", "7"), d_resid = c(1, 1), n_taus = Tt),
             "no Sigma_0 for chromosome")

# =====================================================================
# 7. Array conversion and mashr inputs
# =====================================================================
upper_to_array <- gmoq_internal(".param_cov_upper_to_array")
upper_names    <- gmoq_internal(".param_cov_upper_names")
arr <- upper_to_array(pc_cpp, K = K_par, param_names = colnames(W))
gmoq_close("T0.COV.upper_to_array",
           "Unpacking the triangle reproduces the symmetric K x K matrix per SNP.",
           as.numeric(arr[, , 2]), as.numeric(A %*% Sig2 %*% t(A)), 1e-14)
gmoq_check("T0.COV.upper_names",
           "Parameter covariance columns are named cov_<i>_<j> in row-major order.",
           identical(upper_names(c("beta_mu", "beta_sigma")),
                     c("cov_beta_mu_beta_mu", "cov_beta_mu_beta_sigma",
                       "cov_beta_sigma_beta_sigma")),
           criterion = "documented naming scheme")

mh <- mashr_inputs_from_cov(A %*% b, pc_cpp, mode = "array")
gmoq_close("T0.COV.mashr_shat",
           "mashr Shat is the per-SNP parameter standard error.",
           as.numeric(mh$Shat), as.numeric(se_ref), 1e-12)
gmoq_check("T0.COV.mashr_v_is_correlation",
           "mashr V is a correlation matrix: unit diagonal per SNP.",
           max(abs(apply(mh$V, 3, function(M) diag(M)) - 1)) < 1e-12,
           criterion = "max |diag(V) - 1| < 1e-12")

gmoq_finish()
