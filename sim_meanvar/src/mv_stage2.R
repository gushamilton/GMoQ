## Stage 2 driver for the mean-variance (dispersion) Monte Carlo.
##
## For each replicate: build the pooled cov_data from the single Sigma_0 that
## replicate's phenotype produced, build the vQTL W from the SHIPPED
## make_weights_vqtl(), run the SHIPPED param_gwas() with return_cov = "upper"
## for both tau_cov_estimator settings, then apply the shipped
## estimate_rav() / residualise_dispersion() over ALL 200 variants.
##
## ~/.Rprofile carries a bare .libPaths() that silently overrides R_LIBS_USER,
## so this must be run with Rscript --vanilla and set .libPaths() itself.
.libPaths(c("/user/work/fh6520/gls_coverage/lib", "/user/work/fh6520/r_packages"))
stopifnot("fungwasStage2" %in% rownames(installed.packages()))
suppressPackageStartupMessages(library(fungwasStage2))

GUARDS <- "/user/work/fh6520/wt_pooled-guards/sim_meanvar/src/dispersion-guards.R"
if (!file.exists(GUARDS)) stop("dispersion-guards.R not present yet: ", GUARDS)
source(GUARDS)
stopifnot(is.function(estimate_rav), is.function(residualise_dispersion))

args   <- commandArgs(trailingOnly = TRUE)
indir  <- args[1]
tag    <- args[2]                     # e.g. raw_chunk0000
outdir <- args[3]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, paste0(tag, ".rds"))
if (file.exists(outfile)) { cat("exists:", outfile, "\n"); quit(status = 0) }

TAUS <- seq(0.02, 0.98, length.out = 45); Tt <- 45L

meta <- paste(readLines(file.path(indir, paste0(tag, ".json")), warn = FALSE), collapse = "")
gi <- function(k) as.integer(sub(paste0(".*\"", k, "\": *([0-9]+).*"), "\\1", meta))
R <- gi("reps"); M <- gi("M"); N <- gi("N")
cls <- as.integer(strsplit(gsub(".*\"cls\": *\\[([^]]*)\\].*", "\\1", meta), ", *")[[1]])
stopifnot(length(cls) == M)

rd <- function(k, dims) {
  n <- prod(dims)
  v <- readBin(file.path(indir, paste0(tag, ".", k, ".bin")), "double", n = n,
               size = 8, endian = "little")
  stopifnot(length(v) == n)
  array(aperm(array(v, rev(dims))), dims)   # C order -> R order
}
B    <- rd("beta",      c(R, M, Tt))
BT   <- rd("btrue",     c(R, M, Tt))
D    <- rd("dresid",    c(R, M))
S0   <- rd("sigma0",    c(R, Tt, Tt))
S0P  <- rd("sigma0pop", c(R, Tt, Tt))
QT   <- rd("qtau",      c(R, Tt))
FH   <- rd("fhat",      c(R, Tt))
FP   <- rd("fpop",      c(R, Tt))
MAF  <- rd("maf",       c(R, M))
DEL  <- rd("delta",     c(R, M))
GAM  <- rd("gamma",     c(R, M))
MOM  <- rd("mom",       c(R, 5L))

ests <- c("ols", "gls")
blk  <- rep(1:20, each = M / 20)          # 20 blocks for the r_av jackknife SE

nrow_out <- R * M
mk <- function() matrix(NA_real_, nrow_out, length(ests),
                        dimnames = list(NULL, ests))
out <- list(beta_mu = mk(), se_mu = mk(), beta_s2 = mk(), se_s2 = mk(),
            cov_ms = mk(), true_mu = mk(), true_s2 = mk(),
            d_res = mk(), se_d = mk(), se_d_rav = mk(),
            true_d = mk(), r_av = mk(), se_rav = mk())
rep_id <- rep(seq_len(R), each = M)
var_id <- rep(seq_len(M), times = R)
verify_err <- numeric(0)

for (r in seq_len(R)) {
  S0r <- matrix(S0[r, , ], Tt, Tt)
  S0pr <- matrix(S0P[r, , ], Tt, Tt)
  d <- D[r, ]
  Bm <- matrix(B[r, , ], M, Tt)
  Btm <- matrix(BT[r, , ], M, Tt)

  Sl <- list(S0r); names(Sl) <- "1"
  cd <- fungwasStage2:::.pooled_cov_data(Sl, chr = rep(1L, M), d_resid = d, n_taus = Tt)
  cov_data <- list(cov_vec = cd$cov_vec, offsets = cd$offsets)

  SE_tau <- outer(sqrt(diag(S0r)), 1 / sqrt(d))          # Tt x M
  Qsl <- t(Bm); colnames(Qsl) <- paste0("v", seq_len(M))
  rownames(Qsl) <- paste0("tau", TAUS); rownames(SE_tau) <- rownames(Qsl)
  colnames(SE_tau) <- colnames(Qsl)

  W <- make_weights_vqtl(TAUS, QT[r, ], mu = MOM[r, 1], sd = MOM[r, 2])
  stage1 <- list(taus = TAUS, q_tau = QT[r, ], Q_slope = Qsl, SE_tau = SE_tau)

  A <- solve(crossprod(W), t(W))
  Si <- solve(S0pr)
  Agls <- solve(t(W) %*% Si %*% W, t(W) %*% Si)
  true_ols <- A %*% t(Btm)                 # 2 x M
  true_gls <- Agls %*% t(Btm)

  for (est in ests) {
    fit <- param_gwas(stage1, transform = "custom_W", transform_args = list(W = W),
                      se_mode = "tau_cov", cov_data = cov_data,
                      tau_cov_estimator = est, nu = N,
                      nu_source = "pooled_sigma0_n_samples", return_cov = "upper")
    if (r == 1L && est == "ols") {
      verify_err <- c(verify_err, max(abs(fit$params - A %*% Qsl)))
    }
    cn <- colnames(fit$param_cov_upper)
    stopifnot("cov_beta_mu_beta_sigma2" %in% cn)
    cm <- as.numeric(fit$param_cov_upper[, "cov_beta_mu_beta_sigma2"])
    bmu <- as.numeric(fit$params["beta_mu", ])
    bs2 <- as.numeric(fit$params["beta_sigma2", ])
    smu <- as.numeric(fit$SE_params["beta_mu", ])
    ss2 <- as.numeric(fit$SE_params["beta_sigma2", ])

    rv <- estimate_rav(bmu, bs2, smu, ss2, cov_mu_sigma2 = cm, block = blk)
    rr <- residualise_dispersion(bmu, bs2, smu, ss2, cov_mu_sigma2 = cm,
                                 r_av = rv$r_av, se_rav = 0)
    rr2 <- residualise_dispersion(bmu, bs2, smu, ss2, cov_mu_sigma2 = cm,
                                  r_av = rv$r_av,
                                  se_rav = if (is.finite(rv$se)) rv$se else 0)
    ix <- (r - 1L) * M + seq_len(M)
    tt <- if (est == "ols") true_ols else true_gls
    out$beta_mu[ix, est] <- bmu; out$se_mu[ix, est] <- smu
    out$beta_s2[ix, est] <- bs2; out$se_s2[ix, est] <- ss2
    out$cov_ms[ix, est]  <- cm
    out$true_mu[ix, est] <- tt[1, ]; out$true_s2[ix, est] <- tt[2, ]
    out$d_res[ix, est]   <- rr$beta_sigma2_resid
    out$se_d[ix, est]    <- rr$se_sigma2_resid
    out$se_d_rav[ix, est] <- rr2$se_sigma2_resid
    out$true_d[ix, est]  <- tt[2, ] - rv$r_av * tt[1, ]
    out$r_av[ix, est]    <- rv$r_av
    out$se_rav[ix, est]  <- rv$se
  }
  if (r %% 5 == 0) cat("  rep", r, "/", R, "\n")
}

saveRDS(list(tag = tag, N = N, M = M, R = R, taus = TAUS, cls = cls,
             rep_id = rep_id, var_id = var_id,
             maf = as.numeric(t(MAF)), delta = as.numeric(t(DEL)),
             gamma = as.numeric(t(GAM)),
             mom = MOM, fhat = FH, fpop = FP, out = out,
             verify_ols_max_err = verify_err),
        outfile)
cat("WROTE", outfile, " verify_ols_max_err =", verify_err, "\n")
