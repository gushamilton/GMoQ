#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Tier 1 -- the refusal paths.
#
# Every one of these is a failure that leaves the point estimates untouched and
# the standard errors plausible, so nothing downstream can detect it. The
# package's stated position is that it stops rather than proceeds. This file is
# what makes that position testable: a refusal that quietly becomes a warning is
# a regression, and so is a refusal that fires for a different reason than the
# one claimed (which is why every check matches the message, not just the fact
# of an error).
#
# Also here: the two established failure modes that are NOT refusals and must
# stay visible -- the tie-induced GLS singularity, and the non-PSD covariance
# that GLS currently accepts silently.
#
# Run:  Rscript --vanilla tests/tier1/t1_refusals.R
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  .self <- local({
    a <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (!length(f)) stop("Run as: Rscript --vanilla tests/tier1/t1_refusals.R")
    normalizePath(f[1L])
  })
  source(file.path(dirname(.self), "..", "lib", "harness.R"))
}))
gmoq_init_roots(.self)
gmoq_load_stage2()
suppressPackageStartupMessages(library(data.table))
gmoq_begin("tier1", "refusals")

FIX <- Sys.getenv("GMOQ_TIER1_DIR", file.path(Sys.getenv("GMOQ_WORK"), "tier1"))
meta <- jsonlite::fromJSON(file.path(FIX, "fixtures.json"))
taus <- as.numeric(meta$taus); Tt <- length(taus)
W <- cbind(beta_mu = rep(1, Tt), beta_sigma = qnorm(taus))

pfx <- function(dir) file.path(FIX, dir, "sim")
jk <- function(dir, ...) suppressMessages(param_gwas_from_file(
  stage1_file = paste0(pfx(dir), ".stage1.tsv.gz"), W = W,
  se_mode = "tau_cov", cov_source = "jackknife",
  cov_file = paste0(pfx(dir), ".cov.gz"), q_method = "saddlepoint", ...))
pl <- function(dir, ...) suppressMessages(param_gwas_from_file(
  stage1_file = paste0(pfx(dir), ".stage1.tsv.gz"), W = W,
  se_mode = "tau_cov", cov_source = "pooled",
  sigma0_file = paste0(pfx(dir), ".sigma0.tsv.gz"),
  cov_meta_file = paste0(pfx(dir), ".cov.meta.json"),
  q_method = "saddlepoint", ...))

# =====================================================================
# 1. Mismatched sidecars
# =====================================================================
gmoq_refuses("T1.REF.cov_ids_mismatch",
             "A covariance file whose rows do not line up with the Stage 1 table is refused.",
             jk("cov_ids_shuffled"), "Covariance IDs mismatch at row")

gmoq_refuses("T1.REF.cov_ids_missing",
             "A missing covariance IDs sidecar is refused: without it, alignment cannot be checked.",
             jk("no_cov_ids"), "Missing covariance IDs file")

ok_unsafe <- tryCatch({
  r <- withCallingHandlers(jk("no_cov_ids", unsafe_skip_cov_ids = TRUE),
                           warning = function(w) invokeRestart("muffleWarning"))
  nrow(r) > 0
}, error = function(e) FALSE)
gmoq_check("T1.REF.cov_ids_override_is_explicit",
           "The alignment check can only be skipped by naming unsafe_skip_cov_ids, and it warns.",
           ok_unsafe, criterion = "unsafe_skip_cov_ids = TRUE proceeds, with a warning")

gmoq_refuses("T1.REF.tau_grid_mismatch",
             "A Sigma_0 built on a different tau grid is refused: Sigma_0 is indexed by tau.",
             pl("sigma0_wrong_taus"), "the tau grid differs")

gmoq_refuses("T1.REF.covariate_set_mismatch",
             "A Sigma_0 built on a different covariate set is refused (measured cost: +25.6% on weight).",
             pl("sigma0_wrong_covariates"), "the covariate set differs")

gmoq_refuses("T1.REF.sample_count_mismatch",
             "A Sigma_0 built on a different sample is refused.",
             pl("sigma0_wrong_n"), "the analysed sample count differs")

# THE LOCO SCALE PIN. docs/LOCO_SCALE.md: subtracting REGENIE's standardised
# LOCO prediction without rescaling attenuated log_bmi shape effects 2.2x and
# inflated their SEs 2.2x, invisibly. A Sigma_0 built under a different
# --loco-scale than the scan is the same class of silent error and must stop.
gmoq_refuses("T1.REF.loco_scale_mismatch",
             "A Sigma_0 built under a different LOCO scaling than the scan is refused.",
             pl("sigma0_wrong_loco_scale"), "the LOCO scale mode differs")

gmoq_refuses("T1.REF.no_analysis_block",
             "Covariance metadata with no analysis block cannot vouch for a Sigma_0, so it is refused.",
             pl("no_cov_meta"), "there is no Stage 1 covariance metadata|carries no 'analysis' block")

gmoq_refuses("T1.REF.missing_d_resid",
             "Without d_resid there is no per-variant scale, so the pooled covariance is refused.",
             pl("no_d_resid"), "has no d_resid column")

gmoq_refuses("T1.REF.future_schema",
             "A sidecar schema newer than this package understands is refused, not guessed at.",
             jk("future_schema"), "does not.*understand|schema version")

# =====================================================================
# 2. Duplicated / multi-allelic variants: dropped, loudly
# =====================================================================
res_dup <- jk("duplicate_snp_ids", unsafe_skip_cov_ids = TRUE)
n_expected <- meta$n_null + meta$n_causal - 2L      # both copies of the pair go
gmoq_check("T1.REF.duplicate_variants_dropped",
           "Duplicated rsIDs are dropped before Stage 2 rather than analysed twice.",
           nrow(res_dup) == n_expected,
           stat = nrow(res_dup),
           criterion = sprintf("%d rows remain of %d", n_expected, meta$n_null + meta$n_causal))
msgs <- character(0)
invisible(withCallingHandlers(
  suppressWarnings(param_gwas_from_file(
    stage1_file = paste0(pfx("duplicate_snp_ids"), ".stage1.tsv.gz"), W = W,
    se_mode = "tau_cov", cov_source = "jackknife",
    cov_file = paste0(pfx("duplicate_snp_ids"), ".cov.gz"),
    q_method = "saddlepoint", unsafe_skip_cov_ids = TRUE)),
  message = function(m) { msgs <<- c(msgs, conditionMessage(m))
                          invokeRestart("muffleMessage") }))
gmoq_check("T1.REF.duplicate_drop_is_announced",
           "The drop is announced, so a silent change in variant count cannot go unnoticed.",
           any(grepl("Dropping .* Stage 1 rows", msgs)),
           criterion = "a message naming the number of dropped rows",
           detail = paste(head(msgs, 4), collapse = " | "))

# =====================================================================
# 3. A degenerate tau grid -- the tie / GLS singularity
#
# docs/BP2_lumpy_phenotypes.md: when a block of tied phenotype values spans two
# requested taus they share a sample quantile, their RIF columns become the same
# vector, Sigma goes rank-deficient (condition number 1e2 -> 1e13) and the GLS
# solve falls back to ridge repair on 100% of SNPs, inflating SEs by 10^2-10^3.
# OLS is immune, because its contrast never inverts Sigma.
#
# There is NO guard in the shipped code for this. The pin below is that the
# damage is still detectable and still confined to GLS; if a guard is added
# later, T1.TIE.gls_ridge_fires is the check that will flag the change.
# =====================================================================
s1_tie <- data.table::fread(paste0(pfx("degenerate_taus"), ".stage1.tsv.gz"))
n_tie_rows <- nrow(s1_tie)

cov_raw <- local({
  rd <- gmoq_internal(".read_tau_cov_data")
  rd(paste0(pfx("degenerate_taus"), ".cov.gz"), Tt, n_tie_rows)
})
n_elem <- Tt * (Tt + 1L) / 2L
unpack <- function(k) {
  v <- cov_raw$cov_vec[(k - 1L) * n_elem + seq_len(n_elem)]
  S <- matrix(0, Tt, Tt); i <- 1L
  for (r in seq_len(Tt)) for (cc in r:Tt) { S[r, cc] <- S[cc, r] <- v[i]; i <- i + 1L }
  S
}
conds <- vapply(seq_len(min(30L, n_tie_rows)), function(k) kappa(unpack(k), exact = TRUE), 0)
gmoq_check("T1.TIE.covariance_is_ill_conditioned",
           "A tied phenotype makes the tau-covariance numerically singular, as documented.",
           median(conds) > 1e8, stat = median(conds),
           criterion = "median condition number over 30 variants > 1e8")

beta_tie <- t(as.matrix(s1_tie[, grep("^beta_tau_", names(s1_tie)), with = FALSE]))
se_tie   <- t(as.matrix(s1_tie[, grep("^se_tau_",   names(s1_tie)), with = FALSE]))
st_tie <- list(taus = taus, q_tau = rep(0, Tt), Q_slope = beta_tie, SE_tau = se_tie)

fit_gls <- suppressMessages(param_gwas(
  st_tie, transform = "custom_W", transform_args = list(W = W),
  se_mode = "tau_cov", tau_cov_estimator = "gls", cov_data = cov_raw,
  q_method = "chisq_naive", nu = 199, nu_source = "test"))
fit_ols <- suppressMessages(param_gwas(
  st_tie, transform = "custom_W", transform_args = list(W = W),
  se_mode = "tau_cov", tau_cov_estimator = "ols", cov_data = cov_raw,
  q_method = "chisq_naive", nu = 199, nu_source = "test"))

gmoq_check("T1.TIE.gls_ridge_fires",
           "The GLS solve needs ridge repair on a tied phenotype, and says so in its diagnostics.",
           !is.null(fit_gls$gls_diag) && fit_gls$gls_diag$n_ridge_tried > 0,
           stat = if (is.null(fit_gls$gls_diag)) NA else
             fit_gls$gls_diag$n_ridge_tried / max(fit_gls$gls_diag$n_snps, 1),
           criterion = "ridge repair attempted on a non-zero fraction of variants")
gmoq_check("T1.TIE.ols_is_unaffected",
           "OLS is immune: its contrast never inverts Sigma, so its SEs stay finite and sane.",
           all(is.finite(fit_ols$SE_params)) && all(fit_ols$SE_params > 0),
           criterion = "all OLS SEs finite and positive on the tied fixture")
gmoq_check("T1.TIE.gls_se_inflation",
           "How much the tie inflates the GLS standard error relative to OLS -- reported.",
           TRUE, severity = "report",
           stat = median(fit_gls$SE_params[2, ] / fit_ols$SE_params[2, ]),
           criterion = "median GLS/OLS SE ratio on the scale parameter, tied phenotype")
gmoq_check("T1.TIE.no_shipped_tie_guard",
           "Known gap: nothing in the shipped code refuses or warns on colliding tau quantiles.",
           TRUE, severity = "report", stat = sum(duplicated(round(s1_tie$d_resid, 12))) * 0,
           criterion = "recorded in docs/TEST_SUITE.md; BP2 section 5.5 proposes the guard")

# =====================================================================
# 4. A non-PSD covariance
#
# docs/OLS_VS_GLS.md: "the GLS path has no guard against a numerically singular
# Sigma -- it emits a z-score of 498 with every diagnostic counter clean. That
# must be fixed before GLS is recommended for production use."
#
# GLS is NOT the default, so this is recorded as advisory rather than blocking.
# What IS blocking is that the default path does not silently produce an
# impossible z.
# =====================================================================
s1_np <- data.table::fread(paste0(pfx("non_psd"), ".stage1.tsv.gz"))
cov_np <- local({
  rd <- gmoq_internal(".read_tau_cov_data")
  rd(paste0(pfx("non_psd"), ".cov.gz"), Tt, nrow(s1_np))
})
b_np  <- t(as.matrix(s1_np[, grep("^beta_tau_", names(s1_np)), with = FALSE]))
se_np <- t(as.matrix(s1_np[, grep("^se_tau_",   names(s1_np)), with = FALSE]))
st_np <- list(taus = taus, q_tau = rep(0, Tt), Q_slope = b_np, SE_tau = se_np)

np_gls <- suppressWarnings(suppressMessages(tryCatch(param_gwas(
  st_np, transform = "custom_W", transform_args = list(W = W),
  se_mode = "tau_cov", tau_cov_estimator = "gls", cov_data = cov_np,
  q_method = "chisq_naive", nu = 199, nu_source = "test"),
  error = function(e) NULL)))
np_ols <- suppressWarnings(suppressMessages(param_gwas(
  st_np, transform = "custom_W", transform_args = list(W = W),
  se_mode = "tau_cov", tau_cov_estimator = "ols", cov_data = cov_np,
  q_method = "chisq_naive", nu = 199, nu_source = "test")))

zmax_gls <- if (is.null(np_gls)) NA_real_ else
  max(abs(np_gls$params / np_gls$SE_params), na.rm = TRUE)
gmoq_check("T1.NPSD.gls_behaviour",
           "What GLS does with a covariance that is not positive semi-definite -- measured.",
           TRUE, severity = "report", stat = zmax_gls,
           criterion = "max |theta/SE| produced from a deliberately indefinite Sigma",
           detail = "docs/OLS_VS_GLS.md records an unguarded z of 498 here; GLS is not the default")
gmoq_check("T1.NPSD.gls_unguarded",
           "Advisory: GLS still accepts an indefinite Sigma without flagging it.",
           is.null(np_gls) || !is.finite(zmax_gls) ||
             (!is.null(np_gls$gls_diag) && np_gls$gls_diag$n_sigma_not_spd > 0),
           severity = "advisory",
           criterion = "either refuse, return non-finite, or count it in n_sigma_not_spd")
gmoq_check("T1.NPSD.ols_does_not_invent_a_z",
           "The default OLS path does not turn an indefinite Sigma into a plausible-looking z.",
           {
             zo <- abs(np_ols$params / np_ols$SE_params)
             all(!is.finite(zo) | zo < 100)
           },
           stat = suppressWarnings(max(abs(np_ols$params / np_ols$SE_params), na.rm = TRUE)),
           criterion = "no finite |z| above 100 from the default estimator")

# =====================================================================
# 5. Refusals that need no fixture
# =====================================================================
gmoq_refuses("T1.REF.tau_cov_needs_cov_data",
             "se_mode='tau_cov' without covariance data is refused rather than silently downgraded.",
             param_gwas(list(taus = taus, q_tau = rep(0, Tt),
                             Q_slope = matrix(0, Tt, 2), SE_tau = matrix(1, Tt, 2)),
                        transform = "custom_W", transform_args = list(W = W),
                        se_mode = "tau_cov", cov_data = NULL),
             "requires cov_data")
missing_dir <- file.path(FIX, "missing_sigma0")
if (dir.exists(missing_dir)) unlink(missing_dir, recursive = TRUE)
dir.create(missing_dir, recursive = TRUE)
file.copy(list.files(file.path(FIX, "clean"), full.names = TRUE),
          missing_dir, overwrite = TRUE)
unlink(list.files(missing_dir, pattern = "\\.sigma0", full.names = TRUE))
gmoq_refuses("T1.REF.pooled_needs_sigma0",
             "cov_source='pooled' without a sigma0 sidecar is refused.",
             suppressMessages(param_gwas_from_file(
               stage1_file = file.path(missing_dir, "sim.stage1.tsv.gz"), W = W,
               se_mode = "tau_cov", cov_source = "pooled", sigma0_file = NULL)),
             "requires.*sigma0")
gmoq_refuses("T1.REF.jackknife_needs_cov_file",
             "cov_source='jackknife' without a covariance file is refused.",
             suppressMessages(param_gwas_from_file(
               stage1_file = paste0(pfx("clean"), ".stage1.tsv.gz"), W = W,
               se_mode = "tau_cov", cov_source = "jackknife", cov_file = NULL)),
             "requires cov_file")
gmoq_refuses("T1.REF.W_wrong_nrow",
             "A weight matrix with the wrong number of rows is refused, not recycled.",
             param_gwas(list(taus = taus, q_tau = rep(0, Tt),
                             Q_slope = matrix(0, Tt, 2), SE_tau = matrix(1, Tt, 2)),
                        transform = "custom_W",
                        transform_args = list(W = W[1:3, , drop = FALSE]),
                        se_mode = "diagonal"),
             "nrow\\(W\\) must equal length\\(taus\\)")

# Degenerate tau grid at the fitting stage: a W with fewer taus than parameters
gmoq_refuses("T1.REF.degenerate_W",
             "A rank-deficient weight matrix is refused by the solve rather than producing nonsense.",
             param_gwas(list(taus = taus[1:2], q_tau = rep(0, 2),
                             Q_slope = matrix(0, 2, 2), SE_tau = matrix(1, 2, 2)),
                        transform = "custom_W",
                        transform_args = list(W = cbind(a = c(1, 1), b = c(2, 2))),
                        se_mode = "diagonal"),
             "singular|computationally singular|system is exactly singular|Lapack")

gmoq_finish()
