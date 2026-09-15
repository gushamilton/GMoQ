#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Tier 1 -- the eight model families, and both halves of the tournament claim.
#
# G-MoQ's point is that it can tell distributional effects apart. A tournament
# that always picks the truth is only half a test: the other half is that the
# models which CANNOT represent the truth are actually rejected. Both are here:
#
#   (a) RECOVERY   -- the parsimonious true model wins on BIC, and its own
#                     lack-of-fit statistic Q is null-sized.
#   (b) REJECTION  -- every model whose column space does not contain the truth
#                     is rejected by Q at the designed effect size.
#
# Also checked per family: the fitted theta under the true model recovers the
# exact population projection A b(tau).
#
# The truth map (which model is true, which are merely over-parameterised, and
# which are wrong) is tests/lib/models.R, so it is stated once and used by
# tiers 1 and 2 alike.
#
# Run:  Rscript --vanilla tests/tier1/t1_families.R
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  .self <- local({
    a <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (!length(f)) stop("Run as: Rscript --vanilla tests/tier1/t1_families.R")
    normalizePath(f[1L])
  })
  source(file.path(dirname(.self), "..", "lib", "harness.R"))
}))
gmoq_init_roots(.self)
gmoq_load_stage2()
source(file.path(GMOQ_SUITE_ROOT, "lib", "models.R"))
suppressPackageStartupMessages(library(data.table))
gmoq_begin("tier1", "families")

FIX <- Sys.getenv("GMOQ_TIER1_DIR", file.path(Sys.getenv("GMOQ_WORK"), "tier1"))
idx <- jsonlite::fromJSON(file.path(FIX, "families.json"), simplifyVector = FALSE)
taus <- as.numeric(unlist(idx$taus)); Tt <- length(taus)

read_cov <- gmoq_internal(".read_tau_cov_data")

fit_one <- function(stage1_dt, cov_data, W, nu) {
  bmat <- t(as.matrix(stage1_dt[, grep("^beta_tau_", names(stage1_dt)), with = FALSE]))
  smat <- t(as.matrix(stage1_dt[, grep("^se_tau_",   names(stage1_dt)), with = FALSE]))
  st <- list(taus = taus, q_tau = rep(0, Tt), Q_slope = bmat, SE_tau = smat)
  suppressMessages(param_gwas(st, transform = "custom_W", transform_args = list(W = W),
                              se_mode = "tau_cov", tau_cov_estimator = "ols",
                              cov_data = cov_data, q_method = "saddlepoint",
                              nu = nu, nu_source = "tier1_fixed"))
}

summary_rows <- list()

for (fam in names(idx$families)) {
  info <- idx$families[[fam]]
  d <- file.path(FIX, info$dir)
  s1 <- data.table::fread(file.path(d, "sim.stage1.tsv.gz"))
  cd <- read_cov(file.path(d, "sim.cov.gz"), Tt, nrow(s1))
  b_true <- as.matrix(data.table::fread(file.path(d, "b_true.tsv"), header = FALSE))
  causal <- unlist(info$causal)
  nulls  <- setdiff(seq_len(nrow(s1)), causal)

  mixarg <- NULL
  if (!is.null(info$mix)) {
    mixarg <- list(q_tau = as.numeric(unlist(info$q_pop)),
                   p1 = info$mix$p1, mu1 = info$mix$mu1, sd1 = info$mix$sd1,
                   mu2 = info$mix$mu2, sd2 = info$mix$sd2)
  }
  models <- gmoq_model_set(taus, sigma = info$y_sd, mix = mixarg)
  map <- GMOQ_TRUTH_MAP[[fam]]
  wrong_shape <- if (is.null(map$wrong_shape)) character(0) else map$wrong_shape
  cand <- intersect(c(map$true, map$spanning, map$wrong, wrong_shape), names(models))

  fits <- lapply(models[cand], function(W) fit_one(s1, cd, W, nu = idx$blocks - 1))
  names(fits) <- cand

  # --- (0) the SUITE'S OWN truth machinery, checked against Stage 1 --------
  # b_true is computed here from the conditional CDF, not by simulation. If it
  # disagreed with what the shipped Stage 1 actually estimates, every coverage
  # number in tiers 1-3 would be measuring the wrong thing. This is the check
  # that makes the rest of the suite mean anything.
  if (length(causal)) {
    bh <- t(as.matrix(s1[, grep("^beta_tau_", names(s1)), with = FALSE]))
    sh <- t(as.matrix(s1[, grep("^se_tau_",   names(s1)), with = FALSE]))
    zt <- (bh[, causal, drop = FALSE] - t(b_true)[, causal, drop = FALSE]) /
      sh[, causal, drop = FALSE]
    gmoq_check(sprintf("T1.SIM.%s.exact_estimand_matches_stage1", fam),
               sprintf("The exact population estimand b(tau) computed for the %s family agrees, per tau, with what the shipped Stage 1 estimates.", fam),
               max(abs(zt)) < 5, stat = max(abs(zt)),
               criterion = "max |b_hat(tau) - b_true(tau)| / se_tau < 5 over causal variants x taus")
  }

  # --- (a) recovery of theta under the true model --------------------------
  Wt <- models[[map$true]]
  At <- solve(crossprod(Wt), t(Wt))
  th_pop <- At %*% t(b_true)
  ft <- fits[[map$true]]
  z <- (ft$params - th_pop) / ft$SE_params
  if (length(causal)) {
    gmoq_check(sprintf("T1.FAM.%s.recovers_theta_pop", fam),
               sprintf("Under the true model (%s) the %s family's exact population theta is recovered.",
                       map$true, fam),
               max(abs(z[, causal, drop = FALSE])) < 5.5,
               stat = max(abs(z[, causal, drop = FALSE])),
               criterion = "max |theta_hat - theta_pop| / se < 5.5 over the causal variants")
  }

  # --- (a) the true model's Q is null-sized --------------------------------
  q_true <- mean(ft$Q_pval[causal] < 0.05, na.rm = TRUE)
  gmoq_check(sprintf("T1.FAM.%s.true_model_fits", fam),
             sprintf("The true model (%s) is NOT rejected by Q on the causal variants.",
                     map$true),
             is.nan(q_true) || q_true < 0.30, stat = q_true,
             criterion = "Q rejection at 0.05 below 0.30 for the true model",
             detail = sprintf("n_causal = %d", length(causal)))

  # --- (b) the wrong models are rejected -----------------------------------
  if (length(map$wrong)) {
    rej <- vapply(map$wrong, function(m)
      mean(fits[[m]]$Q_pval[causal] < 0.05, na.rm = TRUE), numeric(1))
    gmoq_check(sprintf("T1.FAM.%s.wrong_models_rejected", fam),
               sprintf("Every model that cannot represent a %s effect is rejected by Q.", fam),
               all(rej > 0.5), stat = min(rej),
               criterion = "Q rejection at 0.05 above 0.50 for every non-spanning model",
               detail = paste(sprintf("%s=%.2f", names(rej), rej), collapse = " "))
  }

  # --- (b') the models that omit only a SHAPE dimension ---------------------
  # Measured rather than asserted. Q rejects a model missing a location or scale
  # dimension every time; it does not reject one missing only a shape dimension,
  # at any effect size where the linearisation still holds. That is a real limit
  # on what the tournament can establish and it belongs in the release report,
  # not in a pass criterion invented to hide it.
  if (length(wrong_shape)) {
    rejs <- vapply(wrong_shape, function(m)
      mean(fits[[m]]$Q_pval[causal] < 0.05, na.rm = TRUE), numeric(1))
    gmoq_check(sprintf("T1.FAM.%s.shape_only_models_rejected", fam),
               sprintf("How often Q rejects a model that omits ONLY the %s dimension -- reported.", fam),
               TRUE, severity = "report", stat = min(rejs),
               criterion = "Q rejection at 0.05 for models missing only the shape column",
               detail = paste(sprintf("%s=%.2f", names(rejs), rejs), collapse = " "))
  }

  # --- (a) selection by BIC ------------------------------------------------
  # AIC = Q + 2k and BIC = Q + k log N, as the package's tournament computes
  # them. docs/KNOWN_ISSUES.md section 4 measures the raw-Q criterion as
  # recovering the truth only 44 (AIC) / 59 (BIC) per cent of the time, so the
  # criterion here is a floor set from that, not an aspiration.
  Nn <- idx$N
  bic <- vapply(cand, function(m)
    fits[[m]]$Q + ncol(models[[m]]) * log(Nn), numeric(nrow(s1)))
  pick <- cand[apply(bic[causal, , drop = FALSE], 1, which.min)]
  acc <- mean(pick %in% c(map$true))
  acc_span <- mean(pick %in% map$spanning)
  # For the shape families the preceding block shows Q cannot separate them, so
  # BIC cannot either, and a pass criterion on selection would be dishonest.
  gmoq_check(sprintf("T1.FAM.%s.selection_is_spanning", fam),
             sprintf("BIC selects a model that CAN represent a %s effect.", fam),
             if (length(wrong_shape)) TRUE else acc_span > 0.8,
             severity = if (length(wrong_shape)) "report" else "blocking",
             stat = acc_span,
             criterion = if (length(wrong_shape))
               "reported only: Q cannot separate shape models, so neither can BIC"
             else "a spanning model selected for over 80 per cent of causal variants",
             detail = paste(sort(unique(pick)), collapse = ", "))
  gmoq_check(sprintf("T1.FAM.%s.selection_accuracy", fam),
             sprintf("How often BIC picks exactly the parsimonious true model for %s -- reported.", fam),
             TRUE, severity = "report", stat = acc,
             criterion = "docs/KNOWN_ISSUES.md section 4 measures 0.59 for BIC on an additive truth")

  # --- Q is calibrated on the null variants under every model --------------
  qn <- vapply(cand, function(m) mean(fits[[m]]$Q_pval[nulls] < 0.05, na.rm = TRUE), numeric(1))
  gmoq_check(sprintf("T1.FAM.%s.Q_null_calibrated", fam),
             "Q does not reject on null variants under any candidate model.",
             all(qn < 0.25), stat = max(qn),
             criterion = "null Q rejection at 0.05 below 0.25 under every candidate",
             detail = paste(sprintf("%s=%.2f", names(qn), qn), collapse = " "))

  summary_rows[[fam]] <- data.frame(family = fam, true_model = map$true,
                                    sel_true = acc, sel_spanning = acc_span,
                                    q_true = q_true,
                                    q_wrong_min = if (length(map$wrong))
                                      min(vapply(map$wrong, function(m)
                                        mean(fits[[m]]$Q_pval[causal] < 0.05, na.rm = TRUE),
                                        numeric(1))) else NA_real_)
}

out <- Sys.getenv("GMOQ_TEST_OUT", "")
if (nzchar(out)) {
  utils::write.table(do.call(rbind, summary_rows),
                     file.path(out, "tier1_families.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
}

gmoq_finish()
