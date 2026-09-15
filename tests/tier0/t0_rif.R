#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Tier 0 -- RIF construction against the analytic influence function.
#
# Everything here has a closed form. No simulation tolerance is tuned: each
# assertion is either exact in floating point or has a stated analytic bound.
#
# Run:  Rscript --vanilla tests/tier0/t0_rif.R
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  .self <- local({
    a <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (!length(f)) stop("Run as: Rscript --vanilla tests/tier0/t0_rif.R")
    normalizePath(f[1L])
  })
  source(file.path(dirname(.self), "..", "lib", "harness.R"))
}))
gmoq_init_roots(.self)
gmoq_load_stage2()
gmoq_begin("tier0", "rif")

build_rif <- gmoq_internal(".build_rif_matrix")
kde_f     <- gmoq_internal(".kde_f_at_q")

set.seed(20260914)

# --- 1. the definition itself ----------------------------------------------
# RIF_i(tau) = q_tau + (tau - 1{Y_i <= q_tau}) / f(q_tau), with q_tau the
# type-8 sample quantile and f the estimator's own density. Rebuild it by hand
# from the pieces the function returns and demand exact agreement.

y <- rnorm(4000, mean = 2, sd = 1.5)
taus <- c(0.02, 0.1, 0.25, 0.5, 0.75, 0.9, 0.98)
r <- build_rif(y, taus)

q_hand <- as.numeric(stats::quantile(y, probs = taus, type = 8))
gmoq_close("T0.RIF.quantile_type8",
           "q_tau is the type-8 sample quantile, the estimator documented in the method.",
           r$q_tau, q_hand, 0)

R_hand <- outer(rep(1, length(y)), r$q_tau) +
  t((taus - t(outer(y, r$q_tau, "<=") * 1)) / r$fhat_tau)
gmoq_close("T0.RIF.definition",
           "The RIF matrix equals q_tau + (tau - 1{Y<=q_tau})/f(q_tau) exactly.",
           as.numeric(r$R), as.numeric(R_hand), 1e-12)

# --- 2. the influence function is a two-valued step -------------------------
# IF(y; q_tau) = (tau - 1{y <= q_tau}) / f(q_tau) takes exactly two values,
# tau/f and (tau - 1)/f, a distance 1/f apart. Anything else means the
# indicator, the density or the recentering has been altered.
ok_two <- TRUE; ok_gap <- TRUE
for (t in seq_along(taus)) {
  v <- sort(unique(round(r$R[, t] - r$q_tau[t], 12)))
  ok_two <- ok_two && length(v) == 2L
  ok_gap <- ok_gap && abs(diff(v) - 1 / r$fhat_tau[t]) < 1e-9
}
gmoq_check("T0.RIF.two_valued",
           "The influence function is the two-valued step the analytic form predicts.",
           ok_two && ok_gap,
           criterion = "exactly two values per tau, gap = 1/f(q_tau)")

# --- 3. the recentering is exact at the empirical CDF -----------------------
# mean_i RIF_i(tau) = q_tau + (tau - Fhat(q_tau)) / f(q_tau). This is the
# "recentered" in RIF: the functional itself is recovered when Fhat(q) = tau.
Fhat <- vapply(seq_along(taus), function(t) mean(y <= r$q_tau[t]), numeric(1))
gmoq_close("T0.RIF.recentering",
           "The column mean is q_tau + (tau - Fhat(q_tau))/f, so RIF recovers the functional.",
           colMeans(r$R), r$q_tau + (taus - Fhat) / r$fhat_tau, 1e-10)

# --- 4. location-scale equivariance ----------------------------------------
# q_tau(aY + b) = a q_tau(Y) + b and f(a q + b; aY + b) = f(q; Y)/a, so the RIF
# must satisfy RIF(aY + b) = a RIF(Y) + b exactly. This is the property that
# makes b/se invariant to the density estimate (docs/BP2 section 1); if it ever
# stops holding, the bandwidth rule has acquired a scale dependence.
a <- 3.25; b <- -7.5
r2 <- build_rif(a * y + b, taus)
gmoq_close("T0.RIF.equivariance",
           "RIF(aY+b) = a RIF(Y) + b exactly, which is why b/se is free of the density estimate.",
           as.numeric(r2$R), as.numeric(a * r$R + b), 1e-9)
gmoq_close("T0.RIF.density_equivariance",
           "f_hat(aY+b) = f_hat(Y)/a exactly.",
           r2$fhat_tau * a, r$fhat_tau, 1e-12)

# --- 5. the kernel density estimator, against its own closed form ----------
# .kde_f_at_q is a Gaussian KDE at Silverman's bw.nrd0. Written out by hand.
yq <- r$q_tau[4]
bw <- stats::bw.nrd0(y)
gmoq_close("T0.RIF.kde_closed_form",
           "The density estimate is the Gaussian KDE at bw.nrd0 the code claims.",
           kde_f(y, yq), mean(dnorm((y - yq) / bw)) / bw, 1e-14)

# NOTE, recorded rather than asserted: the R path uses stats::bw.nrd0
# (0.9 min(sd, IQR/1.349) N^-1/5) while the shipped Stage 1 uses
# 1.06 sd N^-1/5 (stage1/fungwas_stage1/core.py:308). They are different
# bandwidths. Tier 1 measures the consequence; by section 4 above it cancels
# from every per-tau Wald test but not from Q.
h_r  <- stats::bw.nrd0(y)
h_py <- 1.06 * stats::sd(y) * length(y)^(-0.2)
gmoq_check("T0.RIF.bandwidth_rules_differ",
           "Recorded, not a pass/fail: the R quantile_gwas() bandwidth is not the Stage 1 one.",
           TRUE, stat = h_r / h_py,
           criterion = "reported ratio bw.nrd0 / Silverman-1.06",
           detail = "R uses bw.nrd0; Stage 1 core.py uses 1.06*sd*N^-0.2")

# --- 6. the density floor is honoured --------------------------------------
r3 <- build_rif(y, taus, density_floor = 1e6)
gmoq_check("T0.RIF.density_floor",
           "density_floor is a lower bound on f_hat, so the RIF cannot divide by ~0.",
           all(r3$fhat_tau == 1e6),
           criterion = "all f_hat == floor when the floor dominates")

# --- 7. population targets --------------------------------------------------
# E[RIF(tau)] = q_tau and Var[RIF(tau)] = tau(1-tau)/f(q_tau)^2 for the true
# quantile and density. Checked on a large standard normal sample against the
# ANALYTIC q and f, with a Monte Carlo tolerance stated from the variance.
set.seed(11)
n <- 400000
yz <- rnorm(n)
rz <- build_rif(yz, taus)
se_mean <- sqrt(taus * (1 - taus)) / rz$fhat_tau / sqrt(n)
gmoq_check("T0.RIF.population_mean",
           "The RIF is unbiased for the quantile functional: E[RIF(tau)] = q_tau.",
           all(abs(colMeans(rz$R) - stats::qnorm(taus)) < 5 * se_mean),
           stat = max(abs(colMeans(rz$R) - stats::qnorm(taus)) / se_mean),
           criterion = "|mean - qnorm(tau)| < 5 analytic SE, all taus")

# Var[RIF(tau)] = Fhat(q)(1 - Fhat(q)) / f_hat(q)^2 is an EXACT sample identity
# for the estimator's own quantile and density, so it is checked to machine
# precision rather than against a Monte Carlo tolerance.
Fz <- vapply(seq_along(taus), function(t) mean(yz <= rz$q_tau[t]), numeric(1))
var_id <- (n / (n - 1)) * Fz * (1 - Fz) / rz$fhat_tau^2
gmoq_close("T0.RIF.variance_identity",
           "Var[RIF(tau)] = Fhat(1-Fhat)/f_hat^2 exactly, the influence-function variance.",
           apply(rz$R, 2, stats::var), var_id, 1e-8)

# --- 8. the cross-tau covariance is the Brownian bridge --------------------
# Cov(RIF(s), RIF(t)) = (Fhat(min) - Fhat_s Fhat_t) / (f_s f_t) exactly in
# sample, and (min(s,t) - st)/(f_s f_t) in population. This is the structure
# docs/KNOWN_ISSUES.md section 3 relies on, and it is why Q is a correlated
# quadratic form rather than a chi-square.
S_obs <- stats::cov(rz$R)
Fmin <- outer(seq_along(taus), seq_along(taus),
              Vectorize(function(i, j) mean(yz <= min(rz$q_tau[i], rz$q_tau[j]))))
S_id <- (n / (n - 1)) * (Fmin - outer(Fz, Fz)) / outer(rz$fhat_tau, rz$fhat_tau)
gmoq_close("T0.RIF.bridge_identity",
           "The cross-tau covariance is (Fhat(min) - Fhat_s Fhat_t)/(f_s f_t) exactly in sample.",
           as.numeric(S_obs), as.numeric(S_id), 1e-8)

S_pop <- (outer(taus, taus, pmin) - outer(taus, taus)) /
  outer(dnorm(stats::qnorm(taus)), dnorm(stats::qnorm(taus)))
rel <- max(abs(S_obs / S_pop - 1))
gmoq_check("T0.RIF.brownian_bridge_population",
           "The sample covariance converges on the analytic Brownian-bridge form.",
           rel < 0.10, stat = rel,
           criterion = "max relative deviation < 0.10 at N = 4e5",
           detail = "the residual is KDE bias in f_hat at the extreme taus, measured next")

# --- 9. how big the density error itself is, at this N ---------------------
# Reported, not asserted. The bandwidth is Silverman's, and the estimate is
# biased at the tails and at any density curvature; docs/BP2_lumpy_phenotypes.md
# measures up to 2.0x on a bimodal trough. It cancels from b/se (section 4) but
# not from Q.
fd <- max(abs(rz$fhat_tau / dnorm(stats::qnorm(taus)) - 1))
gmoq_check("T0.RIF.kde_bias_at_N4e5",
           "Size of the KDE's own density error on a standard normal at N = 4e5 -- reported.",
           TRUE, stat = fd, severity = "report",
           criterion = "max |f_hat/f_true - 1| over the tau grid")

gmoq_finish()
