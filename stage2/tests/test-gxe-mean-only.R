# ============================================================
# Test: GxE mean-only simulation (E unmodeled)
# Goal: see if GxE can masquerade as shifted log-normal
#
# Model:
#   Y_i = mu0 + (G_i * bG) + (G_i * E_i) * bGE + eps_i
#   eps_i ~ N(0, sd_eps^2)
# ------------------------------------------------------------
# Notes:
# - E is unmodeled in Stage 1.
# - Compare additive-only vs shifted-only vs both (shift+meanlog).
# ============================================================

suppressPackageStartupMessages({
  library(fungwasStage2)
  library(ggplot2)
})

set.seed(123)

# -------------------------
# 1) Simulate genotypes
# -------------------------
N <- 70000
P <- 200
maf <- runif(P, 0.05, 0.5)
G <- matrix(rbinom(N * P, 2, rep(maf, each = N)), nrow = N, ncol = P)
colnames(G) <- paste0("SNP", seq_len(P))

# -------------------------
# 2) Environment and effects
# -------------------------
E <- rnorm(N)
mu0 <- 0.0
sd_eps <- 1.0

# Main SNP effects (additive)
# Rough calibration to BMI-like tau-slope magnitudes
bG_sd <- 0.09
bGE_sd <- 0.28
bG <- rnorm(P, 0, bG_sd)
# GxE effects on mean (half nonzero)
bGE <- c(rnorm(P/2, 0, bGE_sd), rep(0, P/2))

# Mean model: Y = mu0 + G*bG + (G*E)*bGE + eps
GxE <- G * E  # elementwise, each column scaled by E
mu <- as.numeric(mu0 + G %*% bG + GxE %*% bGE)
Y <- mu + rnorm(N, mean = 0, sd = sd_eps)

# -------------------------
# 3) Stage 1: quantile GWAS (E omitted)
# -------------------------
taus <- seq(0.05, 0.95, 0.05)
stage1 <- quantile_gwas(
  Y, G,
  taus = taus,
  benchmark = FALSE,
  verbose = FALSE
)

# -------------------------
# 3b) Multiplicative log-normal model (comparison)
#   Y_mult = shift0_m + exp(beta0 + G*bM) * U0
# -------------------------
shift0_m <- min(Y) - 0.1 * sd(Y)
if (any(Y <= shift0_m)) {
  shift0_m <- min(Y) - 1.0
}
beta0_m <- 0.0
muU_m <- 0.5
sdU_m <- 0.6

# Use same SNP effects as bGE for direct slope-shape comparison
bM <- bGE
U0_m <- rlnorm(N, meanlog = muU_m, sdlog = sdU_m)
Y_mult <- shift0_m + exp(beta0_m + G %*% bM) * U0_m

stage1_mult <- quantile_gwas(
  Y_mult, G,
  taus = taus,
  benchmark = FALSE,
  verbose = FALSE
)

# -------------------------
# 4) Stage 2: three models
#   A) Additive-only (flat W)
#   B) "Shifted" (meanlog-only; curved tau pattern)
#   C) Both (shift + meanlog)
# -------------------------

# A) Additive-only
W_add <- matrix(1, nrow = length(taus), ncol = 1)
rownames(W_add) <- paste0("tau", taus)
colnames(W_add) <- "beta_add"

# Baseline params for shifted log-normal approximation
# Ensure Y - shift0 > 0
shift0 <- min(Y) - 0.1 * sd(Y)
if (any(Y <= shift0)) {
  shift0 <- min(Y) - 1.0
}
logY <- log(Y - shift0)
meanlog0 <- mean(logY)
sdlog0 <- sd(logY)

# CDF/PDF for shifted log-normal
shifted_cdf <- function(y, params) {
  a <- params$shift
  out <- numeric(length(y))
  ok <- (y > a)
  out[ok] <- plnorm(y[ok] - a, meanlog = params$meanlog, sdlog = params$sdlog)
  out
}

shifted_pdf <- function(y, params) {
  a <- params$shift
  out <- numeric(length(y))
  ok <- (y > a)
  out[ok] <- dlnorm(y[ok] - a, meanlog = params$meanlog, sdlog = params$sdlog)
  out
}

# Gradients of CDF
# dF/da = -f0(y-a)
grad_shift <- function(y, params) {
  a <- params$shift
  out <- numeric(length(y))
  ok <- (y > a)
  out[ok] <- -dlnorm(y[ok] - a, meanlog = params$meanlog, sdlog = params$sdlog)
  out
}

# dF/d(meanlog)
grad_meanlog <- function(y, params) {
  a <- params$shift
  out <- numeric(length(y))
  ok <- (y > a)
  yy <- y[ok] - a
  z <- (log(yy) - params$meanlog) / params$sdlog
  out[ok] <- -dnorm(z) / params$sdlog
  out
}

params0 <- list(shift = shift0, meanlog = meanlog0, sdlog = sdlog0)

# B) "Shifted" (meanlog-only; shift fixed)
W_shift <- make_weights_generic(
  taus = taus,
  q_tau = stage1$q_tau,
  dist_cdf = shifted_cdf,
  dist_pdf = shifted_pdf,
  params = params0,
  grad_funcs = list(beta_meanlog = grad_meanlog)
)

# C) Both (shift + meanlog)
W_both <- make_weights_generic(
  taus = taus,
  q_tau = stage1$q_tau,
  dist_cdf = shifted_cdf,
  dist_pdf = shifted_pdf,
  params = params0,
  grad_funcs = list(beta_shift = grad_shift, beta_meanlog = grad_meanlog)
)

fit_add <- param_gwas(
  stage1,
  transform = "custom_W",
  transform_args = list(W = W_add),
  se_mode = "plugin_cor"
)

fit_shift <- param_gwas(
  stage1,
  transform = "custom_W",
  transform_args = list(W = W_shift),
  se_mode = "plugin_cor"
)

fit_both <- param_gwas(
  stage1,
  transform = "custom_W",
  transform_args = list(W = W_both),
  se_mode = "plugin_cor"
)

# -------------------------
# 5) Summaries
# -------------------------
cat("\n--- GxE mean-only sim (E omitted) ---\n")
cat("N=", N, " P=", P, "\n", sep = "")
cat("Model: Y = mu0 + G*bG + (G*E)*bGE + eps\n")
cat("bG_sd=", bG_sd, " bGE_sd=", bGE_sd, " sd_eps=", sd_eps, "\n", sep = "")
cat("Mean Q (add):   ", round(mean(fit_add$Q), 3), "\n")
cat("Mean Q (shift): ", round(mean(fit_shift$Q), 3), "\n")
cat("Mean Q (both):  ", round(mean(fit_both$Q), 3), "\n")

# Compare Q differences per SNP
qdiff_shift_add <- mean(fit_shift$Q - fit_add$Q)
qdiff_both_add  <- mean(fit_both$Q - fit_add$Q)
qdiff_both_shift <- mean(fit_both$Q - fit_shift$Q)

cat("Mean Q diff (shift - add): ", round(qdiff_shift_add, 3), "\n")
cat("Mean Q diff (both  - add): ", round(qdiff_both_add, 3), "\n")
cat("Mean Q diff (both  - shift): ", round(qdiff_both_shift, 3), "\n")

# Correlate additive truth with estimated additive effect (sanity)
add_hat <- as.numeric(t(fit_add$params)[, "beta_add"])
cat("Cor(true bG, beta_add_hat): ", round(cor(bG, add_hat), 3), "\n")

# -------------------------
# 6) Plot selected SNP profiles
# -------------------------
rank_abs <- function(x) order(abs(x), decreasing = TRUE)
idx_add <- rank_abs(bG)[1:2]
idx_gxe <- rank_abs(bGE)[1:2]
idx_both <- rank_abs(bG + bGE)[1:2]
# SNP with both large main and interaction effects
idx_joint <- order(abs(bG) * abs(bGE), decreasing = TRUE)[1]
sel <- unique(c(idx_add, idx_gxe, idx_both, idx_joint))

# Print "both" parameters for the joint SNP
params_both_joint <- t(fit_both$params)[idx_joint, c("beta_shift", "beta_meanlog")]
cat("Joint SNP (", colnames(G)[idx_joint], ") both params: ",
    "beta_shift=", round(params_both_joint[["beta_shift"]], 4),
    " beta_meanlog=", round(params_both_joint[["beta_meanlog"]], 4),
    "\n", sep = "")

# Per-SNP fit diagnostics for joint SNP
q_add_j <- fit_add$Q[idx_joint]
q_shift_j <- fit_shift$Q[idx_joint]
q_both_j <- fit_both$Q[idx_joint]
cat("Joint SNP Q: add=", round(q_add_j, 3),
    " shift=", round(q_shift_j, 3),
    " both=", round(q_both_j, 3), "\n", sep = "")

tau_df <- data.frame(tau = taus, row = seq_along(taus))
out_dir <- file.path("stage2", "tests", "out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

plot_one <- function(j) {
  slope <- stage1$Q_slope[, j]
  pred_add <- as.numeric(W_add %*% t(fit_add$params)[j, "beta_add"])
  pred_shift <- as.numeric(W_shift %*% t(fit_shift$params)[j, "beta_meanlog"])
  pred_both <- as.numeric(W_both %*% as.numeric(t(fit_both$params)[j, c("beta_shift", "beta_meanlog")]))

  df <- data.frame(
    tau = taus,
    slope = slope,
    add = pred_add,
    shift = pred_shift,
    both = pred_both
  )

  ggplot(df, aes(tau, slope)) +
    geom_hline(yintercept = 0, linewidth = 0.5, color = "gray40", linetype = "dashed") +
    geom_point(size = 1.8, color = "black") +
    geom_line(aes(y = add, color = "add"), linewidth = 0.9) +
    geom_line(aes(y = shift, color = "shift"), linewidth = 0.9) +
    geom_line(aes(y = both, color = "both"), linewidth = 0.9, linetype = "dashed") +
    scale_color_manual(values = c(add = "#d62728", shift = "#1f77b4", both = "#2ca02c")) +
    labs(
      title = paste0("SNP ", colnames(G)[j]),
      subtitle = paste0("bG=", round(bG[j], 3), "  bGE=", round(bGE[j], 3)),
      x = "Quantile (tau)",
      y = "Tau-slope"
    ) +
    theme_minimal()
}

for (j in sel) {
  p <- plot_one(j)
  ggsave(
    filename = file.path(out_dir, paste0("gxe_profile_", colnames(G)[j], ".png")),
    plot = p,
    width = 7.2,
    height = 4.2,
    dpi = 130
  )
}

cat("Saved plots to: ", out_dir, "\n", sep = "")

# -------------------------
# 7) Direct comparison plot: GxE vs multiplicative log-normal
# -------------------------
j_cmp <- idx_joint
df_cmp <- data.frame(
  tau = taus,
  slope_gxe = stage1$Q_slope[, j_cmp],
  slope_mult = stage1_mult$Q_slope[, j_cmp]
)

p_cmp <- ggplot(df_cmp, aes(tau)) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "gray40", linetype = "dashed") +
  geom_point(aes(y = slope_gxe, color = "GxE"), size = 1.8) +
  geom_line(aes(y = slope_gxe, color = "GxE"), linewidth = 0.9) +
  geom_point(aes(y = slope_mult, color = "Multiplicative"), size = 1.8) +
  geom_line(aes(y = slope_mult, color = "Multiplicative"), linewidth = 0.9, linetype = "dashed") +
  scale_color_manual(values = c(GxE = "#1f77b4", Multiplicative = "#ff7f0e")) +
  labs(
    title = paste0("Tau-slope comparison for ", colnames(G)[j_cmp]),
    subtitle = paste0("GxE: bG=", round(bG[j_cmp], 3),
                      " bGE=", round(bGE[j_cmp], 3),
                      " | Mult: bM=", round(bM[j_cmp], 3)),
    x = "Quantile (tau)",
    y = "Tau-slope"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(out_dir, paste0("compare_gxe_vs_mult_", colnames(G)[j_cmp], ".png")),
  plot = p_cmp,
  width = 7.2,
  height = 4.2,
  dpi = 130
)
