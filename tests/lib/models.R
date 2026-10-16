# ---------------------------------------------------------------------------
# tests/lib/models.R -- the model catalogue the suite fits and tries to tell apart
#
# G-MoQ's claim is that it can distinguish distributional effects. The suite is
# only a test of that claim if the candidate models it fits actually differ in
# shape, and if for every simulated truth we know which candidate spans it.
# This file is the single place where both are written down; the Python side
# (tests/lib/gmoq_sim.py) generates data under the matching names.
#
# THE DESIGN. Every column is a derivative of a quantile function with respect
# to one parameter, which is what the Stage 2 identity requires:
#
#     W_j(tau) = dq_tau/dtheta_j = -(dF/dtheta_j)(q_tau) / f(q_tau)
#
# The baseline is a SINH-ARCSINH family,
#
#     q_tau = mu + sigma * sinh( (asinh(z_tau) + e) / d ),   z = Phi^-1(tau)
#
# with e the skewness parameter and d the tail-weight parameter. Writing
# u = (asinh(z) + e)/d, the four columns are closed form and hand-checkable:
#
#     mu     : 1
#     sigma  : sinh(u)
#     skew   : sigma cosh(u) / d
#     tail   : -sigma cosh(u) u
#
# A Cornish-Fisher cubic was used first and abandoned: its quantile function
# stops being monotone once the shape coefficients grow, which caps the largest
# simulable kurtosis effect at |Z| ~ 1.5 at N = 50,000 -- a limit of the
# simulation that would have been indistinguishable in the results from a limit
# of the method. sinh-arcsinh is monotone for every e and every d > 0.
#
# These are used rather than make_weights_vqtl()'s [1, z/(2 sd)] because the
# suite needs shape columns too, and mixing parameterisations across models
# would make the tournament a comparison of scalings rather than of shapes.
# make_weights_vqtl() is tested directly in tier 0 against its own closed form,
# and the "LS_pkg" model below fits it unchanged so the shipped builder is
# exercised end to end as well.
#
# For the mixture families the shipped builder make_weights_normal_mixture_vqtl()
# is used as-is; its columns are validated against finite differences of the
# mixture quantile function in tier 0.
# ---------------------------------------------------------------------------

#' Cornish-Fisher weight columns for a Gaussian baseline.
#'
#' @param taus Tau grid.
#' @param sigma Baseline scale of the phenotype.
#' @param e0,d0 Baseline sinh-arcsinh skewness and log tail weight.
#' @return Named list of numeric vectors, one per parameter.
gmoq_cf_columns <- function(taus, sigma, e0 = 0, d0 = 0) {
  z <- stats::qnorm(taus)
  d <- exp(d0)
  u <- (asinh(z) + e0) / d
  list(
    beta_mu    = rep(1, length(taus)),
    beta_sigma = sinh(u),
    beta_skew  = sigma * cosh(u) / d,
    beta_kurt  = -sigma * cosh(u) * u
  )
}

#' The candidate model set fitted in the tournament.
#'
#' @param taus Tau grid.
#' @param sigma Baseline phenotype sd.
#' @param mix Optional list(p1, mu1, sd1, mu2, sd2) and q_tau, enabling the
#'   mixture candidates. Omit for the Gaussian-only tournament.
#' @return Named list of T x K weight matrices.
gmoq_model_set <- function(taus, sigma, mix = NULL, e0 = 0, d0 = 0) {
  cf <- gmoq_cf_columns(taus, sigma, e0 = e0, d0 = d0)
  mk <- function(...) {
    nm <- c(...)
    M <- do.call(cbind, cf[nm])
    colnames(M) <- nm
    rownames(M) <- paste0("tau", taus)
    storage.mode(M) <- "double"
    M
  }
  out <- list(
    ADD   = mk("beta_mu"),                                        # pure location
    SCL   = mk("beta_sigma"),                                     # pure scale
    LS    = mk("beta_mu", "beta_sigma"),                          # location + scale
    LSS   = mk("beta_mu", "beta_sigma", "beta_skew"),             # + skew
    LSK   = mk("beta_mu", "beta_sigma", "beta_kurt"),             # + kurtosis
    LSSK  = mk("beta_mu", "beta_sigma", "beta_skew", "beta_kurt")
  )
  # The shipped vQTL builder, fitted unchanged, so the production code path is
  # in the tournament and not only in the unit tests.
  out$LS_pkg <- make_weights_vqtl(taus, q_tau = stats::qnorm(taus, sd = sigma),
                                  mu = 0, sd = sigma)
  if (!is.null(mix)) {
    Wm <- make_weights_normal_mixture_vqtl(
      taus = taus, q_tau = mix$q_tau, p1 = mix$p1,
      mu1 = mix$mu1, sd1 = mix$sd1, mu2 = mix$mu2, sd2 = mix$sd2,
      density = "model"   # the tournament's registered weights (pre-0.6.0 default)
    )
    out$MIX_W   <- Wm[, "gamma", drop = FALSE]        # SNP moves the mixing weight
    out$MIX_MU1 <- Wm[, "beta_mu1", drop = FALSE]     # SNP moves component 1's mean
    out$MIX_LS  <- Wm[, c("gamma", "beta_mu1", "beta_mu2"), drop = FALSE]
    out$MIX_ALL <- Wm
  }
  out
}

#' Which candidate models span each simulated truth.
#'
#' `true` is the parsimonious model a correct tournament should select.
#' `spanning` also contains it and its supersets: those are not wrong, they are
#' over-parameterised, and their Q is a valid null. `wrong` is the set whose
#' column space does NOT contain the truth; a lack-of-fit statistic that cannot
#' reject those has not demonstrated anything.
GMOQ_TRUTH_MAP <- list(
  null        = list(true = "ADD",     spanning = c("ADD","SCL","LS","LSS","LSK","LSSK","LS_pkg"),
                     wrong = character(0)),
  location    = list(true = "ADD",     spanning = c("ADD","LS","LSS","LSK","LSSK","LS_pkg"),
                     wrong = c("SCL")),
  scale       = list(true = "SCL",     spanning = c("SCL","LS","LSS","LSK","LSSK","LS_pkg"),
                     wrong = c("ADD")),
  locscale    = list(true = "LS",      spanning = c("LS","LSS","LSK","LSSK","LS_pkg"),
                     wrong = c("ADD","SCL")),
  # For the two shape families the models that omit only the shape column are
  # listed separately as `wrong_shape`. Measured (tier 1): the lack-of-fit
  # statistic rejects a model missing a LOCATION or SCALE dimension every time,
  # and does NOT reject one missing only a shape dimension. That is a finding
  # about Q, not a bug, and it is reported rather than asserted.
  skew        = list(true = "LSS",     spanning = c("LSS","LSSK"),
                     wrong = c("ADD","SCL"),
                     wrong_shape = c("LS","LSK","LS_pkg")),
  kurtosis    = list(true = "LSK",     spanning = c("LSK","LSSK"),
                     wrong = c("ADD","SCL"),
                     wrong_shape = c("LS","LSS","LS_pkg")),
  mix_weight  = list(true = "MIX_W",   spanning = c("MIX_W","MIX_LS","MIX_ALL"),
                     wrong = c("ADD","SCL","LS","MIX_MU1")),
  mix_mean    = list(true = "MIX_MU1", spanning = c("MIX_MU1","MIX_LS","MIX_ALL"),
                     wrong = c("ADD","SCL","MIX_W"))
)

#' Names of the eight simulated truths the suite covers.
GMOQ_FAMILIES <- names(GMOQ_TRUTH_MAP)
