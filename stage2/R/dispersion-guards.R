## dispersion-guards.R
##
## Two independent threats turn a variant with a large MEAN effect into a spurious
## dispersion hit.  They have different shapes and need different tools.
##
##  (a) MEAN-VARIANCE COUPLING.  If the analysed scale has a mean-variance
##      relationship -- raw BMI, raw weight, anything right-skewed and untransformed --
##      then a mean effect mechanically drags the variance with it.  The induced
##      dispersion effect is LINEAR in the mean effect with a slope that is a property
##      of the trait, not of the variant.  That is exactly the structure Young et al.
##      2018 (Nat Genet 50:1608) remove genome-wide with
##          d_l = alpha_vl - r_av * alpha_l
##      and it is what `estimate_rav()` / `residualise_dispersion()` below implement.
##      Measured on a lognormal trait it is the dominant artefact by an order of
##      magnitude; on a log- or rank-inverse-normal-transformed trait it is ZERO, so
##      the correction is a no-op there and must not be applied blindly.
##
##  (b) GENOTYPE SKEWNESS x DENSITY CURVATURE.  A pure location effect delta at a
##      variant whose genotype distribution is skewed leaves a second-order residue
##          b(tau) = delta - (delta^2/2) * gskew * f0'(q)/(f0(q) sqrt(Var g)) + O(delta^3)
##      and for a Normal baseline f0'/f0 = -(q - m)/v is exactly proportional to the
##      dispersion column `make_weights_vqtl()` builds, so the residue lands on
##      `beta_sigma2` in full.  It is QUADRATIC in the mean effect, is per-variant in
##      the genotype skewness, vanishes identically at MAF 0.5, and survives every
##      transformation of the phenotype.  A single genome-wide slope cannot absorb it,
##      so `residualise_dispersion()` is the wrong tool.
##
##      This is BIAS, not a variance error, so it can be SUBTRACTED rather than only
##      flagged.  `correct_dispersion_skew()` does that.  Because the residue is not
##      merely proportional to the dispersion weight column but is exactly
##          delta^2 * (mu3/Var g) * W[, "beta_sigma2"],
##      and because every Stage 2 projection A = (W'S^-1 W)^-1 W'S^-1 satisfies AW = I,
##      the induced bias in `beta_sigma2` is exactly delta^2 * (mu3/Var g) and the
##      induced bias in `beta_mu` is exactly zero -- for ols, gls, dwls and plugin_cor
##      alike, on any tau grid, with the phenotype SD cancelling out entirely.
##      `dispersion_skew_guard()` remains, for the cases where the correction cannot be
##      applied: it needs only the MAGNITUDE of the genotype skewness, while the
##      correction needs its SIGN and hence a signed effect-allele frequency.
##
## The two are orthogonal: (a) is linear in the mean effect and common across variants,
## (b) is quadratic and variant-specific.  Applying (a) does not touch (b) and vice versa.

#' Genotype-skewness guard for a dispersion finding
#'
#' Flags dispersion findings that may be the second-order residue of a large mean
#' effect at a skewed genotype rather than a real dispersion effect.
#'
#' The statistic is
#' \deqn{G = c \cdot \mathrm{gskew} \cdot Z_\mu^2 / \sqrt{N}}
#' with \eqn{\mathrm{gskew} = (1 - 2p)/\sqrt{2p(1-p)}} the genotype skewness and
#' \eqn{c = 0.6855}. \code{G} estimates the \eqn{|Z|} that \code{beta_sigma2} will
#' show when the true dispersion effect is zero, so it is directly comparable with
#' the observed dispersion \eqn{|Z|}.
#'
#' The constant was calibrated on an exact sweep with a Normal baseline over
#' MAF 0.005-0.5 and N 50,000-461,320, and is stable to about 1 per cent for
#' \eqn{Z_\mu \le 12} and 3 per cent to \eqn{Z_\mu \le 22}. On a right-skewed
#' baseline (standardised chi-squared) it over-predicts the artefact, i.e. it is
#' conservative there. It has not been validated on heavy-tailed baselines.
#'
#' @details
#' \strong{Supplying the allele frequency.} Stage 1 does not emit an allele
#' frequency, so if \code{maf} is not given it is recovered from the per-variant
#' \code{d_resid} (\eqn{\sum_i g_{\mathrm{resid},i}^2}, the \code{n_eff} vector
#' Stage 2 already threads through) together with \code{N}, using
#' \eqn{v = d_{\mathrm{resid}}/N \approx 2p(1-p)} and hence
#' \eqn{\mathrm{gskew} = \sqrt{(1 - 2v)/v}}. That recovery gives the
#' \emph{magnitude} of the genotype skewness but not its sign, because
#' \eqn{2p(1-p)} cannot distinguish \eqn{p} from \eqn{1-p}. This is why the guard
#' flags rather than corrects: \code{G} is used against a threshold on its
#' magnitude, for which the sign is not needed.
#'
#' @param z_mu Numeric vector of estimated location Z statistics (\code{beta_mu / se}).
#' @param N Sample size. Required unless \code{maf} is supplied.
#' @param maf Optional numeric vector of allele frequencies. If supplied, used directly.
#' @param n_eff Optional numeric vector of per-variant \code{d_resid}. Used with
#'   \code{N} to recover the genotype skewness when \code{maf} is absent.
#' @param c_guard Calibration constant; defaults to 0.6855.
#'
#' @return A numeric vector of \code{G} values, the same length as \code{z_mu}.
#'   \code{NA} where the inputs do not permit a value.
#'
#' @references Young AI, Wauthier FL, Donnelly P (2018) Identifying loci affecting
#'   trait variability and detecting interactions in genome-wide association
#'   studies. Nat Genet 50:1608-1614.
#'
#' @examples
#' # a large mean effect at a rare, strongly skewed genotype
#' dispersion_skew_guard(z_mu = c(8, 16, 30), N = 461320, maf = 0.01)
#' @export
dispersion_skew_guard <- function(z_mu, N = NULL, maf = NULL, n_eff = NULL,
                                  c_guard = 0.6855) {
  z_mu <- as.numeric(z_mu)
  if (!is.null(maf)) {
    p <- as.numeric(maf)
    if (length(p) == 1L) p <- rep(p, length(z_mu))
    if (length(p) != length(z_mu)) stop("maf must be length 1 or length(z_mu).")
    bad <- !is.finite(p) | p <= 0 | p >= 1
    gskew <- rep(NA_real_, length(p))
    gskew[!bad] <- (1 - 2 * p[!bad]) / sqrt(2 * p[!bad] * (1 - p[!bad]))
  } else {
    if (is.null(n_eff) || is.null(N)) {
      stop("Supply either maf, or both n_eff (Stage 1 d_resid) and N.")
    }
    v <- as.numeric(n_eff) / N
    if (length(v) == 1L) v <- rep(v, length(z_mu))
    if (length(v) != length(z_mu)) stop("n_eff must be length 1 or length(z_mu).")
    # v ~= 2p(1-p) in (0, 0.5]; gskew magnitude = sqrt((1 - 2v)/v). Sign is not
    # recoverable from v, so the guard is a magnitude and is used as such.
    gskew <- rep(NA_real_, length(v))
    ok <- is.finite(v) & v > 0 & v < 0.5
    gskew[ok] <- sqrt((1 - 2 * v[ok]) / v[ok])
    gskew[is.finite(v) & v >= 0.5] <- 0    # MAF 0.5: genotype skewness is exactly zero
  }
  if (is.null(N) || !is.finite(N) || N <= 0) stop("N must be a positive sample size.")
  c_guard * gskew * z_mu^2 / sqrt(N)
}


#' Estimate the genome-wide mean-variance coupling r_av
#'
#' Estimates the single slope \eqn{r_{av}} in \eqn{\alpha_{vl} = r_{av}\alpha_l + d_l},
#' the genome-wide regression of dispersion effects on mean effects that Young et al.
#' (2018) use to strip the part of a variance effect that is a mechanical consequence
#' of the mean effect on a mean-variance-coupled scale.
#'
#' @details
#' The estimator is a method of moments corrected for two distinct contaminations
#' that an ordinary regression of \eqn{\hat\alpha_v} on \eqn{\hat\alpha} suffers:
#' \deqn{\hat r_{av} = \frac{\sum_l w_l (\hat\alpha_l \hat\alpha_{vl} - c_l)}
#'                          {\sum_l w_l (\hat\alpha_l^2 - s_{al}^2)}}
#' \itemize{
#'   \item \eqn{-c_l}, the \emph{within-variant sampling covariance} between the two
#'     parameter estimates, removes coupling that is an artefact of the joint fit
#'     rather than a property of the trait. Stage 2 already reports \eqn{c_l} in
#'     \code{param_cov_upper}, so this costs nothing. Omitting it biases
#'     \eqn{\hat r_{av}} whenever the fitted \code{W} does not make the two
#'     parameters orthogonal.
#'   \item \eqn{-s_{al}^2} in the denominator removes regression dilution. Without it
#'     \eqn{\hat r_{av}} is attenuated by \eqn{\mathrm{Var}(\alpha)/(\mathrm{Var}(\alpha)
#'     + \overline{s_a^2})}, which for a scan dominated by null variants is severe.
#' }
#' The estimator is consistent under the assumption that genuine dispersion effects
#' \eqn{d_l} are uncorrelated with mean effects \eqn{\alpha_l}. That assumption is
#' Young et al.'s and is not testable from these data; it is what makes genuine
#' vQTLs add noise to \eqn{\hat r_{av}} rather than bias.
#'
#' \strong{Which variants to use.} Pass an LD-pruned set. Do not pass only the
#' genome-wide-significant variants: \eqn{r_{av}} would then be estimated from the
#' very hits whose status is in question, and the correction would be circular.
#' The denominator correction makes it safe to include null variants, so an
#' LD-pruned genome-wide set is the intended input.
#'
#' @param beta_mu,beta_sigma2 Numeric vectors of parameter estimates.
#' @param se_mu,se_sigma2 Numeric vectors of their standard errors.
#' @param cov_mu_sigma2 Numeric vector of the within-variant sampling covariance
#'   between the two, from \code{param_cov_upper}. If \code{NULL}, treated as zero
#'   and a warning is issued, since that is an assumption rather than a default.
#' @param weights Optional non-negative weights. Default equal weights.
#' @param block Optional grouping vector (chromosome, or LD block) used for a
#'   delete-one-block jackknife standard error. At least 10 distinct blocks are
#'   needed for a usable SE.
#' @param trim Fraction trimmed from each tail of the per-variant contributions
#'   before summing, guarding against a handful of extreme variants dominating the
#'   slope. Default 0.01. Set to 0 for the untrimmed estimator.
#'
#' @return A list with \code{r_av}, \code{se} (jackknife, or \code{NA}),
#'   \code{n_used}, \code{n_blocks}, \code{denominator}, and \code{trim}.
#'
#' @seealso \code{\link{residualise_dispersion}}, \code{\link{dispersion_skew_guard}}
#' @export
estimate_rav <- function(beta_mu, beta_sigma2, se_mu, se_sigma2,
                         cov_mu_sigma2 = NULL, weights = NULL,
                         block = NULL, trim = 0.01) {
  n <- length(beta_mu)
  stopifnot(length(beta_sigma2) == n, length(se_mu) == n, length(se_sigma2) == n)
  if (is.null(cov_mu_sigma2)) {
    warning("cov_mu_sigma2 not supplied; assuming the within-variant sampling ",
            "covariance between beta_mu and beta_sigma2 is zero. That is an ",
            "assumption about the fitted W, not a safe default -- pass ",
            "param_cov_upper if you have it.")
    cov_mu_sigma2 <- rep(0, n)
  }
  if (is.null(weights)) weights <- rep(1, n)
  ok <- is.finite(beta_mu) & is.finite(beta_sigma2) & is.finite(se_mu) &
        is.finite(se_sigma2) & is.finite(cov_mu_sigma2) & is.finite(weights) &
        weights >= 0 & se_mu > 0 & se_sigma2 > 0
  if (sum(ok) < 20L) stop("Fewer than 20 usable variants for estimating r_av.")

  beta_mu <- beta_mu[ok]; beta_sigma2 <- beta_sigma2[ok]
  se_mu <- se_mu[ok]; cov_mu_sigma2 <- cov_mu_sigma2[ok]; weights <- weights[ok]
  num <- weights * (beta_mu * beta_sigma2 - cov_mu_sigma2)
  den <- weights * (beta_mu^2 - se_mu^2)
  blk <- if (is.null(block)) NULL else block[ok]

  keep <- rep(TRUE, length(num))
  if (trim > 0) {
    # trim on the per-variant ratio contribution, not on either sum alone
    lim <- stats::quantile(num, c(trim, 1 - trim), na.rm = TRUE)
    keep <- num >= lim[1] & num <= lim[2]
  }
  num <- num[keep]; den <- den[keep]
  if (!is.null(blk)) blk <- blk[keep]
  r_av <- sum(num) / sum(den)

  se <- NA_real_; nb <- NA_integer_
  if (!is.null(blk)) {
    ub <- unique(blk); nb <- length(ub)
    if (nb >= 10L) {
      # delete-one-block jackknife; sums are additive so each leave-one-out is O(1)
      sn <- sum(num); sd_ <- sum(den)
      bn <- vapply(ub, function(u) sum(num[blk == u]), 0.0)
      bd <- vapply(ub, function(u) sum(den[blk == u]), 0.0)
      ps <- (sn - bn) / (sd_ - bd)
      se <- sqrt((nb - 1) / nb * sum((ps - mean(ps))^2))
    }
  }
  list(r_av = r_av, se = se, n_used = length(num), n_blocks = nb,
       denominator = sum(den), trim = trim)
}


#' Residualise a dispersion effect on the mean effect
#'
#' Applies Young et al.'s genome-wide correction \eqn{d_l = \alpha_{vl} -
#' r_{av}\alpha_l} to per-variant Stage 2 estimates, and returns the corrected
#' statistic \strong{alongside} the uncorrected one rather than in place of it.
#'
#' @details
#' \strong{The standard error.} The correction is applied to the point estimate and
#' to the standard error, which is
#' \deqn{\mathrm{Var}(d_l) = s_{vl}^2 + r_{av}^2 s_{al}^2 - 2 r_{av} c_l
#'                           + \mathrm{Var}(\hat r_{av})\,\alpha_l^2 .}
#' The first three terms are the per-variant delta method and use the within-variant
#' covariance \eqn{c_l} that Stage 2 already reports. The fourth propagates the
#' uncertainty in \eqn{\hat r_{av}} itself; it is negligible for null variants and
#' matters precisely for the large-\eqn{|\alpha|} variants the correction exists to
#' adjudicate, so it is included by default. Set \code{se_rav = 0} to reproduce the
#' uncorrected-\eqn{r_{av}} version used in the original paper.
#'
#' Because \eqn{c_l} is near zero for a well-specified Normal \code{W}, the corrected
#' SE is generally \emph{larger} than the uncorrected one. The correction therefore
#' costs power, and on a scale with no mean-variance coupling (a log- or
#' rank-inverse-normal-transformed trait, where \eqn{r_{av} \approx 0}) it costs
#' power for no benefit. That is the reason this is an extra column and not a
#' replacement.
#'
#' \strong{What it does not fix.} The correction removes only the component of a
#' dispersion effect that is linear in the mean effect with a slope common to all
#' variants. It does nothing to the genotype-skewness residue, which is quadratic in
#' the mean effect and per-variant in the genotype skewness. Use
#' \code{\link{dispersion_skew_guard}} for that.
#'
#' @param beta_mu,beta_sigma2 Numeric vectors of parameter estimates.
#' @param se_mu,se_sigma2 Numeric vectors of their standard errors.
#' @param cov_mu_sigma2 Within-variant sampling covariance, from \code{param_cov_upper}.
#' @param r_av The genome-wide coupling, e.g. from \code{\link{estimate_rav}}.
#' @param se_rav Standard error of \code{r_av}. Default 0, i.e. treat it as known.
#' @param df Optional degrees of freedom for a t reference; Normal if \code{NULL}.
#'
#' @return A data frame with columns \code{beta_sigma2_resid}, \code{se_sigma2_resid},
#'   \code{z_sigma2_resid}, \code{p_sigma2_resid}.
#'
#' @seealso \code{\link{estimate_rav}}, \code{\link{dispersion_skew_guard}}
#' @export
residualise_dispersion <- function(beta_mu, beta_sigma2, se_mu, se_sigma2,
                                   cov_mu_sigma2 = NULL, r_av, se_rav = 0,
                                   df = NULL) {
  n <- length(beta_mu)
  stopifnot(length(beta_sigma2) == n, length(se_mu) == n, length(se_sigma2) == n)
  if (!is.numeric(r_av) || length(r_av) != 1L || !is.finite(r_av)) {
    stop("r_av must be a single finite number.")
  }
  if (is.null(cov_mu_sigma2)) cov_mu_sigma2 <- rep(0, n)

  est <- beta_sigma2 - r_av * beta_mu
  v   <- se_sigma2^2 + r_av^2 * se_mu^2 - 2 * r_av * cov_mu_sigma2 +
         se_rav^2 * beta_mu^2
  v[!is.finite(v) | v < 0] <- NA_real_
  se <- sqrt(v)
  z  <- est / se
  p  <- if (is.null(df)) 2 * stats::pnorm(-abs(z)) else 2 * stats::pt(-abs(z), df = df)
  data.frame(beta_sigma2_resid = est, se_sigma2_resid = se,
             z_sigma2_resid = z, p_sigma2_resid = p)
}


#' Add dispersion guards to a Stage 2 parameter fit
#'
#' Convenience wrapper that takes the output of \code{\link{param_gwas}} together with
#' its \code{param_cov_upper}, and returns a data frame carrying the uncorrected
#' dispersion statistic, the mean-variance-residualised one, and the
#' genotype-skewness guard, side by side.
#'
#' @param fit Output of \code{param_gwas()} run with \code{return_cov = "upper"} and
#'   a weight matrix carrying columns named \code{beta_mu} and \code{beta_sigma2}
#'   (i.e. \code{make_weights_vqtl}).
#' @param r_av,se_rav Genome-wide coupling and its SE; see \code{\link{estimate_rav}}.
#'   If \code{r_av} is \code{NULL} the residualised columns are omitted.
#' @param N Sample size, for the skewness guard.
#' @param maf,n_eff Allele frequency, or per-variant \code{d_resid}; see
#'   \code{\link{dispersion_skew_guard}}. If both are \code{NULL} the guard is omitted.
#' @param eaf Signed effect-allele frequency. When supplied, the
#'   genotype-skewness-corrected dispersion columns are added; see
#'   \code{\link{correct_dispersion_skew}}. This must be oriented to the allele
#'   \code{beta_mu} counts -- \code{maf} is not a substitute.
#' @param debias Passed to \code{\link{correct_dispersion_skew}}.
#' @param mu_name,sigma2_name Parameter row names to use.
#'
#' @return A data frame with one row per variant.
#' @export
dispersion_report <- function(fit, r_av = NULL, se_rav = 0, N = NULL,
                              maf = NULL, n_eff = NULL, eaf = NULL, debias = TRUE,
                              mu_name = "beta_mu", sigma2_name = "beta_sigma2") {
  if (is.null(fit$params) || is.null(fit$SE_params)) {
    stop("fit must carry params and SE_params from param_gwas().")
  }
  pn <- rownames(fit$params)
  if (!all(c(mu_name, sigma2_name) %in% pn)) {
    stop("fit must contain parameters '", mu_name, "' and '", sigma2_name, "'. Found: ",
         paste(pn, collapse = ", "))
  }
  a  <- fit$params[mu_name, ];     v  <- fit$params[sigma2_name, ]
  sa <- fit$SE_params[mu_name, ];  sv <- fit$SE_params[sigma2_name, ]

  cl <- NULL
  if (!is.null(fit$param_cov_upper)) {
    nm <- paste0("cov_", mu_name, "_", sigma2_name)
    alt <- paste0("cov_", sigma2_name, "_", mu_name)
    cn <- colnames(fit$param_cov_upper)
    if (nm %in% cn)       cl <- as.numeric(fit$param_cov_upper[, nm])
    else if (alt %in% cn) cl <- as.numeric(fit$param_cov_upper[, alt])
  }
  if (is.null(cl)) {
    warning("param_cov_upper did not carry the ", mu_name, "/", sigma2_name,
            " covariance; assuming zero. Re-run param_gwas() with return_cov='upper'.")
    cl <- rep(0, length(a))
  }

  out <- data.frame(variant = colnames(fit$params),
                    beta_mu = a, se_mu = sa, z_mu = a / sa,
                    beta_sigma2 = v, se_sigma2 = sv, z_sigma2 = v / sv,
                    cov_mu_sigma2 = cl,
                    p_sigma2 = 2 * stats::pnorm(-abs(v / sv)),
                    stringsAsFactors = FALSE)

  if (!is.null(r_av)) {
    out <- cbind(out, residualise_dispersion(a, v, sa, sv, cl, r_av, se_rav))
  }
  if (!is.null(eaf)) {
    out <- cbind(out, correct_dispersion_skew(a, v, sa, sv, cl, eaf = eaf,
                                              debias = debias))
  }
  if (!is.null(maf) || !is.null(n_eff)) {
    out$G_skew_guard <- dispersion_skew_guard(out$z_mu, N = N, maf = maf, n_eff = n_eff)
    out$flag_skew <- !is.na(out$G_skew_guard) & abs(out$G_skew_guard) >= 0.2
  }
  rownames(out) <- NULL
  out
}


#' Predicted genotype-skewness bias in the dispersion parameter
#'
#' Returns the bias that a pure location effect induces in \code{beta_sigma2},
#' i.e. the quantity \code{\link{correct_dispersion_skew}} subtracts.
#'
#' @details
#' \strong{The result.} For \eqn{Y = \mu + \delta g + \varepsilon} with
#' \eqn{\varepsilon} independent of \eqn{g}, the RIF-OLS estimand at level
#' \eqn{\tau} is not \eqn{\delta} but
#' \deqn{b(\tau) = h'(\bar g) + \tfrac{1}{2} h''(\bar g)\,\mu_3/\mathrm{Var}(g)
#'               = \delta + \frac{\delta^2}{2\sigma}\,\frac{\mu_3}{\mathrm{Var}(g)}\,
#'                 z_\tau + O(\delta^3),}
#' \eqn{z_\tau = \Phi^{-1}(\tau)}, \eqn{\sigma} the phenotype SD,
#' \eqn{\mu_3 = E[(g - \bar g)^3]}. For \eqn{g \sim \mathrm{Binom}(2, p)},
#' \eqn{\mu_3/\mathrm{Var}(g) = 1 - 2p}.
#'
#' \strong{Why this lands entirely on one coefficient.}
#' \code{\link{make_weights_vqtl}} builds the dispersion column as
#' \eqn{W_{\sigma^2}(\tau) = z_\tau/(2\sigma)}. The bias above is therefore
#' \emph{exactly} \eqn{\delta^2 (\mu_3/\mathrm{Var}\,g)} times that column --- not
#' merely proportional to it, the same vector. Stage 2 forms
#' \eqn{\hat\theta = A b} with \eqn{A = (W'S^{-1}W)^{-1}W'S^{-1}} for some
#' weighting \eqn{S}, and every such \eqn{A} satisfies \eqn{AW = I}. Hence
#' \deqn{A \cdot \mathrm{bias} = \delta^2\,\frac{\mu_3}{\mathrm{Var}(g)}\;e_2 ,}
#' so the induced bias in \eqn{\hat\beta_{\sigma^2}} is
#' \eqn{\delta^2 \mu_3/\mathrm{Var}(g)} exactly, the induced bias in
#' \eqn{\hat\beta_\mu} is exactly zero, and both statements hold for
#' \code{ols}, \code{gls}, \code{dwls} and \code{plugin_cor} alike, for any
#' \eqn{\tau} grid. Two consequences worth stating plainly:
#' \itemize{
#'   \item \eqn{\sigma} \strong{cancels}. It appears in the bias and again in the
#'     weight column, so the correction needs no phenotype SD and no new Stage 1
#'     output. (\code{w_sd_ratio} exists only for the unusual case where the
#'     \code{sd} passed to \code{make_weights_vqtl} was not the phenotype SD.)
#'   \item Correcting \eqn{b(\tau)} before the projection and correcting the
#'     coefficient after it are \emph{identical}, not merely similar. Correcting
#'     afterwards is chosen because it is O(1) per variant and needs no refit.
#' }
#'
#' \strong{Order of the expansion.} The derivation is second order in
#' \eqn{\delta}. The third-order term is even in \eqn{z_\tau}, so it is
#' orthogonal to the dispersion column and lands on \eqn{\beta_\mu} and on the
#' fit residual instead. The first term that survives into
#' \eqn{\beta_{\sigma^2}} is therefore \eqn{O(\delta^4)}: measured against an
#' exact population calculation the relative error of this correction is
#' \eqn{\approx 0.065\,(\delta/\sigma)^2}, i.e. 0.07\% at 0.1 SD per allele and
#' 0.3\% at 0.2 SD per allele. It is the \emph{residual} bias after correction
#' that matters, and it is smaller than the correction by that factor.
#'
#' @param beta_mu Numeric vector of estimated location coefficients.
#' @param se_mu Numeric vector of their standard errors. Required when
#'   \code{debias = TRUE}.
#' @param eaf Frequency of the allele that \code{beta_mu} counts --- the
#'   \strong{effect} allele frequency, not the folded MAF. The sign of
#'   \eqn{1 - 2p} is the sign of the bias, so a folded MAF is not enough; see
#'   \code{\link{correct_dispersion_skew}}. Length 1 or \code{length(beta_mu)}.
#' @param mu3_over_var Optional direct value of \eqn{\mu_3/\mathrm{Var}(g)},
#'   overriding \code{eaf}. Use this for dosages, or for any genotype law where
#'   the Hardy-Weinberg identity \eqn{\mu_3/\mathrm{Var}(g) = 1 - 2p} does not
#'   hold.
#' @param debias Subtract \eqn{\mathrm{Var}(\hat\delta)} from \eqn{\hat\delta^2}.
#'   \eqn{E[\hat\delta^2] = \delta^2 + \mathrm{se}_\mu^2}, so the raw square
#'   over-states the bias by one squared standard error; subtracting it makes
#'   the correction unbiased. The estimate can then go negative, which is
#'   correct behaviour for an unbiased estimator of a non-negative quantity and
#'   is not truncated. Default \code{TRUE}.
#' @param w_sd_ratio Ratio (\code{sd} passed to \code{make_weights_vqtl}) /
#'   (phenotype SD). Leave at 1 unless \code{W} was deliberately built on a
#'   different scale.
#'
#' @return Numeric vector, the predicted bias in \code{beta_sigma2}.
#' @seealso \code{\link{correct_dispersion_skew}}, \code{\link{dispersion_skew_guard}}
#' @export
dispersion_skew_bias <- function(beta_mu, se_mu = NULL, eaf = NULL,
                                 mu3_over_var = NULL, debias = TRUE,
                                 w_sd_ratio = 1) {
  beta_mu <- as.numeric(beta_mu)
  n <- length(beta_mu)
  k <- .skew_k(eaf, mu3_over_var, n)
  if (isTRUE(debias)) {
    if (is.null(se_mu)) {
      stop("se_mu is required when debias = TRUE (the correction subtracts ",
           "Var(beta_mu_hat) from beta_mu_hat^2). Pass se_mu, or set debias = FALSE.")
    }
    se_mu <- as.numeric(se_mu)
    if (length(se_mu) == 1L) se_mu <- rep(se_mu, n)
    if (length(se_mu) != n) stop("se_mu must be length 1 or length(beta_mu).")
    delta2 <- beta_mu^2 - se_mu^2
  } else {
    delta2 <- beta_mu^2
  }
  as.numeric(w_sd_ratio) * k * delta2
}

## mu3/Var(g), from an effect-allele frequency under Hardy-Weinberg, or given.
.skew_k <- function(eaf, mu3_over_var, n) {
  if (!is.null(mu3_over_var)) {
    k <- as.numeric(mu3_over_var)
    if (length(k) == 1L) k <- rep(k, n)
    if (length(k) != n) stop("mu3_over_var must be length 1 or length(beta_mu).")
    return(k)
  }
  if (is.null(eaf)) {
    stop("Supply either eaf (the EFFECT allele frequency, signed) or mu3_over_var. ",
         "Note that the folded MAF is not sufficient: 1 - 2p changes sign under ",
         "p -> 1 - p, so an unoriented frequency would correct in the wrong ",
         "direction half the time. dispersion_skew_guard() needs only the ",
         "magnitude and so can fall back on d_resid; this correction cannot.")
  }
  p <- as.numeric(eaf)
  if (length(p) == 1L) p <- rep(p, n)
  if (length(p) != n) stop("eaf must be length 1 or length(beta_mu).")
  k <- rep(NA_real_, n)
  ok <- is.finite(p) & p > 0 & p < 1
  k[ok] <- 1 - 2 * p[ok]
  k
}


#' Correct the dispersion parameter for genotype-skewness bias
#'
#' Removes from \code{beta_sigma2} the second-order residue that a pure location
#' effect leaves behind when the genotype distribution is skewed, and widens the
#' standard error to pay for the estimation of that residue.
#'
#' @details
#' The corrected estimate is
#' \deqn{\tilde\beta_{\sigma^2} = \hat\beta_{\sigma^2}
#'        - k\,(\hat\beta_\mu^2 - \mathrm{se}_\mu^2),\qquad k = \mu_3/\mathrm{Var}(g) = 1 - 2p,}
#' see \code{\link{dispersion_skew_bias}} for why this is the right quantity and
#' why it is exact for \code{ols} and \code{gls} alike.
#'
#' \strong{Why \eqn{\hat\delta} from the same fit, and why no iteration.}
#' The bias is exactly collinear with the dispersion weight column, so it is
#' absorbed by \eqn{\hat\beta_{\sigma^2}} \emph{in full} and
#' \eqn{\hat\beta_\mu} carries none of it (verified exactly: the induced bias in
#' \eqn{\hat\beta_\mu} is below \eqn{10^{-4}} of the bias in
#' \eqn{\hat\beta_{\sigma^2}} at 0.2 SD per allele). \eqn{\hat\beta_\mu} is
#' therefore an unbiased input to its own correction, there is no feedback loop,
#' and iterating would change nothing. What \eqn{\hat\beta_\mu} \emph{is} is
#' noisy and correlated with \eqn{\hat\beta_{\sigma^2}}; both are handled in the
#' variance below, and the upward bias of the square is handled by
#' \code{debias}.
#'
#' \strong{The standard error must widen, and it does.} The correction is a
#' function of an estimate, so
#' \deqn{\mathrm{Var}(\tilde\beta_{\sigma^2}) = \mathrm{se}_{\sigma^2}^2
#'    + 4k^2\hat\beta_\mu^2\,\mathrm{se}_\mu^2 + 2k^2\,\mathrm{se}_\mu^4
#'    - 4k\,\hat\beta_\mu\,c_l ,}
#' where \eqn{c_l} is the within-variant sampling covariance from
#' \code{param_cov_upper}. The middle two terms are the exact variance of
#' \eqn{k\hat\beta_\mu^2} under joint normality (\eqn{4\delta^2 s^2 + 2s^4},
#' with \eqn{\hat\beta_\mu^2} in place of \eqn{\delta^2}); the last is
#' \eqn{-2\mathrm{Cov}(\hat\beta_{\sigma^2}, k\hat\beta_\mu^2)} using
#' \eqn{\mathrm{Cov}(X, Y^2) = 2E[Y]\mathrm{Cov}(X,Y)}. The expression is
#' non-negative for any admissible \eqn{c_l}, since
#' \eqn{|4k\hat\beta_\mu c_l| \le \mathrm{se}_{\sigma^2}^2 + 4k^2\hat\beta_\mu^2
#' \mathrm{se}_\mu^2} by Cauchy-Schwarz and AM-GM.
#'
#' Ignoring \eqn{c_l} is not safe: it is the term that keeps the corrected
#' statistic calibrated when \code{gls} makes the two parameters strongly
#' correlated. Pass \code{param_cov_upper}.
#'
#' \strong{When this is invalid.}
#' \itemize{
#'   \item \strong{It needs a signed effect-allele frequency.} \eqn{1 - 2p} is
#'     odd about \eqn{p = 1/2}; correcting with an unoriented MAF would double
#'     the bias half the time. Stage 1 does not emit a frequency, so
#'     \code{eaf} must come from the genotype source.
#'   \item \strong{It is specific to \code{make_weights_vqtl}} --- a
#'     location/scale \code{W} on a Normal baseline, where
#'     \eqn{f_0'/f_0 \propto z_\tau} makes the residue exactly collinear with
#'     one column. Under a mixture or other \code{W} the residue is not a single
#'     column and this correction does not apply.
#'   \item \strong{The baseline should be Normal.} On a non-Normal baseline the
#'     residue is \eqn{-(\delta^2/2)(\mu_3/\mathrm{Var} g) f_0'(q)/f_0(q)},
#'     which is no longer proportional to \eqn{z_\tau}; the correction then
#'     removes only the part of it that the \code{W} column can see. Measured
#'     behaviour is reported in \code{docs/DISPERSION_SKEW_CORRECTION.md}.
#'   \item It is second order in \eqn{\delta}; see
#'     \code{\link{dispersion_skew_bias}} for the residual.
#' }
#'
#' @param beta_mu Numeric vector of estimated location coefficients; its length
#'   sets the number of variants.
#' @param beta_sigma2 Estimated dispersion coefficients, length 1 or
#'   \code{length(beta_mu)}.
#' @param se_mu,se_sigma2 Their standard errors, length 1 or \code{length(beta_mu)}.
#' @param cov_mu_sigma2 Within-variant sampling covariance between the two, from
#'   \code{param_cov_upper}. If \code{NULL} it is taken as zero and a warning is
#'   issued.
#' @param eaf,mu3_over_var,debias,w_sd_ratio See \code{\link{dispersion_skew_bias}}.
#' @param df Optional degrees of freedom for a t reference; Normal if \code{NULL}.
#'
#' @return A data frame with \code{bias_skew} (the amount subtracted),
#'   \code{beta_sigma2_skewcorr}, \code{se_sigma2_skewcorr},
#'   \code{z_sigma2_skewcorr}, \code{p_sigma2_skewcorr}.
#'
#' @examples
#' # a 0.2 SD/allele location effect at EAF 0.01, no true dispersion effect
#' correct_dispersion_skew(beta_mu = 0.2, beta_sigma2 = 0.0392,
#'                         se_mu = 0.01, se_sigma2 = 0.012,
#'                         cov_mu_sigma2 = 0, eaf = 0.01)
#' @seealso \code{\link{dispersion_skew_bias}}, \code{\link{dispersion_skew_guard}},
#'   \code{\link{residualise_dispersion}}
#' @export
correct_dispersion_skew <- function(beta_mu, beta_sigma2, se_mu, se_sigma2,
                                    cov_mu_sigma2 = NULL, eaf = NULL,
                                    mu3_over_var = NULL, debias = TRUE,
                                    w_sd_ratio = 1, df = NULL) {
  beta_mu <- as.numeric(beta_mu); n <- length(beta_mu)
  rc <- function(x, nm) {
    x <- as.numeric(x)
    if (length(x) == 1L) x <- rep(x, n)
    if (length(x) != n) stop(nm, " must be length 1 or length(beta_mu).")
    x
  }
  beta_sigma2 <- rc(beta_sigma2, "beta_sigma2")
  se_mu <- rc(se_mu, "se_mu"); se_sigma2 <- rc(se_sigma2, "se_sigma2")
  if (is.null(cov_mu_sigma2)) {
    warning("cov_mu_sigma2 not supplied; assuming the within-variant sampling ",
            "covariance between beta_mu and beta_sigma2 is zero. Under GLS it ",
            "is not, and the corrected SE will be wrong -- pass param_cov_upper.")
    cov_mu_sigma2 <- rep(0, n)
  }
  cov_mu_sigma2 <- rc(cov_mu_sigma2, "cov_mu_sigma2")

  k <- .skew_k(eaf, mu3_over_var, n) * as.numeric(w_sd_ratio)
  bias <- dispersion_skew_bias(beta_mu, se_mu = se_mu, mu3_over_var = k,
                               debias = debias, w_sd_ratio = 1)
  est <- beta_sigma2 - bias
  v <- se_sigma2^2 + 4 * k^2 * beta_mu^2 * se_mu^2 + 2 * k^2 * se_mu^4 -
       4 * k * beta_mu * cov_mu_sigma2
  v[!is.finite(v) | v < 0] <- NA_real_
  se <- sqrt(v)
  z <- est / se
  p <- if (is.null(df)) 2 * stats::pnorm(-abs(z)) else 2 * stats::pt(-abs(z), df = df)
  data.frame(bias_skew = bias, beta_sigma2_skewcorr = est,
             se_sigma2_skewcorr = se, z_sigma2_skewcorr = z,
             p_sigma2_skewcorr = p)
}
