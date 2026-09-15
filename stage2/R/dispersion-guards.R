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
##
## ---------------------------------------------------------------------------
## WHY THE CORRECTION IS APPLIED AFTER THE PROJECTION, NOT BEFORE
##
## Correcting b(tau) BEFORE the projection is right for every W; correcting the
## coefficient AFTER is right only for make_weights_vqtl's design.  On that
## design the two are provably identical (asserted in test-dispersion-guards.R),
## and the before-projection route is no more expensive than it looks: the bias
## vector is a per-variant scalar delta^2 * k times the FIXED vector
## u = qnorm(taus)/(2*sigma_Y), so A %*% u is computed once and the per-variant
## cost is still O(1).  Cost is not the reason to prefer the shipped choice.
##
## The reason is that before-projection does not actually generalise:
##
##  1. sigma_Y stops cancelling.  After the projection it cancels exactly,
##     because the bias IS a multiple of the W column that carries sigma_Y in
##     its denominator.  Before the projection you must build
##     u = qnorm(taus)/(2*sigma_Y) numerically, so the phenotype SD becomes a
##     REQUIRED input that Stage 1 does not emit and Stage 2's report does not
##     thread through.  Getting it wrong scales the correction linearly, with no
##     symptom -- the same failure mode as a wrong w_sd_ratio, made mandatory.
##
##  2. delta stops being identifiable.  The bias is delta^2 * k where delta is
##     the per-allele location shift in phenotype units.  Under
##     make_weights_vqtl, W[,"beta_mu"] = 1, so AW = I makes beta_mu EXACTLY
##     that delta.  Under a general W there is no such guarantee: for the
##     5-parameter mixture, "beta_mu1" is a component mean, not the shift of the
##     whole distribution, which is precisely why that design over-subtracts by
##     2.4x today.  Shipping before-projection would move the silent failure
##     from "wrong W" to "wrong delta", not remove it.
##
## So the failure class is not closed by moving the correction earlier; it is
## closed by refusing to correct on a design the derivation does not cover.
## check_vqtl_design() does that, by checking the NUMBERS in W rather than its
## column names, and it is a total function: it either verifies the exact
## identity the proof needs or it refuses.  Under that gate the after-projection
## form is exact, so a general-W path would be dead code for every input that
## passes and would still be unsafe for every input that does not.
##
## Generalising properly needs both prerequisites above -- a phenotype SD
## threaded from Stage 1, and a declared location contrast for a general W.
## That is a redesign, and it is not what this change is.
## ---------------------------------------------------------------------------

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
#' @section What this guard is blind to:
#' \code{G} is proportional to the genotype skewness, which is \strong{exactly
#' zero at MAF 0.50}. The guard is therefore identically zero there, for any
#' \eqn{Z_\mu} however large. That is correct for the mechanism it was built for
#' --- the genotype-skewness residue really does vanish at MAF 0.50 --- but a
#' strongly skewed \emph{phenotype} produces dispersion inflation at MAF 0.50 as
#' well (0.163 at nominal 0.05 and \eqn{4.5\times10^{-6}} at
#' \eqn{5\times10^{-8}} on a \eqn{\chi^2_4} baseline with \eqn{Z_\mu = 60}).
#' That third mechanism lives exactly where this guard reads zero and where
#' \code{\link{correct_dispersion_skew}} is identically zero by construction, so
#' neither tool sees it. Nothing in this package addresses it; transform the
#' phenotype (\code{docs/phenotype_preparation.md}).
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
    # Trim the per-variant ratio contribution, not the numerator alone.  The
    # estimator is a ratio of two corrected sums, so a large numerator is not
    # necessarily influential when its denominator is large as well.
    ratio <- num / den
    finite <- is.finite(ratio)
    if (sum(finite) < 20L) {
      stop("Fewer than 20 finite per-variant ratio contributions remain after ",
           "the r_av trim denominator check.")
    }
    lim <- stats::quantile(ratio[finite], c(trim, 1 - trim), na.rm = TRUE)
    keep <- finite & ratio >= lim[1] & ratio <= lim[2]
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
#' @section Three ways to get this wrong, and what stops them:
#' \strong{1. A folded MAF in \code{eaf}.} \code{maf} is the folded magnitude the
#' guard uses; \code{eaf} must be the \emph{signed} frequency of the allele
#' \code{beta_mu} counts. Taking both from the same column makes the correction
#' add the bias rather than remove it --- 0.916 to \strong{1.000} at
#' \eqn{5\times10^{-8}} for MAF 0.01, \eqn{Z_\mu = 30}, against
#' \eqn{1.5\times10^{-9}} done right. Passing the same object, or element-wise
#' identical values, to \code{maf} and \code{eaf} is refused; so is an \code{eaf}
#' column that lies below 0.5 at every variant; and \code{eaf_orientation} must
#' be set deliberately.
#'
#' \strong{2. A \code{W} that is merely \emph{named} \code{beta_sigma2}.} The
#' correction is exact for \code{\link{make_weights_vqtl}}'s numbers, not for its
#' column names. By default this function verifies \code{fit$W} elementwise with
#' \code{\link{check_vqtl_design}} and refuses otherwise, which needs
#' \code{W_sd}.
#'
#' \strong{3. Reading the corrected column without its caveats.} They are in
#' \code{\link{correct_dispersion_skew}}: a ~26\% power loss at MAF 0.01 /
#' \eqn{Z_\mu = 30}, a usable skewness envelope of about 0.65 under GLS (about
#' 1.0 under OLS), and \code{w_sd_ratio} having to match the \code{sd} that built
#' \code{W}.
#'
#' @param fit Output of \code{param_gwas()} run with \code{return_cov = "upper"} and
#'   a weight matrix carrying columns named \code{beta_mu} and \code{beta_sigma2}
#'   (i.e. \code{make_weights_vqtl}). Its \code{W} element is what the design
#'   check verifies.
#' @param r_av,se_rav Genome-wide coupling and its SE; see \code{\link{estimate_rav}}.
#'   If \code{r_av} is \code{NULL} the residualised columns are omitted.
#' @param N Sample size, for the skewness guard.
#' @param maf,n_eff Allele frequency, or per-variant \code{d_resid}; see
#'   \code{\link{dispersion_skew_guard}}. If both are \code{NULL} the guard is
#'   omitted. \code{maf} is \emph{folded} and is not a substitute for \code{eaf}.
#' @param eaf Signed effect-allele frequency, oriented to the allele
#'   \code{beta_mu} counts. When supplied, the genotype-skewness-corrected
#'   dispersion columns are added; see \code{\link{correct_dispersion_skew}}.
#' @param eaf_orientation Required whenever \code{eaf} is supplied. Set to
#'   \code{"effect"} to state deliberately that \code{eaf} is the effect-allele
#'   frequency. There is no default: Stage 1 emits no frequency column, so the
#'   one you have almost certainly came from a source reporting folded MAF, and
#'   the sign of \eqn{k = 1 - 2\,\mathrm{eaf}} is what the whole correction turns
#'   on.
#' @param mu3_over_var Optional direct \eqn{\mu_3/\mathrm{Var}(g)}, overriding
#'   \code{eaf}. Use this for dosages, for non-HWE genotype laws, or when you
#'   have the signed skewness but not a signed frequency.
#' @param W_sd The \code{sd} argument that was passed to
#'   \code{\link{make_weights_vqtl}} when the fit's \code{W} was built. Required
#'   when \code{eaf} or \code{mu3_over_var} is supplied and
#'   \code{check_design = TRUE}.
#' @param w_sd_ratio \code{W_sd} divided by the phenotype SD. Leave at 1 for the
#'   usual \code{make_weights_vqtl(taus, q_tau, mu = mean(y), sd = sd(y))}; see
#'   \code{\link{dispersion_skew_bias}}.
#' @param taus Optional tau grid for the design check; recovered from
#'   \code{rownames(fit$W)} when absent.
#' @param check_design Verify \code{fit$W} numerically before correcting.
#'   Default \code{TRUE}. Setting it \code{FALSE} is a deliberate act and warns
#'   prominently.
#' @param debias Passed to \code{\link{correct_dispersion_skew}}.
#' @param df Optional degrees of freedom for the parameter reference. When
#'   omitted, the calibrated reference carried by \code{fit} is used. If the
#'   fit has no reference distribution and \code{df} is not supplied, p-values
#'   are returned as \code{NA} rather than silently using a Normal reference.
#'   The same df is used for the derived residual and skewness-corrected
#'   statistics; their extra estimated-correction uncertainty is represented in
#'   their SEs, while this shared t reference is an approximation.
#' @param mu_name,sigma2_name Parameter row names to use.
#'
#' @return A data frame with one row per variant. The uncorrected p-value reuses
#'   the calibrated reference carried by \code{fit}; when no supported reference
#'   is present it is \code{NA}. When the skewness correction is applied,
#'   \code{attr(result, "na_counts")} carries the count of variants left
#'   uncorrected and why.
#' @export
dispersion_report <- function(fit, r_av = NULL, se_rav = 0, N = NULL,
                              maf = NULL, n_eff = NULL, eaf = NULL,
                              eaf_orientation = NULL, mu3_over_var = NULL,
                              W_sd = NULL, w_sd_ratio = 1, taus = NULL,
                              check_design = TRUE, debias = TRUE,
                              mu_name = "beta_mu", sigma2_name = "beta_sigma2",
                              df = NULL) {
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

  ## The fit already carries the calibrated reference used by param_gwas().
  ## Re-use its p-values when no override is requested, and thread its degrees
  ## of freedom into every derived dispersion statistic.  A report assembled
  ## from an old fit with no reference must remain explicit: pnorm() here would
  ## recreate the silent tail error this wrapper is intended to prevent.
  fit_df <- NULL
  if (is.null(df) && !is.null(fit$param_ref_df)) {
    ref_df <- as.matrix(fit$param_ref_df)
    if (!is.null(rownames(ref_df)) && sigma2_name %in% rownames(ref_df)) {
      fit_df <- as.numeric(ref_df[sigma2_name, ])
    }
  }
  report_df <- if (!is.null(df)) as.numeric(df) else fit_df
  if (!is.null(report_df)) {
    if (length(report_df) == 1L) report_df <- rep(report_df, length(v))
    if (length(report_df) != length(v)) {
      stop("df must be length 1 or the number of fitted variants.")
    }
  } else if (is.null(df) && identical(fit$param_ref_dist, "normal")) {
    ## Explicit normal-reference fits are retained for compatibility, but are
    ## represented as df = Inf so downstream helpers use the same reference.
    report_df <- rep(Inf, length(v))
  }

  p_sigma2 <- rep(NA_real_, length(v))
  if (!is.null(df)) {
    p_sigma2 <- 2 * stats::pt(-abs(v / sv), df = report_df)
  } else if (!is.null(fit$P_params)) {
    pp <- as.matrix(fit$P_params)
    if (!is.null(rownames(pp)) && sigma2_name %in% rownames(pp)) {
      p_sigma2 <- as.numeric(pp[sigma2_name, ])
    }
  } else if (!is.null(report_df)) {
    p_sigma2 <- 2 * stats::pt(-abs(v / sv), df = report_df)
  }

  out <- data.frame(variant = colnames(fit$params),
                    beta_mu = a, se_mu = sa, z_mu = a / sa,
                    beta_sigma2 = v, se_sigma2 = sv, z_sigma2 = v / sv,
                    cov_mu_sigma2 = cl,
                    p_sigma2 = p_sigma2,
                    stringsAsFactors = FALSE)

  if (!is.null(r_av)) {
    rr <- residualise_dispersion(a, v, sa, sv, cl, r_av, se_rav,
                                 df = report_df)
    if (is.null(report_df)) rr$p_sigma2_resid[] <- NA_real_
    out <- cbind(out, rr)
  }

  nac <- NULL
  if (!is.null(eaf) || !is.null(mu3_over_var)) {
    if (!is.null(eaf)) {
      .refuse_maf_as_eaf(eaf, maf, length(a))
      if (is.null(eaf_orientation)) {
        stop("dispersion_report(): `eaf` was supplied but `eaf_orientation` was ",
             "not. Set eaf_orientation = \"effect\" to state deliberately that ",
             "`eaf` is the frequency of the allele beta_mu counts. There is no ",
             "default because Stage 1 emits no frequency column, so the one you ",
             "have very likely came from a source reporting a FOLDED minor-allele ",
             "frequency -- and ", .WRONG_WAY, " If you have only a folded MAF, ",
             "pass mu3_over_var = the signed mu3/Var(g), or use ",
             "dispersion_skew_guard(), which needs only the magnitude.",
             call. = FALSE)
      }
    }
    ## GUARDRAIL 2: the correction is exact for make_weights_vqtl's NUMBERS.
    if (isTRUE(check_design)) {
      check_vqtl_design(fit$W, sd = W_sd, taus = taus,
                        mu_name = mu_name, sigma2_name = sigma2_name)
    } else {
      warning("dispersion_report(): check_design = FALSE, so the design was NOT ",
              "verified. The genotype-skewness correction is exact only for ",
              "make_weights_vqtl's numeric design; a W merely NAMED '",
              sigma2_name, "' is silently wrong. The same Normal model ",
              "parameterised in the SD rather than the variance over-subtracts ",
              "by a factor of 2*sd and, at sd = 1.3, drives a truly null estimate from +1.51e-02 to ",
              "-2.41e-02 -- 1.6x the original bias with the sign flipped; the ",
              "5-parameter mixture design over-subtracts by 2.4x.",
              call. = FALSE)
    }
    cs  <- correct_dispersion_skew(a, v, sa, sv, cl, eaf = eaf,
                                   eaf_orientation = eaf_orientation,
                                   mu3_over_var = mu3_over_var,
                                   debias = debias, w_sd_ratio = w_sd_ratio,
                                   df = report_df)
    if (is.null(report_df)) cs$p_sigma2_skewcorr[] <- NA_real_
    nac <- attr(cs, "na_counts")
    out <- cbind(out, cs)
  }

  if (!is.null(maf) || !is.null(n_eff)) {
    out$G_skew_guard <- dispersion_skew_guard(out$z_mu, N = N, maf = maf, n_eff = n_eff)
    out$flag_skew <- !is.na(out$G_skew_guard) & abs(out$G_skew_guard) >= 0.2
  }
  rownames(out) <- NULL
  if (!is.null(nac)) attr(out, "na_counts") <- nac
  out
}


## GUARDRAIL 1, at the wrapper that takes BOTH frequencies: refuse when the
## folded `maf` has been handed to `eaf` as well.  This is the single most
## likely way to get the correction pointed the wrong way.
.refuse_maf_as_eaf <- function(eaf, maf, n) {
  if (is.null(maf)) return(invisible(NULL))
  same_obj <- identical(eaf, maf)
  same_val <- FALSE
  if (!same_obj) {
    em <- suppressWarnings(as.numeric(eaf))
    mm <- suppressWarnings(as.numeric(maf))
    if (length(em) && length(mm)) {
      L  <- max(length(em), length(mm), n)
      em <- rep_len(em, L); mm <- rep_len(mm, L)
      same_val <- identical(is.na(em), is.na(mm)) &&
                  isTRUE(all.equal(em, mm, tolerance = 0))
    }
  }
  if (!(same_obj || same_val)) return(invisible(NULL))
  stop("dispersion_report(): `eaf` and `maf` are ",
       if (same_obj) "the same object" else
         "element-wise identical after recycling",
       ". They are different things. `maf` is the FOLDED minor-allele frequency, ",
       "and the skewness guard needs only its magnitude. `eaf` must be the ",
       "SIGNED frequency of the allele beta_mu counts, and ", .WRONG_WAY,
       " Supply a properly oriented effect-allele frequency, or pass ",
       "mu3_over_var = the signed mu3/Var(g). If the effect allele genuinely is ",
       "the minor allele at every variant, say so explicitly by passing that ",
       "column as `eaf` with eaf_orientation = \"effect\" and leaving `maf` ",
       "unset (the guard will fall back on n_eff).", call. = FALSE)
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
#'     \code{sd} passed to \code{make_weights_vqtl} was not the phenotype SD. It
#'     must match that \code{sd}: leaving it at 1 when \code{W} was built on
#'     another scale scales the whole correction by exactly that ratio, with no
#'     symptom.)
#'   \item Correcting \eqn{b(\tau)} before the projection and correcting the
#'     coefficient after it are \emph{identical} \strong{on this design}, not
#'     merely similar, because the bias vector lies exactly in the column space
#'     of \eqn{W}. Correcting afterwards is chosen because it is O(1) per variant
#'     and needs no refit. Outside \code{make_weights_vqtl} the equality fails,
#'     which is why \code{\link{check_vqtl_design}} gates on the numeric identity
#'     of \eqn{W} rather than on its column names.
#' }
#'
#' \strong{Order of the expansion, and the condition on the tau grid.} The
#' derivation is second order in \eqn{\delta}. The third-order term is even in
#' \eqn{z_\tau}, so \emph{on a tau grid symmetric about 0.5} it is orthogonal to
#' the dispersion column and lands on \eqn{\beta_\mu} and on the fit residual
#' instead. The first term that survives into \eqn{\beta_{\sigma^2}} is then
#' \eqn{O(\delta^4)}: measured against an exact population calculation the
#' relative error of this correction is \eqn{\approx 0.065\,(\delta/\sigma)^2},
#' i.e. 0.07\% at 0.1 SD per allele and 0.3\% at 0.2 SD per allele.
#'
#' The symmetry is a real condition, not a technicality. On
#' \code{seq(0.02, 0.70, length.out = 45)} the orthogonality is lost: the
#' exponent drops from 4.00 to 3.01 and the relative error grows about 70-fold,
#' with a sign that flips depending on which tail is truncated.
#' \code{\link{check_vqtl_design}} warns when the grid is asymmetric.
#'
#' @param beta_mu Numeric vector of estimated location coefficients.
#' @param se_mu Numeric vector of their standard errors. Required when
#'   \code{debias = TRUE}.
#' @param eaf Frequency of the allele that \code{beta_mu} counts --- the
#'   \strong{effect} allele frequency, not the folded MAF. The sign of
#'   \eqn{1 - 2p} is the sign of the bias, so a folded MAF is not merely
#'   insufficient, it corrects in the wrong direction and is worse than no
#'   correction; see \code{\link{correct_dispersion_skew}}. A column lying below
#'   0.5 at every variant is refused as a folded-MAF signature. Length 1 or
#'   \code{length(beta_mu)}.
#' @param eaf_orientation Set to \code{"effect"} to assert deliberately that
#'   \code{eaf} is oriented to the allele \code{beta_mu} counts. \code{NULL}
#'   (default) leaves the assertion unmade: a conclusive folded-MAF signature is
#'   then refused rather than warned about, and an inconclusive one warns.
#'   \code{"folded"} is accepted as an explicit statement and always refused ---
#'   use \code{\link{dispersion_skew_guard}}, which needs only the magnitude.
#' @param n_min_orient Number of distinct \code{eaf} values below 0.5 at which
#'   the folded-MAF signature becomes conclusive. Default 10.
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
                                 eaf_orientation = NULL,
                                 mu3_over_var = NULL, debias = TRUE,
                                 w_sd_ratio = 1, n_min_orient = 10L) {
  beta_mu <- as.numeric(beta_mu)
  n <- length(beta_mu)
  k <- .skew_k(eaf, mu3_over_var, n, eaf_orientation = eaf_orientation,
               n_min_orient = n_min_orient, where = "dispersion_skew_bias()")
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

## ---------------------------------------------------------------------------
## GUARDRAIL 1: the EAF orientation.
##
## k = 1 - 2*eaf is SIGNED.  Hand this a FOLDED minor-allele frequency at a
## variant whose effect allele is the major allele and k flips sign, so the
## correction ADDS the bias instead of removing it.  Measured: a variant at
## |Z| = 3.27 goes to |Z| = 6.20, p from 1.1e-03 to 5.6e-10 -- a sub-threshold
## artefact promoted to a genome-wide-significant false vQTL.  At MAF 0.01 with
## Z_mu 30 the rate at 5e-8 goes from 9.16e-01 uncorrected to 1.000 with a folded
## k, against 1.54e-09 with the correctly oriented one: the wrong orientation is
## strictly worse than not correcting at all.
##
## It is the MOST LIKELY user error, because Stage 1 emits no frequency column,
## so the frequency almost certainly comes from a source reporting folded MAF.
## Three things stand in the way:
##   * a signature check -- a genuine EAF column straddles 0.5 over any real set
##     of variants, a folded one cannot;
##   * `eaf_orientation`, which the caller must set deliberately;
##   * `dispersion_report()` refusing outright when `eaf` and `maf` are the same
##     values.
## ---------------------------------------------------------------------------

.EAF_ORIENT_EFFECT <- c("effect", "effect_allele", "signed")
.EAF_ORIENT_FOLDED <- c("folded", "maf", "minor", "minor_allele", "unknown")

## Does a frequency vector carry the signature of a FOLDED minor-allele column?
##   "ok"           -- some value is >= 0.5, so it straddles: cannot be folded
##   "folded"       -- every value below 0.5, over enough distinct values to be sure
##   "unverifiable" -- every value below 0.5, but too few distinct values to tell
.eaf_orientation_signature <- function(p, n_min = 10L) {
  fin <- p[is.finite(p)]
  if (!length(fin)) return("ok")
  if (any(fin >= 0.5)) return("ok")
  if (length(unique(fin)) >= n_min) return("folded")
  "unverifiable"
}

.WRONG_WAY <- paste0(
  "k = 1 - 2*eaf is SIGNED, so an unoriented frequency makes this function ADD ",
  "the bias rather than remove it. Measured: a variant at |Z| = 3.27 becomes ",
  "|Z| = 6.20, p 1.1e-03 -> 5.6e-10 -- a sub-threshold artefact promoted to a ",
  "genome-wide-significant false vQTL. At MAF 0.01 / Z_mu 30 the rate at 5e-8 ",
  "goes to 1.000, against 9.2e-01 with no correction at all and 1.5e-09 with ",
  "the right orientation.")

## Validate the caller's orientation claim and the data's own signature.
.check_eaf_orientation <- function(p, eaf_orientation, n_min_orient = 10L,
                                   where = "correct_dispersion_skew()") {
  asserted <- FALSE
  if (!is.null(eaf_orientation)) {
    o <- tolower(as.character(eaf_orientation)[1])
    if (o %in% .EAF_ORIENT_FOLDED) {
      stop(where, ": eaf_orientation = '", o, "' states that the supplied ",
           "frequency is NOT oriented to the allele beta_mu counts. This ",
           "correction cannot use such a column. ", .WRONG_WAY, " Supply the ",
           "effect-allele frequency, or pass mu3_over_var = mu3/Var(g) directly. ",
           "dispersion_skew_guard() needs only the magnitude and does accept a ",
           "folded MAF.", call. = FALSE)
    }
    if (!(o %in% .EAF_ORIENT_EFFECT)) {
      stop(where, ": eaf_orientation must be one of '",
           paste(c(.EAF_ORIENT_EFFECT, .EAF_ORIENT_FOLDED), collapse = "', '"),
           "', or NULL. Got '", o, "'.", call. = FALSE)
    }
    asserted <- TRUE
  }

  sig <- .eaf_orientation_signature(p, n_min_orient)
  nfin <- sum(is.finite(p))
  if (identical(sig, "folded")) {
    m <- paste0(
      where, ": every one of the ", nfin, " supplied `eaf` values lies below ",
      "0.5. That is the signature of a FOLDED minor-allele frequency, not of an ",
      "effect-allele frequency -- a genuine EAF column straddles 0.5 over any ",
      "real set of variants. ", .WRONG_WAY)
    if (asserted) {
      warning(m, " You set eaf_orientation = 'effect', so this is a warning ",
              "rather than a refusal; check the column before you believe the ",
              "corrected p-values.", call. = FALSE)
    } else {
      stop(m, " Supply the frequency of the allele beta_mu counts, or pass ",
           "mu3_over_var directly. If this really is an effect-allele frequency ",
           "that happens to be below 0.5 at every variant, set ",
           "eaf_orientation = 'effect' to assert that and proceed with a warning.",
           call. = FALSE)
    }
  } else if (identical(sig, "unverifiable") && !asserted) {
    warning(where, ": the supplied `eaf` is below 0.5 everywhere but has too ",
            "few distinct values (", length(unique(p[is.finite(p)])), ") for the ",
            "folded-MAF check to be conclusive, so the orientation could NOT be ",
            "verified. ", .WRONG_WAY, " Set eaf_orientation = 'effect' to confirm ",
            "the frequency is oriented to the allele beta_mu counts and silence ",
            "this warning.", call. = FALSE)
  }
  invisible(sig)
}

## mu3/Var(g), from an effect-allele frequency under Hardy-Weinberg, or given.
.skew_k <- function(eaf, mu3_over_var, n, eaf_orientation = NULL,
                    n_min_orient = 10L, where = "correct_dispersion_skew()") {
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
  .check_eaf_orientation(p, eaf_orientation, n_min_orient = n_min_orient,
                         where = where)
  k <- rep(NA_real_, n)
  ok <- is.finite(p) & p > 0 & p < 1
  k[ok] <- 1 - 2 * p[ok]
  k
}


## ---------------------------------------------------------------------------
## GUARDRAIL 2: gate on the numeric identity of W, not on its column names.
##
## The correction is exact because the bias vector IS delta^2 k times
## `make_weights_vqtl()`'s dispersion column -- not merely proportional to it,
## the same vector -- so AW = I sends all of it onto beta_sigma2 and none onto
## beta_mu.  That argument is about the NUMBERS in W, and a check on column
## NAMES does not test it.  Two designs that pass a name check and are wrong:
##
##  * the same Normal model built with `make_weights_generic` parameterised in
##    the SD rather than the variance, whose column is z_tau/sd rather than
##    z_tau/(2 sd).  The column ratio is 2*sd, so the correction over-subtracts
##    by that factor: at sd = 1.3 it drives the estimate from +1.51e-02
##    (truth 0) to -2.41e-02, 1.6x the original bias with the sign flipped.  A shape check alone does NOT catch this -- the column is still
##    proportional to z_tau -- which is why the gate needs the `sd` that built W.
##  * the 5-parameter Normal-mixture design with mu_name = "beta_mu1", which
##    over-subtracts by 2.4x.
##
## `check_vqtl_design()` therefore compares W elementwise against
## qnorm(taus)/(2*sd) and against the constant 1, and refuses otherwise.
## ---------------------------------------------------------------------------

#' Verify that a weight matrix really is \code{make_weights_vqtl}'s design
#'
#' The genotype-skewness correction is exact for one specific design and silently
#' wrong for others that carry the same column names. This checks the numeric
#' identity rather than the names, and is the gate
#' \code{\link{dispersion_report}} runs before it corrects anything.
#'
#' @details
#' Two conditions are required, and both are checked elementwise:
#' \deqn{W[, \mathrm{mu\_name}] = 1, \qquad
#'       W[, \mathrm{sigma2\_name}] = \Phi^{-1}(\tau)/(2\,\mathrm{sd}).}
#' The second cannot be reduced to a shape check. A Normal location/scale design
#' parameterised in the \emph{standard deviation} rather than the variance gives
#' \eqn{\Phi^{-1}(\tau)/\mathrm{sd}}, which is still proportional to
#' \eqn{\Phi^{-1}(\tau)} and still passes any check that does not know the scale;
#' correcting on it over-subtracts by a factor of \eqn{2\,\mathrm{sd}} and, at
#' \eqn{\mathrm{sd} = 1.3}, drives a truly null estimate from
#' \eqn{+1.51\times10^{-2}} to \eqn{-2.41\times10^{-2}}. This is
#' why \code{sd} is a required argument: pass the \emph{same} \code{sd} you
#' passed to \code{\link{make_weights_vqtl}}.
#'
#' Extra columns are permitted. \eqn{AW = I} holds for the whole design, so a
#' bias vector equal to a scalar times one verified column still projects onto
#' that column's coefficient alone.
#'
#' \strong{The tau grid.} The \eqn{O(\delta^3)} term is even in
#' \eqn{\Phi^{-1}(\tau)} and therefore orthogonal to the dispersion column
#' \emph{only on a tau grid symmetric about 0.5}. On an asymmetric grid it leaks
#' in: on \code{seq(0.02, 0.70, length.out = 45)} the order of the residual in
#' \eqn{\delta} falls from 4.00 to 3.01 and its relative size grows about
#' 70-fold, with a sign that depends on which tail was truncated. An asymmetric
#' grid is warned about, not refused.
#'
#' @param W The weight matrix the fit used --- \code{param_gwas()} returns it as
#'   \code{fit$W}.
#' @param sd The \code{sd} argument that was passed to
#'   \code{\link{make_weights_vqtl}} when \code{W} was built. Required: the check
#'   is not identifiable without it.
#' @param taus Quantile levels, length \code{nrow(W)}. If \code{NULL} they are
#'   recovered from \code{rownames(W)}, which \code{make_weights_vqtl} sets to
#'   \code{paste0("tau", taus)}.
#' @param mu_name,sigma2_name Column names of the location and dispersion columns.
#' @param tol Tolerance. \code{W[, mu_name]} must be within \code{tol} of 1 in
#'   absolute terms; \code{W[, sigma2_name]} within \code{tol} of the target
#'   relative to the largest target entry.
#' @param warn_asymmetric_taus Warn when the tau grid is not symmetric about 0.5.
#'
#' @return Invisibly, a list with \code{ok}, \code{taus}, \code{sd},
#'   \code{dev_mu}, \code{dev_sigma2}, \code{sd_implied} and
#'   \code{taus_symmetric}. Stops with a diagnostic message if the design does
#'   not match.
#'
#' @examples
#' taus <- seq(0.05, 0.95, length.out = 21)
#' W <- make_weights_vqtl(taus, q_tau = 1.3 * qnorm(taus), mu = 0, sd = 1.3)
#' check_vqtl_design(W, sd = 1.3)
#' @seealso \code{\link{correct_dispersion_skew}}, \code{\link{dispersion_report}}
#' @export
check_vqtl_design <- function(W, sd, taus = NULL,
                              mu_name = "beta_mu", sigma2_name = "beta_sigma2",
                              tol = 1e-8, warn_asymmetric_taus = TRUE) {
  if (is.null(W)) {
    stop("check_vqtl_design(): W is NULL, so the design cannot be verified. ",
         "param_gwas() returns the design it fitted as fit$W; pass that. ",
         "The genotype-skewness correction is exact only for ",
         "make_weights_vqtl's design and is silently wrong for others carrying ",
         "the same column names.", call. = FALSE)
  }
  W <- as.matrix(W)
  if (!is.numeric(W)) stop("check_vqtl_design(): W must be a numeric matrix.",
                           call. = FALSE)
  cn <- colnames(W)
  if (is.null(cn)) {
    stop("check_vqtl_design(): W has no column names, so '", mu_name, "' and '",
         sigma2_name, "' cannot be located.", call. = FALSE)
  }
  if (!all(c(mu_name, sigma2_name) %in% cn)) {
    stop("check_vqtl_design(): W must carry columns '", mu_name, "' and '",
         sigma2_name, "'. Found: ", paste(cn, collapse = ", "), call. = FALSE)
  }
  if (missing(sd) || is.null(sd) || !is.numeric(sd) || length(sd) != 1L ||
      !is.finite(sd) || sd <= 0) {
    stop("check_vqtl_design(): `sd` must be a single positive number -- the ",
         "same `sd` you passed to make_weights_vqtl() when you built W. ",
         "Without it the check is not identifiable: a design parameterised in ",
         "the SD rather than the variance gives qnorm(taus)/sd, which is still ",
         "proportional to qnorm(taus) and would pass any scale-free check while ",
         "making the correction over-subtract by a factor of 2*sd.",
         call. = FALSE)
  }

  if (is.null(taus)) {
    rn <- rownames(W)
    if (!is.null(rn) && all(grepl("^tau", rn))) {
      taus <- suppressWarnings(as.numeric(sub("^tau", "", rn)))
    }
  }
  if (is.null(taus) || length(taus) != nrow(W) || anyNA(taus)) {
    stop("check_vqtl_design(): the tau grid could not be recovered from ",
         "rownames(W) (make_weights_vqtl sets them to paste0(\"tau\", taus)). ",
         "Pass `taus` explicitly.", call. = FALSE)
  }
  taus <- as.numeric(taus)
  if (any(taus <= 0 | taus >= 1)) {
    stop("check_vqtl_design(): taus must all lie strictly inside (0, 1).",
         call. = FALSE)
  }

  z      <- stats::qnorm(taus)
  target <- z / (2 * sd)
  wmu    <- as.numeric(W[, mu_name])
  wsg    <- as.numeric(W[, sigma2_name])

  dev_mu    <- max(abs(wmu - 1))
  scale_sg  <- max(abs(target), .Machine$double.eps)
  dev_sigma <- max(abs(wsg - target)) / scale_sg

  ## best-fitting proportional column, for a message that says what went wrong
  ssz        <- sum(z * z)
  cfit       <- if (ssz > 0) sum(wsg * z) / ssz else NA_real_
  prop_dev   <- if (is.finite(cfit) && abs(cfit) > 0) {
                  max(abs(wsg - cfit * z)) / max(abs(cfit * z), .Machine$double.eps)
                } else NA_real_
  sd_implied <- if (is.finite(cfit) && cfit > 0) 1 / (2 * cfit) else NA_real_

  if (dev_mu > tol) {
    stop("check_vqtl_design(): W[, \"", mu_name, "\"] is not the constant 1 ",
         "(max |W - 1| = ", format(dev_mu, digits = 3), ", tol = ", tol, "). ",
         "The genotype-skewness correction assumes beta_mu is the coefficient ",
         "on a column of ones, i.e. that it is the per-allele location shift ",
         "delta in phenotype units. It is not, so delta cannot be read off this ",
         "fit and the correction does not apply. Use make_weights_vqtl(), or ",
         "pass mu3_over_var and correct b(tau) yourself before the projection.",
         call. = FALSE)
  }

  if (dev_sigma > tol) {
    diag_msg <- if (is.finite(prop_dev) && prop_dev <= tol) {
      paste0("It IS proportional to qnorm(taus), but implies sd = ",
             format(sd_implied, digits = 6), " rather than the sd = ",
             format(sd, digits = 6), " you passed (ratio ",
             format(sd / sd_implied, digits = 6), "). Either `sd` is wrong, or ",
             "W was parameterised in the standard deviation rather than the ",
             "variance -- d q_tau / d sd = qnorm(tau) against ",
             "d q_tau / d sigma^2 = qnorm(tau)/(2 sd) -- whose column ratio is ",
             "exactly 2*sd, and which makes the correction over-subtract by that ",
             "same factor: at sd = 1.3 it drives a truly null estimate from ",
             "+1.51e-02 to -2.41e-02, 1.6x the original bias with the sign ",
             "flipped.")
    } else {
      paste0("It is not even proportional to qnorm(taus) (best-fit ",
             "proportional column still leaves a relative deviation of ",
             format(prop_dev, digits = 3), "). The bias residue is only ",
             "collinear with the dispersion column on a Normal location/scale ",
             "design; under a mixture or other W it is not a single column and ",
             "this correction does not apply at all -- the 5-parameter Normal ",
             "mixture with mu_name = \"beta_mu1\" over-subtracts by 2.4x.")
    }
    stop("check_vqtl_design(): W[, \"", sigma2_name, "\"] is not ",
         "qnorm(taus)/(2*sd) (max relative deviation ",
         format(dev_sigma, digits = 3), ", tol = ", tol, "). ", diag_msg,
         call. = FALSE)
  }

  st      <- sort(taus)
  sym_dev <- max(abs(st - rev(1 - st)))
  symmetric <- sym_dev <= 1e-8
  if (!symmetric && isTRUE(warn_asymmetric_taus)) {
    warning("check_vqtl_design(): the tau grid is not symmetric about 0.5 (max ",
            "deviation ", format(sym_dev, digits = 3), "). The O(delta^3) term ",
            "is even in qnorm(tau) and hence orthogonal to the dispersion ",
            "column ONLY on a symmetric grid; on an asymmetric one it leaks ",
            "into beta_sigma2. Measured on seq(0.02, 0.70, length.out = 45) the ",
            "order of the residual in delta falls from 4.00 to 3.01 and its ",
            "relative size grows about 70-fold, with a sign that depends on ",
            "which tail is truncated. The correction is still applied, but its ",
            "residual is larger than the O(delta^4) figure quoted for a ",
            "symmetric grid.", call. = FALSE)
  }

  invisible(list(ok = TRUE, taus = taus, sd = sd, dev_mu = dev_mu,
                 dev_sigma2 = dev_sigma, sd_implied = sd_implied,
                 taus_symmetric = symmetric))
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
#' \strong{When this is invalid. Read this before using the corrected column.}
#' \itemize{
#'   \item \strong{It needs a signed effect-allele frequency, and getting that
#'     wrong is worse than not correcting.} \eqn{k = 1 - 2p} is odd about
#'     \eqn{p = 1/2}, so a \emph{folded} MAF at a variant whose effect allele is
#'     the major allele makes this function \strong{add} the bias instead of
#'     removing it: a variant at \eqn{|Z| = 3.27} goes to \eqn{|Z| = 6.20},
#'     \eqn{p} from \eqn{1.1\times10^{-3}} to \eqn{5.6\times10^{-10}}. At
#'     MAF 0.01 with \eqn{Z_\mu = 30} the rate at \eqn{5\times10^{-8}} goes to
#'     \strong{1.000} with a folded \eqn{k}, against 0.916 with no correction at
#'     all and \eqn{1.5\times10^{-9}} with the right one. Stage 1 emits no
#'     frequency column, so \code{eaf} must come from the genotype source and is
#'     very likely to arrive folded. A column that lies below 0.5 at every
#'     variant is refused; set \code{eaf_orientation = "effect"} to assert the
#'     orientation deliberately.
#'   \item \strong{It is specific to \code{make_weights_vqtl}}, and the check is
#'     on the \emph{numbers} in \code{W}, not its column names. The same Normal
#'     model built with \code{\link{make_weights_generic}} parameterised in the
#'     SD rather than the variance has a dispersion column
#'     \eqn{z_\tau/\mathrm{sd}} instead of \eqn{z_\tau/(2\,\mathrm{sd})}; the
#'     correction then over-subtracts by a factor of \eqn{2\,\mathrm{sd}} and,
#'     at \eqn{\mathrm{sd} = 1.3}, drives a truly null estimate from
#'     \eqn{+1.51\times10^{-2}} to \eqn{-2.41\times10^{-2}}, i.e. 1.6x the
#'     original bias with the sign flipped. The 5-parameter Normal-mixture
#'     design with \code{mu_name = "beta_mu1"} over-subtracts by 2.4x. Pass
#'     \code{W} and \code{W_sd} here, or use \code{\link{dispersion_report}},
#'     which gates on \code{\link{check_vqtl_design}} by default.
#'   \item \strong{\code{w_sd_ratio} must match the \code{sd} that built
#'     \code{W}.} It is (the \code{sd} argument passed to
#'     \code{\link{make_weights_vqtl}}) divided by (the phenotype SD). With the
#'     usual \code{make_weights_vqtl(taus, q_tau, mu = mean(y), sd = sd(y))} it
#'     is 1 and the phenotype SD cancels. If \code{W} was deliberately built on
#'     another scale and \code{w_sd_ratio} is left at 1, the correction is wrong
#'     by exactly that ratio --- linearly, with no symptom.
#'   \item \strong{The baseline must be close to Normal, and how close depends on
#'     the estimator.} These two are coupled and the choice of estimator moves
#'     the limit. Corrected rates at \eqn{5\times10^{-8}}, MAF 0.01,
#'     \eqn{Z_\mu = 30}, \strong{GLS}: skewness 0.63 gives
#'     \eqn{3.8\times10^{-10}}; 0.71 gives \eqn{3.6\times10^{-9}}; 0.82 gives
#'     \eqn{5.0\times10^{-7}}; 0.89 gives \eqn{4.8\times10^{-5}}; 0.95 gives
#'     \strong{1.00}. The usable envelope under GLS is therefore about
#'     \strong{|skewness| 0.65}, not the 1.4 the uncorrected bracket suggests.
#'     Under OLS it runs out to about 1.0. Rank-inverse-normal or log transform
#'     the phenotype (see \code{docs/phenotype_preparation.md}) and this is not a
#'     live concern.
#'   \item \strong{The failure is not summarised by skewness alone, nor by
#'     kurtosis alone.} At matched skewness 0.80 a skew-normal baseline (excess
#'     kurtosis 0.65) leaves \eqn{2.0\times10^{-4}} under GLS while
#'     \eqn{\chi^2_{12}} (excess kurtosis 1.00) leaves
#'     \eqn{5.0\times10^{-7}} --- 390x apart, with the \emph{lower}-kurtosis
#'     distribution the worse of the two. Symmetric \eqn{t_6}, skewness exactly
#'     0, still leaves a 62x residual. Neither moment is the determining
#'     property; what matters is how far \eqn{f_0'/f_0} departs from
#'     proportionality to \eqn{z_\tau}, which no single moment captures.
#'   \item \strong{The tau grid should be symmetric about 0.5.} The
#'     \eqn{O(\delta^3)} term is orthogonal to the dispersion column only then;
#'     see \code{\link{check_vqtl_design}} for the measured cost of an
#'     asymmetric grid.
#'   \item It is second order in \eqn{\delta}; see
#'     \code{\link{dispersion_skew_bias}} for the residual.
#' }
#'
#' \strong{What the correction costs.} Essentially nothing on null variants and
#' nothing measurable on a genuine dispersion effect with no location effect
#' (power at \eqn{5\times10^{-8}} unchanged in the fifth decimal). The one real
#' cost is at a rare variant carrying a colossal location effect, where the
#' correction has genuine variance of its own: \strong{at MAF 0.01 with
#' \eqn{Z_\mu = 30} a marginal dispersion effect loses about 26\% of its honest
#' power} (0.372 against an honest 0.50 at \eqn{5\times10^{-8}}; the SE widens by
#' a factor 1.094). At MAF 0.01 with \eqn{Z_\mu = 16} the loss is about 12\%
#' (0.562 against 0.50 honest, SE ratio 1.025 --- here the corrected estimator is
#' still slightly over-powered). That worst cell is the same cell where the
#' \emph{uncorrected} rate is 0.916, and \code{\link{dispersion_skew_guard}}
#' flags it with \eqn{G = 6.3}.
#'
#' \strong{Missing corrections are counted, not hidden.} An \code{eaf} of 0, 1 or
#' \code{NA} gives no \eqn{k}; a \code{cov_mu_sigma2} taken from a different fit
#' than the standard errors violates Cauchy-Schwarz and yields a negative
#' delta-method variance (at \eqn{|c| = 1.01\,\mathrm{se}_\mu
#' \mathrm{se}_{\sigma^2}} that is 15,316 variants per 4M). Both produce
#' \code{NA}, both are counted, a warning names the counts, and the counts are
#' attached as \code{attr(result, "na_counts")}.
#'
#' @param beta_mu Numeric vector of estimated location coefficients; its length
#'   sets the number of variants.
#' @param beta_sigma2 Estimated dispersion coefficients, length 1 or
#'   \code{length(beta_mu)}.
#' @param se_mu,se_sigma2 Their standard errors, length 1 or \code{length(beta_mu)}.
#' @param cov_mu_sigma2 Within-variant sampling covariance between the two, from
#'   \code{param_cov_upper}. If \code{NULL} it is taken as zero and a warning is
#'   issued.
#' @param eaf,eaf_orientation,mu3_over_var,debias,w_sd_ratio See
#'   \code{\link{dispersion_skew_bias}}.
#' @param W Optional weight matrix the fit used (\code{fit$W}). When supplied it
#'   is verified against \code{\link{check_vqtl_design}} and the call refuses if
#'   the design is not \code{make_weights_vqtl}'s. This function takes
#'   coefficients rather than a fit, so it \emph{cannot} verify the design unless
#'   you pass it; \code{\link{dispersion_report}} does so by default.
#' @param W_sd The \code{sd} argument that built \code{W}. Required when
#'   \code{W} is supplied.
#' @param taus Optional tau grid for the \code{W} check; recovered from
#'   \code{rownames(W)} when absent.
#' @param n_min_orient Number of distinct \code{eaf} values below 0.5 at which
#'   the folded-MAF signature becomes conclusive. Default 10.
#' @param df Optional degrees of freedom for a t reference; Normal if \code{NULL}.
#'
#' @return A data frame with \code{bias_skew} (the amount subtracted),
#'   \code{beta_sigma2_skewcorr}, \code{se_sigma2_skewcorr},
#'   \code{z_sigma2_skewcorr}, \code{p_sigma2_skewcorr}, carrying an
#'   \code{"na_counts"} attribute.
#'
#' @examples
#' # a 0.2 SD/allele location effect at EAF 0.01, no true dispersion effect.
#' # eaf_orientation asserts that 0.01 is the frequency of the allele beta_mu
#' # counts, not a folded MAF -- with one variant the check cannot tell.
#' correct_dispersion_skew(beta_mu = 0.2, beta_sigma2 = 0.0392,
#'                         se_mu = 0.01, se_sigma2 = 0.012,
#'                         cov_mu_sigma2 = 0, eaf = 0.01,
#'                         eaf_orientation = "effect")
#' @seealso \code{\link{dispersion_skew_bias}}, \code{\link{dispersion_skew_guard}},
#'   \code{\link{check_vqtl_design}}, \code{\link{residualise_dispersion}}
#' @export
correct_dispersion_skew <- function(beta_mu, beta_sigma2, se_mu, se_sigma2,
                                    cov_mu_sigma2 = NULL, eaf = NULL,
                                    eaf_orientation = NULL,
                                    mu3_over_var = NULL, debias = TRUE,
                                    w_sd_ratio = 1, W = NULL, W_sd = NULL,
                                    taus = NULL, n_min_orient = 10L, df = NULL) {
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

  ## GUARDRAIL 2, when the caller can supply the design.
  if (!is.null(W)) check_vqtl_design(W, sd = W_sd, taus = taus)

  ## GUARDRAIL 1 lives inside .skew_k().
  k <- .skew_k(eaf, mu3_over_var, n, eaf_orientation = eaf_orientation,
               n_min_orient = n_min_orient,
               where = "correct_dispersion_skew()") * as.numeric(w_sd_ratio)
  bias <- dispersion_skew_bias(beta_mu, se_mu = se_mu, mu3_over_var = k,
                               debias = debias, w_sd_ratio = 1)
  est <- beta_sigma2 - bias
  v <- se_sigma2^2 + 4 * k^2 * beta_mu^2 * se_mu^2 + 2 * k^2 * se_mu^4 -
       4 * k * beta_mu * cov_mu_sigma2

  ## Count what is about to become NA rather than returning it silently.
  nac <- .count_skewcorr_na(n = n, eaf = eaf, mu3_over_var = mu3_over_var, k = k,
                            v = v, se_mu = se_mu, se_sigma2 = se_sigma2,
                            cov_mu_sigma2 = cov_mu_sigma2)

  v[!is.finite(v) | v < 0] <- NA_real_
  se <- sqrt(v)
  z <- est / se
  p <- if (is.null(df)) 2 * stats::pnorm(-abs(z)) else 2 * stats::pt(-abs(z), df = df)
  out <- data.frame(bias_skew = bias, beta_sigma2_skewcorr = est,
                    se_sigma2_skewcorr = se, z_sigma2_skewcorr = z,
                    p_sigma2_skewcorr = p)
  nac$n_na <- sum(!is.finite(est) | !is.finite(se))
  attr(out, "na_counts") <- nac
  .warn_skewcorr_na(nac, where = "correct_dispersion_skew()")
  out
}


## ---------------------------------------------------------------------------
## Count and report NA corrections rather than returning them silently.
## Two distinct causes, and they mean different things:
##   * eaf of 0, 1 or NA -- no k, so no correction for that variant;
##   * a negative delta-method variance, which requires
##     |cov_mu_sigma2| > se_mu * se_sigma2.  That breaches Cauchy-Schwarz and
##     CANNOT happen when the covariance and the standard errors come from the
##     same param_gwas() fit, so it is diagnostic of a mismatched covariance
##     file.  At |c| = 1.01 * se_mu * se_sigma2 it is 15,316 variants per 4M.
## ---------------------------------------------------------------------------
.count_skewcorr_na <- function(n, eaf, mu3_over_var, k, v, se_mu, se_sigma2,
                               cov_mu_sigma2) {
  n_eaf_na <- 0L; n_eaf_bound <- 0L
  if (is.null(mu3_over_var) && !is.null(eaf)) {
    p <- rep_len(as.numeric(eaf), n)
    n_eaf_na    <- sum(is.na(p))
    n_eaf_bound <- sum(!is.na(p) & (p <= 0 | p >= 1))
  }
  cs_lim <- se_mu * se_sigma2
  list(n = n,
       n_k_unusable = sum(!is.finite(k)),
       n_eaf_na = n_eaf_na,
       n_eaf_out_of_range = n_eaf_bound,
       n_var_negative = sum(is.finite(v) & v < 0),
       n_var_nonfinite = sum(!is.finite(v)),
       n_cov_cauchy_schwarz = sum(is.finite(cov_mu_sigma2) & is.finite(cs_lim) &
                                  abs(cov_mu_sigma2) > cs_lim * (1 + 1e-12)),
       n_na = NA_integer_)
}

.warn_skewcorr_na <- function(nac, where) {
  if (!is.finite(nac$n_na) || nac$n_na == 0L) return(invisible(NULL))
  msg <- paste0(
    where, ": ", nac$n_na, " of ", nac$n, " variants (",
    format(100 * nac$n_na / nac$n, digits = 3), "%) have NO corrected estimate. ",
    "Causes: ", nac$n_k_unusable, " with an unusable eaf / mu3_over_var (",
    nac$n_eaf_na, " NA, ", nac$n_eaf_out_of_range, " outside (0,1)); ",
    nac$n_var_negative, " with a negative delta-method variance.")
  if (nac$n_cov_cauchy_schwarz > 0L) {
    msg <- paste0(
      msg, " ", nac$n_cov_cauchy_schwarz, " variants have ",
      "|cov_mu_sigma2| > se_mu * se_sigma2, which violates Cauchy-Schwarz and ",
      "cannot happen when the covariance and the standard errors come from the ",
      "SAME param_gwas() fit -- check that param_cov_upper belongs to these SEs.")
  }
  warning(msg, " Counts are on attr(<result>, \"na_counts\").", call. = FALSE)
  invisible(NULL)
}
