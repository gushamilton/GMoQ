#' Build tau-by-K weights for a 2-component Normal mixture
#'
#' Constructs a weight matrix \eqn{W} that maps RIF tau-slopes \eqn{b(\tau)}
#' into component-parameter effects under a two-Normal mixture baseline.
#'
#' Let \eqn{Y \sim p_1 N(\mu_1, \sigma_1^2) + (1-p_1) N(\mu_2, \sigma_2^2)}.
#' For baseline unconditional quantiles \eqn{q_\tau} and mixture pdf
#' \eqn{f(q_\tau)}, the default columns of \eqn{W} correspond to component means:
#' \deqn{W_1(\tau) = \frac{p_1 f_1(q_\tau)}{f(q_\tau)}, \quad
#'       W_2(\tau) = \frac{(1-p_1) f_2(q_\tau)}{f(q_\tau)}}
#'
#' If \code{include_membership = TRUE}, a first column is added for membership
#' (log-odds) perturbation:
#' \deqn{W_\gamma(\tau) = - \frac{p_1 (1-p_1)}{f(q_\tau)}
#'       \left\{F_1(q_\tau) - F_2(q_\tau)\right\}}
#'
#' @param taus Numeric vector of quantile levels.
#' @param q_tau Numeric vector of baseline quantiles at \code{taus} (\code{type=8} recommended).
#' @param p1 Proportion of class 1 in (0,1); class 2 is \code{1 - p1}.
#' @param mu1,sd1 Mean and standard deviation of component 1.
#' @param mu2,sd2 Mean and standard deviation of component 2.
#' @param include_membership Logical; include a first column \code{"gamma"} if TRUE.
#' @param tiny Small positive floor for stabilizing divisions (default: mixture pdf clamped to \code{tiny}).
#' @param density Which density divides each column; see
#'   \code{docs/MIXTURE_WEIGHT_DENSITY.md}. The default (\code{NULL}) is
#'   \code{"hybrid"} when \code{include_membership = TRUE} and \code{"model"}
#'   otherwise (since 0.6.0; before that it was always \code{"model"}).
#'   \code{"model"}: every column is divided by the fitted mixture pdf. In a
#'   well-separated mixture this attenuates \code{gamma} (by about 30\% in the
#'   1000 Genomes example) and leaks it into the component locations.
#'   \code{"hybrid"}: the component columns stay the component shares
#'   \eqn{P(k \mid q_\tau)} (they sum to 1) and the \code{"gamma"} column is
#'   divided by the Stage 1 density \code{stage1_density} instead, because the
#'   Stage 1 coefficient carries that density's scaling.
#'   \code{"stage1"}: every column is divided by \code{stage1_density}; this
#'   breaks \eqn{W_1 + W_2 = 1} and warns (sensitivity only).
#' @param stage1_density Required for \code{"hybrid"} (so by default whenever
#'   \code{include_membership = TRUE}) and \code{"stage1"}: the
#'   Stage 1 density \eqn{\hat f(q_\tau)} at \code{taus}. A numeric vector,
#'   a list with \code{fhat} (and optionally \code{taus}, \code{q_tau}, which
#'   are checked), or the path of a Stage 1 density JSON
#'   (\code{\link{read_stage1_density}}). Compute it from the phenotype with
#'   \code{\link{stage1_kde_density}}.
#' @param identity_tol Tolerance for the \eqn{W_1 + W_2 = 1} check under
#'   \code{"hybrid"}.
#'
#' @return A numeric matrix (\code{T x K}) with column names:
#'   \itemize{
#'     \item If \code{include_membership = FALSE}: \code{c("beta_1","beta_2")}
#'     \item Else: \code{c("gamma","beta_1","beta_2")}
#'   }
#'   Unless \code{density = "model"} the matrix carries attributes
#'   \code{density_mode}, \code{density_weight_mode},
#'   \code{stage1_density_columns}, \code{stage1_density} and
#'   \code{model_density}; read them with \code{\link{weight_density_mode}}.
#'   A \code{"model"} W carries none, so it is identical to earlier versions.
#'
#' @examples
#' taus <- seq(0.10, 0.90, by = 0.05)
#' y <- rnorm(2000, 2, 1)
#' q_tau <- as.numeric(quantile(y, taus, type = 8))
#' W <- make_weights_normal_mixture(
#'   taus, q_tau,
#'   p1 = 0.5, mu1 = 1.2, sd1 = 0.45,
#'   mu2 = 3.0, sd2 = 0.7,
#'   include_membership = TRUE,
#'   stage1_density = stage1_kde_density(y, q_tau)   # default density = "hybrid"
#' )
#' dim(W); colnames(W)
#'
#' @export
make_weights_normal_mixture <- function(
  taus,
  q_tau,
  p1,
  mu1, sd1,
  mu2, sd2,
  include_membership = FALSE,
  tiny = 1e-12,
  density = NULL,
  stage1_density = NULL,
  identity_tol = 1e-8
) {
  stopifnot(length(taus) == length(q_tau))
  density <- .resolve_density_mode(density, has_membership = isTRUE(include_membership))
  Tt <- length(taus)

  F1 <- pnorm(q_tau, mean = mu1, sd = sd1)
  f1 <- dnorm(q_tau, mean = mu1, sd = sd1)
  F2 <- pnorm(q_tau, mean = mu2, sd = sd2)
  f2 <- dnorm(q_tau, mean = mu2, sd = sd2)

  f  <- pmax(p1 * f1 + (1 - p1) * f2, tiny)
  # Density for the location columns (f_loc) and for the membership column
  # (f_gam). Under the default both are the model pdf, so the arithmetic is
  # exactly that of earlier versions.
  f_loc <- f; f_gam <- f; s1 <- NULL
  if (density != "model") {
    s1 <- .resolve_stage1_density(stage1_density, taus, q_tau, density)
    fh <- pmax(s1$fhat, tiny)
    f_gam <- fh
    if (density == "stage1") f_loc <- fh
  }

  W1 <- p1 * f1 / f_loc
  W2 <- (1 - p1) * f2 / f_loc

  if (isTRUE(include_membership)) {
    Wp <- -(p1 * (1 - p1) / f_gam) * (F1 - F2)
    W  <- cbind(Wp, W1, W2)
    colnames(W) <- c("gamma", "beta_1", "beta_2")
  } else {
    W  <- cbind(W1, W2)
    colnames(W) <- c("beta_1", "beta_2")
  }
  rownames(W) <- paste0("tau", taus)
  storage.mode(W) <- "double"
  if (density != "model") {
    .check_weight_identities(W, list(c("beta_1", "beta_2")), density, identity_tol)
    cols <- if (density == "stage1") colnames(W) else
      intersect("gamma", colnames(W))
    W <- .tag_weight_density(W, density, s1$fhat, cols, f, s1$source)
  }
  W
}


#' Build tau-by-2 weights for Normal mean/variance perturbations
#'
#' Constructs a weight matrix W that maps RIF tau-slopes b(tau) into
#' effects on the mean and variance of a Normal baseline phenotype.
#'
#' Baseline: Y ~ N(mu, sigma^2).
#' For tau-quantiles \eqn{q_\tau = \mu + \sigma z_\tau}, with
#' \eqn{z_\tau = \Phi^{-1}(\tau)}:
#' \deqn{ W_mu(tau) = 1, \quad W_sigma2(tau) = z_tau / (2 * sigma). }
#'
#' @param taus Numeric vector of quantile levels (length T).
#' @param q_tau Numeric vector of baseline quantiles at taus (length T, type=8
#'   recommended). Only its length is checked: the Normal weights depend on
#'   \code{taus} and \code{sd} alone. Kept for a common builder interface.
#' @param mu Baseline mean of Y. Unused, for the same reason; kept for the
#'   interface.
#' @param sd Baseline standard deviation of Y.
#'
#' @return A T x 2 numeric matrix W with columns:
#'   \itemize{
#'     \item \code{"beta_mu"} effect of SNP on the mean
#'     \item \code{"beta_sigma2"} effect of SNP on the variance
#'   }
#' @examples
#' taus <- seq(0.1, 0.9, by = 0.2)
#' y <- rnorm(2000, mean = 2, sd = 1.5)
#' q_tau <- as.numeric(quantile(y, taus, type = 8))
#' W <- make_weights_vqtl(taus, q_tau, mu = mean(y), sd = sd(y))
#' W
#' @export
make_weights_vqtl <- function(
  taus,
  q_tau,
  mu,
  sd
) {
  stopifnot(length(taus) == length(q_tau))
  if (sd <= 0) stop("sd must be positive.")

  z_tau <- qnorm(taus)  # baseline standard normal quantiles

  W_mu     <- rep(1, length(taus))
  W_sigma2 <- z_tau / (2 * sd)

  W <- cbind(W_mu, W_sigma2)
  colnames(W) <- c("beta_mu", "beta_sigma2")
  rownames(W) <- paste0("tau", taus)
  storage.mode(W) <- "double"
  W
}

#' Build tau-by-K weights for a two-component Normal mixture with mean, variance, and membership effects
#'
#' Constructs a weight matrix \code{W} that maps RIF tau-slopes \eqn{b(tau)}
#' into SNP effects on the parameters of a two-component Normal mixture:
#' \deqn{Y ~ p1 * N(mu1, sd1^2) + (1 - p1) * N(mu2, sd2^2).}
#'
#' Supported SNP effect targets:
#' \itemize{
#'   \item \code{gamma}: log-odds of membership in component 1 (class proportion).
#'   \item \code{beta_mu1}, \code{beta_mu2}: effects on the component means.
#'   \item \code{beta_sigma1}, \code{beta_sigma2}: effects on the component standard deviations.
#' }
#'
#' @details
#' For baseline quantiles \eqn{q_tau} and mixture pdf
#' \deqn{f(q) = p1 f1(q) + (1 - p1) f2(q),}
#' where \eqn{f1, f2} are the component pdfs and \eqn{F1, F2} their CDFs, the weights are:
#'
#' - Membership:
#'   \deqn{W_gamma(tau) = - (p1 * (1 - p1) / f(q_tau)) * (F1(q_tau) - F2(q_tau))}
#'
#' - Means:
#'   \deqn{W_mu1(tau) = (p1 * f1(q_tau)) / f(q_tau), \quad
#'         W_mu2(tau) = ((1 - p1) * f2(q_tau)) / f(q_tau)}
#'
#' - Standard deviations:
#'   \deqn{W_sigma1(tau) = (p1 * f1(q_tau) * (q_tau - mu1) / sd1) / f(q_tau),}
#'   \deqn{W_sigma2(tau) = ((1 - p1) * f2(q_tau) * (q_tau - mu2) / sd2) / f(q_tau).}
#'
#' These arise from the derivative of the mixture CDF with respect to each parameter,
#' mapped into the quantile effect scale via the implicit function theorem.
#'
#' @param taus Numeric vector of quantile levels (length T).
#' @param q_tau Numeric vector of baseline quantiles of Y at \code{taus} (length T, \code{type=8} recommended).
#' @param p1 Proportion of component 1 (in (0,1)); component 2 proportion = 1 - p1.
#' @param mu1,sd1 Mean and sd of component 1.
#' @param mu2,sd2 Mean and sd of component 2.
#' @param tiny Small positive floor to stabilize divisions (mixture pdf clamped to at least \code{tiny}).
#' @param density,stage1_density,identity_tol As in
#'   \code{\link{make_weights_normal_mixture}}. This W always has a membership
#'   column, so the default is \code{"hybrid"} and \code{stage1_density} is
#'   required unless \code{density = "model"}. Under \code{"hybrid"} only the
#'   \code{gamma} column is divided by the Stage 1 density; the mean and sd
#'   columns (component share times a location or scale score) keep the model
#'   pdf, so \code{beta_mu1 + beta_mu2 = 1} still holds.
#'
#' @return A T x 5 numeric matrix \code{W} with columns:
#'   \code{c("gamma","beta_mu1","beta_mu2","beta_sigma1","beta_sigma2")}.
#'   A non-"model" \code{density} adds the attributes described in
#'   \code{\link{make_weights_normal_mixture}}.
#'
#' @examples
#' taus <- seq(0.1, 0.9, 0.1)
#' set.seed(1)
#' y <- c(rnorm(500, 1.2, 0.6), rnorm(500, 3.0, 0.9))
#' q_tau <- as.numeric(quantile(y, taus, type = 8))
#' W <- make_weights_normal_mixture_vqtl(
#'   taus, q_tau,
#'   p1 = 0.5, mu1 = 1.2, sd1 = 0.6,
#'   mu2 = 3.0, sd2 = 0.9,
#'   stage1_density = stage1_kde_density(y, q_tau)
#' )
#' head(W)
#' @export
make_weights_normal_mixture_vqtl <- function(
  taus,
  q_tau,
  p1,
  mu1, sd1,
  mu2, sd2,
  tiny = 1e-12,
  density = NULL,
  stage1_density = NULL,
  identity_tol = 1e-8
) {
  stopifnot(length(taus) == length(q_tau))
  density <- .resolve_density_mode(density, has_membership = TRUE)

  # Component pdfs and cdfs
  f1 <- dnorm(q_tau, mean = mu1, sd = sd1)
  F1 <- pnorm(q_tau, mean = mu1, sd = sd1)
  f2 <- dnorm(q_tau, mean = mu2, sd = sd2)
  F2 <- pnorm(q_tau, mean = mu2, sd = sd2)

  # Mixture density (guarded)
  f <- pmax(p1 * f1 + (1 - p1) * f2, tiny)
  f_loc <- f; f_gam <- f; s1 <- NULL
  if (density != "model") {
    s1 <- .resolve_stage1_density(stage1_density, taus, q_tau, density)
    fh <- pmax(s1$fhat, tiny)
    f_gam <- fh
    if (density == "stage1") f_loc <- fh
  }

  # Weights
  W_gamma    <- -(p1 * (1 - p1) / f_gam) * (F1 - F2)
  W_mu1      <- (p1 * f1) / f_loc
  W_mu2      <- ((1 - p1) * f2) / f_loc
  W_sigma1   <- (p1 * f1 * (q_tau - mu1) / sd1) / f_loc
  W_sigma2   <- ((1 - p1) * f2 * (q_tau - mu2) / sd2) / f_loc

  W <- cbind(W_gamma, W_mu1, W_mu2, W_sigma1, W_sigma2)
  colnames(W) <- c("gamma","beta_mu1","beta_mu2","beta_sigma1","beta_sigma2")
  rownames(W) <- paste0("tau", taus)
  storage.mode(W) <- "double"
  if (density != "model") {
    .check_weight_identities(W, list(c("beta_mu1", "beta_mu2")), density, identity_tol)
    cols <- if (density == "stage1") colnames(W) else "gamma"
    W <- .tag_weight_density(W, density, s1$fhat, cols, f, s1$source)
  }
  W
}


#' Construct a finite-difference gradient function
#'
#' Builds a function that computes the numerical gradient of a distribution's CDF
#' with respect to a single parameter using central differences.
#'
#' @param param Name of the parameter (string).
#' @param dist_cdf Function (y, params) -> numeric. CDF of the distribution.
#' @param eps Finite-difference step size (default 1e-6).
#'
#' @return A function with signature (y, params) that returns
#'   dF(y; params) / d(param) for the given parameter.
#'
#' @examples
#' dist_cdf <- function(y, params) pgamma(y, shape = params$shape, scale = params$scale)
#' grad_shape <- make_fd_grad("shape", dist_cdf)
#' grad_shape(1.5, list(shape = 2, scale = 1))
#'
#' @export
make_fd_grad <- function(param, dist_cdf, eps = 1e-6) {
  force(param); force(dist_cdf); force(eps)
  function(y, params) {
    params_up   <- params; params_up[[param]]   <- params_up[[param]] + eps
    params_down <- params; params_down[[param]] <- params_down[[param]] - eps
    (dist_cdf(y, params_up) - dist_cdf(y, params_down)) / (2 * eps)
  }
}

#' Generic weight builder for parametric GWAS
#'
#' Constructs a tau-by-K weight matrix \eqn{W} for an arbitrary parametric family,
#' mapping RIF tau-slopes \eqn{b(tau)} into parameter effects.
#'
#' @param taus Numeric vector of quantile levels (length T).
#' @param q_tau Numeric vector of baseline quantiles at \code{taus} (length T).
#' @param dist_cdf Function (y, params) -> numeric. CDF of the distribution. Not
#'   called by the builder itself (the gradients carry the CDF); kept so the
#'   family is fully specified and for \code{make_fd_grad()}-style gradients.
#' @param dist_pdf Function (y, params) -> numeric. PDF of the distribution.
#' @param params Named list of baseline parameters (e.g., \code{list(mu = 0, sd = 1)}).
#' @param grad_funcs Named list of gradient functions. Each element is a function
#'   (y, params) -> numeric giving derivative of F(y; params) wrt that parameter.
#'   Can be analytic, or constructed with \code{make_fd_grad()}.
#' @param tiny Small positive floor to stabilize divisions (\eqn{f(q_tau)} clamped).
#' @param density Which density divides each column (see
#'   \code{docs/MIXTURE_WEIGHT_DENSITY.md}). The default (\code{NULL}) is
#'   \code{"hybrid"} when a gradient is named in \code{stage1_columns} (a
#'   membership column) and \code{"model"} otherwise. \code{"model"} divides
#'   every column by \code{dist_pdf}; \code{"hybrid"}
#'   divides the columns named in \code{stage1_columns} (the membership
#'   derivative of a mixture) by the Stage 1 density and every other column by
#'   \code{dist_pdf}; \code{"stage1"} divides every column by the Stage 1
#'   density and warns (sensitivity only).
#' @param stage1_density Stage 1 density at \code{taus}, required unless
#'   \code{density = "model"}; see \code{\link{make_weights_normal_mixture}}.
#' @param stage1_columns Columns divided by the Stage 1 density under
#'   \code{"hybrid"} (default \code{"gamma"}). At least one must be present,
#'   unless \code{character(0)} is passed to state that the family has no
#'   membership column (then \code{"hybrid"} equals \code{"model"}).
#' @param unit_sum Optional identity check: a character vector, or a list of
#'   them, naming groups of columns whose row sums must equal 1 (for example
#'   \code{c("beta_left", "beta_right")} for component locations, or
#'   \code{"additive"} for a shift column, which must be identically 1). A
#'   violation beyond \code{identity_tol} stops under \code{"model"} and
#'   \code{"hybrid"}; under \code{"stage1"} the deviations are reported in a
#'   warning. Not checked when \code{NULL} (default).
#' @param identity_tol Tolerance for \code{unit_sum}; the default \code{1e-6}
#'   allows for finite-difference gradients (\code{\link{make_fd_grad}}).
#'
#' @return A T x K numeric matrix of weights with columns named by parameter.
#'   A non-"model" \code{density} adds the attributes described in
#'   \code{\link{make_weights_normal_mixture}}.
#'
#' @details
#' The generic approach uses the identity:
#' \deqn{W_j(tau) = - (dF/dtheta_j)(q_tau; params) / f(q_tau; params),}
#' where \eqn{F} is the CDF and \eqn{f} is the PDF.
#'
#' @examples
#' taus <- seq(0.1, 0.9, 0.1)
#' y <- rnorm(5000, 2, 1)
#' q_tau <- as.numeric(quantile(y, taus, type = 8))
#'
#' dist_cdf <- function(y, params) pnorm(y, mean = params$mu, sd = params$sd)
#' dist_pdf <- function(y, params) dnorm(y, mean = params$mu, sd = params$sd)
#' grad_mu  <- function(y, params) -dnorm((y - params$mu)/params$sd)/params$sd
#' grad_sd  <- function(y, params) { z <- (y - params$mu)/params$sd; -(z*dnorm(z))/params$sd }
#'
#' W <- make_weights_generic(
#'   taus, q_tau, dist_cdf, dist_pdf,
#'   params = list(mu = 2, sd = 1),
#'   grad_funcs = list(beta_mu = grad_mu, beta_sd = grad_sd)
#' )
#'
#' @export
make_weights_generic <- function(
  taus, q_tau, dist_cdf, dist_pdf,
  params, grad_funcs,
  tiny = 1e-12,
  density = NULL,
  stage1_density = NULL,
  stage1_columns = "gamma",
  unit_sum = NULL,
  identity_tol = 1e-6
) {
  stopifnot(length(taus) == length(q_tau))
  density <- .resolve_density_mode(
    density, has_membership = any(names(grad_funcs) %in% stage1_columns))
  Tt <- length(taus)
  f  <- pmax(dist_pdf(q_tau, params), tiny)

  s1 <- NULL
  use_s1 <- rep(FALSE, length(grad_funcs))
  if (density != "model") {
    s1 <- .resolve_stage1_density(stage1_density, taus, q_tau, density)
    fh <- pmax(s1$fhat, tiny)
    if (density == "stage1") {
      use_s1[] <- TRUE
    } else {
      use_s1 <- names(grad_funcs) %in% stage1_columns
      # character(0) states explicitly that the family has no membership
      # column, in which case "hybrid" equals "model" (still recorded).
      if (!any(use_s1) && length(stage1_columns)) {
        stop("density = \"hybrid\" needs at least one membership column: none of ",
             "`stage1_columns` (", paste(stage1_columns, collapse = ", "),
             ") is among the gradient names (", paste(names(grad_funcs), collapse = ", "),
             ").", call. = FALSE)
      }
    }
  }

  W <- matrix(NA_real_, nrow = Tt, ncol = length(grad_funcs))
  colnames(W) <- names(grad_funcs)
  rownames(W) <- paste0("tau", taus)

  for (j in seq_along(grad_funcs)) {
    gfun <- grad_funcs[[j]]
    dF   <- vapply(q_tau, function(q) gfun(q, params), numeric(1))
    W[, j] <- if (use_s1[j]) -dF / fh else -dF / f
  }
  storage.mode(W) <- "double"
  if (density != "model" || !is.null(unit_sum)) {
    .check_weight_identities(W, unit_sum, density, identity_tol)
  }
  if (density != "model") {
    W <- .tag_weight_density(W, density, s1$fhat, colnames(W)[use_s1], f, s1$source)
  }
  W
}
