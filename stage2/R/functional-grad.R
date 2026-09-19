# ===========================================================================
# Functional gradients of a fitted baseline distribution.
#
# A weight matrix W maps Stage 2 parameters theta onto quantile effects:
#   W_j(tau) = d q_tau / d theta_j.
# It does NOT provide a route to a functional of the distribution (a mean, a
# variance, P(Y > c)).  Averaging the tau grid -- colMeans(W) %*% theta or an
# unweighted projection -- is only valid when the quantile map is smooth and
# well conditioned.  When two mixture components barely overlap the gamma
# column develops a near-singular spike and the grid average is materially
# wrong (see docs/FUNCTIONAL_GRADIENTS.md for the worked example).
#
# The information needed for the exact answer is already in the fitted object:
# the mixture stores its components (weight/xi/omega/alpha) and the lognormal
# families store params0 (shift/meanlog/sdlog).  These helpers route that
# information into a finite-difference gradient of a user-supplied functional
# f of the baseline parameters, taken on the parametric family itself.
# ===========================================================================

## The gradient only needs W for the Stage 2 parameter names and their order;
## the W values never enter it.  Accept a bare matrix or a weight-object list.
.functional_weight_matrix <- function(object) {
  W <- if (is.matrix(object)) object else object$W
  if (is.null(W) || !is.matrix(W)) {
    stop("`object` must be a weight matrix or a list with a matrix `W`.",
         call. = FALSE)
  }
  if (is.null(colnames(W)) || any(!nzchar(colnames(W)))) {
    stop("The weight matrix must have non-empty column names: they name the ",
         "Stage 2 parameters and fix the order of the gradient.", call. = FALSE)
  }
  W
}

.functional_align_theta <- function(theta, nms) {
  if (is.null(theta)) return(stats::setNames(rep(0, length(nms)), nms))
  if (!is.null(names(theta))) {
    missing <- setdiff(nms, names(theta))
    if (length(missing)) {
      stop("`theta` is missing Stage 2 parameter(s): ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    out <- as.numeric(theta[nms])
    names(out) <- nms
    return(out)
  }
  if (length(theta) != length(nms)) {
    stop("Unnamed `theta` has length ", length(theta), "; expected ", length(nms),
         " (the number of W columns).", call. = FALSE)
  }
  out <- as.numeric(theta)
  names(out) <- nms
  out
}

## Build the function theta -> baseline distribution parameters for a weight
## object.  The mapping is the parameterisation the package supplies; the
## semantics of every supported column are tabulated in
## docs/FUNCTIONAL_GRADIENTS.md.
.functional_param_map <- function(object) {
  if (is.list(object) && is.function(object$param_map)) return(object$param_map)
  W <- .functional_weight_matrix(object)
  nms <- colnames(W)
  components <- if (is.list(object)) object$components else NULL
  params0 <- if (is.list(object)) object$params0 else NULL

  ## --- two-component (skew-)normal mixture --------------------------------
  ## columns: gamma, beta_left, beta_right.  `gamma` is the log-odds of
  ## membership in the LOWER (first) component, so p_lower = plogis(gamma).
  if (!is.null(components)) {
    comp <- as.data.frame(components)
    need <- c("component", "weight", "xi", "omega", "alpha")
    if (!all(need %in% names(comp))) {
      stop("`components` must have columns: ", paste(need, collapse = ", "),
           call. = FALSE)
    }
    comp <- comp[order(match(comp$component, c("lower", "upper"))), , drop = FALSE]
    if (nrow(comp) != 2L) {
      stop("`components` must have exactly two rows (lower, upper).", call. = FALSE)
    }
    base_w <- as.numeric(comp$weight)
    base_xi <- as.numeric(comp$xi)
    base_omega <- as.numeric(comp$omega)
    base_alpha <- as.numeric(comp$alpha)
    return(function(theta) {
      w_lower <- base_w[1L]
      if ("gamma" %in% nms) {
        w_lower <- stats::plogis(stats::qlogis(base_w[1L]) + theta[["gamma"]])
      }
      xi <- base_xi
      if ("beta_left" %in% nms)  xi[1L] <- base_xi[1L] + theta[["beta_left"]]
      if ("beta_right" %in% nms) xi[2L] <- base_xi[2L] + theta[["beta_right"]]
      list(weight = c(w_lower, 1 - w_lower), xi = xi,
           omega = base_omega, alpha = base_alpha)
    })
  }

  ## --- shifted lognormal, and both log-scale vQTL parameterisations -------
  ## columns: additive (shift), multiplicative (meanlog);
  ##          vqtl: b1_mu (meanlog), g1_sigma2 (sdlog^2);
  ##          log-scale vQTL: meanlog, sigma2 (sdlog^2).
  if (!is.null(params0) && all(c("shift", "meanlog", "sdlog") %in% names(params0))) {
    p0 <- params0
    return(function(theta) {
      p <- p0
      if ("additive" %in% nms)       p$shift   <- p0$shift   + theta[["additive"]]
      if ("multiplicative" %in% nms) p$meanlog <- p0$meanlog + theta[["multiplicative"]]
      if ("meanlog" %in% nms)        p$meanlog <- p0$meanlog + theta[["meanlog"]]
      if ("b1_mu" %in% nms)          p$meanlog <- p0$meanlog + theta[["b1_mu"]]
      if ("sigma2" %in% nms)         p$sdlog   <- sqrt(p0$sdlog^2 + theta[["sigma2"]])
      if ("g1_sigma2" %in% nms)      p$sdlog   <- sqrt(p0$sdlog^2 + theta[["g1_sigma2"]])
      p
    })
  }

  ## --- Normal mean / variance vQTL ---------------------------------------
  ## columns: beta_mu (mean), beta_sigma2 (variance).  Baseline params must be
  ## supplied with the object, since make_weights_vqtl() returns a bare matrix.
  if (all(c("beta_mu", "beta_sigma2") %in% nms)) {
    if (is.null(params0) || is.null(params0$mu) ||
        (is.null(params0$sigma2) && is.null(params0$sd))) {
      stop("A Normal mean/variance weight object needs ",
           "`params0 = list(mu = , sigma2 = )` (or sd).", call. = FALSE)
    }
    mu0 <- as.numeric(params0$mu)[1L]
    s20 <- if (!is.null(params0$sigma2)) as.numeric(params0$sigma2)[1L] else
      as.numeric(params0$sd)[1L]^2
    return(function(theta) {
      mu <- mu0 + theta[["beta_mu"]]
      s2 <- s20 + theta[["beta_sigma2"]]
      list(mu = mu, sd = sqrt(s2), sigma2 = s2)
    })
  }

  stop("Unsupported weight object: no `components`, no ",
       "`shift`/`meanlog`/`sdlog` params0, and no Normal mean/variance ",
       "columns. Supply `object$param_map` or see docs/FUNCTIONAL_GRADIENTS.md.",
       call. = FALSE)
}

#' Gradient of a functional of the fitted baseline distribution
#'
#' Computes the gradient of a smooth functional \eqn{f} of a fitted baseline
#' distribution with respect to the Stage 2 parameters \eqn{\theta}, the
#' columns of the weight matrix \eqn{W}.  The gradient is taken by central
#' finite differences on the parametric family itself, never by averaging the
#' \eqn{\tau} grid; see \code{docs/FUNCTIONAL_GRADIENTS.md} for why that
#' distinction matters and for the parameter semantics of each weight family.
#'
#' @param object A weight object: a list with a matrix \code{W} and the
#'   baseline parameters, or a bare weight matrix when the family can be
#'   inferred.  The supported families and the baseline-parameter list passed
#'   to \code{f} are:
#'   \describe{
#'     \item{Skew-normal mixture}{an object with \code{components} (columns
#'       \code{component}, \code{weight}, \code{xi}, \code{omega},
#'       \code{alpha}) and \eqn{W} columns \code{gamma}, \code{beta_left},
#'       \code{beta_right}; \code{f} receives \code{list(weight, xi, omega,
#'       alpha)} as length-two vectors.}
#'     \item{Shifted lognormal}{an object with \code{params0 =
#'       list(shift, meanlog, sdlog)} and columns \code{additive} and/or
#'       \code{multiplicative}; \code{f} receives the \code{params0} list.}
#'     \item{Log-scale vQTL}{an object with the same \code{params0} and columns
#'       \code{meanlog}/\code{sigma2} or \code{b1_mu}/\code{g1_sigma2};
#'       \code{f} receives the \code{params0} list.}
#'     \item{Normal mean/variance vQTL}{\code{params0 = list(mu, sigma2)} and
#'       columns \code{beta_mu}, \code{beta_sigma2}; \code{f} receives
#'       \code{list(mu, sd, sigma2)}.}
#'   }
#'   An object may instead carry \code{param_map}, a function
#'   \code{function(theta) -> params}, which overrides the built-in families.
#' @param f A function of the baseline parameters returning a single number
#'   (for example the distribution mean, variance, or a tail probability).
#' @param theta Optional named vector of Stage 2 parameters at which to take
#'   the gradient.  Defaults to all zeros, i.e. the fitted baseline.  Names
#'   must match \code{colnames(object$W)}; an unnamed vector is used in the W
#'   column order.
#' @param eps Central finite-difference step (default \code{1e-6}).
#'
#' @return A named numeric vector, names \code{colnames(object$W)}, giving
#'   \eqn{d f / d \theta_j} at \code{theta} in the same units and order as the
#'   Stage 2 parameters.
#'
#' @seealso \code{\link{functional_estimate}}, \code{\link{make_fd_grad}}
#'
#' @examples
#' comp <- data.frame(
#'   component = c("lower", "upper"),
#'   weight = c(0.5, 0.5), xi = c(0.3, 2.7),
#'   omega = c(0.55, 0.38), alpha = c(-0.93, -1.03)
#' )
#' obj <- list(W = matrix(0, 3, 3,
#'                        dimnames = list(NULL, c("gamma", "beta_left", "beta_right"))),
#'             components = comp)
#' f_mean <- function(p) {
#'   d <- p$alpha / sqrt(1 + p$alpha^2)
#'   E <- p$xi + p$omega * d * sqrt(2 / pi)
#'   p$weight[1] * E[1] + p$weight[2] * E[2]
#' }
#' functional_grad(obj, f_mean)
#'
#' @export
functional_grad <- function(object, f, theta = NULL, eps = 1e-6) {
  if (!is.function(f)) {
    stop("`f` must be a function of the baseline parameters.", call. = FALSE)
  }
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("`eps` must be a single positive finite number.", call. = FALSE)
  }
  W <- .functional_weight_matrix(object)
  nms <- colnames(W)
  theta <- .functional_align_theta(theta, nms)
  map <- .functional_param_map(object)

  g <- vapply(nms, function(nm) {
    up <- theta; up[[nm]] <- up[[nm]] + eps
    dn <- theta; dn[[nm]] <- dn[[nm]] - eps
    fu <- f(map(up)); fd <- f(map(dn))
    if (length(fu) != 1L || length(fd) != 1L) {
      stop("`f` must return a single number.", call. = FALSE)
    }
    (as.numeric(fu) - as.numeric(fd)) / (2 * eps)
  }, numeric(1))
  names(g) <- nms
  g
}

#' Functional estimate and delta-method standard error
#'
#' Evaluates a functional \code{f} at a Stage 2 parameter vector \code{theta}
#' and combines the gradient from \code{\link{functional_grad}} with the
#' parameter covariance \code{cov_theta} to give a delta-method standard error.
#'
#' @param object A weight object, as for \code{\link{functional_grad}}.
#' @param f A function of the baseline parameters returning a single number.
#' @param theta Named (or positional) vector of Stage 2 parameter estimates.
#' @param cov_theta Optional covariance matrix of \code{theta}, with rows and
#'   columns in the order of \code{colnames(object$W)}.  When \code{NULL} the
#'   returned \code{se} is \code{NA}.
#' @param eps Central finite-difference step (default \code{1e-6}).
#'
#' @return A list with \code{estimate} (the functional at \code{theta}),
#'   \code{se} (the delta-method standard error, or \code{NA} when
#'   \code{cov_theta} is \code{NULL}) and \code{grad} (the gradient at
#'   \code{theta}, see \code{\link{functional_grad}}).
#'
#' @seealso \code{\link{functional_grad}}
#'
#' @examples
#' comp <- data.frame(
#'   component = c("lower", "upper"),
#'   weight = c(0.5, 0.5), xi = c(0.3, 2.7),
#'   omega = c(0.55, 0.38), alpha = c(-0.93, -1.03)
#' )
#' obj <- list(W = matrix(0, 3, 3,
#'                        dimnames = list(NULL, c("gamma", "beta_left", "beta_right"))),
#'             components = comp)
#' f_mean <- function(p) {
#'   d <- p$alpha / sqrt(1 + p$alpha^2)
#'   E <- p$xi + p$omega * d * sqrt(2 / pi)
#'   p$weight[1] * E[1] + p$weight[2] * E[2]
#' }
#' functional_estimate(obj, f_mean,
#'                     theta = c(gamma = 0.1, beta_left = 0, beta_right = 0),
#'                     cov_theta = diag(c(0.01, 0.02, 0.03)))
#'
#' @export
functional_estimate <- function(object, f, theta, cov_theta = NULL, eps = 1e-6) {
  W <- .functional_weight_matrix(object)
  nms <- colnames(W)
  theta <- .functional_align_theta(theta, nms)
  map <- .functional_param_map(object)
  estimate <- as.numeric(f(map(theta)))
  if (length(estimate) != 1L) {
    stop("`f` must return a single number.", call. = FALSE)
  }
  grad <- functional_grad(object, f, theta, eps)
  se <- NA_real_
  if (!is.null(cov_theta)) {
    cov_theta <- as.matrix(cov_theta)
    if (!identical(dim(cov_theta), c(length(nms), length(nms)))) {
      stop("`cov_theta` must be ", length(nms), " x ", length(nms), ".",
           call. = FALSE)
    }
    v <- as.numeric(crossprod(grad, cov_theta %*% grad))
    if (!is.finite(v) || v < 0) {
      stop("The delta-method variance is not finite and non-negative.",
           call. = FALSE)
    }
    se <- sqrt(v)
  }
  list(estimate = estimate, se = se, grad = grad)
}