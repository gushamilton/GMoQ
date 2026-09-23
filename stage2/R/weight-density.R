# ===========================================================================
# Which density divides each column of a mixture weight matrix.
#
# A weight column is W_j(tau) = -(dF/dtheta_j)(q_tau) / f(q_tau).  Stage 1
# returns b(tau) = -(dF/dtheta)(q_tau) / fhat(q_tau), where fhat is the
# Stage 1 kernel density (Silverman Gaussian KDE at the empirical quantile), so
# the coefficient carries fhat's scaling.  Three choices are offered:
#
#   "model"   every column divided by the fitted parametric density (the
#             historical behaviour and the default).  Component-share columns
#             are P(k | q) and sum to 1; the additive column of a shifted
#             family is exactly 1.
#   "hybrid"  location (component-share) columns as in "model"; the membership
#             column(s) (-dF/dgamma, a CDF-level quantity) divided by the
#             Stage 1 density fhat -- the density the RIF actually used.  The
#             identities of "model" still hold.
#   "stage1"  every column divided by fhat.  Breaks both identities; retained
#             as a sensitivity analysis only, and warns.
#
# See docs/MIXTURE_WEIGHT_DENSITY.md.  The chosen mode is recorded as
# attributes on the returned matrix (only when it is not the default, so the
# default output is unchanged); weight_density_mode() reads it back and the
# Stage 2 file wrappers can refuse a mismatch via `require_density_mode`.
# ===========================================================================

#' @keywords internal
WEIGHT_DENSITY_MODES <- c("model", "hybrid", "stage1")

#' Provenance labels for the weight density modes
#'
#' The strings the paper pipeline records as \code{density_weight_mode}.
#' @keywords internal
WEIGHT_DENSITY_LABELS <- c(
  model  = "parametric_model_density",
  hybrid = "hybrid_stage1_membership_density",
  stage1 = "stage1_rif_density"
)

## Map a mode or a provenance label onto the canonical mode, or NA.
.canonical_density_mode <- function(x) {
  x <- as.character(x)
  if (length(x) != 1L || is.na(x)) return(NA_character_)
  if (x %in% WEIGHT_DENSITY_MODES) return(x)
  hit <- names(WEIGHT_DENSITY_LABELS)[WEIGHT_DENSITY_LABELS == x]
  if (length(hit) == 1L) hit else NA_character_
}

#' Stage 1 kernel density at the baseline quantiles
#'
#' Reproduces the density Stage 1 uses to build the RIF
#' (\code{fungwas_stage1.core.compute_rif}): a Gaussian kernel with the
#' Silverman bandwidth \eqn{h = 1.06\,\mathrm{sd}(Y)\,N^{-1/5}}, where
#' \eqn{\mathrm{sd}} is the population (divide-by-\eqn{N}) standard
#' deviation, evaluated at each \code{q_tau} and floored at
#' \code{density_floor}. Pass the same phenotype vector Stage 1 used (after any
#' LOCO subtraction), and the Stage 1 \code{q_tau} (median-unbiased, R
#' \code{type = 8}).
#'
#' @param y Numeric phenotype vector used by Stage 1.
#' @param q_tau Numeric vector of baseline quantiles.
#' @param method \code{"silverman_gaussian"} (the Stage 1 default) or
#'   \code{"fourth_order_gaussian"} (kernel \eqn{(1.5 - 0.5u^2)\phi(u)}, same
#'   bandwidth).
#' @param density_floor Positive floor, as in Stage 1 (default \code{1e-8}).
#' @return Numeric vector of densities, one per \code{q_tau}.
#' @examples
#' set.seed(1)
#' y <- c(rnorm(500, 0, 1), rnorm(500, 5, 1))
#' taus <- seq(0.1, 0.9, 0.1)
#' q_tau <- as.numeric(quantile(y, taus, type = 8))
#' stage1_kde_density(y, q_tau)
#' @export
stage1_kde_density <- function(y, q_tau,
                               method = c("silverman_gaussian", "fourth_order_gaussian"),
                               density_floor = 1e-8) {
  method <- match.arg(method)
  y <- as.numeric(y)
  if (length(y) < 2L || any(!is.finite(y))) {
    stop("`y` must be a finite numeric vector with at least two values.", call. = FALSE)
  }
  if (!is.finite(density_floor) || density_floor <= 0) {
    stop("`density_floor` must be finite and positive.", call. = FALSE)
  }
  n <- length(y)
  h <- 1.06 * sqrt(mean((y - mean(y))^2)) * n^(-0.2)
  if (!is.finite(h) || h <= 0) stop("`y` has zero spread; no bandwidth.", call. = FALSE)
  norm <- sqrt(2 * pi)
  out <- vapply(as.numeric(q_tau), function(q) {
    u <- (y - q) / h
    k <- exp(-0.5 * u^2)
    if (method == "fourth_order_gaussian") k <- k * (1.5 - 0.5 * u^2)
    mean(k) / (h * norm)
  }, numeric(1))
  pmax(out, density_floor)
}

#' Read the Stage 1 density at the taus from a JSON sidecar
#'
#' Reads a Stage 1 RIF nuisance JSON (as written by the paper pipeline next to
#' each Stage 1 scan, \code{rif_nuisance.json}): an object with numeric arrays
#' \code{fhat} (required), \code{taus} and \code{q_tau} (optional but
#' recommended; used to check the grid), and optionally scalar \code{n}.
#' Uses \pkg{jsonlite} when installed, otherwise the package's base-R JSON
#' readers.
#'
#' @param path Path to the JSON file.
#' @return A list with \code{fhat}, \code{taus}, \code{q_tau} (either may be
#'   \code{NULL}), \code{n} and \code{source}.
#' @export
read_stage1_density <- function(path) {
  if (!is.character(path) || length(path) != 1L || !file.exists(path)) {
    stop("Stage 1 density file not found: ", paste(path, collapse = ", "), call. = FALSE)
  }
  num <- function(v) if (is.null(v) || !length(v)) NULL else as.numeric(unlist(v))
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    m <- jsonlite::fromJSON(path, simplifyVector = TRUE)
    fhat <- num(m$fhat); taus <- num(m$taus); q_tau <- num(m$q_tau); n <- num(m$n)
  } else {
    txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
    fhat <- num(.json_number_array(txt, "fhat"))
    taus <- num(.json_number_array(txt, "taus"))
    q_tau <- num(.json_number_array(txt, "q_tau"))
    n <- .json_number(txt, "n"); if (is.na(n)) n <- NULL
  }
  if (is.null(fhat)) stop("Stage 1 density file ", path, " has no `fhat` array.", call. = FALSE)
  list(fhat = fhat, taus = taus, q_tau = q_tau, n = n,
       source = normalizePath(path, mustWork = TRUE))
}

## Resolve `stage1_density` (vector, list or path) to a checked fhat vector on
## exactly this tau / q_tau grid.
.resolve_stage1_density <- function(stage1_density, taus, q_tau, density) {
  if (is.null(stage1_density)) {
    stop("density = \"", density, "\" needs `stage1_density`: the Stage 1 density ",
         "fhat at q_tau (a numeric vector, a list with `fhat`, or the path of a ",
         "Stage 1 density JSON; see stage1_kde_density() and read_stage1_density()).",
         call. = FALSE)
  }
  src <- NULL
  if (is.character(stage1_density)) stage1_density <- read_stage1_density(stage1_density)
  if (is.list(stage1_density)) {
    obj <- stage1_density
    if (is.null(obj$fhat)) stop("`stage1_density` list has no `fhat`.", call. = FALSE)
    if (!is.null(obj$taus)) {
      if (length(obj$taus) != length(taus) ||
          max(abs(as.numeric(obj$taus) - taus)) > 5e-4) {
        stop("`stage1_density` was computed on a different tau grid.", call. = FALSE)
      }
    }
    if (!is.null(obj$q_tau)) {
      qs <- as.numeric(obj$q_tau)
      scale <- max(1, diff(range(q_tau)))
      if (length(qs) != length(q_tau) || max(abs(qs - q_tau)) > 1e-6 * scale) {
        stop("`stage1_density` was computed at different q_tau than the weights; ",
             "use the Stage 1 q_tau for both.", call. = FALSE)
      }
    }
    src <- obj$source
    fhat <- as.numeric(obj$fhat)
  } else {
    fhat <- as.numeric(stage1_density)
  }
  if (length(fhat) != length(taus)) {
    stop("`stage1_density` has length ", length(fhat), "; expected ", length(taus),
         " (one Stage 1 density per tau).", call. = FALSE)
  }
  if (any(!is.finite(fhat)) || any(fhat <= 0)) {
    stop("`stage1_density` must be finite and strictly positive.", call. = FALSE)
  }
  list(fhat = fhat, source = src)
}

## Check that each group of columns sums to 1 row-wise (a one-column group
## means the column is identically 1).
.check_weight_identities <- function(W, unit_sum, density, tol) {
  if (is.null(unit_sum) || !length(unit_sum)) {
    if (identical(density, "stage1")) {
      warning("density = \"stage1\" divides every column by the Stage 1 density: ",
              "component-share columns no longer sum to 1 and an additive column is ",
              "no longer 1. Use it as a sensitivity analysis only ",
              "(docs/MIXTURE_WEIGHT_DENSITY.md).", call. = FALSE)
    }
    return(invisible(NULL))
  }
  if (!is.list(unit_sum)) unit_sum <- list(unit_sum)
  devs <- vapply(unit_sum, function(g) {
    miss <- setdiff(g, colnames(W))
    if (length(miss)) stop("`unit_sum` names column(s) not in W: ",
                           paste(miss, collapse = ", "), call. = FALSE)
    max(abs(rowSums(W[, g, drop = FALSE]) - 1))
  }, numeric(1))
  labs <- vapply(unit_sum, paste, character(1), collapse = " + ")
  if (identical(density, "stage1")) {
    warning("density = \"stage1\" breaks the weight identities (max |sum - 1|: ",
            paste(sprintf("%s %.3g", labs, devs), collapse = "; "),
            "). Use it as a sensitivity analysis only (docs/MIXTURE_WEIGHT_DENSITY.md).",
            call. = FALSE)
  } else if (any(devs > tol)) {
    bad <- devs > tol
    stop("Weight identity violated under density = \"", density, "\": ",
         paste(sprintf("max |%s - 1| = %.3g", labs[bad], devs[bad]), collapse = "; "),
         " (tolerance ", tol, ").", call. = FALSE)
  }
  invisible(devs)
}

## Record the mode on a non-default weight matrix.
.tag_weight_density <- function(W, density, fhat, columns, f_model, source = NULL) {
  attr(W, "density_mode") <- density
  attr(W, "density_weight_mode") <- unname(WEIGHT_DENSITY_LABELS[[density]])
  attr(W, "stage1_density_columns") <- columns
  attr(W, "stage1_density") <- as.numeric(fhat)
  if (!is.null(f_model)) attr(W, "model_density") <- as.numeric(f_model)
  if (!is.null(source)) attr(W, "stage1_density_source") <- source
  W
}

#' Which density a weight matrix was built with
#'
#' Reads the density mode recorded by the mixture and generic weight builders
#' (\code{density = "model"}, \code{"hybrid"} or \code{"stage1"}). A matrix
#' without the attribute was built with the default parametric model density
#' (every builder before this option existed, and the default now), so it
#' reports \code{"model"}. A weight-object list is also accepted: its
#' \code{W} attribute is used first, then \code{provenance$density_weight_mode}.
#'
#' @param W A weight matrix, or a list with a matrix \code{W}.
#' @param label If \code{TRUE}, return the provenance label
#'   (\code{"parametric_model_density"}, \code{"hybrid_stage1_membership_density"},
#'   \code{"stage1_rif_density"}) instead of the mode.
#' @return A single string.
#' @examples
#' taus <- seq(0.1, 0.9, 0.1)
#' W <- make_weights_normal_mixture(taus, qnorm(taus), 0.5, -1, 1, 1, 1)
#' weight_density_mode(W)
#' @export
weight_density_mode <- function(W, label = FALSE) {
  mode <- NA_character_
  M <- if (is.list(W) && !is.null(W$W)) W$W else W
  a <- attr(M, "density_mode", exact = TRUE)
  if (!is.null(a)) mode <- .canonical_density_mode(a)
  if (is.na(mode) && is.list(W) && !is.null(W$provenance$density_weight_mode)) {
    mode <- .canonical_density_mode(W$provenance$density_weight_mode)
    if (is.na(mode)) {
      stop("Unrecognised density_weight_mode in provenance: ",
           W$provenance$density_weight_mode, call. = FALSE)
    }
  }
  if (is.na(mode)) {
    if (!is.null(a)) stop("Unrecognised density_mode attribute: ", a, call. = FALSE)
    mode <- "model"
  }
  if (isTRUE(label)) unname(WEIGHT_DENSITY_LABELS[[mode]]) else mode
}

## Stage 2 provenance gate used by the file wrappers.
.check_require_density_mode <- function(W, require_density_mode, what = "W") {
  if (is.null(require_density_mode)) return(invisible(NULL))
  want <- .canonical_density_mode(require_density_mode)
  if (is.na(want)) {
    stop("`require_density_mode` must be one of ",
         paste(c(WEIGHT_DENSITY_MODES, WEIGHT_DENSITY_LABELS), collapse = ", "),
         call. = FALSE)
  }
  got <- weight_density_mode(W)
  if (!identical(got, want)) {
    stop(what, " was built with density = \"", got, "\" but require_density_mode = \"",
         want, "\"; refusing mismatched weights.", call. = FALSE)
  }
  invisible(got)
}
