# ---------------------------------------------------------------------------
# Reference distribution for Stage 2 parameter tests
#
# Stage 1 estimates the T x T tau-covariance Sigma by a delete-one-block
# jackknife with B blocks.  Stage 2 projects it onto model parameters: for a
# fixed contrast a, SE^2 = a' Sigma_jk a.  Because Sigma_jk is ESTIMATED and
# not known,
#
#     a' Sigma_jk a / a' Sigma a  ~  chi2_nu / nu
#
# for a Wishart-like Sigma_jk on nu degrees of freedom, and therefore
#
#     z = theta / SE  ~  t_nu     exactly, not approximately.
#
# Reading z against a normal reference is what the package used to force
# downstream code to do -- it emitted no parameter p-value at all -- and it
# inflates Type-I error by 17x at p = 5e-8 when B = 65: a nominal 5e-8
# threshold actually delivers about 8.7e-7.  At p = 0.05 the same error is
# only 1.09x, which is why it went unnoticed.  See docs/KNOWN_ISSUES.md.
#
# nu comes from Stage 1 via <prefix>.cov.meta.json.  Covariance files written
# before that field existed carry no block count and no nu; for those the
# package refuses to emit a p-value rather than quietly reverting to the
# normal reference, because a silent fallback is exactly how this defect
# survived in the first place.
# ---------------------------------------------------------------------------

#' Path of the Stage 1 covariance metadata sidecar for a .cov.gz file
#'
#' @param cov_file Path to a Stage 1 \code{.cov.gz} file.
#' @return Path to the matching \code{.cov.meta.json}.
#' @keywords internal
.cov_meta_file <- function(cov_file) {
  sub("\\.cov\\.gz$", ".cov.meta.json", cov_file)
}

#' Pull one scalar out of a JSON document with base R
#'
#' Only ever applied to \code{.cov.meta.json}, which Stage 1 writes itself and
#' whose keys are unique across the document. Used when \pkg{jsonlite} is not
#' installed, so that reading the degrees of freedom never depends on a
#' package Stage 2 does not otherwise need.
#'
#' @param txt Whole JSON document as a single string.
#' @param pattern Regular expression with one capture group.
#' @return The captured string, or \code{NA_character_}.
#' @keywords internal
.json_capture <- function(txt, pattern) {
  m <- regexec(pattern, txt)
  g <- regmatches(txt, m)[[1]]
  if (length(g) < 2L) return(NA_character_)
  g[2L]
}

#' Extract a numeric scalar from a JSON document
#'
#' @param txt Whole JSON document as a single string.
#' @param key Key name.
#' @return Numeric scalar, or \code{NA_real_}.
#' @keywords internal
.json_number <- function(txt, key) {
  v <- .json_capture(
    txt,
    paste0("\"", key, "\"[[:space:]]*:[[:space:]]*(-?[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?)")
  )
  if (is.na(v)) return(NA_real_)
  suppressWarnings(as.numeric(v))
}

#' Extract a string scalar from a JSON document
#'
#' @param txt Whole JSON document as a single string.
#' @param key Key name.
#' @return Character scalar, or \code{NA_character_}.
#' @keywords internal
.json_string <- function(txt, key) {
  .json_capture(txt, paste0("\"", key, "\"[[:space:]]*:[[:space:]]*\"([^\"]*)\""))
}

#' Read the Stage 1 covariance metadata sidecar
#'
#' @param path Path to a \code{.cov.meta.json} file.
#' @return A list, or \code{NULL} if the file does not exist.
#' @keywords internal
.read_cov_meta <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    return(jsonlite::fromJSON(path, simplifyVector = TRUE))
  }
  txt <- paste(readLines(path, warn = FALSE), collapse = " ")
  list(
    cov_meta_version = .json_number(txt, "cov_meta_version"),
    jackknife = list(
      n_blocks  = .json_number(txt, "n_blocks"),
      nu        = .json_number(txt, "nu"),
      nu_source = .json_string(txt, "nu_source")
    )
  )
}

#' Resolve the degrees of freedom of the Stage 1 covariance
#'
#' Resolution order, matching \code{fungwas_stage1.cov_meta.resolve_nu}: an
#' explicitly recorded \code{jackknife$nu} wins; otherwise
#' \code{jackknife$n_blocks - 1}; otherwise the file predates the field and
#' \code{nu} is not knowable.
#'
#' Storing \code{nu} separately from \code{n_blocks} is deliberate. A separate
#' workstream is measuring whether the effective degrees of freedom are
#' materially lower than \code{B - 1}; when it reports, a better number can be
#' written into the same field with a different \code{nu_source} and picked up
#' here with no format change.
#'
#' @param meta List as returned by \code{.read_cov_meta()}, or \code{NULL}.
#' @return List with \code{nu} and \code{nu_source}.
#' @keywords internal
.resolve_nu <- function(meta) {
  none <- list(nu = NA_real_, nu_source = NA_character_)
  if (is.null(meta) || !is.list(meta)) return(none)
  jk <- meta[["jackknife"]]
  if (is.null(jk) || !is.list(jk)) return(none)

  nu <- suppressWarnings(as.numeric(jk[["nu"]]))
  nu <- if (length(nu)) nu[1L] else NA_real_
  if (is.finite(nu) && nu > 0) {
    src <- as.character(jk[["nu_source"]])
    src <- if (length(src)) src[1L] else NA_character_
    if (is.na(src) || !nzchar(src)) src <- "blocks_minus_one"
    return(list(nu = nu, nu_source = src))
  }

  nb <- suppressWarnings(as.numeric(jk[["n_blocks"]]))
  nb <- if (length(nb)) nb[1L] else NA_real_
  if (is.finite(nb) && nb >= 2) {
    return(list(nu = nb - 1, nu_source = "blocks_minus_one"))
  }
  none
}

#' Message used whenever nu could not be determined
#'
#' @param where Human-readable description of the covariance in question.
#' @return A character scalar.
#' @keywords internal
.unknown_nu_message <- function(where = "this covariance") {
  paste0(
    "Degrees of freedom of the Stage 1 covariance are unknown for ", where, ". ",
    "theta/SE is t distributed, not normal, because Sigma is estimated by a ",
    "block jackknife; with B = 65 blocks a normal reference inflates Type-I ",
    "error by 17x at p = 5e-8. Stage 1 now records the block count in ",
    "<prefix>.cov.meta.json, but this covariance predates that field. ",
    "Re-run Stage 1, or pass nu explicitly (nu = n_blocks - 1)."
  )
}

#' Choose the reference distribution for theta/SE
#'
#' @param nu Degrees of freedom, or \code{NA}/\code{NULL} if unknown.
#' @param nu_source Provenance string recorded alongside \code{nu}.
#' @param unknown_nu What to do when \code{nu} is unknown. \code{"omit"}
#'   (default) emits no p-value and warns; \code{"normal"} emits the old,
#'   wrong-in-the-tail normal p-value and warns, recording
#'   \code{ref_dist = "normal"} in the output so a reader can see it;
#'   \code{"error"} stops.
#' @param where Human-readable description used in messages.
#' @return List with \code{dist} ("t", "normal" or "none"), \code{df} and
#'   \code{nu_source}.
#' @keywords internal
.param_reference <- function(nu = NA_real_, nu_source = NA_character_,
                             unknown_nu = c("omit", "normal", "error"),
                             where = "this covariance") {
  unknown_nu <- match.arg(unknown_nu)
  nu <- suppressWarnings(as.numeric(nu))
  nu <- if (length(nu)) nu[1L] else NA_real_

  if (is.finite(nu) && nu > 0) {
    src <- as.character(nu_source)
    src <- if (length(src)) src[1L] else NA_character_
    if (is.na(src) || !nzchar(src)) src <- "supplied"
    return(list(dist = "t", df = nu, nu_source = src))
  }

  msg <- .unknown_nu_message(where)
  if (unknown_nu == "error") {
    stop(msg, call. = FALSE)
  }
  if (unknown_nu == "normal") {
    warning(msg,
            " Emitting NORMAL-reference p-values because unknown_nu = \"normal\";",
            " these are anti-conservative in the tail.",
            " The output records ref_dist = \"normal\".", call. = FALSE)
    return(list(dist = "normal", df = Inf, nu_source = NA_character_))
  }
  warning(msg,
          " No parameter p-value will be emitted (unknown_nu = \"omit\").",
          " Do not recompute one with pnorm(): that is the defect this",
          " refusal exists to prevent.", call. = FALSE)
  list(dist = "none", df = NA_real_, nu_source = NA_character_)
}

#' Two-sided p-values for parameter estimates
#'
#' Exact for the fixed-contrast projection (\code{tau_cov_estimator = "ols"}),
#' where the contrast does not depend on the estimated covariance. For GLS the
#' weight matrix is itself a function of \code{Sigma_jk}, so \code{t_nu} is a
#' very good approximation there rather than an identity.
#'
#' @param params K x P matrix of estimates.
#' @param SE_params K x P matrix of standard errors.
#' @param ref List from \code{.param_reference()}.
#' @return K x P matrix of p-values, or \code{NULL} when no reference applies.
#' @keywords internal
.param_pvalues <- function(params, SE_params, ref) {
  if (is.null(ref) || identical(ref$dist, "none")) return(NULL)
  params <- as.matrix(params)
  SE_params <- as.matrix(SE_params)
  z <- params / SE_params
  p <- if (identical(ref$dist, "t")) {
    2 * stats::pt(-abs(z), df = ref$df)
  } else {
    2 * stats::pnorm(-abs(z))
  }
  p[!is.finite(z)] <- NA_real_
  dim(p) <- dim(params)
  dimnames(p) <- dimnames(params)
  p
}
