# ---------------------------------------------------------------------------
# Reference distribution for Stage 2 parameter tests
#
# Sigma is ESTIMATED, not known: Stage 1 gets it from a delete-one-block
# jackknife. For a fixed contrast a and a Wishart-like Sigma_jk on nu degrees
# of freedom,
#
#     a' Sigma_jk a / a' Sigma a  ~  chi2_nu / nu      =>   theta/SE ~ t_nu
#
# Reading that against a normal is what the package used to force downstream
# code to do, and it inflates Type-I error by 17x at p = 5e-8 when B = 65.
#
# HOW MANY DEGREES OF FREEDOM. Not B - 1, and not a constant. Measured
# (fungwas_paper/docs/BP1_effective_df.md, 2026-09-14):
#
#   * nu at N = 50,000 is 62.4 / 445.9 / 801.8 / 1351.0 at B = 65/500/1000/2000
#   * at B = 1000 it runs 530 (MAF 0.01) to 979 (MAF 0.50)
#   * at B = 1000 it is 566 for a scale-type contrast against 802 for a
#     location-type one, and that gap widens with B
#
# The stable object is the law, fitted over 28 cells to a median error of 0.8%:
#
#     1/nu = 1/(B - 1) + c / N_eff        N_eff = D^2 / sum_i g_ri^4
#     c = 0.88 (location-type contrast), 2.62 (scale-type)
#
# The first term is the nominal jackknife df; the second is a Satterthwaite
# penalty for unequal block contributions (block_ids are drawn independently
# per individual, so block totals D_k are multinomial). For a genotype the
# fourth central moment equals the second, so N_eff collapses to
# D = sum_i g_resid_i^2, which Stage 1 emits per SNP as the d_resid column.
#
# WHY c IS NOT COMPUTED ANALYTICALLY HERE. It is derivable. The per-individual
# score for contrast a is u(v) = sum_t (a_t/f_t) (tau_t - 1{v <= tau_t}), a step
# function of the phenotype rank v, so its moments are exact finite sums and
#
#     c = (1 - 1/B) + gamma_e / 2,   gamma_e = excess kurtosis of u(v)
#
# For a normal phenotype on the production grid (45 taus, 0.02 to 0.98) this
# gives u proportional to a clipped z for a flat contrast and to a clipped
# z^2 for the q-centred one, whence gamma_e = -0.506 and +1.58 and therefore
# c = 0.75 and 1.79. The measured values are 0.88 and 2.62. The derivation
# gets the mechanism and the ordering right and the magnitude wrong, by 15%
# and 32%, in the ANTI-CONSERVATIVE direction (too small a c means too large
# a nu means too small a p-value); the report independently notes that the
# second-order kurtosis expansion under-predicts, by 2.4x per-tau. It also
# needs f_tau, which Stage 1 does not emit. So the measured constants are used
# and the derivation is recorded here rather than shipped. See the package
# report for the full working.
#
# Contrast files predating d_resid, and covariances predating the block count,
# are handled explicitly rather than by a silent fallback -- a silent fallback
# is exactly how the normal-reference defect survived.
# ---------------------------------------------------------------------------

#' Satterthwaite coefficient for a flat (location-type) contrast over taus
#' @keywords internal
C_LOCATION_DEFAULT <- 0.88

#' Satterthwaite coefficient for a tail-loaded (scale-type) contrast over taus
#' @keywords internal
C_SCALE_DEFAULT <- 2.62

#' N_eff assumed when the per-SNP value is unavailable
#'
#' The production value at MAF 0.01, so conservative at every commoner
#' variant. Reproduces the report's recommended flat values (nu = 911 for a
#' location contrast, 776 for a scale one at B = 1000) and, unlike a bare
#' constant, generalises to other block counts.
#' @keywords internal
FALLBACK_N_EFF_DEFAULT <- 9134

#' Stage 1 TSV column carrying the per-SNP N_eff
#' @keywords internal
N_EFF_COLUMN <- "d_resid"

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
      n_blocks          = .json_number(txt, "n_blocks"),
      n_blocks_nonempty = .json_number(txt, "n_blocks_nonempty"),
      nu                = .json_number(txt, "nu"),
      nu_source         = .json_string(txt, "nu_source"),
      df_model = list(
        n_eff_column   = .json_string(txt, "n_eff_column"),
        c_location     = .json_number(txt, "c_location"),
        c_scale        = .json_number(txt, "c_scale"),
        fallback_n_eff = .json_number(txt, "fallback_n_eff")
      )
    )
  )
}

#' First finite positive scalar among its arguments
#' @keywords internal
.first_positive <- function(...) {
  for (v in list(...)) {
    v <- suppressWarnings(as.numeric(v))
    if (length(v) && is.finite(v[1L]) && v[1L] > 0) return(v[1L])
  }
  NA_real_
}

#' Number of jackknife blocks that carried a replicate
#'
#' Prefers \code{n_blocks_nonempty}: deleting an empty block reproduces the
#' full-sample estimate and contributes no replicate.
#'
#' @param meta List from \code{.read_cov_meta()}, or \code{NULL}.
#' @return Numeric block count, or \code{NA_real_}.
#' @keywords internal
.resolve_block_count <- function(meta) {
  if (is.null(meta) || !is.list(meta)) return(NA_real_)
  jk <- meta[["jackknife"]]
  if (is.null(jk) || !is.list(jk)) return(NA_real_)
  b <- .first_positive(jk[["n_blocks_nonempty"]], jk[["n_blocks"]])
  if (is.finite(b) && b >= 2) b else NA_real_
}

#' Degrees-of-freedom law recorded by Stage 1, with package defaults filled in
#'
#' @param meta List from \code{.read_cov_meta()}, or \code{NULL}.
#' @return List with \code{c_location}, \code{c_scale}, \code{fallback_n_eff}
#'   and \code{n_eff_column}.
#' @keywords internal
.resolve_df_model <- function(meta) {
  model <- NULL
  if (is.list(meta)) {
    jk <- meta[["jackknife"]]
    if (is.list(jk)) model <- jk[["df_model"]]
  }
  if (!is.list(model)) model <- list()
  n_eff_col <- as.character(model[["n_eff_column"]])[1L]
  if (is.null(n_eff_col) || is.na(n_eff_col) || !nzchar(n_eff_col)) n_eff_col <- N_EFF_COLUMN
  list(
    c_location     = .first_positive(model[["c_location"]], C_LOCATION_DEFAULT),
    c_scale        = .first_positive(model[["c_scale"]], C_SCALE_DEFAULT),
    fallback_n_eff = .first_positive(model[["fallback_n_eff"]], FALLBACK_N_EFF_DEFAULT),
    n_eff_column   = n_eff_col
  )
}

#' Read a constant nu deliberately recorded in the metadata
#'
#' There is no \code{B - 1} fallback here on purpose: the effective degrees of
#' freedom are measurably below \code{B - 1} and vary by SNP and by parameter,
#' so deriving a constant from the block count would manufacture a number that
#' measurement contradicts.
#'
#' @param meta List from \code{.read_cov_meta()}, or \code{NULL}.
#' @return List with \code{nu} and \code{nu_source}.
#' @keywords internal
.resolve_nu <- function(meta) {
  none <- list(nu = NA_real_, nu_source = NA_character_)
  if (is.null(meta) || !is.list(meta)) return(none)
  jk <- meta[["jackknife"]]
  if (is.null(jk) || !is.list(jk)) return(none)
  nu <- .first_positive(jk[["nu"]])
  if (!is.finite(nu)) return(none)
  src <- as.character(jk[["nu_source"]])[1L]
  if (is.null(src) || is.na(src) || !nzchar(src)) src <- "unrecorded"
  list(nu = nu, nu_source = src)
}

#' How flat a contrast is across taus
#'
#' 1 for a contrast that weights every tau equally (a location-type contrast,
#' e.g. the first row of \code{(W'W)^-1 W'} when \code{W} has an intercept
#' column), 0 for one that sums to zero (a scale-type contrast, e.g. the row
#' for a mean-centred \code{q_tau} column).
#'
#' @param a Numeric contrast vector over taus.
#' @return Numeric in [0, 1], or \code{NA_real_}.
#' @keywords internal
.contrast_flatness <- function(a) {
  a <- as.numeric(a)
  n <- length(a)
  if (!n) return(NA_real_)
  nrm <- sqrt(sum(a * a))
  if (!is.finite(nrm) || nrm <= 0) return(NA_real_)
  abs(sum(a)) / (sqrt(n) * nrm)
}

#' Satterthwaite coefficient c for a contrast
#'
#' The two values were measured for a location-type and a scale-type contrast
#' under \code{W = [1, q_tau - mean(q_tau)]} on a 45-tau grid. A contrast that
#' is not clearly flat gets the scale value, which is the conservative choice:
#' larger \code{c} means smaller \code{nu} means a heavier-tailed reference and
#' a larger p-value.
#'
#' @param a Numeric contrast vector over taus.
#' @param model List from \code{.resolve_df_model()}.
#' @param flat_threshold Flatness at or above which a contrast counts as
#'   location-type.
#' @return Numeric coefficient.
#' @keywords internal
.contrast_c <- function(a, model = .resolve_df_model(NULL), flat_threshold = 0.9) {
  f <- .contrast_flatness(a)
  if (is.finite(f) && f >= flat_threshold) model$c_location else model$c_scale
}

#' Effective degrees of freedom of the jackknife covariance
#'
#' \code{1/nu = 1/(B - 1) + c/N_eff}. Vectorised over \code{n_eff}.
#'
#' @param n_blocks Number of populated jackknife blocks.
#' @param n_eff Effective sample size of the SNP (Stage 1's \code{d_resid}).
#' @param c Satterthwaite coefficient for the contrast.
#' @return Numeric vector of degrees of freedom.
#' @keywords internal
.effective_nu <- function(n_blocks, n_eff, c) {
  n_blocks <- as.numeric(n_blocks)[1L]
  inv <- 1 / (n_blocks - 1) + as.numeric(c) / as.numeric(n_eff)
  nu <- 1 / inv
  nu[!is.finite(nu) | nu <= 0] <- NA_real_
  nu
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
    "<prefix>.cov.meta.json and the per-SNP N_eff in the d_resid column, but ",
    "this covariance predates both. Re-run Stage 1, or pass nu explicitly."
  )
}

#' Resolve the reference distribution and its degrees of freedom
#'
#' Ladder, best first:
#' \enumerate{
#'   \item an explicit \code{nu} from the caller;
#'   \item the block count and a per-SNP \code{n_eff}: the measured law, one
#'     \code{nu} per SNP per parameter;
#'   \item the block count alone: the same law at \code{fallback_n_eff}, with
#'     a warning;
#'   \item a constant \code{nu} recorded in the metadata, with a warning that a
#'     constant is not defensible at \code{B >= 500};
#'   \item nothing, in which case \code{unknown_nu} decides.
#' }
#'
#' @param A Mapping matrix (K x T); its rows are the contrasts over taus.
#' @param P Number of SNPs.
#' @param nu Explicit degrees of freedom, or \code{NULL}.
#' @param nu_source Provenance for \code{nu}.
#' @param n_blocks Number of populated jackknife blocks, or \code{NULL}.
#' @param n_eff Per-SNP effective sample size (length \code{P}), or \code{NULL}.
#' @param df_model List from \code{.resolve_df_model()}.
#' @param df_coef Optional override for \code{c}: one number, or a vector named
#'   by parameter.
#' @param unknown_nu One of "omit", "normal", "error".
#' @param where Human-readable description used in messages.
#' @return List with \code{dist}, \code{df} (K x P matrix or \code{NULL}),
#'   \code{source} and \code{c}.
#' @keywords internal
.param_reference <- function(A, P,
                             nu = NULL, nu_source = NULL,
                             n_blocks = NULL, n_eff = NULL,
                             df_model = NULL, df_coef = NULL,
                             unknown_nu = c("omit", "normal", "error"),
                             where = "this covariance") {
  unknown_nu <- match.arg(unknown_nu)
  A <- as.matrix(A)
  K <- nrow(A)
  param_names <- rownames(A)
  if (is.null(df_model)) df_model <- .resolve_df_model(NULL)

  as_df_matrix <- function(v) {
    m <- matrix(as.numeric(v), nrow = K, ncol = P)
    rownames(m) <- param_names
    m
  }

  # --- 1. explicit nu from the caller ------------------------------------
  nu1 <- if (is.null(nu)) NA_real_ else .first_positive(nu)
  if (is.finite(nu1)) {
    src <- as.character(nu_source)[1L]
    if (is.null(src) || is.na(src) || !nzchar(src)) src <- "user_supplied"
    return(list(dist = "t", df = as_df_matrix(nu1), source = src, c = NULL))
  }

  # --- c per parameter, needed by rungs 2 and 3 --------------------------
  c_vec <- vapply(seq_len(K), function(k) .contrast_c(A[k, ], df_model), numeric(1))
  names(c_vec) <- param_names
  if (!is.null(df_coef)) {
    if (length(df_coef) == 1L && is.null(names(df_coef))) {
      c_vec[] <- as.numeric(df_coef)
    } else {
      hit <- intersect(names(df_coef), names(c_vec))
      if (!length(hit)) {
        stop("df_coef must be a single number or a vector named by parameter; ",
             "names given: ", paste(names(df_coef), collapse = ", "), call. = FALSE)
      }
      c_vec[hit] <- as.numeric(df_coef[hit])
    }
  }

  B <- if (is.null(n_blocks)) NA_real_ else .first_positive(n_blocks)
  if (is.finite(B) && B >= 2) {
    n_eff_v <- if (is.null(n_eff)) rep(NA_real_, P) else as.numeric(n_eff)
    if (length(n_eff_v) != P) n_eff_v <- rep(n_eff_v, length.out = P)
    good <- is.finite(n_eff_v) & n_eff_v > 0

    # --- 2. the measured law, per SNP -----------------------------------
    if (any(good)) {
      if (!all(good)) {
        warning(sum(!good), " of ", P, " SNPs have no usable ",
                df_model$n_eff_column, "; using fallback_n_eff = ",
                df_model$fallback_n_eff, " for those.", call. = FALSE)
      }
      n_eff_v[!good] <- df_model$fallback_n_eff
      df <- matrix(NA_real_, nrow = K, ncol = P, dimnames = list(param_names, NULL))
      for (k in seq_len(K)) df[k, ] <- .effective_nu(B, n_eff_v, c_vec[k])
      return(list(dist = "t", df = df, source = "per_snp_satterthwaite", c = c_vec))
    }

    # --- 3. the same law at a fallback N_eff ----------------------------
    warning("No per-SNP ", df_model$n_eff_column, " available for ", where,
            ". Falling back to a flat effective df computed at N_eff = ",
            df_model$fallback_n_eff, " (production MAF 0.01, so conservative at ",
            "commoner variants). Re-run Stage 1 to get the per-SNP value; a flat ",
            "df costs up to 2.4% of nominal Type-I error for a location ",
            "parameter and 6.5% for a scale one.", call. = FALSE)
    df <- matrix(NA_real_, nrow = K, ncol = P, dimnames = list(param_names, NULL))
    for (k in seq_len(K)) {
      df[k, ] <- .effective_nu(B, df_model$fallback_n_eff, c_vec[k])
    }
    return(list(dist = "t", df = df, source = "flat_fallback_n_eff", c = c_vec))
  }

  # --- 5. nothing usable -------------------------------------------------
  msg <- .unknown_nu_message(where)
  if (unknown_nu == "error") stop(msg, call. = FALSE)
  if (unknown_nu == "normal") {
    warning(msg,
            " Emitting NORMAL-reference p-values because unknown_nu = \"normal\";",
            " these are anti-conservative in the tail.",
            " The output records ref_dist = \"normal\".", call. = FALSE)
    return(list(dist = "normal", df = as_df_matrix(Inf), source = NA_character_, c = NULL))
  }
  warning(msg,
          " No parameter p-value will be emitted (unknown_nu = \"omit\").",
          " Do not recompute one with pnorm(): that is the defect this",
          " refusal exists to prevent.", call. = FALSE)
  list(dist = "none", df = NULL, source = NA_character_, c = NULL)
}

#' Two-sided p-values for parameter estimates
#'
#' Exact for the fixed-contrast projection (\code{tau_cov_estimator = "ols"}),
#' where the contrast does not depend on the estimated covariance. For GLS the
#' weight matrix is itself a function of \code{Sigma_jk}, so \code{t_nu} there
#' is a very good approximation rather than an identity, and the contrast used
#' to classify \code{c} is the OLS one.
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
    df <- ref$df
    if (length(df) != length(z)) df <- matrix(as.numeric(df), nrow = nrow(z), ncol = ncol(z))
    2 * stats::pt(-abs(as.numeric(z)), df = as.numeric(df))
  } else {
    2 * stats::pnorm(-abs(as.numeric(z)))
  }
  p[!is.finite(as.numeric(z))] <- NA_real_
  dim(p) <- dim(params)
  dimnames(p) <- dimnames(params)
  p
}
