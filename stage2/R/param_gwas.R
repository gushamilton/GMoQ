#' Model tournament for parametric GWAS
#'
#' Fits multiple parametric weight matrices to Stage 1 results in one pass
#' and returns per-SNP AIC/BIC values for model selection.
#'
#' @param stage1 Output list from quantile_gwas() containing Q_slope and SE_tau.
#' @param W_list List of weight matrices (each T x K).
#' @param N Sample size used for BIC.
#' @param n_threads Number of OpenMP threads.
#' @param se_mode One of \code{"tau_cov"} or \code{"dwls"} (or \code{"diagonal"}).
#'   \code{"tau_cov"} uses per-SNP covariance (recommended when available).
#' @param cov_data Optional list with \code{cov_vec} and \code{offsets} for
#'   \code{se_mode = "tau_cov"}.
#' @param residual_loss Which residual loss the score uses. \code{"auto"}
#'   (default) preserves the historical behaviour: \code{se_mode = "tau_cov"}
#'   fits GLS and scores \code{r' Sigma^-1 r + 2k}, while \code{"dwls"} /
#'   \code{"diagonal"} score the diagonal loss. \code{"full_cov"} fits OLS
#'   parameters and scores the residuals with the per-SNP full covariance,
#'   \code{Q_full = r' Sigma^-1 r}, returning the coherent OLS score
#'   \code{Q_full + 2 rank(W)} and \code{Q_full + rank(W) log(N)}. It requires
#'   \code{se_mode = "tau_cov"} and \code{cov_data}.
#' @param return_params Include the fitted OLS parameter matrices in the
#'   \code{residual_loss = "full_cov"} return value. Default \code{FALSE}
#'   because a genome-wide scan does not need to hold them.
#'
#' @return With \code{residual_loss = "auto"}, a matrix with AIC/BIC columns per
#'   model, exactly as before. With \code{residual_loss = "full_cov"}, a list
#'   with \code{Q_full} (P x n_models), \code{rank} (length n_models),
#'   \code{AIC} and \code{BIC} (both P x n_models), \code{reason} (P x
#'   n_models character), \code{rank_deficient} (logical per model),
#'   \code{model_names}, and (when requested) \code{params}.
#'
#' @export
param_gwas_multi <- function(stage1, W_list, N, n_threads = 1L,
                             se_mode = c("tau_cov", "dwls", "diagonal"),
                             cov_data = NULL,
                             residual_loss = c("auto", "full_cov"),
                             return_params = FALSE) {
  se_mode <- match.arg(se_mode)
  residual_loss <- match.arg(residual_loss)
  if (is.null(stage1$Q_slope) || is.null(stage1$SE_tau)) {
    stop("stage1 must include Q_slope and SE_tau.")
  }
  if (!is.list(W_list) || length(W_list) == 0) {
    stop("W_list must be a non-empty list of weight matrices.")
  }

  if (residual_loss == "full_cov") {
    if (se_mode != "tau_cov") {
      stop("residual_loss = 'full_cov' requires se_mode = 'tau_cov'.")
    }
    if (is.null(cov_data) || is.null(cov_data$cov_vec) ||
        is.null(cov_data$offsets)) {
      stop("residual_loss = 'full_cov' requires cov_data with cov_vec and offsets.")
    }
    return(.ols_full_cov_multi(stage1, W_list, N, n_threads, cov_data,
                               return_params = return_params))
  }

  beta_stage1 <- t(stage1$Q_slope)
  se_stage1 <- t(stage1$SE_tau)

  if (se_mode == "tau_cov") {
    if (is.null(cov_data) || is.null(cov_data$cov_vec) || is.null(cov_data$offsets)) {
      stop("se_mode = 'tau_cov' requires cov_data with cov_vec and offsets.")
    }
    result <- fit_multi_models_tau_cov(
      beta_stage1, cov_data$cov_vec, cov_data$offsets, W_list, N, n_threads,
      cov_scales = cov_data$cov_scales
    )
  } else {
    result <- fit_multi_models(beta_stage1, se_stage1, W_list, N, n_threads)
  }

  model_names <- names(W_list)
  if (is.null(model_names) || any(model_names == "")) {
    model_names <- paste0("model", seq_along(W_list))
  }
  colnames(result) <- as.vector(rbind(
    paste0("AIC_", model_names),
    paste0("BIC_", model_names)
  ))

  result
}

## Numeric reason codes returned by fit_multi_models_ols_full_cov(). Kept here
## rather than in C++ so the labels are visible and testable.
.ols_full_cov_reason_codes <- c(
  "0" = "ok",
  "1" = "beta_nonfinite",
  "2" = "offset_bad",
  "3" = "cov_scale_invalid",
  "4" = "cov_nonfinite",
  "5" = "sigma_not_pd",
  "6" = "sigma_ill_conditioned",
  "7" = "sigma_solve_fail",
  "8" = "theta_solve_fail",
  "9" = "q_invalid"
)

## OLS parameter fit + full-covariance residual score. Internal: the public
## entry point is param_gwas_multi(residual_loss = "full_cov").
.ols_full_cov_multi <- function(stage1, W_list, N, n_threads, cov_data,
                                return_params = FALSE) {
  beta_stage1 <- t(stage1$Q_slope)
  model_names <- names(W_list)
  if (is.null(model_names) || any(model_names == "")) {
    model_names <- paste0("model", seq_along(W_list))
  }
  snp_names <- colnames(stage1$Q_slope)

  res <- fit_multi_models_ols_full_cov(
    beta_stage1 = beta_stage1,
    cov_vec = cov_data$cov_vec,
    offsets = cov_data$offsets,
    W_list = W_list,
    N = N,
    n_threads = as.integer(n_threads),
    cov_scales = cov_data$cov_scales
  )

  dimnames(res$Q_full) <- list(snp_names, model_names)
  dimnames(res$AIC) <- list(snp_names, model_names)
  dimnames(res$BIC) <- list(snp_names, model_names)

  rank <- res$rank
  names(rank) <- model_names
  rank_deficient <- res$rank_deficient
  names(rank_deficient) <- model_names

  reason <- matrix(
    .ols_full_cov_reason_codes[as.character(res$reason_code)],
    nrow = nrow(res$reason_code), ncol = ncol(res$reason_code),
    dimnames = list(snp_names, model_names)
  )

  out <- list(
    Q_full = res$Q_full,
    rank = rank,
    AIC = res$AIC,
    BIC = res$BIC,
    reason = reason,
    rank_deficient = rank_deficient,
    model_names = model_names,
    N = N
  )

  if (isTRUE(return_params)) {
    params <- vector("list", length(W_list))
    names(params) <- model_names
    for (m in seq_along(W_list)) {
      k_m <- ncol(W_list[[m]])
      p <- matrix(as.numeric(res$params[[m]]), nrow = k_m)
      if (!is.null(colnames(W_list[[m]]))) rownames(p) <- colnames(W_list[[m]])
      colnames(p) <- snp_names
      params[[m]] <- p
    }
    out$params <- params
  }

  out
}
