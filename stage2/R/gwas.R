#' Quantile-level RIF GWAS
#'
#' Runs a genome-wide association scan at user-specified quantile levels (taus)
#' using Recentered Influence Functions (RIF). Supports two covariate modes:
#' - "residualize": partial out C from SNPs (but never from Y).
#' - "include": include C directly in the regression alongside each SNP.
#'
#' @param Y Numeric vector (length N), phenotype.
#' @param G Numeric matrix (N x P) of SNP dosages or genotypes.
#' @param taus Numeric vector of quantile levels.
#' @param C Optional N x K covariate matrix.
#' @param covar_mode Either "residualize" (default) or "include".
#'   - "residualize": partial out C from G before regression.
#'   - "include": include C directly in regression.
#' @param density_floor Positive scalar; lower bound for estimated densities.
#' @param benchmark Logical; if TRUE, include timing info.
#' @param verbose Logical; print progress messages.
#'
#' @return A list with taus, q_tau, fhat_tau, Q_slope, SE_tau, timing.
#' 
#' @export
quantile_gwas <- function(
  Y, G,
  taus = seq(0.10, 0.90, 0.05),
  C = NULL,
  covar_mode = c("residualize", "include"),
  density_floor = 1e-8,
  benchmark = TRUE,
  verbose = TRUE
) {
  covar_mode <- match.arg(covar_mode)
  Y <- as.numeric(Y)
  G <- as.matrix(G)
  N <- length(Y); P <- ncol(G); Tt <- length(taus)
  if (nrow(G) != N) stop("nrow(G) must equal length(Y).")

  t_all0 <- proc.time()[["elapsed"]]

  # (1) Build RIF on raw Y
  if (verbose) message("Building RIF matrix on raw Y...")
  t0 <- proc.time()[["elapsed"]]
  rif <- .build_rif_matrix(Y, taus, density_floor = density_floor)
  RIF     <- rif$R
  q_tau   <- rif$q_tau
  fhat    <- rif$fhat_tau
  R_mean  <- colMeans(RIF)
  t1 <- proc.time()[["elapsed"]]
  time_rif <- t1 - t0

  # (2) SNP regressions
  if (verbose) message("Computing per-SNP tau-slopes...")
  t2 <- proc.time()[["elapsed"]]
  Q_slope <- matrix(NA_real_, nrow = Tt, ncol = P,
                    dimnames = list(paste0("tau", taus), colnames(G)))
  SE_tau  <- matrix(NA_real_, nrow = Tt, ncol = P,
                    dimnames = list(paste0("tau", taus), colnames(G)))

  if (covar_mode == "residualize") {
    # Partial out C from G (never from Y)
    if (!is.null(C)) {
      if (verbose) message("Residualizing SNPs on covariates...")
      for (j in seq_len(P)) {
        G[, j] <- .residualize_on_C(G[, j], C)
      }
    }
    # Univariate regression of RIF on each SNP
    for (j in seq_len(P)) {
      x <- G[, j]
      x_c <- x - mean(x)
      Sxx <- sum(x_c * x_c)
      if (!(is.finite(Sxx) && Sxx > 0)) next

      b <- as.numeric(crossprod(x_c, RIF) / Sxx)
      a <- R_mean - b * mean(x)
      term1 <- matrix(rep(a, each = N), nrow = N, ncol = Tt)
      term2 <- outer(x, b)
      E     <- RIF - term1 - term2

      SSE    <- colSums(E * E)
      sigma2 <- SSE / max(N - 2L, 1L)
      se_b   <- sqrt(sigma2 / Sxx)

      Q_slope[, j] <- b
      SE_tau[, j]  <- se_b
    }

  } else if (covar_mode == "include") {
    stop("covar_mode = 'include' is still in development")

    # Joint regression RIF ~ SNP + C
    if (is.null(C)) stop("Must provide C if covar_mode = 'include'.")
    C <- as.matrix(C)

    for (j in seq_len(P)) {
      X <- cbind(G[, j], C)  # SNP + covariates
      fit <- lm.fit(X, RIF)  # multi-response regression (N x T outcome)
      b <- fit$coefficients[1, ]      # SNP slope
      Q_slope[, j] <- b

      # compute SEs for SNP slope (per tau)
      resid <- RIF - X %*% fit$coefficients
      sigma2 <- colSums(resid^2) / (N - ncol(X))
      XtXinv <- chol2inv(chol(crossprod(X)))
      se_b <- sqrt(sigma2 * XtXinv[1, 1])
      SE_tau[, j] <- se_b
    }
  }

  t3 <- proc.time()[["elapsed"]]
  time_gwas <- t3 - t2

  t_all1 <- proc.time()[["elapsed"]]
  timing <- list(rif_build_sec = time_rif,
                 gwas_loop_sec = time_gwas,
                 total_sec     = t_all1 - t_all0)

  out <- list(
    taus     = taus,
    q_tau    = q_tau,
    fhat_tau = fhat,
    Q_slope  = Q_slope,
    SE_tau   = SE_tau
  )
  if (benchmark) out$timing <- timing
  out
}


#' Parametric GWAS mapping
#'
#' Maps quantile-level GWAS results into effects on user-specified
#' parametric system parameters via a weight matrix \code{W}.
#'
#' @param stage1 Output list from \code{quantile_gwas()}.
#' @param transform One of:
#'   \itemize{
#'     \item \code{"custom_W"} – supply a weight matrix \code{W} in \code{transform_args}.
#'     \item \code{"two_normal"} – construct weights for a two-component normal mixture.
#'   }
#' @param transform_args Arguments for the chosen transform.
#' @param se_mode Standard error mode:
#'   \itemize{
#'     \item \code{"tau_cov"} – use per-SNP jackknife covariance from Stage 1 (most accurate).
#'     \item \code{"diagonal"} – assume independence across taus (fast fallback only).
#'     \item \code{"plugin_cor"} – plugin correlation with near-PD repair.
#'     \item \code{"dwls"} – diagonal weighted least squares (Q test calibrated).
#'   }
#' @param plugin_R Optional T x T correlation matrix for \code{"plugin_cor"}.
#' @param cov_data For \code{se_mode = "tau_cov"}: either a list with \code{cov_vec}
#'   (numeric vector of covariance data) and \code{offsets} (byte offsets per SNP),
#'   optionally with \code{cov_scales} (one positive multiplier per SNP for a
#'   shared pooled payload),
#'   \code{se_mode = "tau_cov"} requires covariance data; use
#'   \code{se_mode = "diagonal"} explicitly only for fallback/debug runs.
#' @param tau_cov_estimator Estimator used when \code{se_mode = "tau_cov"} and
#'   covariance data are available. Defaults to \code{"ols"}: the fixed
#'   unweighted projection with covariance-propagated SEs.
#'
#'   \code{"gls"} uses the per-SNP tau-by-tau covariance as the weight matrix,
#'   and is the efficient estimator when the fitted model is correct. Measured
#'   on 28 known-truth cells at 20,000 replicates (see \code{docs/OLS_VS_GLS.md}),
#'   GLS estimates are 17-24 percent less variable than OLS on the additive,
#'   sigma2 and kurtosis parameters, at the same coverage: both give 0.9499
#'   under the null against a 0.95 nominal.
#'
#'   OLS is the default because it degrades more gracefully when the fitted
#'   model is wrong. Where lack of fit exceeds three per-tau standard errors,
#'   GLS coverage falls to 0.9359 against OLS's 0.9448 and its bias is four
#'   times larger (0.128 against 0.033 empirical SDs); it also appears most
#'   efficient there, because the estimates cluster tightly around the wrong
#'   answer. Those cells are concentrated at rare MAF and skewed phenotypes
#'   carrying a scale effect, which is where a variance-QTL claim would live.
#'
#'   Prefer \code{"gls"} for analyses whose model fit has been checked --
#'   \code{Q} measures exactly the separating quantity per variant.
#' @param gls_ridge Relative diagonal ridge added only when a per-SNP covariance
#'   matrix cannot be inverted for GLS.
#' @param n_threads Number of OpenMP threads for tau-covariance estimation.
#' @param return_cov Optional output of parameter covariance derived from Stage 1
#'   tau covariance. One of \code{"none"} (default), \code{"upper"} (upper-triangle
#'   matrix), \code{"full"} (K x K x P array), or \code{"mashr"} (Bhat/Shat/V list).
#' @param q_method Reference for the lack-of-fit statistic: calibrated
#'   \code{"saddlepoint"}, \code{"liu"}, or \code{"satterthwaite"} on the OLS
#'   path; \code{"hotelling"} on the GLS path; or explicit
#'   \code{"chisq_naive"}, which returns \code{Q} without \code{Q_pval}.
#' @param nu Explicit degrees of freedom, overriding everything else.
#'   \code{theta/SE} is \code{t_nu} distributed, not normal, because
#'   \code{Sigma} is estimated: for a fixed contrast \code{a},
#'   \code{a' Sigma_jk a / a' Sigma a ~ chi2_nu/nu}. Normally leave this
#'   \code{NULL} and supply \code{n_blocks} and \code{n_eff} instead, from
#'   which a per-SNP, per-parameter \code{nu} is computed; a single \code{nu}
#'   is not defensible, since the measured value ranges 530 to 979 across MAF
#'   at \code{n_blocks = 1000} and differs by 236 between a scale and a
#'   location contrast.
#' @param nu_source Free-text provenance for \code{nu}, carried into the output
#'   so a result can be traced back to how its degrees of freedom were obtained.
#' @param n_blocks Number of populated jackknife blocks behind \code{cov_data}.
#'   \code{param_gwas_from_file()} reads it from \code{.cov.meta.json}.
#' @param n_eff Per-SNP effective sample size, length \code{ncol(Q_slope)}.
#'   This is Stage 1's \code{d_resid} column, \code{sum_i g_resid_i^2}.
#'   Together with \code{n_blocks} it gives the measured law
#'   \code{1/nu = 1/(n_blocks - 1) + c/n_eff}.
#' @param df_model Optional list overriding the coefficients of that law
#'   (\code{c_location}, \code{c_scale}, \code{fallback_n_eff}); normally read
#'   from the Stage 1 metadata.
#' @param df_coef Optional override for \code{c} itself: one number, or a
#'   vector named by parameter. Use it when \code{W} is not a location/scale
#'   design, since \code{c} was calibrated for one that is.
#' @param unknown_nu What to do when \code{nu} is unknown. \code{"omit"}
#'   (default) emits no parameter p-value and warns; \code{"normal"} emits the
#'   old normal-reference p-value and warns, recording
#'   \code{param_ref_dist = "normal"}; \code{"error"} stops. There is
#'   deliberately no silent normal fallback: it inflates Type-I error 17x at
#'   p = 5e-8 when \code{n_blocks = 65}, and being silent is why that went
#'   unnoticed.
#'
#' @return A list with the following elements:
#' \item{\code{W}}{Weight matrix used in the parametric mapping.}
#' \item{\code{A}}{Mapping matrix (tau-by-parameter).}
#' \item{\code{params}}{Estimated parameter values.}
#' \item{\code{SE_params}}{Delta-method standard errors of parameter estimates.}
#' \item{\code{P_params}}{Two-sided parameter p-values, \code{2 * pt(-|theta/SE|, nu)},
#'   or \code{NULL} when no reference distribution could be established.}
#' \item{\code{param_ref_dist}}{Reference distribution actually used for those
#'   p-values: "t", "normal" or "none".}
#' \item{\code{param_ref_df}}{K x P matrix of degrees of freedom, one per
#'   parameter per SNP.}
#' \item{\code{param_ref_source}}{How those degrees of freedom were obtained.}
#' \item{\code{param_ref_c}}{Satterthwaite coefficient used per parameter.}
#' \item{\code{Q}}{Vector of Q statistics per SNP (if computed).}
#' \item{\code{df}}{Degrees of freedom for Q statistic.}
#' \item{\code{Q_pval}}{Calibrated lack-of-fit p-value. On the OLS path the
#'   reference is the weighted chi-square mixture selected by \code{q_method};
#'   on the GLS path it is Hotelling's \code{F_{T-K, nu-(T-K)+1}}, which is
#'   exact. Absent when \code{q_method = "chisq_naive"}, and absent on the GLS
#'   path when \code{nu} is unknown and \code{unknown_nu = "omit"}.}
#' \item{\code{Q_df_eff}, \code{Q_scale}}{Chosen so that
#'   \code{Q_df_eff * Q_scale} is the null mean of \code{Q}. The calibrated
#'   analogue of the raw \code{Q / (T - K)} ratio is therefore
#'   \code{Q / (Q_df_eff * Q_scale)}, which has null mean 1 on both paths.}
#' \item{\code{Q_ref_dist}, \code{Q_ref_df}, \code{Q_ref_source}}{GLS path
#'   only. Which reference produced \code{Q_pval}, the per-SNP \code{nu} it
#'   used, and how that \code{nu} was obtained.}
#' 
#' @export
param_gwas <- function(
  stage1,
  transform = c("custom_W", "two_normal"),
  transform_args = list(),
  se_mode = c("tau_cov", "diagonal", "plugin_cor", "dwls"),
  plugin_R = NULL,
  cov_data = NULL,
  tau_cov_estimator = c("ols", "gls"),
  gls_ridge = 1e-8,
  n_threads = 1L,
  return_cov = c("none", "upper", "full", "mashr"),
  q_method = c("saddlepoint", "chisq_naive", "liu", "satterthwaite",
               "hotelling"),
  nu = NULL,
  nu_source = NULL,
  unknown_nu = c("omit", "normal", "error"),
  n_blocks = NULL,
  n_eff = NULL,
  df_model = NULL,
  df_coef = NULL
) {
  transform <- match.arg(transform)
  se_mode   <- match.arg(se_mode)
  tau_cov_estimator <- match.arg(tau_cov_estimator)
  return_cov <- match.arg(return_cov)
  q_method <- match.arg(q_method)
  unknown_nu <- match.arg(unknown_nu)

  # The two estimators compute different statistics, so their references are
  # not interchangeable and the package refuses rather than silently
  # substituting one. See the derivation at the foot of R/tail-reference.R.
  gls_path <- (se_mode == "tau_cov" && tau_cov_estimator == "gls")
  if (q_method == "hotelling" && !gls_path) {
    stop("q_method = \"hotelling\" applies only to ",
         "se_mode = \"tau_cov\" with tau_cov_estimator = \"gls\". The OLS Q is ",
         "r' D^-1 r with D diagonal, whose null law is a weighted sum of ",
         "chi-squares with a per-SNP spectrum, not a Hotelling T^2. Use ",
         "q_method = \"saddlepoint\" there.", call. = FALSE)
  }
  if (gls_path && q_method %in% c("liu", "satterthwaite")) {
    stop("q_method = \"", q_method, "\" does not apply to ",
         "tau_cov_estimator = \"gls\". Those methods match the first two or ",
         "four moments of the eigenvalue spectrum of ",
         "D^-1/2 M Sigma M' D^-1/2, which is what the OLS Q needs. The GLS Q ",
         "weights by Sigma^-1, so that spectrum is degenerate -- every ",
         "eigenvalue is 1 -- and moment matching to it returns chi2(T - K), ",
         "which ignores that Sigma is ESTIMATED and is anti-conservative at ",
         "every finite nu (12% rejection at nominal 5% when nu = 500, 99% when ",
         "nu = 60). The exact reference is available in closed form: use ",
         "q_method = \"hotelling\" (or the default \"saddlepoint\", which ",
         "resolves to it on this path).", call. = FALSE)
  }

  taus    <- stage1$taus
  q_tau   <- stage1$q_tau
  Q_slope <- stage1$Q_slope
  SE_tau  <- stage1$SE_tau
  Tt <- length(taus)
  P  <- ncol(Q_slope)

  # (1) Construct W
  if (transform == "custom_W") {
    if (is.null(transform_args$W)) stop("Provide transform_args$W.")
    W <- as.matrix(transform_args$W)
    if (nrow(W) != Tt) stop("nrow(W) must equal length(taus).")
  } else {
    req <- c("p1","mu1","sd1","mu2","sd2","include_membership")
    missing_args <- setdiff(req, names(transform_args))
    if (length(missing_args))
      stop("Missing args for two_normal: ", paste(missing_args, collapse=", "))
    W <- make_weights_normal_mixture(
      taus = taus, q_tau = q_tau,
      p1 = transform_args$p1,
      mu1 = transform_args$mu1, sd1 = transform_args$sd1,
      mu2 = transform_args$mu2, sd2 = transform_args$sd2,
      include_membership = isTRUE(transform_args$include_membership)
    )
  }
  K <- ncol(W)
  A <- solve(crossprod(W), t(W))

  # (2) Map tau-slopes -> parameters
  params <- A %*% Q_slope
  rownames(params) <- colnames(W)
  colnames(params) <- colnames(Q_slope)

  # (3) SEs + Q statistics
  SE_params <- matrix(NA_real_, nrow = K, ncol = P,
                      dimnames = list(colnames(W), colnames(Q_slope)))
  Q_stat <- rep(NA_real_, P)

  Q_stat <- rep(NA_real_, P)

if (se_mode == "diagonal") {
  A2 <- A * A
  Vb <- SE_tau * SE_tau
  SE_params <- sqrt(A2 %*% Vb)
  
  # Q statistic (naive, assumes independence across taus)
  for (j in seq_len(P)) {
    resid_j <- Q_slope[, j] - W %*% params[, j]
    var_j   <- SE_tau[, j]^2
    Q_stat[j] <- sum((resid_j^2) / pmax(var_j, 1e-12))
  }

} else if (se_mode == "plugin_cor") {
  if (is.null(plugin_R)) {
    R_tau <- stats::cor(t(Q_slope), use = "pairwise.complete.obs")
  } else {
    R_tau <- as.matrix(plugin_R)
  }
  R_tau <- .nearPD_eig(R_tau, eps = 1e-8)

  for (j in seq_len(P)) {
    d <- SE_tau[, j]
    Dj <- diag(d, nrow = Tt, ncol = Tt)
    Sigma_j <- Dj %*% R_tau %*% Dj
    Vtheta  <- A %*% Sigma_j %*% t(A)
    Vtheta  <- .nearPD_eig(Vtheta, eps = 1e-12)
    SE_params[, j] <- sqrt(pmax(diag(Vtheta), 0))
    
    # Q statistic with full covariance
    resid_j <- Q_slope[, j] - W %*% params[, j]
    Q_stat[j] <- as.numeric(t(resid_j) %*% solve(Sigma_j, resid_j))
  }

} else if (se_mode == "dwls") {
  for (j in seq_len(P)) {
    var_j <- SE_tau[, j]^2
    Wgt   <- diag(1 / pmax(var_j, 1e-12), Tt)
    # GLS estimator (diagonal weight = DWLS)
    A_dwls <- solve(t(W) %*% Wgt %*% W, t(W) %*% Wgt)
    theta_j <- A_dwls %*% Q_slope[, j]
    params[, j] <- theta_j
    # SEs
    cov_theta <- solve(t(W) %*% Wgt %*% W)
    SE_params[, j] <- sqrt(diag(cov_theta))
    # Q statistic (DWLS)
    resid_j <- Q_slope[, j] - W %*% theta_j
    Q_stat[j] <- sum((resid_j^2) / var_j)
  }

} else if (se_mode == "tau_cov") {
  # Use per-SNP jackknife covariance from Stage 1
  if (is.null(cov_data)) {
    stop("se_mode = 'tau_cov' requires cov_data. Use se_mode = 'diagonal' explicitly only for fallback/debug runs.")
  } else if (tau_cov_estimator == "gls") {
    gls_fit <- compute_gls_param_gwas(
      beta_stage1 = Q_slope,
      cov_vec = cov_data$cov_vec,
      offsets = cov_data$offsets,
      W = W,
      return_cov = return_cov != "none",
      gls_ridge = gls_ridge,
      n_threads = n_threads,
      cov_scales = cov_data$cov_scales
    )
    params <- gls_fit$params
    rownames(params) <- colnames(W)
    colnames(params) <- colnames(Q_slope)
    SE_params <- gls_fit$se
    rownames(SE_params) <- colnames(W)
    colnames(SE_params) <- colnames(Q_slope)
    Q_stat <- gls_fit$Q
    if (!is.null(gls_fit$diag)) {
      d <- gls_fit$diag
      message(sprintf(
        paste0("GLS tau_cov diagnostics: %.0f SNPs, %.0f solved (%.4f%%); ",
               "ridge repair attempted %.0f (%.4f%%), rescued %.0f, ",
               "Sigma still not SPD %.0f; skipped: beta non-finite %.0f, ",
               "bad offset %.0f, cov non-finite %.0f, Sigma solve fail %.0f, ",
               "WtSinvW not SPD %.0f, theta fail %.0f, Vtheta fail %.0f, Q invalid %.0f"),
        d$n_snps, d$n_ok, 100 * d$n_ok / max(d$n_snps, 1),
        d$n_ridge_tried, 100 * d$n_ridge_tried / max(d$n_snps, 1),
        d$n_ridge_rescued, d$n_sigma_not_spd,
        d$n_beta_nonfinite, d$n_offset_bad, d$n_cov_nonfinite,
        d$n_solve_fail, d$n_xtx_not_spd, d$n_theta_fail,
        d$n_vtheta_fail, d$n_q_bad
      ))
      gls_diag <- d
    }

    # Reference for the GLS lack-of-fit statistic.
    #
    # Q_gls = (b - W theta)' S^-1 (b - W theta) is Hotelling's T^2 on the
    # m = T - K residual dimensions and nu degrees of freedom, EXACTLY:
    #
    #     Q_gls ~ (nu m / (nu - m + 1)) F_{m, nu - m + 1}
    #
    # not chi2(T - K), and not the OLS eigenvalue spectrum. The full
    # derivation and the simulation that confirms it are in
    # R/tail-reference.R. "saddlepoint" is the argument default, and on this
    # path the quantity it approximates is available exactly, so it resolves
    # to the exact reference rather than erroring.
    if (q_method != "chisq_naive") {
      if (q_method == "saddlepoint") {
        message("GLS path: q_method = \"saddlepoint\" resolves to the exact ",
                "Hotelling F reference for Q. The GLS quadratic form has a ",
                "degenerate eigenvalue spectrum, so there is nothing for a ",
                "saddlepoint to approximate. Recorded as ",
                "Q_ref_dist = \"hotelling_f\".")
      }
      q_ref <- .q_reference_gls(
        P = P,
        nu = nu, nu_source = nu_source,
        n_blocks = n_blocks, n_eff = n_eff,
        df_model = df_model, df_coef = df_coef,
        unknown_nu = unknown_nu,
        where = "this call to param_gwas(tau_cov_estimator = \"gls\")"
      )
      q_cal <- .q_pvalue_gls(Q = Q_stat, m = Tt - K, ref = q_ref)
    }
  } else {
    # Use C++ helper for fast computation
    se_mat <- compute_calibrated_se(
      cov_vec = cov_data$cov_vec,
      offsets = cov_data$offsets,
      A = A,
      K = Tt,
      cov_scales = cov_data$cov_scales
    )
    SE_params <- t(se_mat)
    rownames(SE_params) <- colnames(W)
    colnames(SE_params) <- colnames(Q_slope)
  }
  
  if (is.null(cov_data) || tau_cov_estimator == "ols") {
    # Q statistic (using diagonal approximation)
    for (j in seq_len(P)) {
      resid_j <- Q_slope[, j] - W %*% params[, j]
      var_j   <- SE_tau[, j]^2
      Q_stat[j] <- sum((resid_j^2) / pmax(var_j, 1e-12))
    }
    # The STATISTIC above is unchanged. What q_method changes is only the
    # reference distribution used to turn it into a p-value: the residuals
    # are strongly correlated across taus, so chisq(T - k) is wrong.
    if (q_method != "chisq_naive" && !is.null(cov_data)) {
      q_cal <- .q_pvalue_calibrated(
        Q = Q_stat, cov_data = cov_data, SE_tau = SE_tau, W = W,
        method = q_method, n_threads = n_threads
      )
    }
  }
}

  param_cov_upper <- NULL
  param_cov_full <- NULL
  mashr_inputs <- NULL
  if (return_cov != "none") {
    if (se_mode != "tau_cov" || is.null(cov_data)) {
      warning("return_cov requires se_mode = 'tau_cov' with cov_data; returning NULL")
    } else if (tau_cov_estimator == "gls") {
      param_cov_upper <- gls_fit$param_cov_upper
      if (!is.null(colnames(W))) {
        colnames(param_cov_upper) <- .param_cov_upper_names(colnames(W))
      }
      rownames(param_cov_upper) <- colnames(Q_slope)
      if (return_cov == "full") {
        param_cov_full <- .param_cov_upper_to_array(
          param_cov_upper,
          K = nrow(params),
          param_names = colnames(W)
        )
      } else if (return_cov == "mashr") {
        mashr_inputs <- mashr_inputs_from_cov(params, param_cov_upper, mode = "array")
      }
    } else {
      param_cov_upper <- compute_param_cov(
        cov_vec = cov_data$cov_vec,
        offsets = cov_data$offsets,
        A = A,
        K = Tt,
        cov_scales = cov_data$cov_scales
      )
      if (!is.null(colnames(W))) {
        colnames(param_cov_upper) <- .param_cov_upper_names(colnames(W))
      }
      rownames(param_cov_upper) <- colnames(Q_slope)
      if (return_cov == "full") {
        param_cov_full <- .param_cov_upper_to_array(
          param_cov_upper,
          K = nrow(params),
          param_names = colnames(W)
        )
      } else if (return_cov == "mashr") {
        mashr_inputs <- mashr_inputs_from_cov(params, param_cov_upper, mode = "array")
      }
    }
  }

  # Wald p-values for the parameters. Sigma is ESTIMATED (block jackknife), so
  # theta/SE is t_nu and not normal; see R/tail-reference.R for the identity and
  # docs/KNOWN_ISSUES.md for the 17x tail cost of getting it wrong.
  #
  # Q and Q_pval are deliberately untouched. The saddlepoint reference for Q
  # conditions on Sigma as known and inherits the same defect, but the fix there
  # is not a pt() substitution and is tracked separately.
  param_ref <- .param_reference(
    A = A, P = P,
    nu = nu, nu_source = nu_source,
    n_blocks = n_blocks, n_eff = n_eff,
    df_model = df_model, df_coef = df_coef,
    unknown_nu = unknown_nu,
    where = "this call to param_gwas()"
  )
  P_params <- .param_pvalues(params, SE_params, param_ref)

  out <- list(
    W         = W,
    A         = A,
    params    = params,
    SE_params = SE_params,
    P_params  = P_params,
    param_ref_dist   = param_ref$dist,
    param_ref_df     = param_ref$df,
    param_ref_source = param_ref$source,
    param_ref_c      = param_ref$c,
    Q         = Q_stat,
    df        = Tt - K,
    q_method  = q_method
  )
  if (exists("q_cal", inherits = FALSE) && !is.null(q_cal)) {
    out$Q_pval   <- q_cal$pval
    out$Q_df_eff <- q_cal$df_eff
    out$Q_scale  <- q_cal$scale
    # Provenance, GLS path only: the OLS output is left byte-for-byte as it was.
    if (!is.null(q_cal$ref_dist)) {
      out$Q_ref_dist   <- q_cal$ref_dist
      out$Q_ref_df     <- q_cal$nu
      out$Q_ref_source <- q_cal$ref_source
    }
  }
  if (return_cov == "upper") out$param_cov_upper <- param_cov_upper
  if (return_cov == "full") out$param_cov <- param_cov_full
  if (return_cov == "mashr") out$mashr <- mashr_inputs
  if (se_mode == "plugin_cor") out$R_tau <- R_tau
  if (se_mode == "tau_cov") out$tau_cov_estimator <- tau_cov_estimator
  if (exists("gls_diag", inherits = FALSE)) out$gls_diag <- gls_diag
  out
}

.validate_cov_ids <- function(stage1_dt, cov_ids_file) {
  cov_ids_dt <- data.table::fread(cov_ids_file)
  if (!("snp_id" %in% names(cov_ids_dt))) {
    stop("Covariance IDs file missing snp_id column: ", cov_ids_file)
  }
  if (nrow(cov_ids_dt) != nrow(stage1_dt)) {
    stop("Covariance IDs rows (", nrow(cov_ids_dt),
         ") do not match Stage 1 rows (", nrow(stage1_dt), ").")
  }
  key_cols <- intersect(
    c("snp_id", "chr", "bp", "effect_allele", "other_allele"),
    intersect(names(cov_ids_dt), names(stage1_dt))
  )
  if (!length(key_cols)) {
    key_cols <- "snp_id"
  }
  for (col in key_cols) {
    mism <- which(as.character(cov_ids_dt[[col]]) != as.character(stage1_dt[[col]]))
    if (length(mism) > 0) {
      i <- mism[1]
      stop("Covariance IDs mismatch at row ", i, " column ", col,
           ": cov_ids=", cov_ids_dt[[col]][i],
           " stage1=", stage1_dt[[col]][i])
    }
  }
  invisible(TRUE)
}

.filter_multiallelic_stage1 <- function(stage1_dt, drop_multiallelic = TRUE) {
  keep_idx <- seq_len(nrow(stage1_dt))
  if (!isTRUE(drop_multiallelic) || !("snp_id" %in% names(stage1_dt))) {
    return(list(dt = stage1_dt, keep_idx = keep_idx, dropped = data.table::data.table()))
  }

  dup_snp_mask <- duplicated(stage1_dt$snp_id) | duplicated(stage1_dt$snp_id, fromLast = TRUE)
  dup_pos_mask <- rep(FALSE, nrow(stage1_dt))
  if (all(c("chr", "bp") %in% names(stage1_dt))) {
    pos_dt <- stage1_dt[, c("chr", "bp"), with = FALSE]
    dup_pos_mask <- duplicated(pos_dt) | duplicated(pos_dt, fromLast = TRUE)
  }

  drop_mask <- dup_snp_mask | dup_pos_mask
  if (!any(drop_mask)) {
    return(list(dt = stage1_dt, keep_idx = keep_idx, dropped = data.table::data.table()))
  }

  dup_ids <- unique(stage1_dt$snp_id[dup_snp_mask])
  dropped_cols <- intersect(
    c("snp_id", "chr", "bp", "effect_allele", "other_allele"),
    names(stage1_dt)
  )
  dropped <- stage1_dt[drop_mask, dropped_cols, with = FALSE]
  keep_idx <- which(!drop_mask)
  message(
    "Dropping ", sum(drop_mask), " Stage 1 rows from ",
    length(dup_ids), " duplicated rsIDs and ",
    sum(dup_pos_mask), " duplicated-position rows before Stage 2."
  )
  if (!length(keep_idx)) {
    stop("All Stage 1 rows were duplicated/multi-allelic; no rows remain for Stage 2.")
  }
  list(dt = stage1_dt[keep_idx, ], keep_idx = keep_idx, dropped = dropped)
}

.cov_scale_file <- function(cov_file) {
  sub("\\.cov\\.gz$", ".cov.scale.gz", cov_file)
}

.read_tau_cov_data <- function(cov_file, n_taus, n_snps) {
  n_cov_elements <- as.integer(n_taus * (n_taus + 1) / 2)
  scale_file <- .cov_scale_file(cov_file)

  if (file.exists(scale_file)) {
    con_cov <- gzfile(cov_file, "rb")
    on.exit(close(con_cov), add = TRUE)
    cov_q <- readBin(
      con_cov,
      what = "integer",
      n = n_cov_elements * n_snps,
      size = 1,
      signed = TRUE,
      endian = "little"
    )
    n_expect <- n_cov_elements * n_snps
    if (length(cov_q) < n_expect) {
      warning(
        "Int8 covariance shorter than expected in ", cov_file, ": expected ",
        n_expect, " bytes, got ", length(cov_q), ". Padding with NA.",
        call. = FALSE
      )
      cov_q <- c(cov_q, rep(NA_integer_, n_expect - length(cov_q)))
    } else if (length(cov_q) > n_expect) {
      warning(
        "Int8 covariance longer than expected in ", cov_file, ": expected ",
        n_expect, " bytes, got ", length(cov_q), ". Truncating extra bytes.",
        call. = FALSE
      )
      cov_q <- cov_q[seq_len(n_expect)]
    }

    con_scale <- gzfile(scale_file, "rb")
    on.exit(close(con_scale), add = TRUE)
    scales <- readBin(
      con_scale,
      what = "numeric",
      n = n_snps,
      size = 4,
      endian = "little"
    )
    if (length(scales) < n_snps) {
      warning(
        "Scale file shorter than expected in ", scale_file, ": expected ",
        n_snps, " records, got ", length(scales), ". Padding with NA.",
        call. = FALSE
      )
      scales <- c(scales, rep(NA_real_, n_snps - length(scales)))
    } else if (length(scales) > n_snps) {
      warning(
        "Scale file longer than expected in ", scale_file, ": expected ",
        n_snps, " records, got ", length(scales), ". Truncating extra records.",
        call. = FALSE
      )
      scales <- scales[seq_len(n_snps)]
    }

    cov_vec <- as.numeric(cov_q) * rep(scales, each = n_cov_elements)
    offsets <- seq(0, by = n_cov_elements * 4, length.out = n_snps)
    return(list(cov_vec = cov_vec, offsets = offsets))
  }

  con <- gzfile(cov_file, "rb")
  on.exit(close(con), add = TRUE)
  cov_vec <- readBin(
    con,
    what = "numeric",
    n = n_cov_elements * n_snps,
    size = 4,
    endian = "little"
  )
  n_expect <- n_cov_elements * n_snps
  if (length(cov_vec) < n_expect) {
    warning(
      "Float32 covariance shorter than expected in ", cov_file, ": expected ",
      n_expect, " values, got ", length(cov_vec), ". Padding with NA.",
      call. = FALSE
    )
    cov_vec <- c(cov_vec, rep(NA_real_, n_expect - length(cov_vec)))
  } else if (length(cov_vec) > n_expect) {
    warning(
      "Float32 covariance longer than expected in ", cov_file, ": expected ",
      n_expect, " values, got ", length(cov_vec), ". Truncating extra values.",
      call. = FALSE
    )
    cov_vec <- cov_vec[seq_len(n_expect)]
  }
  offsets <- seq(0, by = n_cov_elements * 4, length.out = n_snps)
  list(cov_vec = cov_vec, offsets = offsets)
}

## Stage 1 writes the pooled sidecar beside the stage1 TSV. Inferring that
## conventional path keeps the pooled default usable without making callers
## repeat a filename; an absent sidecar still fails explicitly rather than
## silently falling back to the per-SNP covariance.
.infer_sigma0_file <- function(stage1_file) {
  path <- as.character(stage1_file)[1L]
  candidate <- sub("\\.stage1\\.tsv\\.gz$", ".sigma0.tsv.gz", path)
  if (identical(candidate, path)) {
    candidate <- sub("\\.stage1\\.tsv$", ".sigma0.tsv.gz", path)
  }
  if (identical(candidate, path)) candidate <- paste0(path, ".sigma0.tsv.gz")
  if (file.exists(candidate)) candidate else NULL
}

#' Parametric GWAS from Stage 1 files
#'
#' Convenience wrapper for \code{param_gwas()} that reads Stage 1 output files
#' directly, including covariance data for \code{se_mode = "tau_cov"}.
#'
#' @param stage1_file Path to Stage 1 TSV file (.stage1.tsv.gz)
#' @param W Weight matrix (T x K)
#' @param se_mode SE computation mode (see \code{param_gwas})
#' @param cov_file Path to covariance file (.cov.gz) for tau_cov mode. If a
#'   sidecar file with suffix .cov.ids.tsv.gz exists, it is used to validate
#'   row alignment against the Stage 1 file. If .cov.scale.gz is present, the
#'   covariance payload is read as int8 and dequantized using per-SNP scales.
#' @param unsafe_skip_cov_ids If TRUE, allow tau_cov to proceed without a
#'   covariate IDs sidecar (not recommended).
#' @param rtau_file Path to R_tau correlation matrix (.Rtau.tsv) for plugin_cor mode
#' @param return_cov Optional output of parameter covariance derived from Stage 1
#'   tau covariance. One of \code{"none"} (default), \code{"upper"}, \code{"full"},
#'   or \code{"mashr"}.
#' @param param_cov_file Optional output file for parameter covariance upper triangles
#'   (binary float32, same layout as Stage 1 .cov.gz).
#' @param cov_source Which Stage 1 covariance to use. \code{"pooled"}
#'   (default) reads the trait-level
#'   \code{Sigma_0} sidecar from \code{sigma0_file} and forms
#'   \code{Sigma_v = Sigma_0 / d_resid} per SNP, which reproduces an exact
#'   per-individual HC0 sandwich to a median diagonal ratio of 1.00000 and a
#'   median cross-tau correlation deviation of 0.0053 on real chr22 data, against
#'   the jackknife's 0.725-1.307 5-95\% spread on the variance. The pooled path
#'   refuses rather than proceeds if the recorded covariate set, tau grid or
#'   sample count does not match the scan. Set \code{cov_source = "jackknife"}
#'   explicitly to read the per-SNP delete-one-block covariance from
#'   \code{cov_file}; that path remains fully supported for analyses where a
#'   lead locus may violate the pooled working-independence assumption.
#' @param sigma0_file Path to the Stage 1 \code{.sigma0.tsv.gz} sidecar, required
#'   when \code{cov_source = "pooled"}. If omitted, the wrapper looks for a
#'   sibling sidecar by replacing \code{.stage1.tsv.gz} with
#'   \code{.sigma0.tsv.gz}; pass this argument explicitly when the sidecar has
#'   another name. Its \code{.sigma0.meta.json} and the scan's
#'   \code{.cov.meta.json} must both be present: the second is what the first
#'   is checked against.
#' @param require_covariates Optional character vector of covariate column names
#'   the \code{Sigma_0} sidecar must have been built on. Supply it when you know
#'   what the scan used and want the mismatch caught here rather than trusted.
#' @param require_n_samples Optional analysed sample count the \code{Sigma_0}
#'   sidecar must match.
#' @param drop_multiallelic If TRUE, drop duplicated rsIDs or positions before
#'   fitting.
#' @param tau_cov_estimator Estimator for tau-covariance mode: \code{"ols"}
#'   (default) or \code{"gls"}.
#' @param gls_ridge Relative diagonal ridge used only when GLS inversion needs
#'   repair.
#' @param n_threads Number of threads for covariance calculations.
#' @param q_method Reference method for the lack-of-fit statistic. An explicitly
#'   naive or unsupported reference omits \code{Q_pval}.
#' @param cov_meta_file Path to the Stage 1 \code{.cov.meta.json} sidecar that
#'   records the jackknife block count. Defaults to the \code{.cov.meta.json}
#'   next to \code{cov_file}. This is where the degrees of freedom for the
#'   \code{t} reference come from.
#' @param nu Override for the degrees of freedom, applied to every SNP and
#'   every parameter. Normally leave this \code{NULL}: the per-SNP law is
#'   better, and a constant is not defensible at \code{n_blocks >= 500}.
#' @param df_coef Optional override for the Satterthwaite coefficient \code{c}:
#'   one number, or a vector named by parameter. The shipped values (0.88 for a
#'   location-type contrast, 2.62 for a scale-type one) were calibrated for
#'   \code{W = [1, q_tau - mean(q_tau)]}; a very different \code{W} needs its
#'   own.
#' @param unknown_nu What to do when \code{nu} cannot be determined:
#'   \code{"omit"} (default) emits no \code{<param>_p} columns and warns,
#'   \code{"normal"} emits normal-reference p-values and warns, \code{"error"}
#'   stops. Whatever happens, the returned table carries \code{ref_dist} and
#'   \code{ref_df_source} columns, and a \code{<param>_df} column per fitted
#'   parameter, so a reader of the results can tell exactly which reference
#'   distribution produced each p-value.
#'
#' @return A data.table with SNP info and parameter estimates
#' @export
param_gwas_from_file <- function(
  stage1_file,
  W,
  se_mode = c("tau_cov", "diagonal", "plugin_cor", "dwls"),
  cov_file = NULL,
  cov_source = c("pooled", "jackknife"),
  sigma0_file = NULL,
  require_covariates = NULL,
  require_n_samples = NULL,
  rtau_file = NULL,
  return_cov = c("none", "upper", "full", "mashr"),
  param_cov_file = NULL,
  unsafe_skip_cov_ids = FALSE,
  drop_multiallelic = TRUE,
  tau_cov_estimator = c("ols", "gls"),
  gls_ridge = 1e-8,
  n_threads = 1L,
  q_method = c("saddlepoint", "chisq_naive", "liu", "satterthwaite",
               "hotelling"),
  cov_meta_file = NULL,
  nu = NULL,
  unknown_nu = c("omit", "normal", "error"),
  df_coef = NULL
) {
  se_mode <- match.arg(se_mode)
  cov_source <- match.arg(cov_source)
  return_cov <- match.arg(return_cov)
  tau_cov_estimator <- match.arg(tau_cov_estimator)
  q_method <- match.arg(q_method)
  unknown_nu <- match.arg(unknown_nu)
  
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table package required for param_gwas_from_file")
  }
  
  # Read Stage 1
  message("Reading Stage 1 file: ", stage1_file)
  stage1_dt <- data.table::fread(stage1_file)
  data.table::setDT(stage1_dt)
  stage1_n_raw <- nrow(stage1_dt)
  stage1_dt_raw <- data.table::copy(stage1_dt)
  stage1_filter <- .filter_multiallelic_stage1(stage1_dt, drop_multiallelic)
  stage1_dt <- stage1_filter$dt
  
  # Extract beta/SE columns
  beta_cols <- grep("^beta_tau_", names(stage1_dt), value = TRUE)
  se_cols <- grep("^se_tau_", names(stage1_dt), value = TRUE)
  
  taus <- as.numeric(sub("^beta_tau_", "", beta_cols))
  
  Q_slope <- t(as.matrix(stage1_dt[, beta_cols, with = FALSE]))
  SE_tau <- t(as.matrix(stage1_dt[, se_cols, with = FALSE]))
  
  # Build stage1 list
  stage1 <- list(
    taus = taus,
    q_tau = rep(0, length(taus)),  # Not used for param_gwas
    Q_slope = Q_slope,
    SE_tau = SE_tau
  )
  
  # Load covariance if needed
  cov_data <- NULL
  pooled <- NULL
  if (se_mode == "tau_cov" && cov_source == "pooled") {
    # The pooled path deliberately does NOT read cov_file: the whole point is
    # that the 229-289 bytes/SNP payload does not have to exist. It does read
    # the scan's .cov.meta.json, because that is what the sidecar is checked
    # against, and it reads d_resid out of the Stage 1 TSV for the per-SNP scale.
    if (is.null(sigma0_file)) {
      sigma0_file <- .infer_sigma0_file(stage1_file)
      if (is.null(sigma0_file)) {
        stop("cov_source = 'pooled' requires the Stage 1 .sigma0.tsv.gz ",
             "sidecar. Pass sigma0_file explicitly or place the sidecar next ",
             "to the .stage1.tsv.gz file; no silent jackknife fallback is used.",
             call. = FALSE)
      }
    }
    message("Loading pooled covariance: ", sigma0_file)
    d_resid_col <- if (N_EFF_COLUMN %in% names(stage1_dt)) {
      as.numeric(stage1_dt[[N_EFF_COLUMN]])
    } else NULL
    pooled <- .pooled_tau_cov(
      sigma0_file = sigma0_file,
      chr = stage1_dt$chr,
      d_resid = d_resid_col,
      taus = taus,
      cov_meta_file = cov_meta_file,
      require_covariates = require_covariates,
      require_n_samples = require_n_samples
    )
    cov_data <- pooled$cov_data
    message("Pooled covariance: Sigma_0 on chromosome(s) ",
            paste(pooled$meta$chromosomes, collapse = ", "), ", ",
            pooled$meta$n_samples, " samples, ",
            length(pooled$meta$covariates), " covariates (",
            paste(pooled$meta$covariates, collapse = ", "), ").")
  }
  if (se_mode == "tau_cov" && cov_source == "jackknife" && is.null(cov_file)) {
    stop("se_mode = 'tau_cov' requires cov_file. Use se_mode = 'diagonal' explicitly only for fallback/debug runs.")
  }
  if (se_mode == "tau_cov" && cov_source == "jackknife" && !is.null(cov_file)) {
    message("Loading covariance file: ", cov_file)

    # Validate row alignment if cov IDs sidecar exists
    cov_ids_file <- sub("\\.cov\\.gz$", ".cov.ids.tsv.gz", cov_file)
    if (!file.exists(cov_ids_file)) {
      if (!isTRUE(unsafe_skip_cov_ids)) {
        stop("Missing covariance IDs file: ", cov_ids_file,
             ". This file is required to validate row alignment. ",
             "Set unsafe_skip_cov_ids=TRUE to proceed without validation.")
      }
      warning("Missing covariance IDs file: ", cov_ids_file,
              ". Proceeding without validation because unsafe_skip_cov_ids=TRUE.",
              call. = FALSE)
    } else {
      message("Validating covariance SNP IDs: ", cov_ids_file)
      .validate_cov_ids(stage1_dt_raw, cov_ids_file)
    }
    
    cov_data <- .read_tau_cov_data(
      cov_file = cov_file,
      n_taus = length(taus),
      n_snps = stage1_n_raw
    )
    cov_data$offsets <- cov_data$offsets[stage1_filter$keep_idx]
  }
  
  # Load R_tau if needed
  plugin_R <- NULL
  if (se_mode == "plugin_cor" && !is.null(rtau_file)) {
    message("Loading R_tau file: ", rtau_file)
    rtau_dt <- data.table::fread(rtau_file)
    plugin_R <- as.matrix(rtau_dt)
  }
  
  # Degrees of freedom for the t reference on theta/SE.
  #
  #   1. an explicit nu from the caller wins;
  #   2. otherwise the block count from .cov.meta.json plus the per-SNP n_eff
  #      from the d_resid column give one nu per SNP per parameter;
  #   3. otherwise the block count alone gives a flat fallback (warned);
  #   4. otherwise a constant nu recorded in the metadata (warned);
  #   5. otherwise unknown_nu decides.
  #
  # nu = n_blocks - 1 is NOT among these. It is measurably wrong: 801.8 against
  # 999 at B = 1000, N = 50,000.
  nu_value <- NULL
  nu_prov <- NULL
  n_blocks_value <- NULL
  df_model_value <- .resolve_df_model(NULL)
  if (!is.null(nu)) {
    nu_value <- as.numeric(nu)[1]
    nu_prov <- "user_supplied"
  } else if (!is.null(pooled)) {
    # The jackknife's t reference exists because Sigma is estimated from B
    # replicates, giving nu around 800 at B = 1000 and 62 at B = 65. Sigma_0 is
    # estimated too, but from n individuals, and the same matrix serves every
    # SNP: its Satterthwaite degrees of freedom are O(n / (1 + kurt/2)), which
    # on the production chr22 inputs is tens of thousands. t on that many df is
    # the normal distribution to well past the twelfth significant figure at
    # p = 5e-8, so the per-SNP nu law simply does not apply here -- but a finite
    # nu is still recorded rather than asserting normality, because the package
    # does not assume a reference it cannot name. Pass `nu` to pin a smaller one.
    nu_value <- as.numeric(pooled$nu)[1]
    nu_prov <- "pooled_sigma0_n_samples"
    message("Pooled covariance: reference is t_", format(nu_value, scientific = FALSE),
            " (numerically normal); the per-SNP jackknife df law does not apply.")
  } else if (!is.null(cov_file)) {
    if (is.null(cov_meta_file)) cov_meta_file <- .cov_meta_file(cov_file)
    cov_meta <- .read_cov_meta(cov_meta_file)
    if (is.null(cov_meta)) {
      message("No Stage 1 covariance metadata at ", cov_meta_file,
              " (covariance predates the block-count field).")
    }
    df_model_value <- .resolve_df_model(cov_meta)
    b <- .resolve_block_count(cov_meta)
    if (is.finite(b)) {
      n_blocks_value <- b
      message("Stage 1 jackknife blocks: ", b,
              "; effective df computed per SNP as 1/nu = 1/(B-1) + c/n_eff.")
    } else {
      recorded <- .resolve_nu(cov_meta)
      if (is.finite(recorded$nu)) {
        nu_value <- recorded$nu
        nu_prov <- recorded$nu_source
        warning("Using the constant nu = ", nu_value, " recorded in ",
                cov_meta_file, " (", nu_prov, "). A constant nu is not ",
                "defensible at n_blocks >= 500: the effective df varies by SNP ",
                "and by parameter. Re-run Stage 1 for the per-SNP form.",
                call. = FALSE)
      }
    }
  }

  # Per-SNP effective sample size, if Stage 1 emitted it.
  n_eff_value <- NULL
  n_eff_col <- df_model_value$n_eff_column
  if (n_eff_col %in% names(stage1_dt)) {
    n_eff_value <- as.numeric(stage1_dt[[n_eff_col]])
  } else if (!is.null(n_blocks_value)) {
    message("Stage 1 file has no ", n_eff_col,
            " column; the effective df will fall back to a flat value.")
  }

  # Run param_gwas
  result <- param_gwas(
    stage1 = stage1,
    transform = "custom_W",
    transform_args = list(W = W),
    se_mode = se_mode,
    plugin_R = plugin_R,
    cov_data = cov_data,
    tau_cov_estimator = tau_cov_estimator,
    gls_ridge = gls_ridge,
    n_threads = n_threads,
    return_cov = return_cov,
    q_method = q_method,
    nu = nu_value,
    nu_source = nu_prov,
    unknown_nu = unknown_nu,
    n_blocks = n_blocks_value,
    n_eff = n_eff_value,
    df_model = df_model_value,
    df_coef = df_coef
  )
  
  # Build output data.table
  out_dt <- stage1_dt[, c("snp_id", "chr", "bp", "effect_allele", "other_allele"), with = FALSE]
  
  # Add parameters, SEs and p-values. The p-value is emitted by the package on
  # purpose: downstream code used to derive it as 2 * pnorm(-|beta/se|), which is
  # the wrong reference and is exactly what this change exists to stop.
  for (i in seq_len(nrow(result$params))) {
    param_name <- rownames(result$params)[i]
    out_dt[[param_name]] <- result$params[i, ]
    out_dt[[paste0(param_name, "_se")]] <- result$SE_params[i, ]
    if (!is.null(result$P_params)) {
      out_dt[[paste0(param_name, "_p")]] <- result$P_params[i, ]
      # nu varies by SNP and by parameter, so it is a column per parameter and
      # not a single number for the run.
      out_dt[[paste0(param_name, "_df")]] <- result$param_ref_df[i, ]
    }
  }

  # Which covariance produced these standard errors, so a results file can be
  # told apart from one built the other way without re-running anything.
  out_dt$cov_source <- if (se_mode == "tau_cov") cov_source else NA_character_

  # Which reference distribution produced those p-values, and how its degrees of
  # freedom were arrived at, recorded in the output itself so a reader of a
  # results file can tell without re-running anything.
  out_dt$ref_dist <- result$param_ref_dist
  out_dt$ref_df_source <- if (is.null(result$param_ref_source)) {
    NA_character_
  } else {
    result$param_ref_source
  }
  
  out_dt$Q <- result$Q
  if (!is.null(result$Q_pval)) {
    # Calibrated reference distribution (q_method != "chisq_naive").
    out_dt$Q_pval   <- result$Q_pval
    out_dt$Q_df_eff <- result$Q_df_eff
    out_dt$Q_scale  <- result$Q_scale
    if (!is.null(result$Q_ref_dist)) {
      out_dt$Q_ref_dist   <- result$Q_ref_dist
      out_dt$Q_ref_df     <- result$Q_ref_df
      out_dt$Q_ref_source <- result$Q_ref_source
    }
  }
  ## Do not manufacture a chi-square p-value when param_gwas() deliberately
  ## omitted one (for example q_method = "chisq_naive" or unknown GLS df).
  ## The old fallback made an uncalibrated result indistinguishable from a
  ## calibrated one and was catastrophically anti-conservative for estimated
  ## covariance matrices.

  if (!is.null(param_cov_file) && !is.null(result$param_cov_upper)) {
    .write_param_cov_upper(result$param_cov_upper, param_cov_file)
  }

  if (return_cov == "none" && is.null(param_cov_file)) {
    return(out_dt)
  }
  list(results = out_dt, param_cov_upper = result$param_cov_upper, param_cov = result$param_cov, mashr = result$mashr)
}

#' Convenience wrapper for param_gwas_multi() that reads Stage 1 and covariance
#'
#' @param stage1_file Path to Stage 1 results (.stage1.tsv.gz)
#' @param W_list List of weight matrices (each T x K).
#' @param N Sample size used for BIC.
#' @param se_mode One of \code{"tau_cov"} or \code{"dwls"} (or \code{"diagonal"}).
#' @param cov_file Path to covariance file (.cov.gz) for tau_cov mode. If
#'   .cov.scale.gz is present, the covariance payload is read as int8 and
#'   dequantized using per-SNP scales.
#' @param cov_source \code{"pooled"} (default) or \code{"jackknife"}; see
#'   \code{param_gwas_from_file}. Set it explicitly to \code{"jackknife"}
#'   when supplying a per-SNP covariance file.
#' @param sigma0_file Path to the Stage 1 \code{.sigma0.tsv.gz}, required when
#'   \code{cov_source = "pooled"}. If omitted, the conventional sibling sidecar
#'   next to \code{stage1_file} is used when present.
#' @param unsafe_skip_cov_ids If TRUE, allow tau_cov to proceed without a
#'   covariate IDs sidecar (not recommended).
#' @param drop_multiallelic If TRUE, drop duplicated rsIDs or positions before
#'   fitting.
#' @param n_threads Number of OpenMP threads.
#'
#' @return Matrix with AIC/BIC columns per model.
#'
#' @export
param_gwas_multi_from_file <- function(
  stage1_file,
  W_list,
  N,
  se_mode = c("tau_cov", "dwls", "diagonal"),
  cov_file = NULL,
  cov_source = c("pooled", "jackknife"),
  sigma0_file = NULL,
  unsafe_skip_cov_ids = FALSE,
  drop_multiallelic = TRUE,
  n_threads = 1L
) {
  se_mode <- match.arg(se_mode)
  cov_source <- match.arg(cov_source)
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table package required for param_gwas_multi_from_file")
  }

  message("Reading Stage 1 file: ", stage1_file)
  stage1_dt <- data.table::fread(stage1_file)
  data.table::setDT(stage1_dt)
  # Covariance payloads and their IDs are indexed against the unfiltered Stage 1
  # rows. Keep that table for validation before dropping duplicated/multiallelic
  # records; offsets are subset below using the same keep_idx.
  stage1_dt_raw <- data.table::copy(stage1_dt)
  stage1_n_raw <- nrow(stage1_dt)
  stage1_filter <- .filter_multiallelic_stage1(stage1_dt, drop_multiallelic)
  stage1_dt <- stage1_filter$dt

  beta_cols <- grep("^beta_tau_", names(stage1_dt), value = TRUE)
  se_cols <- grep("^se_tau_", names(stage1_dt), value = TRUE)
  taus <- as.numeric(sub("^beta_tau_", "", beta_cols))

  Q_slope <- t(as.matrix(stage1_dt[, beta_cols, with = FALSE]))
  SE_tau <- t(as.matrix(stage1_dt[, se_cols, with = FALSE]))

  stage1 <- list(
    taus = taus,
    q_tau = rep(0, length(taus)),
    Q_slope = Q_slope,
    SE_tau = SE_tau
  )

  cov_data <- NULL
  if (se_mode == "tau_cov" && cov_source == "pooled") {
    if (is.null(sigma0_file)) {
      sigma0_file <- .infer_sigma0_file(stage1_file)
      if (is.null(sigma0_file)) {
        stop("cov_source = 'pooled' requires the Stage 1 .sigma0.tsv.gz ",
             "sidecar. Pass sigma0_file explicitly or place the sidecar next ",
             "to the .stage1.tsv.gz file; no silent jackknife fallback is used.",
             call. = FALSE)
      }
    }
    cov_data <- .pooled_tau_cov(
      sigma0_file = sigma0_file,
      chr = stage1_dt$chr,
      d_resid = if (N_EFF_COLUMN %in% names(stage1_dt)) as.numeric(stage1_dt[[N_EFF_COLUMN]]) else NULL,
      taus = taus
    )$cov_data
  } else if (se_mode == "tau_cov") {
    if (is.null(cov_file)) {
      stop("se_mode = 'tau_cov' requires cov_file.")
    }
    message("Loading covariance file: ", cov_file)

    cov_ids_file <- sub("\\.cov\\.gz$", ".cov.ids.tsv.gz", cov_file)
    if (!file.exists(cov_ids_file)) {
      if (!isTRUE(unsafe_skip_cov_ids)) {
        stop("Missing covariance IDs file: ", cov_ids_file,
             ". This file is required to validate row alignment. ",
             "Set unsafe_skip_cov_ids=TRUE to proceed without validation.")
      }
      warning("Missing covariance IDs file: ", cov_ids_file,
              ". Proceeding without validation because unsafe_skip_cov_ids=TRUE.",
              call. = FALSE)
    } else {
      message("Validating covariance SNP IDs: ", cov_ids_file)
      .validate_cov_ids(stage1_dt_raw, cov_ids_file)
    }

    cov_data <- .read_tau_cov_data(
      cov_file = cov_file,
      n_taus = length(taus),
      n_snps = stage1_n_raw
    )
    cov_data$offsets <- cov_data$offsets[stage1_filter$keep_idx]
  }

  param_gwas_multi(
    stage1 = stage1,
    W_list = W_list,
    N = N,
    n_threads = n_threads,
    se_mode = se_mode,
    cov_data = cov_data
  )
}

#' Calibrated p-value for the Stage 2 lack-of-fit statistic Q
#'
#' The Q statistic sums squared residuals divided by their own per-tau
#' variances. The 45 residuals are strongly correlated (same individuals,
#' overlapping influence functions), so Q is NOT chi-square on \code{T - k}
#' degrees of freedom: it is a quadratic form in a correlated normal vector,
#'
#'   Q = r' D^-1 r,  r = M b,  M = I - W (W'W)^-1 W',  Var(r) = M Sigma M',
#'
#' whose exact null distribution is \code{sum_i lambda_i chi2_1} with
#' \code{lambda} the eigenvalues of \code{D^-1/2 M Sigma t(M) D^-1/2}.
#'
#' Critically this needs \code{Sigma} only and never \code{Sigma^-1}, so it is
#' valid on the stored int8 covariance, which is accurate in quadratic forms
#' but useless for inversion.
#'
#' @param Q Numeric vector of observed Q statistics.
#' @param cov_data List with \code{cov_vec} and \code{offsets}.
#' @param SE_tau Numeric matrix (T x n_snps) of Stage 1 per-tau SEs.
#' @param W Weight matrix (T x k).
#' @param method One of "satterthwaite", "liu", "saddlepoint".
#' @param n_threads OpenMP threads.
#' @return List with \code{pval}, \code{df_eff}, \code{scale}, \code{n_eig}.
#' @keywords internal
.q_pvalue_calibrated <- function(Q, cov_data, SE_tau, W,
                                 method = c("saddlepoint", "liu", "satterthwaite"),
                                 n_threads = 1L) {
  method <- match.arg(method)
  Tt <- nrow(SE_tau)
  M <- diag(Tt) - W %*% solve(crossprod(W), t(W))
  Dvar <- SE_tau * SE_tau

  cal <- compute_q_calibration(
    cov_vec = cov_data$cov_vec,
    offsets = cov_data$offsets,
    Qstat = as.numeric(Q),
    Dvar = Dvar,
    M = M,
    K = as.integer(Tt),
    n_threads = as.integer(n_threads),
    cov_scales = cov_data$cov_scales
  )

  c1 <- cal[, 1]; c2 <- cal[, 2]; c3 <- cal[, 3]; c4 <- cal[, 4]
  p_saddle <- cal[, 5]; n_eig <- cal[, 6]

  ok <- is.finite(c1) & is.finite(c2) & c1 > 0 & c2 > 0
  n <- length(Q)
  pval <- rep(NA_real_, n)

  # Two-moment Satterthwaite: Q ~ scale * chisq(df_eff).
  # c1 = E[Q] = tr(B) and c2 = Var[Q]/2 = tr(B^2) are both exact for a
  # two-sided eigenvalue spectrum (c2 is a sum of squares, so always > 0).
  # The scale/df match needs c1 > 0, which `ok` already requires.
  scale_a <- rep(NA_real_, n); df_eff <- rep(NA_real_, n)
  scale_a[ok] <- c2[ok] / c1[ok]              # Var/(2E) = 2c2/(2c1)
  df_eff[ok]  <- c1[ok]^2 / c2[ok]            # 2E^2/Var
  p_satt <- rep(NA_real_, n)
  p_satt[ok] <- stats::pchisq(Q[ok] / scale_a[ok], df = df_eff[ok],
                              lower.tail = FALSE)

  # Four-moment (Liu, Tang & Zhang 2009): matches skewness and kurtosis too,
  # which the two-moment version does not, and which dominates the far tail.
  p_liu <- rep(NA_real_, n)
  okl <- ok & is.finite(c3) & is.finite(c4) & c4 > 0
  if (any(okl)) {
    s1 <- c3[okl] / c2[okl]^1.5
    s2 <- c4[okl] / c2[okl]^2
    muQ <- c1[okl]
    sigmaQ <- sqrt(2 * c2[okl])
    disc <- s1^2 - s2
    # The moments c1..c4 are exact for a two-sided spectrum -- the cumulants
    # of sum lam_i chi2_1 are 2^(k-1) (k-1)! c_k whatever the signs. What is
    # sign-sensitive is the FIT: both branches below map the moments onto a
    # (non-central) chi-square, whose skewness is necessarily positive, so
    # they are only meaningful when s1 = c3/c2^1.5 > 0. With enough negative
    # weight c3 can go non-positive, and the matching then silently returns a
    # number from the wrong shape family. Refuse it instead.
    use_nc <- disc > 0 & s1 > 0
    a <- numeric(length(s1)); delta <- numeric(length(s1)); l <- numeric(length(s1))
    a[use_nc] <- 1 / (s1[use_nc] - sqrt(disc[use_nc]))
    delta[use_nc] <- s1[use_nc] * a[use_nc]^3 - a[use_nc]^2
    l[use_nc] <- a[use_nc]^2 - 2 * delta[use_nc]
    a[!use_nc] <- 1 / sqrt(s2[!use_nc])
    delta[!use_nc] <- 0
    l[!use_nc] <- 1 / s2[!use_nc]
    muX <- l + delta
    sigmaX <- sqrt(2) * a
    Qn <- (Q[okl] - muQ) / sigmaQ * sigmaX + muX
    good <- is.finite(Qn) & is.finite(l) & l > 0 & is.finite(delta) & delta >= 0 &
            is.finite(s1) & s1 > 0
    v <- rep(NA_real_, length(Qn))
    v[good] <- stats::pchisq(Qn[good], df = l[good], ncp = delta[good],
                             lower.tail = FALSE)
    p_liu[okl] <- v
  }

  if (method == "satterthwaite") {
    pval <- p_satt
  } else if (method == "liu") {
    pval <- p_liu
  } else {
    # Saddlepoint where it applies (upper tail only); Liu elsewhere, which is
    # the region where p is large and the approximation hardly matters.
    pval <- p_liu
    use <- is.finite(p_saddle) & p_saddle >= 0
    pval[use] <- p_saddle[use]
    still <- !is.finite(pval)
    pval[still] <- p_satt[still]
  }

  list(pval = pval, df_eff = df_eff, scale = scale_a, n_eig = n_eig,
       p_satterthwaite = p_satt, p_liu = p_liu, p_saddlepoint = p_saddle)
}
