#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericMatrix fit_multi_models(const arma::mat& beta_stage1,
                               const arma::mat& se_stage1,
                               const List& W_list,
                               double N,
                               int n_threads = 1) {
  const uword P = beta_stage1.n_rows;
  const uword T = beta_stage1.n_cols;
  const int n_models = W_list.size();

  if (se_stage1.n_rows != P || se_stage1.n_cols != T) {
    stop("se_stage1 must have the same shape as beta_stage1.");
  }
  if (n_models == 0) {
    return NumericMatrix(P, 0);
  }

  std::vector<arma::mat> W_mats;
  std::vector<uword> k_sizes;
  W_mats.reserve(n_models);
  k_sizes.reserve(n_models);

  for (int m = 0; m < n_models; m++) {
    arma::mat W = as<arma::mat>(W_list[m]);
    if (W.n_rows != T) {
      stop("Each W matrix must have T rows.");
    }
    W_mats.push_back(W);
    k_sizes.push_back(W.n_cols);
  }

  NumericMatrix out(P, 2 * n_models);

  #ifdef _OPENMP
  omp_set_num_threads(n_threads);
  #endif

  #pragma omp parallel for schedule(dynamic, 20)
  for (uword i = 0; i < P; i++) {
    arma::rowvec beta = beta_stage1.row(i);
    arma::rowvec se = se_stage1.row(i);
    arma::rowvec weight = 1.0 / arma::square(se);

    for (uword t = 0; t < T; t++) {
      if (!std::isfinite(weight(t)) || !std::isfinite(beta(t))) {
        weight(t) = 0.0;
      }
    }

    arma::vec beta_vec = beta.t();
    arma::vec weight_vec = weight.t();

    for (int m = 0; m < n_models; m++) {
      const arma::mat& W = W_mats[m];
      uword k = k_sizes[m];

      arma::mat W_weighted = W.each_col() % weight_vec;
      arma::mat WtSinvW = W.t() * W_weighted;
      arma::vec WtSinvBeta = W.t() * (weight_vec % beta_vec);

      arma::vec theta;
      bool ok = arma::solve(theta, WtSinvW, WtSinvBeta,
                            arma::solve_opts::fast + arma::solve_opts::likely_sympd);
      if (!ok || !theta.is_finite()) {
        out(i, 2 * m) = NA_REAL;
        out(i, 2 * m + 1) = NA_REAL;
        continue;
      }

      arma::vec resid = beta_vec - W * theta;
      double Q = arma::dot(resid % weight_vec, resid);
      out(i, 2 * m) = Q + 2.0 * k;
      out(i, 2 * m + 1) = Q + k * std::log(N);
    }
  }

  return out;
}

// [[Rcpp::export]]
NumericMatrix fit_multi_models_tau_cov(const arma::mat& beta_stage1,
                                       const NumericVector& cov_vec,
                                       const NumericVector& offsets,
                                       const List& W_list,
                                       double N,
                                       int n_threads = 1,
                                       Nullable<NumericVector> cov_scales = R_NilValue) {
  const uword P = beta_stage1.n_rows;
  const uword T = beta_stage1.n_cols;
  const int n_models = W_list.size();

  if (offsets.size() < P) {
    stop("offsets must have at least one entry per SNP.");
  }
  NumericVector scales;
  if (cov_scales.isNotNull()) scales = NumericVector(cov_scales);
  if (scales.size() != 0 && scales.size() < P) {
    stop("cov_scales must have at least one entry per SNP when supplied.");
  }
  if (n_models == 0) {
    return NumericMatrix(P, 0);
  }

  std::vector<arma::mat> W_mats;
  std::vector<uword> k_sizes;
  W_mats.reserve(n_models);
  k_sizes.reserve(n_models);

  for (int m = 0; m < n_models; m++) {
    arma::mat W = as<arma::mat>(W_list[m]);
    if (W.n_rows != T) {
      stop("Each W matrix must have T rows.");
    }
    W_mats.push_back(W);
    k_sizes.push_back(W.n_cols);
  }

  const size_t n_cov = static_cast<size_t>(T) * (T + 1) / 2;
  const size_t cov_len = static_cast<size_t>(cov_vec.size());

  NumericMatrix out(P, 2 * n_models);

  #ifdef _OPENMP
  omp_set_num_threads(n_threads);
  #endif

  #pragma omp parallel for schedule(dynamic, 20)
  for (uword i = 0; i < P; i++) {
    arma::rowvec beta = beta_stage1.row(i);
    if (!beta.is_finite()) {
      for (int m = 0; m < n_models; m++) {
        out(i, 2 * m) = NA_REAL;
        out(i, 2 * m + 1) = NA_REAL;
      }
      continue;
    }

    size_t start = static_cast<size_t>(offsets[i]) / 4;
    if (start + n_cov > cov_len) {
      for (int m = 0; m < n_models; m++) {
        out(i, 2 * m) = NA_REAL;
        out(i, 2 * m + 1) = NA_REAL;
      }
      continue;
    }

    double scale = 1.0;
    if (scales.size() != 0) {
      scale = scales[i];
      if (!std::isfinite(scale) || scale <= 0.0) {
        for (int m = 0; m < n_models; m++) {
          out(i, 2 * m) = NA_REAL;
          out(i, 2 * m + 1) = NA_REAL;
        }
        continue;
      }
    }

    arma::mat Sigma(T, T, arma::fill::zeros);
    bool cov_ok = true;
    size_t idx = 0;
    for (uword r = 0; r < T; r++) {
      for (uword c = r; c < T; c++) {
        double v = cov_vec[start + idx] * scale;
        if (!std::isfinite(v)) {
          cov_ok = false;
          break;
        }
        Sigma(r, c) = v;
        Sigma(c, r) = v;
        idx++;
      }
      if (!cov_ok) break;
    }

    if (!cov_ok) {
      for (int m = 0; m < n_models; m++) {
        out(i, 2 * m) = NA_REAL;
        out(i, 2 * m + 1) = NA_REAL;
      }
      continue;
    }

    arma::vec beta_vec = beta.t();
    arma::vec SinvBeta;
    bool ok = arma::solve(SinvBeta, Sigma, beta_vec,
                          arma::solve_opts::fast + arma::solve_opts::likely_sympd);
    if (!ok || !SinvBeta.is_finite()) {
      for (int m = 0; m < n_models; m++) {
        out(i, 2 * m) = NA_REAL;
        out(i, 2 * m + 1) = NA_REAL;
      }
      continue;
    }

    for (int m = 0; m < n_models; m++) {
      const arma::mat& W = W_mats[m];
      uword k = k_sizes[m];

      arma::mat SinvW;
      ok = arma::solve(SinvW, Sigma, W,
                       arma::solve_opts::fast + arma::solve_opts::likely_sympd);
      if (!ok || !SinvW.is_finite()) {
        out(i, 2 * m) = NA_REAL;
        out(i, 2 * m + 1) = NA_REAL;
        continue;
      }

      arma::mat WtSinvW = W.t() * SinvW;
      arma::vec WtSinvBeta = W.t() * SinvBeta;

      arma::vec theta;
      ok = arma::solve(theta, WtSinvW, WtSinvBeta,
                       arma::solve_opts::fast + arma::solve_opts::likely_sympd);
      if (!ok || !theta.is_finite()) {
        out(i, 2 * m) = NA_REAL;
        out(i, 2 * m + 1) = NA_REAL;
        continue;
      }

      arma::vec resid = beta_vec - W * theta;
      arma::vec SinvResid = SinvBeta - SinvW * theta;
      double Q = arma::dot(resid, SinvResid);

      out(i, 2 * m) = Q + 2.0 * k;
      out(i, 2 * m + 1) = Q + k * std::log(N);
    }
  }

  return out;
}

// [[Rcpp::export]]
List fit_multi_models_ols_full_cov(const arma::mat& beta_stage1,
                                   const NumericVector& cov_vec,
                                   const NumericVector& offsets,
                                   const List& W_list,
                                   double N,
                                   int n_threads = 1,
                                   Nullable<NumericVector> cov_scales = R_NilValue) {
  const uword P = beta_stage1.n_rows;
  const uword T = beta_stage1.n_cols;
  const int n_models = W_list.size();

  if (offsets.size() < P) {
    stop("offsets must have at least one entry per SNP.");
  }
  NumericVector scales;
  if (cov_scales.isNotNull()) scales = NumericVector(cov_scales);
  if (scales.size() != 0 && scales.size() < P) {
    stop("cov_scales must have at least one entry per SNP when supplied.");
  }
  if (n_models == 0) {
    return List::create(Named("Q_full") = NumericMatrix(P, 0),
                        Named("AIC") = NumericMatrix(P, 0),
                        Named("BIC") = NumericMatrix(P, 0),
                        Named("rank") = IntegerVector(0),
                        Named("rank_deficient") = LogicalVector(0),
                        Named("reason_code") = IntegerMatrix(P, 0),
                        Named("params") = List::create());
  }

  // OLS coefficients are the least-squares fit; the penalty uses the actual
  // rank of W, not its column count, so a duplicated column is visible in the
  // reported rank and is not scored as a full-rank design.
  std::vector<arma::mat> W_mats;
  std::vector<uword> k_sizes;
  std::vector<uword> ranks;
  std::vector<bool> rank_def;
  W_mats.reserve(n_models);
  k_sizes.reserve(n_models);
  ranks.reserve(n_models);
  rank_def.reserve(n_models);

  for (int m = 0; m < n_models; m++) {
    arma::mat W = as<arma::mat>(W_list[m]);
    if (W.n_rows != T) {
      stop("Each W matrix must have T rows.");
    }
    uword r = arma::rank(W);
    W_mats.push_back(W);
    k_sizes.push_back(W.n_cols);
    ranks.push_back(r);
    rank_def.push_back(r < W.n_cols);
  }

  const size_t n_cov = static_cast<size_t>(T) * (T + 1) / 2;
  const size_t cov_len = static_cast<size_t>(cov_vec.size());

  NumericMatrix Q_full(P, n_models);
  NumericMatrix AIC(P, n_models);
  NumericMatrix BIC(P, n_models);
  IntegerMatrix reason_code(P, n_models);
  std::fill(Q_full.begin(), Q_full.end(), NA_REAL);
  std::fill(AIC.begin(), AIC.end(), NA_REAL);
  std::fill(BIC.begin(), BIC.end(), NA_REAL);
  std::fill(reason_code.begin(), reason_code.end(), 0);

  std::vector<arma::mat> params(n_models);
  for (int m = 0; m < n_models; m++) {
    params[m] = arma::mat(k_sizes[m], P, arma::fill::value(NA_REAL));
  }

  IntegerVector rank_out(n_models);
  LogicalVector rank_deficient(n_models);
  for (int m = 0; m < n_models; m++) {
    rank_out[m] = static_cast<int>(ranks[m]);
    rank_deficient[m] = rank_def[m];
  }

  const double logN = std::log(N);

  #ifdef _OPENMP
  omp_set_num_threads(n_threads);
  #endif

  #pragma omp parallel for schedule(dynamic, 20)
  for (uword i = 0; i < P; i++) {
    arma::rowvec beta = beta_stage1.row(i);
    if (!beta.is_finite()) {
      for (int m = 0; m < n_models; m++) reason_code(i, m) = 1;
      continue;
    }

    const double off_i = offsets[i];
    if (!std::isfinite(off_i) || off_i < 0.0) {
      for (int m = 0; m < n_models; m++) reason_code(i, m) = 2;
      continue;
    }
    const size_t start = static_cast<size_t>(off_i) / 4;
    if (start + n_cov > cov_len) {
      for (int m = 0; m < n_models; m++) reason_code(i, m) = 2;
      continue;
    }

    double scale = 1.0;
    if (scales.size() != 0) {
      scale = scales[i];
      if (!std::isfinite(scale) || scale <= 0.0) {
        for (int m = 0; m < n_models; m++) reason_code(i, m) = 3;
        continue;
      }
    }

    arma::mat Sigma(T, T, arma::fill::zeros);
    bool cov_ok = true;
    size_t idx = 0;
    for (uword r = 0; r < T; r++) {
      for (uword c = r; c < T; c++) {
        double v = cov_vec[start + idx] * scale;
        if (!std::isfinite(v)) {
          cov_ok = false;
          break;
        }
        Sigma(r, c) = v;
        Sigma(c, r) = v;
        idx++;
      }
      if (!cov_ok) break;
    }
    if (!cov_ok) {
      for (int m = 0; m < n_models; m++) reason_code(i, m) = 4;
      continue;
    }

    // No silent repair: if Sigma is not positive definite, or is numerically
    // singular, the SNP is refused with a reason rather than ridge-repaired.
    arma::mat chol_tmp;
    if (!arma::chol(chol_tmp, Sigma, "lower")) {
      for (int m = 0; m < n_models; m++) reason_code(i, m) = 5;
      continue;
    }
    double rc = arma::rcond(Sigma);
    if (!std::isfinite(rc) || rc < 1e-12) {
      for (int m = 0; m < n_models; m++) reason_code(i, m) = 6;
      continue;
    }

    arma::vec beta_vec = beta.t();
    arma::vec SinvBeta;
    bool ok = arma::solve(SinvBeta, Sigma, beta_vec,
                          arma::solve_opts::fast + arma::solve_opts::likely_sympd);
    if (!ok || !SinvBeta.is_finite()) {
      for (int m = 0; m < n_models; m++) reason_code(i, m) = 7;
      continue;
    }

    for (int m = 0; m < n_models; m++) {
      const arma::mat& W = W_mats[m];

      arma::mat SinvW;
      ok = arma::solve(SinvW, Sigma, W,
                       arma::solve_opts::fast + arma::solve_opts::likely_sympd);
      if (!ok || !SinvW.is_finite()) {
        reason_code(i, m) = 7;
        continue;
      }

      arma::vec theta;
      if (rank_def[m]) {
        // Minimum-norm least-squares solution: Q_full is invariant to the
        // choice of solution because the residual is the projection residual.
        theta = arma::pinv(W) * beta_vec;
      } else {
        arma::mat WtW = W.t() * W;
        arma::vec Wtb = W.t() * beta_vec;
        ok = arma::solve(theta, WtW, Wtb,
                         arma::solve_opts::fast + arma::solve_opts::likely_sympd);
      }
      if (!theta.is_finite()) {
        reason_code(i, m) = 8;
        continue;
      }

      arma::vec resid = beta_vec - W * theta;
      arma::vec SinvResid = SinvBeta - SinvW * theta;
      double Q = arma::dot(resid, SinvResid);
      if (!std::isfinite(Q)) {
        reason_code(i, m) = 9;
        continue;
      }
      if (Q < 0.0) {
        if (Q > -1e-8) {
          Q = 0.0;
        } else {
          reason_code(i, m) = 9;
          continue;
        }
      }

      Q_full(i, m) = Q;
      AIC(i, m) = Q + 2.0 * static_cast<double>(ranks[m]);
      BIC(i, m) = Q + static_cast<double>(ranks[m]) * logN;
      params[m].col(i) = theta;
    }
  }

  return List::create(
    Named("Q_full") = Q_full,
    Named("AIC") = AIC,
    Named("BIC") = BIC,
    Named("rank") = rank_out,
    Named("rank_deficient") = rank_deficient,
    Named("reason_code") = reason_code,
    Named("params") = wrap(params)
  );
}
