#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

#ifdef _OPENMP
#include <omp.h>
#endif

//' Compute Calibrated Standard Errors
//' 
//' @param cov_vec Numeric vector containing concatenated upper triangles of covariance matrices.
//' @param offsets Integer vector of byte offsets (0-based) into the original file (needs conversion to index).
//' @param A The transformation matrix (n_params x K).
//' @param K Dimension of the covariance matrix (e.g. 17).
//' @return NumericMatrix of standard errors (n_snps x n_params).
// [[Rcpp::export]]
NumericMatrix compute_calibrated_se(NumericVector cov_vec, NumericVector offsets, NumericMatrix A, int K) {
  int n_snps = offsets.size();
  int n_params = A.nrow();
  int n_cov_elements = K * (K + 1) / 2;
  
  arma::mat A_mat(A.begin(), n_params, K, false); // Wrap existing memory
  arma::mat Sigma(K, K, arma::fill::zeros);
  NumericMatrix out_se(n_snps, n_params);
  
  // R uses 32-bit integers for vector indexing, but offsets might be large if file is >2GB.
  // However, cov_vec is in memory, so it fits in R vector (max 2^31-1 elements).
  // The offsets passed here are from the .idx file, which are BYTE offsets.
  // cov_vec is a vector of floats (4 bytes).
  
  for (int i = 0; i < n_snps; i++) {
    double byte_offset = offsets[i];
    
    // Convert byte offset to vector index (0-based)
    // dividing by 4 because floats are 4 bytes
    // R is 1-based, but C++ is 0-based.
    // offsets[i] is start of record in bytes.
    long long start_idx = (long long)(byte_offset / 4);
    
    if (start_idx < 0 || start_idx + n_cov_elements > cov_vec.size()) {
       // Index out of bounds, skip or produce NA
       for (int p=0; p<n_params; p++) out_se(i, p) = NA_REAL;
       continue;
    }
    
    // Fill Sigma upper triangle (row-major order)
    // numpy: Sigma[np.triu_indices(K)] yields row-major order:
    // Row 0: (0,0),(0,1),(0,2),...
    // Row 1: (1,1),(1,2),...
    int idx = 0;
    for (int r = 0; r < K; r++) {
      for (int c = r; c < K; c++) {
        Sigma(r, c) = cov_vec[start_idx + idx];
        idx++;
      }
    }
    // Symmetrize
    for (int r = 0; r < K; r++) {
      for (int c = r + 1; c < K; c++) {
        Sigma(c, r) = Sigma(r, c);
      }
    }
    
    // V_param = A * Sigma * A.t()
    // We only need the diagonal (variances)
    // Optimization: V_j = (A.row(j)) * Sigma * (A.row(j)).t()
    
    for (int p = 0; p < n_params; p++) {
       arma::rowvec a_row = A_mat.row(p);
       double var = arma::as_scalar(a_row * Sigma * a_row.t());
       if (std::isfinite(var) && var > 0.0) {
         out_se(i, p) = std::sqrt(var);
       } else {
         out_se(i, p) = NA_REAL;
       }
    }
  }
  
  return out_se;
}

//' Compute parameter covariance upper triangles
//'
//' @param cov_vec Numeric vector containing concatenated upper triangles of covariance matrices.
//' @param offsets Integer vector of byte offsets (0-based) into the original file (needs conversion to index).
//' @param A The transformation matrix (n_params x K).
//' @param K Dimension of the covariance matrix (e.g. 17).
//' @return NumericMatrix of upper-triangle covariance (n_snps x n_cov_params).
// [[Rcpp::export]]
NumericMatrix compute_param_cov(NumericVector cov_vec, NumericVector offsets, NumericMatrix A, int K) {
  int n_snps = offsets.size();
  int n_params = A.nrow();
  int n_cov_elements = K * (K + 1) / 2;
  int n_cov_params = n_params * (n_params + 1) / 2;

  arma::mat A_mat(A.begin(), n_params, K, false);
  arma::mat Sigma(K, K, arma::fill::zeros);
  NumericMatrix out_cov(n_snps, n_cov_params);

  for (int i = 0; i < n_snps; i++) {
    double byte_offset = offsets[i];
    long long start_idx = (long long)(byte_offset / 4);

    if (start_idx < 0 || start_idx + n_cov_elements > cov_vec.size()) {
      for (int p = 0; p < n_cov_params; p++) out_cov(i, p) = NA_REAL;
      continue;
    }

    int idx = 0;
    for (int r = 0; r < K; r++) {
      for (int c = r; c < K; c++) {
        Sigma(r, c) = cov_vec[start_idx + idx];
        idx++;
      }
    }
    for (int r = 0; r < K; r++) {
      for (int c = r + 1; c < K; c++) {
        Sigma(c, r) = Sigma(r, c);
      }
    }

    arma::mat Vtheta = A_mat * Sigma * A_mat.t();

    int out_idx = 0;
    for (int r = 0; r < n_params; r++) {
      for (int c = r; c < n_params; c++) {
        out_cov(i, out_idx) = Vtheta(r, c);
        out_idx++;
      }
    }
  }

  return out_cov;
}

//' Compute GLS Stage 2 Parameter GWAS
//'
//' Uses the per-SNP tau covariance matrix from Stage 1 to estimate SNP-specific
//' generalized least-squares parameter effects and their covariance.
//'
//' @param beta_stage1 Tau-level SNP effects as a K x n_snps matrix.
//' @param cov_vec Numeric vector containing concatenated upper triangles of covariance matrices.
//' @param offsets Integer vector of byte offsets (0-based) into the original file.
//' @param W Weight/design matrix with K rows and n_params columns.
//' @param return_cov Whether to return parameter covariance upper triangles.
//' @param gls_ridge Relative diagonal ridge added only if the covariance solve fails.
//' @param n_threads Number of OpenMP threads.
//' @return List with params, se, Q, and optionally param_cov_upper.
// [[Rcpp::export]]
List compute_gls_param_gwas(const arma::mat& beta_stage1,
                            NumericVector cov_vec,
                            NumericVector offsets,
                            NumericMatrix W,
                            bool return_cov = false,
                            double gls_ridge = 1e-8,
                            int n_threads = 1) {
  const int K = beta_stage1.n_rows;
  const int n_snps = beta_stage1.n_cols;
  const int n_params = W.ncol();
  const int n_cov_elements = K * (K + 1) / 2;
  const int n_cov_params = n_params * (n_params + 1) / 2;

  if (W.nrow() != K) {
    stop("W must have the same number of rows as beta_stage1.");
  }
  if (offsets.size() < n_snps) {
    stop("offsets must have at least one entry per SNP.");
  }

  arma::mat W_mat(W.begin(), K, n_params, false);
  arma::mat params(n_params, n_snps, arma::fill::value(NA_REAL));
  arma::mat se(n_params, n_snps, arma::fill::value(NA_REAL));
  NumericVector Q(n_snps, NA_REAL);
  NumericMatrix param_cov_upper;
  if (return_cov) {
    param_cov_upper = NumericMatrix(n_snps, n_cov_params);
    std::fill(param_cov_upper.begin(), param_cov_upper.end(), NA_REAL);
  }

  // Diagnostic counters. Purely observational: they do not affect any numeric
  // result. Used to report how often the ridge repair path fires on real data.
  long long d_beta_nonfinite = 0;
  long long d_offset_bad = 0;
  long long d_cov_nonfinite = 0;
  long long d_ridge_tried = 0;
  long long d_ridge_rescued = 0;
  long long d_sigma_not_spd = 0;
  long long d_solve_fail = 0;
  long long d_xtx_not_spd = 0;
  long long d_theta_fail = 0;
  long long d_vtheta_fail = 0;
  long long d_q_bad = 0;
  long long d_ok = 0;

  #ifdef _OPENMP
  omp_set_num_threads(n_threads);
  #endif

  #pragma omp parallel for schedule(dynamic, 20) \
    reduction(+:d_beta_nonfinite,d_offset_bad,d_cov_nonfinite,d_ridge_tried, \
                d_ridge_rescued,d_sigma_not_spd,d_solve_fail,d_xtx_not_spd, \
                d_theta_fail,d_vtheta_fail,d_q_bad,d_ok)
  for (int i = 0; i < n_snps; i++) {
    const arma::vec beta_vec = beta_stage1.col(i);
    if (!beta_vec.is_finite()) { d_beta_nonfinite++; continue; }

    const long long start_idx = (long long)(offsets[i] / 4);
    if (start_idx < 0 || start_idx + n_cov_elements > cov_vec.size()) { d_offset_bad++; continue; }

    arma::mat Sigma(K, K, arma::fill::zeros);
    bool cov_ok = true;
    int idx = 0;
    for (int r = 0; r < K; r++) {
      for (int c = r; c < K; c++) {
        double v = cov_vec[start_idx + idx];
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
    if (!cov_ok) { d_cov_nonfinite++; continue; }

    arma::mat chol_tmp;
    bool sigma_spd = arma::chol(chol_tmp, Sigma, "lower");
    if (!sigma_spd && gls_ridge > 0.0) {
      d_ridge_tried++;
      const arma::vec d = Sigma.diag();
      double diag_scale = arma::mean(d.elem(arma::find_finite(d)));
      if (!std::isfinite(diag_scale) || diag_scale <= 0.0) diag_scale = 1.0;
      Sigma.diag() += gls_ridge * diag_scale;
      sigma_spd = arma::chol(chol_tmp, Sigma, "lower");
      if (sigma_spd) d_ridge_rescued++;
    }
    if (!sigma_spd) { d_sigma_not_spd++; continue; }

    arma::mat rhs(K, n_params + 1);
    rhs.cols(0, n_params - 1) = W_mat;
    rhs.col(n_params) = beta_vec;

    arma::mat SinvRhs;
    bool ok = arma::solve(SinvRhs, Sigma, rhs,
                          arma::solve_opts::fast + arma::solve_opts::likely_sympd);

    if (!ok || !SinvRhs.is_finite()) { d_solve_fail++; continue; }

    arma::mat SinvW = SinvRhs.cols(0, n_params - 1);
    arma::vec SinvBeta = SinvRhs.col(n_params);
    arma::mat WtSinvW = W_mat.t() * SinvW;
    WtSinvW = 0.5 * (WtSinvW + WtSinvW.t());
    if (!arma::chol(chol_tmp, WtSinvW, "lower")) { d_xtx_not_spd++; continue; }

    arma::vec WtSinvBeta = W_mat.t() * SinvBeta;

    arma::vec theta;
    ok = arma::solve(theta, WtSinvW, WtSinvBeta,
                     arma::solve_opts::fast + arma::solve_opts::likely_sympd);
    if (!ok || !theta.is_finite()) { d_theta_fail++; continue; }

    arma::mat Vtheta;
    ok = arma::solve(Vtheta, WtSinvW, arma::eye(n_params, n_params),
                     arma::solve_opts::fast + arma::solve_opts::likely_sympd);
    if (!ok || !Vtheta.is_finite()) { d_vtheta_fail++; continue; }

    arma::vec resid = beta_vec - W_mat * theta;
    arma::vec SinvResid = SinvBeta - SinvW * theta;
    double q_val = arma::dot(resid, SinvResid);
    if (!std::isfinite(q_val)) { d_q_bad++; continue; }
    if (q_val < -1e-8) { d_q_bad++; continue; }
    if (q_val < 0.0) q_val = 0.0;

    d_ok++;
    params.col(i) = theta;
    for (int p = 0; p < n_params; p++) {
      double var = Vtheta(p, p);
      se(p, i) = (std::isfinite(var) && var > 0.0) ? std::sqrt(var) : NA_REAL;
    }
    Q[i] = q_val;

    if (return_cov) {
      int out_idx = 0;
      for (int r = 0; r < n_params; r++) {
        for (int c = r; c < n_params; c++) {
          param_cov_upper(i, out_idx) = Vtheta(r, c);
          out_idx++;
        }
      }
    }
  }

  List out = List::create(
    _["params"] = params,
    _["se"] = se,
    _["Q"] = Q
  );
  out["diag"] = List::create(
    _["n_snps"] = (double) n_snps,
    _["n_ok"] = (double) d_ok,
    _["n_beta_nonfinite"] = (double) d_beta_nonfinite,
    _["n_offset_bad"] = (double) d_offset_bad,
    _["n_cov_nonfinite"] = (double) d_cov_nonfinite,
    _["n_ridge_tried"] = (double) d_ridge_tried,
    _["n_ridge_rescued"] = (double) d_ridge_rescued,
    _["n_sigma_not_spd"] = (double) d_sigma_not_spd,
    _["n_solve_fail"] = (double) d_solve_fail,
    _["n_xtx_not_spd"] = (double) d_xtx_not_spd,
    _["n_theta_fail"] = (double) d_theta_fail,
    _["n_vtheta_fail"] = (double) d_vtheta_fail,
    _["n_q_bad"] = (double) d_q_bad
  );
  if (return_cov) {
    out["param_cov_upper"] = param_cov_upper;
  }
  return out;
}
