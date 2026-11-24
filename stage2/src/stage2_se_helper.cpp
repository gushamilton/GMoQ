#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
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
//' @param cov_scales Optional per-SNP positive multiplier for a shared covariance payload.
//' @param A The transformation matrix (n_params x K).
//' @param K Dimension of the covariance matrix (e.g. 17).
//' @return NumericMatrix of standard errors (n_snps x n_params).
// [[Rcpp::export]]
NumericMatrix compute_calibrated_se(NumericVector cov_vec, NumericVector offsets,
                                    NumericMatrix A, int K,
                                    Nullable<NumericVector> cov_scales = R_NilValue) {
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
  NumericVector scales;
  if (cov_scales.isNotNull()) scales = NumericVector(cov_scales);
  if (scales.size() != 0 && scales.size() < n_snps) {
    stop("cov_scales must have at least one entry per SNP when supplied.");
  }
  
  // Cache base variances diag(A * Sigma_0 * A') when consecutive SNPs share a
  // pooled Sigma_0 record (same start_idx, differing only by scalar cov_scales).
  long long cached_start_idx = -1;
  arma::vec base_var(n_params, arma::fill::value(NA_REAL));

  for (int i = 0; i < n_snps; i++) {
    double byte_offset = offsets[i];

    // Convert byte offset to vector index (0-based)
    // dividing by 4 because floats are 4 bytes
    // R is 1-based, but C++ is 0-based.
    // offsets[i] is start of record in bytes.
    long long start_idx = (long long)(byte_offset / 4);

    double scale = 1.0;
    if (scales.size() != 0) {
      scale = scales[i];
      if (!std::isfinite(scale) || scale <= 0.0) {
        for (int p=0; p<n_params; p++) out_se(i, p) = NA_REAL;
        continue;
      }
    }

    if (start_idx < 0 || start_idx + n_cov_elements > cov_vec.size()) {
       // Index out of bounds, skip or produce NA
       for (int p=0; p<n_params; p++) out_se(i, p) = NA_REAL;
       continue;
    }

    if (start_idx != cached_start_idx) {
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
      for (int p = 0; p < n_params; p++) {
        arma::rowvec a_row = A_mat.row(p);
        base_var[p] = arma::as_scalar(a_row * Sigma * a_row.t());
      }
      cached_start_idx = start_idx;
    }

    for (int p = 0; p < n_params; p++) {
       double var = base_var[p] * scale;
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
//' @param cov_scales Optional per-SNP positive multiplier for a shared covariance payload.
//' @param A The transformation matrix (n_params x K).
//' @param K Dimension of the covariance matrix (e.g. 17).
//' @return NumericMatrix of upper-triangle covariance (n_snps x n_cov_params).
// [[Rcpp::export]]
NumericMatrix compute_param_cov(NumericVector cov_vec, NumericVector offsets,
                                NumericMatrix A, int K,
                                Nullable<NumericVector> cov_scales = R_NilValue) {
  int n_snps = offsets.size();
  int n_params = A.nrow();
  int n_cov_elements = K * (K + 1) / 2;
  int n_cov_params = n_params * (n_params + 1) / 2;

  arma::mat A_mat(A.begin(), n_params, K, false);
  arma::mat Sigma(K, K, arma::fill::zeros);
  NumericMatrix out_cov(n_snps, n_cov_params);

  NumericVector scales;
  if (cov_scales.isNotNull()) scales = NumericVector(cov_scales);
  if (scales.size() != 0 && scales.size() < n_snps) {
    stop("cov_scales must have at least one entry per SNP when supplied.");
  }

  long long cached_start_idx = -1;
  arma::vec base_vtheta_upper(n_cov_params, arma::fill::value(NA_REAL));

  for (int i = 0; i < n_snps; i++) {
    double byte_offset = offsets[i];
    long long start_idx = (long long)(byte_offset / 4);

    double scale = 1.0;
    if (scales.size() != 0) {
      scale = scales[i];
      if (!std::isfinite(scale) || scale <= 0.0) {
        for (int p = 0; p < n_cov_params; p++) out_cov(i, p) = NA_REAL;
        continue;
      }
    }

    if (start_idx < 0 || start_idx + n_cov_elements > cov_vec.size()) {
      for (int p = 0; p < n_cov_params; p++) out_cov(i, p) = NA_REAL;
      continue;
    }

    if (start_idx != cached_start_idx) {
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
          base_vtheta_upper[out_idx++] = Vtheta(r, c);
        }
      }
      cached_start_idx = start_idx;
    }

    for (int out_idx = 0; out_idx < n_cov_params; out_idx++) {
      out_cov(i, out_idx) = base_vtheta_upper[out_idx] * scale;
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
//' @param cov_scales Optional per-SNP positive multiplier for a shared covariance payload.
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
                            int n_threads = 1,
                            Nullable<NumericVector> cov_scales = R_NilValue) {
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
  NumericVector scales;
  if (cov_scales.isNotNull()) scales = NumericVector(cov_scales);
  if (scales.size() != 0 && scales.size() < n_snps) {
    stop("cov_scales must have at least one entry per SNP when supplied.");
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

    double scale = 1.0;
    if (scales.size() != 0) {
      scale = scales[i];
      if (!std::isfinite(scale) || scale <= 0.0) { d_cov_nonfinite++; continue; }
    }

    arma::mat Sigma(K, K, arma::fill::zeros);
    bool cov_ok = true;
    int idx = 0;
    for (int r = 0; r < K; r++) {
      for (int c = r; c < K; c++) {
        double v = cov_vec[start_idx + idx] * scale;
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

// ---------------------------------------------------------------------------
// Calibrated null distribution for the Stage 2 lack-of-fit statistic Q
// ---------------------------------------------------------------------------
//
// Q = r' D^-1 r with r = M b, M = I - W (W'W)^-1 W', D = diag(se_tau^2).
// Under the null b ~ N(W theta, Sigma), so r ~ N(0, M Sigma M') and
//
//     Q ~ sum_i lambda_i chi2_1,
//
// where lambda_i are the eigenvalues of B = D^-1/2 M Sigma M' D^-1/2.
// This requires Sigma only -- never Sigma^-1 -- which is what makes it usable
// on the stored int8 covariance.
//
// Returns, per SNP, the first four "power sums" c_k = sum lambda^k (from which
// R forms the Satterthwaite and Liu approximations) together with a
// Lugannani-Rice saddlepoint tail probability, which is accurate far enough
// into the tail to be used at genome-wide thresholds.
//
// Only pure double arithmetic happens inside the parallel region; every R API
// call is left to the caller.

namespace {

// Upper-tail probability of the standard normal.
inline double norm_sf(double x) {
  return 0.5 * std::erfc(x / std::sqrt(2.0));
}

// Lugannani-Rice saddlepoint upper-tail probability for sum_i lam_i chi2_1.
// Returns a negative value if it is not applicable (caller falls back).
double saddlepoint_sf(double q, const arma::vec& lam) {
  const arma::uword n = lam.n_elem;
  if (n == 0) return -1.0;
  double mu = arma::accu(lam);          // = K'(0)
  double lmax = lam.max();
  if (!(q > mu) || lmax <= 0.0) return -1.0;   // only valid in the upper tail

  // DOMAIN OF CONVERGENCE.  K(t) = -0.5 sum log(1 - 2 lam_i t) is finite only
  // where every 1 - 2 lam_i t > 0, i.e. on the open interval
  //
  //     1/(2 lam_min) < t < 1/(2 lam_max)      (lower bound -inf if lam_min>=0)
  //
  // With a two-sided spectrum -- which is exactly what keeping the negative
  // eigenvalues produces -- the interval is bounded at BOTH ends, and stepping
  // outside it does not crash: 1 - 2 lam_i t merely goes negative for some i
  // and K'/K'' silently return values from the wrong branch.  The bracket must
  // therefore never be allowed to leave the interval.
  //
  // K' is strictly increasing on the interval (K'' = sum 2 lam^2/(1-2 lam t)^2
  // > 0), so a bracket is valid iff K'(lo) <= q <= K'(hi) with both endpoints
  // inside.  We are only ever called with q > mu = K'(0), and t = 0 is in the
  // interval for ANY spectrum, so lo = 0 is a provably safe lower bracket and
  // no downward expansion is needed at all.  The root zeta is in (0, hi_lim).
  const double hi_lim = 0.5 / lmax;                       // > 0, always finite

  auto kprime = [&](double z) {
    double s = 0.0;
    for (arma::uword i = 0; i < n; ++i) s += lam[i] / (1.0 - 2.0 * lam[i] * z);
    return s;
  };

  double lo = 0.0;
  // Upper bracket: K'(t) -> +inf as t -> hi_lim from below, so a point with
  // K'(hi) > q always exists strictly inside.  Walk in from a safe interior
  // start by halving the remaining gap; never set hi >= hi_lim.
  double hi = hi_lim * 0.5;
  int guard = 0;
  while (!(kprime(hi) > q)) {
    double next = hi + 0.5 * (hi_lim - hi);
    if (!(next > hi) || !(next < hi_lim)) return -1.0;    // gap underflowed
    hi = next;
    if (++guard > 200) return -1.0;
  }

  // Safeguarded Newton-Raphson on (lo, hi): K' is strictly increasing and
  // convex for positive eigenvalues; falls back to bisection whenever a
  // Newton step would leave the current bracket.
  double zeta = 0.5 * (lo + hi);
  for (int it = 0; it < 200; ++it) {
    double kp = 0.0, kpp = 0.0;
    for (arma::uword i = 0; i < n; ++i) {
      const double li = lam[i];
      const double inv_den = 1.0 / (1.0 - 2.0 * li * zeta);
      const double term = li * inv_den;
      kp += term;
      kpp += 2.0 * term * term;
    }
    const double diff = kp - q;
    if (diff > 0.0) hi = zeta; else lo = zeta;
    if (hi - lo < 1e-14 * std::max(1.0, std::fabs(hi))) break;
    double cand = (kpp > 0.0) ? (zeta - diff / kpp) : (0.5 * (lo + hi));
    if (!(cand > lo && cand < hi)) cand = 0.5 * (lo + hi);
    if (std::fabs(cand - zeta) < 1e-14 * std::max(1.0, std::fabs(cand))) {
      zeta = cand;
      break;
    }
    zeta = cand;
  }

  double K = 0.0, K2 = 0.0;
  for (arma::uword i = 0; i < n; ++i) {
    double den = 1.0 - 2.0 * lam[i] * zeta;
    if (den <= 0.0) return -1.0;
    K -= 0.5 * std::log(den);
    K2 += 2.0 * lam[i] * lam[i] / (den * den);
  }
  double arg = 2.0 * (zeta * q - K);
  if (!(arg > 0.0) || !(K2 > 0.0)) return -1.0;

  double w = (zeta >= 0 ? 1.0 : -1.0) * std::sqrt(arg);
  double v = zeta * std::sqrt(K2);
  if (w == 0.0 || v <= 0.0) return -1.0;
  double p = norm_sf(w + std::log(v / w) / w);
  if (!std::isfinite(p) || p < 0.0 || p > 1.0) return -1.0;
  return p;
}

}  // namespace

//' Null-distribution quantities for the Stage 2 lack-of-fit statistic
//'
//' Computes, for each SNP, the power sums of the eigenvalues of
//' \code{D^-1/2 M Sigma t(M) D^-1/2} together with a saddlepoint tail
//' probability for the observed \code{Q}. Uses \code{Sigma} only; it never
//' inverts or factorises it, so it is valid on the stored int8 covariance.
//'
//' @param cov_vec Concatenated upper triangles of the per-SNP covariances.
//' @param offsets Byte offsets (0-based) of each SNP's record.
//' @param cov_scales Optional per-SNP positive multiplier for a shared covariance payload.
//' @param Qstat Numeric vector of observed Q statistics (n_snps).
//' @param Dvar Numeric matrix (K x n_snps) of per-tau variances, se_tau^2.
//' @param M Residual projector I - W (W'W)^-1 W' (K x K).
//' @param K Number of taus.
//' @param n_threads OpenMP threads.
//' @return Numeric matrix (n_snps x 6): c1, c2, c3, c4, p_saddle, n_eig.
//'   \code{p_saddle} is negative where the saddlepoint does not apply.
// [[Rcpp::export]]
NumericMatrix compute_q_calibration(NumericVector cov_vec, NumericVector offsets,
                                    NumericVector Qstat, NumericMatrix Dvar,
                                    NumericMatrix M, int K, int n_threads,
                                    Nullable<NumericVector> cov_scales = R_NilValue) {
  const int n_snps = offsets.size();
  const int n_cov_elements = K * (K + 1) / 2;

  NumericVector scales;
  if (cov_scales.isNotNull()) scales = NumericVector(cov_scales);
  if (scales.size() != 0 && scales.size() < n_snps) {
    stop("cov_scales must have at least one entry per SNP when supplied.");
  }

  arma::mat M_mat(M.begin(), K, K, false);
  arma::mat D_mat(Dvar.begin(), K, n_snps, false);
  NumericMatrix out(n_snps, 6);
  const double* cov_ptr = cov_vec.begin();
  const double* off_ptr = offsets.begin();
  const double* q_ptr = Qstat.begin();
  const double* scale_ptr = (scales.size() != 0) ? scales.begin() : nullptr;
  const long long cov_len = static_cast<long long>(cov_vec.size());
  double* out_ptr = out.begin();

#ifdef _OPENMP
  if (n_threads > 1) omp_set_num_threads(n_threads);
#pragma omp parallel
#endif
  {
    // Per-thread 1-entry cache for pooled covariance records: when consecutive
    // SNPs share Sigma_0 (same start_idx) and Dvar[:, i] is proportional to
    // scale * diag(Sigma_0), B = (1/r0) B_0 has eigenvalues lam_0 / r0, so the
    // K x K eigendecomposition runs once per unique Sigma_0 instead of once per SNP.
    long long cached_start = -1;
    bool cached_sigma_ok = false;
    bool cached_diag_pos = false;
    bool cached_eig_ready = false;
    bool cached_eig_ok = false;
    arma::mat Sigma0(K, K, arma::fill::zeros);
    arma::vec diag0(K, arma::fill::zeros);
    arma::vec lam0;
    double c1_0 = 0.0, c2_0 = 0.0, c3_0 = 0.0, c4_0 = 0.0;
    arma::vec inv_sqrt(K);
    arma::mat Sigma(K, K);

#ifdef _OPENMP
#pragma omp for schedule(dynamic, 256)
#endif
    for (int i = 0; i < n_snps; i++) {
      for (int p = 0; p < 6; p++) out_ptr[i + static_cast<size_t>(p) * n_snps] = NA_REAL;

      const long long start_idx = static_cast<long long>(off_ptr[i] / 4.0);
      if (start_idx < 0 || start_idx + n_cov_elements > cov_len) continue;

      double scale = 1.0;
      if (scale_ptr != nullptr) {
        scale = scale_ptr[i];
        if (!std::isfinite(scale) || scale <= 0.0) continue;
      }

      if (start_idx != cached_start) {
        int idx = 0;
        for (int r = 0; r < K; r++) {
          for (int c = r; c < K; c++) {
            const double v = cov_ptr[start_idx + idx++];
            Sigma0(r, c) = v;
            Sigma0(c, r) = v;
          }
        }
        cached_sigma_ok = Sigma0.is_finite();
        cached_diag_pos = cached_sigma_ok;
        for (int t = 0; t < K; t++) {
          diag0[t] = Sigma0(t, t);
          if (!(diag0[t] > 0.0) || !std::isfinite(diag0[t])) cached_diag_pos = false;
        }
        cached_eig_ready = false;
        cached_eig_ok = false;
        cached_start = start_idx;
      }
      if (!cached_sigma_ok) continue;

      const double* d = D_mat.colptr(i);
      bool bad = false;
      bool prop_to_diag0 = cached_diag_pos;
      double r_min = std::numeric_limits<double>::infinity();
      double r_max = 0.0;
      double r_sum = 0.0;
      for (int t = 0; t < K; t++) {
        const double v = d[t];
        if (!(v > 0.0) || !std::isfinite(v)) { bad = true; break; }
        if (prop_to_diag0) {
          const double rt = v / (scale * diag0[t]);
          if (rt < r_min) r_min = rt;
          if (rt > r_max) r_max = rt;
          r_sum += rt;
        }
      }
      if (bad) continue;
      // Match constant ratio within 1e-5 relative (absorbs %.6g TSV formatting of
      // pooled se_tau while always missing jackknife se_tau, which varies by ~2-10%).
      if (prop_to_diag0 && !(r_min > 0.0 && (r_max - r_min) <= 1e-5 * r_max)) {
        prop_to_diag0 = false;
      }

      if (prop_to_diag0) {
        if (!cached_eig_ready) {
          cached_eig_ready = true;
          for (int t = 0; t < K; t++) inv_sqrt[t] = 1.0 / std::sqrt(diag0[t]);
          arma::mat C0 = arma::diagmat(inv_sqrt) * M_mat;
          arma::mat B0 = C0 * Sigma0 * C0.t();
          B0 = 0.5 * (B0 + B0.t());
          arma::vec ev;
          if (arma::eig_sym(ev, B0)) {
            const double ev_tol = 1e-8;
            const double emax = arma::abs(ev).max();
            if (emax > 0.0) {
              lam0 = ev.elem(arma::find(arma::abs(ev) > ev_tol * emax));
              if (lam0.n_elem >= 2) {
                c1_0 = c2_0 = c3_0 = c4_0 = 0.0;
                for (arma::uword j = 0; j < lam0.n_elem; ++j) {
                  const double l = lam0[j], l2 = l * l;
                  c1_0 += l; c2_0 += l2; c3_0 += l2 * l; c4_0 += l2 * l2;
                }
                cached_eig_ok = true;
              }
            }
          }
        }
        if (!cached_eig_ok) continue;
        const double r0 = r_sum / static_cast<double>(K);
        const double inv_r = 1.0 / r0;
        const double inv_r2 = inv_r * inv_r;
        out_ptr[i + static_cast<size_t>(0) * n_snps] = c1_0 * inv_r;
        out_ptr[i + static_cast<size_t>(1) * n_snps] = c2_0 * inv_r2;
        out_ptr[i + static_cast<size_t>(2) * n_snps] = c3_0 * inv_r2 * inv_r;
        out_ptr[i + static_cast<size_t>(3) * n_snps] = c4_0 * inv_r2 * inv_r2;
        out_ptr[i + static_cast<size_t>(4) * n_snps] = saddlepoint_sf(q_ptr[i] * r0, lam0);
        out_ptr[i + static_cast<size_t>(5) * n_snps] = static_cast<double>(lam0.n_elem);
        continue;
      }

      // General per-SNP path (jackknife covariance or non-proportional Dvar).
      Sigma = Sigma0 * scale;
      for (int t = 0; t < K; t++) inv_sqrt[t] = 1.0 / std::sqrt(d[t]);
      arma::mat C = arma::diagmat(inv_sqrt) * M_mat;
      arma::mat B = C * Sigma * C.t();
      B = 0.5 * (B + B.t());

      arma::vec ev;
      if (!arma::eig_sym(ev, B)) continue;

      const double ev_tol = 1e-8;
      double emax = arma::abs(ev).max();          // spectral radius
      if (!(emax > 0.0)) continue;
      arma::vec lam = ev.elem(arma::find(arma::abs(ev) > ev_tol * emax));
      if (lam.n_elem < 2) continue;

      double c1 = 0, c2 = 0, c3 = 0, c4 = 0;
      for (arma::uword j = 0; j < lam.n_elem; ++j) {
        double l = lam[j], l2 = l * l;
        c1 += l; c2 += l2; c3 += l2 * l; c4 += l2 * l2;
      }
      out_ptr[i + static_cast<size_t>(0) * n_snps] = c1;
      out_ptr[i + static_cast<size_t>(1) * n_snps] = c2;
      out_ptr[i + static_cast<size_t>(2) * n_snps] = c3;
      out_ptr[i + static_cast<size_t>(3) * n_snps] = c4;
      out_ptr[i + static_cast<size_t>(4) * n_snps] = saddlepoint_sf(q_ptr[i], lam);
      out_ptr[i + static_cast<size_t>(5) * n_snps] = static_cast<double>(lam.n_elem);
    }
  }

  return out;
}
