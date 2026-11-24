// Fast pooled-only Stage 1 kernels.
//
// Two paths live here:
//
// 1. PooledBinPlan + pooled_bin_fast: exact O(N) bin-prefix kernel.
//    Because RIF_i(tau_t) = c_t - 1{y_i <= q_t} / f_t with c_t = q_t + tau_t/f_t,
//    each individual belongs to a single quantile bin b_i in {0, ..., T} where
//    1{y_i <= q_sorted[t]} == 1 iff b_i <= t.  For any centred genotype g_c = g - g_bar
//    (missing values imputed to g_bar, so g_c = 0 on missing entries),
//
//        g_c^T 1{y <= q_sorted[t]} = sum_{b=0}^{t} sum_{i : b_i = b} g_c[i],
//
//    which computes all T quantile cross-products in ONE O(N) pass over g_c
//    regardless of T, without ever allocating or streaming the N x T RIF matrix.
//    By Frisch-Waugh-Lovell with orthonormal covariate basis Q (Q^T Q = I_K),
//    let w = Q^T 1, e = 1 - Q w (identically 0 whenever 1 lies in col(Q)),
//    QtI = Q^T I, and QtR[:, t] = c_t * w - QtI[:, t] / f_t (precomputed once
//    in O(N K) time).  Then for every SNP:
//
//        QtGc = Q^T g_c
//        D    = ||g_c||^2 - ||QtGc||^2  [ - 2*g_bar*(w . QtGc) + g_bar^2*(N - ||w||^2) ]
//        num  = -g_c^T I[:, t] / f_t - QtGc^T QtR[:, t]  [ + g_bar * (e^T RIF[:, t]) ]
//        beta = num / D
//
//    where the bracketed terms vanish whenever col(Q) spans the intercept (all
//    standard CLI runs).  Sigma_0 = (1/N) (M RIF)^T (M RIF) is likewise computed
//    in O(T^2 K) from the cumulative bin counts and QtI, with zero N x T storage.
//
// 2. pooled_panel_fast: hoisted, register-blocked (AVX2+FMA on x86_64) pooled
//    kernel for callers that already pass a dense (N, T) residualised RIF matrix
//    (core.compute_pooled_panel).  Uses the same FWL covariate hoist and row-tiled
//    sweep so RIF is streamed in L2-resident chunks without a per-(i, p) K-loop.

#ifndef GMOQ_STAGE1_POOLED_FAST_HPP
#define GMOQ_STAGE1_POOLED_FAST_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(__GNUC__) && defined(__x86_64__) && !defined(GMOQ_NO_AVX2_CLONE)
#define GMOQ_POOLED_HAVE_AVX2 1
#include <immintrin.h>
#endif

namespace gmoq_pooled {

namespace py = pybind11;

enum { PBLK = 8, CBLK = 6, RBLK = 2048 };

// ---------------------------------------------------------------------------
// Register-blocked row-chunk cross-product: accm[c, p] += sum_r Mat[r, c] * Gt[r, p]
// Gt is (rlen, 8) row-major; Mat is (rlen, Cpad) row-major with Cpad % 6 == 0.
// ---------------------------------------------------------------------------

static void sweep_rows_base(
    const double* __restrict__ Gt,
    const double* __restrict__ Mat,
    size_t Cpad,
    size_t rlen,
    double* __restrict__ accm
) {
    for (size_t cb = 0; cb < Cpad; cb += CBLK) {
        double* __restrict__ ap = accm + cb * PBLK;
        for (size_t s = 0; s < rlen; s++) {
            const double* __restrict__ g = Gt + s * PBLK;
            const double* __restrict__ r = Mat + s * Cpad + cb;
            for (size_t c = 0; c < CBLK; c++) {
                const double rv = r[c];
                double* __restrict__ ac = ap + c * PBLK;
                for (size_t j = 0; j < PBLK; j++) ac[j] += g[j] * rv;
            }
        }
    }
}

#ifdef GMOQ_POOLED_HAVE_AVX2
#pragma GCC push_options
#pragma GCC target("avx,avx2,fma")
static void sweep_rows_avx2(
    const double* __restrict__ Gt,
    const double* __restrict__ Mat,
    size_t Cpad,
    size_t rlen,
    double* __restrict__ accm
) {
    for (size_t cb = 0; cb < Cpad; cb += CBLK) {
        double* __restrict__ ap = accm + cb * PBLK;
        __m256d a0 = _mm256_loadu_pd(ap +  0), a1 = _mm256_loadu_pd(ap +  4);
        __m256d a2 = _mm256_loadu_pd(ap +  8), a3 = _mm256_loadu_pd(ap + 12);
        __m256d a4 = _mm256_loadu_pd(ap + 16), a5 = _mm256_loadu_pd(ap + 20);
        __m256d a6 = _mm256_loadu_pd(ap + 24), a7 = _mm256_loadu_pd(ap + 28);
        __m256d a8 = _mm256_loadu_pd(ap + 32), a9 = _mm256_loadu_pd(ap + 36);
        __m256d aa = _mm256_loadu_pd(ap + 40), ab = _mm256_loadu_pd(ap + 44);
        for (size_t s = 0; s < rlen; s++) {
            const double* g = Gt + s * PBLK;
            const double* r = Mat + s * Cpad + cb;
            const __m256d g0 = _mm256_loadu_pd(g);
            const __m256d g1 = _mm256_loadu_pd(g + 4);
            __m256d rv;
            rv = _mm256_broadcast_sd(r + 0);
            a0 = _mm256_fmadd_pd(g0, rv, a0); a1 = _mm256_fmadd_pd(g1, rv, a1);
            rv = _mm256_broadcast_sd(r + 1);
            a2 = _mm256_fmadd_pd(g0, rv, a2); a3 = _mm256_fmadd_pd(g1, rv, a3);
            rv = _mm256_broadcast_sd(r + 2);
            a4 = _mm256_fmadd_pd(g0, rv, a4); a5 = _mm256_fmadd_pd(g1, rv, a5);
            rv = _mm256_broadcast_sd(r + 3);
            a6 = _mm256_fmadd_pd(g0, rv, a6); a7 = _mm256_fmadd_pd(g1, rv, a7);
            rv = _mm256_broadcast_sd(r + 4);
            a8 = _mm256_fmadd_pd(g0, rv, a8); a9 = _mm256_fmadd_pd(g1, rv, a9);
            rv = _mm256_broadcast_sd(r + 5);
            aa = _mm256_fmadd_pd(g0, rv, aa); ab = _mm256_fmadd_pd(g1, rv, ab);
        }
        _mm256_storeu_pd(ap +  0, a0); _mm256_storeu_pd(ap +  4, a1);
        _mm256_storeu_pd(ap +  8, a2); _mm256_storeu_pd(ap + 12, a3);
        _mm256_storeu_pd(ap + 16, a4); _mm256_storeu_pd(ap + 20, a5);
        _mm256_storeu_pd(ap + 24, a6); _mm256_storeu_pd(ap + 28, a7);
        _mm256_storeu_pd(ap + 32, a8); _mm256_storeu_pd(ap + 36, a9);
        _mm256_storeu_pd(ap + 40, aa); _mm256_storeu_pd(ap + 44, ab);
    }
}
#pragma GCC pop_options
#endif

inline bool use_avx2() {
#ifdef GMOQ_POOLED_HAVE_AVX2
    static const bool ok = __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
    return ok;
#else
    return false;
#endif
}

// ---------------------------------------------------------------------------
// PooledBinPlan: O(N) bin-prefix plan for --cov-mode pooled-only.
// ---------------------------------------------------------------------------

struct PooledBinPlan {
    size_t N = 0, T = 0, K = 0, Kpad = 0;
    bool has_intercept = false;
    double e_norm2 = 0.0;                  // N - ||Q^T 1||^2

    std::vector<size_t> tau_ord;           // sorted index -> original tau index
    std::vector<double> inv_f_sorted;      // 1 / fhat[tau_ord[st]]
    std::vector<double> c_sorted;          // q[st] + tau[st] / f[st]
    std::vector<uint16_t> bin;             // (N,) in [0, T]: y_i <= q_sorted[st] iff bin[i] <= st
    std::vector<double> bin_cnt;           // (T + 1,) counts per bin
    std::vector<double> cum_cnt;           // (T,) cumulative counts sum_{b<=st} bin_cnt[b]

    std::vector<double> w;                 // (K,) Q^T 1 (Neumaier compensated)
    std::vector<double> QtI_sorted;        // (K, T) Q^T I in sorted-tau order
    std::vector<double> QtR_sorted;        // (K, T) Q^T RIF in sorted-tau order
    std::vector<double> e_dot_R_sorted;    // (T,) e^T RIF in sorted-tau order
    std::vector<double> Q_pad;             // (N, Kpad) row-major (empty if K == 1 && has_intercept)
    std::vector<double> sigma0_mat;        // (T, T) row-major in original tau order

    PooledBinPlan(
        py::array_t<double, py::array::c_style | py::array::forcecast> y_arr,
        py::array_t<double, py::array::c_style | py::array::forcecast> taus_arr,
        py::array_t<double, py::array::c_style | py::array::forcecast> q_tau_arr,
        py::array_t<double, py::array::c_style | py::array::forcecast> fhat_arr,
        py::array_t<double, py::array::c_style | py::array::forcecast> Q_arr,
        int n_threads = 1
    ) {
        auto y_buf = y_arr.request();
        auto taus_buf = taus_arr.request();
        auto q_buf = q_tau_arr.request();
        auto f_buf = fhat_arr.request();
        auto Q_buf = Q_arr.request();

        if (y_buf.ndim != 1 || taus_buf.ndim != 1 || q_buf.ndim != 1 ||
            f_buf.ndim != 1 || Q_buf.ndim != 2) {
            throw std::runtime_error(
                "PooledBinPlan expects y (N,), taus (T,), q_tau (T,), fhat (T,), Q (N, K).");
        }
        N = static_cast<size_t>(y_buf.shape[0]);
        T = static_cast<size_t>(taus_buf.shape[0]);
        K = static_cast<size_t>(Q_buf.shape[1]);
        if (N < 2) throw std::runtime_error("PooledBinPlan requires N >= 2.");
        if (T < 1 || T >= 65535) throw std::runtime_error("PooledBinPlan requires 1 <= T < 65535.");
        if (static_cast<size_t>(q_buf.shape[0]) != T ||
            static_cast<size_t>(f_buf.shape[0]) != T ||
            static_cast<size_t>(Q_buf.shape[0]) != N) {
            throw std::runtime_error("PooledBinPlan dimension mismatch.");
        }

        const double* y_ptr = static_cast<const double*>(y_buf.ptr);
        const double* taus_ptr = static_cast<const double*>(taus_buf.ptr);
        const double* q_ptr = static_cast<const double*>(q_buf.ptr);
        const double* f_ptr = static_cast<const double*>(f_buf.ptr);
        const double* Q_ptr = static_cast<const double*>(Q_buf.ptr);

        // Sort tau indices by q_tau so lower_bound works even if caller passed unsorted taus.
        tau_ord.resize(T);
        std::iota(tau_ord.begin(), tau_ord.end(), 0);
        std::stable_sort(tau_ord.begin(), tau_ord.end(),
                         [q_ptr](size_t a, size_t b) { return q_ptr[a] < q_ptr[b]; });

        std::vector<double> q_sorted(T);
        inv_f_sorted.resize(T);
        c_sorted.resize(T);
        for (size_t st = 0; st < T; st++) {
            const size_t ot = tau_ord[st];
            if (!(f_ptr[ot] > 0.0) || !std::isfinite(f_ptr[ot])) {
                throw std::runtime_error("fhat values must be finite and positive.");
            }
            q_sorted[st] = q_ptr[ot];
            inv_f_sorted[st] = 1.0 / f_ptr[ot];
            c_sorted[st] = q_ptr[ot] + taus_ptr[ot] * inv_f_sorted[st];
        }

        bin.resize(N);
        bin_cnt.assign(T + 1, 0.0);
        cum_cnt.assign(T, 0.0);
        w.assign(K, 0.0);
        QtI_sorted.assign(K * T, 0.0);
        QtR_sorted.assign(K * T, 0.0);
        e_dot_R_sorted.assign(T, 0.0);

#ifdef _OPENMP
        if (n_threads > 0) omp_set_num_threads(n_threads);
#else
        (void)n_threads;
#endif

        {
        py::gil_scoped_release gil_release;

        // 1. Assign quantile bins in parallel.
        #pragma omp parallel for schedule(static)
        for (int64_t i = 0; i < static_cast<int64_t>(N); i++) {
            const double yi = y_ptr[i];
            const size_t b = static_cast<size_t>(
                std::lower_bound(q_sorted.begin(), q_sorted.end(), yi) - q_sorted.begin());
            bin[static_cast<size_t>(i)] = static_cast<uint16_t>(b);
        }

        // 2. Accumulate bin counts, Q_bin = sum_{i in bin b} Q[i,:], and w = Q^T 1 (Neumaier).
        std::vector<double> Q_bin((T + 1) * K, 0.0);
        std::vector<double> wc(K, 0.0);
        for (size_t i = 0; i < N; i++) {
            const size_t b = static_cast<size_t>(bin[i]);
            bin_cnt[b] += 1.0;
            const double* qrow = Q_ptr + i * K;
            double* qb = &Q_bin[b * K];
            for (size_t k = 0; k < K; k++) {
                const double x = qrow[k];
                qb[k] += x;
                const double yv = w[k], t = yv + x;
                wc[k] += (std::abs(yv) >= std::abs(x)) ? ((yv - t) + x) : ((x - t) + yv);
                w[k] = t;
            }
        }
        for (size_t k = 0; k < K; k++) w[k] += wc[k];

        double w2 = 0.0;
        for (size_t k = 0; k < K; k++) w2 += w[k] * w[k];
        e_norm2 = std::max(0.0, static_cast<double>(N) - w2);
        has_intercept = (e_norm2 <= 1e-10 * static_cast<double>(N));
        if (has_intercept) e_norm2 = 0.0;

        // Prefix sums for cum_cnt and QtI_sorted.
        double run_cnt = 0.0;
        std::vector<double> run_q(K, 0.0);
        for (size_t st = 0; st < T; st++) {
            run_cnt += bin_cnt[st];
            cum_cnt[st] = run_cnt;
            const double* qb = &Q_bin[st * K];
            for (size_t k = 0; k < K; k++) {
                run_q[k] += qb[k];
                QtI_sorted[k * T + st] = run_q[k];
            }
        }

        // QtR_sorted and e_dot_R_sorted.
        for (size_t st = 0; st < T; st++) {
            const double c_st = c_sorted[st];
            const double invf = inv_f_sorted[st];
            double w_qtr = 0.0;
            for (size_t k = 0; k < K; k++) {
                const double v = c_st * w[k] - QtI_sorted[k * T + st] * invf;
                QtR_sorted[k * T + st] = v;
                w_qtr += w[k] * v;
            }
            if (!has_intercept) {
                const double one_r = c_st * static_cast<double>(N) - cum_cnt[st] * invf;
                e_dot_R_sorted[st] = one_r - w_qtr;
            }
        }

        // 3. Exact Sigma_0 = (1/N) (M RIF)^T (M RIF) without allocating (N, T).
        sigma0_mat.assign(T * T, 0.0);
        const double invN = 1.0 / static_cast<double>(N);
        for (size_t s = 0; s < T; s++) {
            const size_t os = tau_ord[s];
            const double invfs = inv_f_sorted[s];
            const double cs = c_sorted[s];
            const double ns = cum_cnt[s];
            for (size_t t = s; t < T; t++) {
                const size_t ot = tau_ord[t];
                const double invft = inv_f_sorted[t];
                double val;
                if (has_intercept) {
                    // M * 1 = 0 => (M RIF)^T (M RIF)_{s,t} = (N_min - QtI_s . QtI_t) / (f_s f_t)
                    double qdot = 0.0;
                    for (size_t k = 0; k < K; k++) {
                        qdot += QtI_sorted[k * T + s] * QtI_sorted[k * T + t];
                    }
                    val = (ns - qdot) * (invfs * invft) * invN;
                } else {
                    const double ct = c_sorted[t];
                    const double nt = cum_cnt[t];
                    double rtr = static_cast<double>(N) * cs * ct
                               - cs * nt * invft - ct * ns * invfs
                               + ns * invfs * invft;
                    double qdot = 0.0;
                    for (size_t k = 0; k < K; k++) {
                        qdot += QtR_sorted[k * T + s] * QtR_sorted[k * T + t];
                    }
                    val = (rtr - qdot) * invN;
                }
                sigma0_mat[os * T + ot] = val;
                sigma0_mat[ot * T + os] = val;
            }
        }

        // 4. Store padded Q only when covariate projection is actually needed.
        if (K > 1 || !has_intercept) {
            Kpad = ((K + CBLK - 1) / CBLK) * CBLK;
            Q_pad.assign(N * Kpad, 0.0);
            #pragma omp parallel for schedule(static)
            for (int64_t i = 0; i < static_cast<int64_t>(N); i++) {
                std::memcpy(&Q_pad[static_cast<size_t>(i) * Kpad],
                            Q_ptr + static_cast<size_t>(i) * K,
                            K * sizeof(double));
            }
        } else {
            Kpad = 0;
        }

        }  // gil_release
    }

    py::array_t<double> sigma0() const {
        py::array_t<double> out({T, T});
        std::memcpy(out.request().ptr, sigma0_mat.data(), T * T * sizeof(double));
        return out;
    }
};

inline py::tuple pooled_bin_fast(
    py::array_t<double, py::array::f_style | py::array::forcecast> G_arr,
    const PooledBinPlan& plan,
    double min_mac,
    double min_var,
    int n_threads
) {
    auto G_buf = G_arr.request();
    if (G_buf.ndim != 2) throw std::runtime_error("G must be 2-D.");
    const size_t N = static_cast<size_t>(G_buf.shape[0]);
    const size_t M = static_cast<size_t>(G_buf.shape[1]);
    if (N != plan.N) {
        throw std::runtime_error("G has a different number of rows than PooledBinPlan.");
    }
    const size_t T = plan.T;
    const size_t K = plan.K;

    py::array_t<double, py::array::f_style> betas_arr({M, T});
    py::array_t<double> d_resid_arr(std::vector<size_t>{M});
    double* betas_ptr = static_cast<double*>(betas_arr.request().ptr);
    double* d_ptr = static_cast<double*>(d_resid_arr.request().ptr);

    const double nan_val = std::numeric_limits<double>::quiet_NaN();
    std::fill(betas_ptr, betas_ptr + M * T, nan_val);
    std::fill(d_ptr, d_ptr + M, nan_val);
    if (M == 0) return py::make_tuple(betas_arr, d_resid_arr);

    const double* __restrict__ G_ptr = static_cast<const double*>(G_buf.ptr);
    const uint16_t* __restrict__ bin_ptr = plan.bin.data();
    const double* __restrict__ bin_cnt = plan.bin_cnt.data();
    const size_t* __restrict__ tau_ord = plan.tau_ord.data();
    const double* __restrict__ inv_f = plan.inv_f_sorted.data();
    const double* __restrict__ QtR = plan.QtR_sorted.data();
    const double* __restrict__ w = plan.w.data();
    const double* __restrict__ e_dot_R = plan.e_dot_R_sorted.data();
    const bool has_intercept = plan.has_intercept;
    const bool need_q = (K > 1 || !has_intercept);
    const bool avx2 = use_avx2();

#ifdef _OPENMP
    if (n_threads > 0) omp_set_num_threads(n_threads);
#else
    (void)n_threads;
#endif

    {
    py::gil_scoped_release gil_release;

    std::vector<char> keep(M, 0);
    std::vector<double> gmean_all(M, 0.0);
    std::vector<double> gc2_all(M, 0.0);
    // gc_dot_I_all stored in sorted-tau order only when Pass 2 is needed.
    std::vector<double> gc_dot_I_all(need_q ? (M * T) : 0, 0.0);

    // Pass 1: single sequential scan per SNP for MAC filter, mean, shifted ||g_c||^2,
    // and O(N) bin-prefix sums g_c^T I[:, st] for all T quantiles at once.
    #pragma omp parallel
    {
        std::vector<double> bsum(T + 1);
        std::vector<uint32_t> bmiss(T + 1);

        #pragma omp for schedule(static)
        for (int64_t j = 0; j < static_cast<int64_t>(M); j++) {
            const size_t jj = static_cast<size_t>(j);
            const double* __restrict__ col = G_ptr + jj * N;

            // Find first finite value for shift-stable one-pass variance.
            size_t first = 0;
            while (first < N && !std::isfinite(col[first])) first++;
            if (N - first < 100) continue;
            const double g0 = col[first];

            std::fill(bsum.begin(), bsum.end(), 0.0);
            std::fill(bmiss.begin(), bmiss.end(), 0u);
            for (size_t i = 0; i < first; i++) bmiss[bin_ptr[i]]++;

            size_t n_valid = 0;
            double d_sum = 0.0;
            double d_sum2 = 0.0;
            for (size_t i = first; i < N; i++) {
                const double v = col[i];
                const size_t b = static_cast<size_t>(bin_ptr[i]);
                if (std::isfinite(v)) {
                    const double dv = v - g0;
                    d_sum += dv;
                    d_sum2 += dv * dv;
                    bsum[b] += dv;
                    n_valid++;
                } else {
                    bmiss[b]++;
                }
            }
            if (n_valid < 100) continue;

            const double n_v = static_cast<double>(n_valid);
            const double d_mean = d_sum / n_v;
            const double g_mean = g0 + d_mean;
            const double g_sum = g_mean * n_v;
            const double mac = std::min(g_sum, 2.0 * n_v - g_sum);
            if (mac < min_mac) continue;

            const double gc2 = std::max(0.0, d_sum2 - n_v * d_mean * d_mean);

            if (!need_q) {
                // K == 1 with intercept: Q^T g_c == 0 identically; done in one pass!
                if (!(gc2 >= min_var)) continue;
                const double invD = 1.0 / gc2;
                double prefix = 0.0;
                for (size_t st = 0; st < T; st++) {
                    const double cnt_valid = bin_cnt[st] - static_cast<double>(bmiss[st]);
                    prefix += bsum[st] - cnt_valid * d_mean;
                    betas_ptr[jj + M * tau_ord[st]] = (-prefix * inv_f[st]) * invD;
                }
                d_ptr[jj] = gc2;
            } else {
                keep[jj] = 1;
                gmean_all[jj] = g_mean;
                gc2_all[jj] = gc2;
                double* __restrict__ out_I = &gc_dot_I_all[jj * T];
                double prefix = 0.0;
                for (size_t st = 0; st < T; st++) {
                    const double cnt_valid = bin_cnt[st] - static_cast<double>(bmiss[st]);
                    prefix += bsum[st] - cnt_valid * d_mean;
                    out_I[st] = prefix;
                }
            }
        }
    }

    if (need_q) {
        std::vector<size_t> kept;
        kept.reserve(M);
        for (size_t j = 0; j < M; j++) if (keep[j]) kept.push_back(j);

        const size_t n_kept = kept.size();
        const size_t n_panels = (n_kept + PBLK - 1) / PBLK;
        const size_t Kpad = plan.Kpad;
        const double* __restrict__ Q_pad = plan.Q_pad.data();

        // Pass 2: compute QtGc = Q^T G_c in L2-resident row tiles, parallel over 8-SNP panels.
        #pragma omp parallel
        {
            std::vector<double> gt(RBLK * PBLK, 0.0);
            std::vector<double> accm(Kpad * PBLK, 0.0);

            #pragma omp for schedule(static)
            for (int64_t ip = 0; ip < static_cast<int64_t>(n_panels); ip++) {
                const size_t p0 = static_cast<size_t>(ip) * PBLK;
                const size_t np = std::min<size_t>(PBLK, n_kept - p0);
                size_t cols_p[PBLK];
                double gm_p[PBLK];
                for (size_t p = 0; p < np; p++) {
                    cols_p[p] = kept[p0 + p];
                    gm_p[p] = gmean_all[cols_p[p]];
                }

                std::fill(accm.begin(), accm.end(), 0.0);

                for (size_t i0 = 0; i0 < N; i0 += RBLK) {
                    const size_t rlen = std::min<size_t>(RBLK, N - i0);
                    // Zero padding columns once if np < PBLK.
                    if (np < PBLK) {
                        std::fill(gt.begin(), gt.begin() + rlen * PBLK, 0.0);
                    }
                    for (size_t p = 0; p < np; p++) {
                        const double* __restrict__ col = G_ptr + cols_p[p] * N + i0;
                        const double gm = gm_p[p];
                        double* __restrict__ dst = gt.data() + p;
                        for (size_t r = 0; r < rlen; r++) {
                            const double v = col[r];
                            dst[r * PBLK] = std::isfinite(v) ? (v - gm) : 0.0;
                        }
                    }
                    const double* __restrict__ qchunk = Q_pad + i0 * Kpad;
                    if (avx2) {
#ifdef GMOQ_POOLED_HAVE_AVX2
                        sweep_rows_avx2(gt.data(), qchunk, Kpad, rlen, accm.data());
#endif
                    } else {
                        sweep_rows_base(gt.data(), qchunk, Kpad, rlen, accm.data());
                    }
                }

                // Finalise each SNP in the 8-SNP panel.
                for (size_t p = 0; p < np; p++) {
                    const size_t jj = cols_p[p];
                    double q2 = 0.0;
                    double w_qtg = 0.0;
                    for (size_t k = 0; k < K; k++) {
                        const double qk = accm[k * PBLK + p];
                        q2 += qk * qk;
                        if (!has_intercept) w_qtg += w[k] * qk;
                    }
                    double D = gc2_all[jj] - q2;
                    if (!has_intercept) {
                        const double gm = gm_p[p];
                        D += -2.0 * gm * w_qtg + gm * gm * plan.e_norm2;
                    }
                    if (!(D >= min_var)) continue;

                    const double invD = 1.0 / D;
                    const double* __restrict__ gcI = &gc_dot_I_all[jj * T];
                    const double gm = gm_p[p];
                    for (size_t st = 0; st < T; st++) {
                        double num = -gcI[st] * inv_f[st];
                        for (size_t k = 0; k < K; k++) {
                            num -= accm[k * PBLK + p] * QtR[k * T + st];
                        }
                        if (!has_intercept) {
                            num += gm * e_dot_R[st];
                        }
                        betas_ptr[jj + M * tau_ord[st]] = num * invD;
                    }
                    d_ptr[jj] = D;
                }
            }
        }
    }

    }  // gil_release

    return py::make_tuple(betas_arr, d_resid_arr);
}

// ---------------------------------------------------------------------------
// Hoisted, register-blocked pooled panel kernel for callers passing (G, RIF_resid, Q).
// Replaces the scalar triple loop in stage1_pooled_panel_impl with the same
// FWL covariate hoist + L2 row-tiled AVX2/portable sweep.
// ---------------------------------------------------------------------------

inline py::tuple pooled_panel_fast(
    py::array_t<double, py::array::f_style | py::array::forcecast> G_arr,
    py::array_t<double, py::array::c_style | py::array::forcecast> RIF_resid_arr,
    py::array_t<double, py::array::c_style | py::array::forcecast> Q_arr,
    double min_mac,
    double min_var,
    int n_threads,
    int panel_size
) {
    (void)panel_size;
    auto G_buf = G_arr.request();
    auto RIF_buf = RIF_resid_arr.request();
    auto Q_buf = Q_arr.request();

    if (G_buf.ndim != 2 || RIF_buf.ndim != 2 || Q_buf.ndim != 2) {
        throw std::runtime_error("G, RIF_resid and Q must all be 2-D.");
    }
    const size_t N = static_cast<size_t>(G_buf.shape[0]);
    const size_t M = static_cast<size_t>(G_buf.shape[1]);
    const size_t T = static_cast<size_t>(RIF_buf.shape[1]);
    const size_t K = static_cast<size_t>(Q_buf.shape[1]);
    if (static_cast<size_t>(RIF_buf.shape[0]) != N || static_cast<size_t>(Q_buf.shape[0]) != N) {
        throw std::runtime_error("RIF_resid and Q must have the same number of rows as G.");
    }

    py::array_t<double, py::array::f_style> betas_arr({M, T});
    py::array_t<double> d_resid_arr(std::vector<size_t>{M});
    double* betas_ptr = static_cast<double*>(betas_arr.request().ptr);
    double* d_ptr = static_cast<double*>(d_resid_arr.request().ptr);

    const double nan_val = std::numeric_limits<double>::quiet_NaN();
    std::fill(betas_ptr, betas_ptr + M * T, nan_val);
    std::fill(d_ptr, d_ptr + M, nan_val);
    if (M == 0) return py::make_tuple(betas_arr, d_resid_arr);

    const double* __restrict__ G_ptr = static_cast<const double*>(G_buf.ptr);
    const double* __restrict__ RIF_ptr = static_cast<const double*>(RIF_buf.ptr);
    const double* __restrict__ Q_ptr = static_cast<const double*>(Q_buf.ptr);
    const bool avx2 = use_avx2();

#ifdef _OPENMP
    if (n_threads > 0) omp_set_num_threads(n_threads);
#else
    (void)n_threads;
#endif

    {
    py::gil_scoped_release gil_release;

    // Pass 1: MAC filter, mean, shifted ||g_c||^2.
    std::vector<char> keep(M, 0);
    std::vector<double> gmean_all(M, 0.0);
    std::vector<double> gc2_all(M, 0.0);

    #pragma omp parallel for schedule(static)
    for (int64_t j = 0; j < static_cast<int64_t>(M); j++) {
        const size_t jj = static_cast<size_t>(j);
        const double* __restrict__ col = G_ptr + jj * N;
        size_t first = 0;
        while (first < N && !std::isfinite(col[first])) first++;
        if (N - first < 100) continue;
        const double g0 = col[first];
        size_t n_valid = 0;
        double d_sum = 0.0, d_sum2 = 0.0;
        for (size_t i = first; i < N; i++) {
            const double v = col[i];
            if (std::isfinite(v)) {
                const double dv = v - g0;
                d_sum += dv;
                d_sum2 += dv * dv;
                n_valid++;
            }
        }
        if (n_valid < 100) continue;
        const double n_v = static_cast<double>(n_valid);
        const double d_mean = d_sum / n_v;
        const double g_mean = g0 + d_mean;
        const double g_sum = g_mean * n_v;
        const double mac = std::min(g_sum, 2.0 * n_v - g_sum);
        if (mac < min_mac) continue;
        keep[jj] = 1;
        gmean_all[jj] = g_mean;
        gc2_all[jj] = std::max(0.0, d_sum2 - n_v * d_mean * d_mean);
    }

    std::vector<size_t> kept;
    kept.reserve(M);
    for (size_t j = 0; j < M; j++) if (keep[j]) kept.push_back(j);
    const size_t n_kept = kept.size();
    if (n_kept > 0) {

    // Precompute w = Q^T 1, QtR = Q^T RIF_resid, and e^T RIF_resid.
    // (In standard calls RIF_resid is already orthogonal to Q and 1 in col(Q),
    // so these correction terms are ~0, but computing them once keeps arbitrary
    // callers exact to rounding.)
    std::vector<double> w(K, 0.0), wc(K, 0.0);
    std::vector<double> one_R(T, 0.0);
    std::vector<double> QtR(K * T, 0.0);
    for (size_t i = 0; i < N; i++) {
        const double* qrow = Q_ptr + i * K;
        const double* rrow = RIF_ptr + i * T;
        for (size_t k = 0; k < K; k++) {
            const double x = qrow[k], yv = w[k], t = yv + x;
            wc[k] += (std::abs(yv) >= std::abs(x)) ? ((yv - t) + x) : ((x - t) + yv);
            w[k] = t;
            double* qk = &QtR[k * T];
            for (size_t tt = 0; tt < T; tt++) qk[tt] += x * rrow[tt];
        }
        for (size_t tt = 0; tt < T; tt++) one_R[tt] += rrow[tt];
    }
    for (size_t k = 0; k < K; k++) w[k] += wc[k];
    double w2 = 0.0;
    for (size_t k = 0; k < K; k++) w2 += w[k] * w[k];
    const double e_norm2 = std::max(0.0, static_cast<double>(N) - w2);
    const bool has_intercept = (e_norm2 <= 1e-10 * static_cast<double>(N));
    std::vector<double> e_dot_R(T, 0.0);
    if (!has_intercept) {
        for (size_t tt = 0; tt < T; tt++) {
            double wq = 0.0;
            for (size_t k = 0; k < K; k++) wq += w[k] * QtR[k * T + tt];
            e_dot_R[tt] = one_R[tt] - wq;
        }
    }

    // Build padded row-major [RIF_resid | Q] chunked per row-tile (no full N x (T+K) copy!).
    const size_t C = T + K;
    const size_t Cpad = ((C + CBLK - 1) / CBLK) * CBLK;
    const size_t n_panels = (n_kept + PBLK - 1) / PBLK;

    // Parallel over 8-SNP panels; each thread packs small L2-resident (RBLK x Cpad)
    // and (RBLK x 8) tiles per chunk.
    #pragma omp parallel
    {
        std::vector<double> gt(RBLK * PBLK, 0.0);
        std::vector<double> rq_tile(RBLK * Cpad, 0.0);
        std::vector<double> accm(Cpad * PBLK, 0.0);

        #pragma omp for schedule(static)
        for (int64_t ip = 0; ip < static_cast<int64_t>(n_panels); ip++) {
            const size_t p0 = static_cast<size_t>(ip) * PBLK;
            const size_t np = std::min<size_t>(PBLK, n_kept - p0);
            size_t cols_p[PBLK];
            double gm_p[PBLK];
            for (size_t p = 0; p < np; p++) {
                cols_p[p] = kept[p0 + p];
                gm_p[p] = gmean_all[cols_p[p]];
            }
            std::fill(accm.begin(), accm.end(), 0.0);

            for (size_t i0 = 0; i0 < N; i0 += RBLK) {
                const size_t rlen = std::min<size_t>(RBLK, N - i0);
                if (np < PBLK) std::fill(gt.begin(), gt.begin() + rlen * PBLK, 0.0);
                for (size_t p = 0; p < np; p++) {
                    const double* __restrict__ col = G_ptr + cols_p[p] * N + i0;
                    const double gm = gm_p[p];
                    double* __restrict__ dst = gt.data() + p;
                    for (size_t r = 0; r < rlen; r++) {
                        const double v = col[r];
                        dst[r * PBLK] = std::isfinite(v) ? (v - gm) : 0.0;
                    }
                }
                for (size_t r = 0; r < rlen; r++) {
                    double* dst = &rq_tile[r * Cpad];
                    std::memcpy(dst, RIF_ptr + (i0 + r) * T, T * sizeof(double));
                    std::memcpy(dst + T, Q_ptr + (i0 + r) * K, K * sizeof(double));
                    for (size_t c = C; c < Cpad; c++) dst[c] = 0.0;
                }
                if (avx2) {
#ifdef GMOQ_POOLED_HAVE_AVX2
                    sweep_rows_avx2(gt.data(), rq_tile.data(), Cpad, rlen, accm.data());
#endif
                } else {
                    sweep_rows_base(gt.data(), rq_tile.data(), Cpad, rlen, accm.data());
                }
            }

            for (size_t p = 0; p < np; p++) {
                const size_t jj = cols_p[p];
                double q2 = 0.0, w_qtg = 0.0;
                for (size_t k = 0; k < K; k++) {
                    const double qk = accm[(T + k) * PBLK + p];
                    q2 += qk * qk;
                    if (!has_intercept) w_qtg += w[k] * qk;
                }
                double D = gc2_all[jj] - q2;
                if (!has_intercept) {
                    const double gm = gm_p[p];
                    D += -2.0 * gm * w_qtg + gm * gm * e_norm2;
                }
                if (!(D >= min_var)) continue;
                const double invD = 1.0 / D;
                const double gm = gm_p[p];
                for (size_t tt = 0; tt < T; tt++) {
                    double num = accm[tt * PBLK + p];
                    for (size_t k = 0; k < K; k++) {
                        num -= accm[(T + k) * PBLK + p] * QtR[k * T + tt];
                    }
                    if (!has_intercept) num += gm * e_dot_R[tt];
                    betas_ptr[jj + M * tt] = num * invD;
                }
                d_ptr[jj] = D;
            }
        }
    }

    }  // if (n_kept > 0)

    }  // gil_release

    return py::make_tuple(betas_arr, d_resid_arr);
}

}  // namespace gmoq_pooled

#endif  // GMOQ_STAGE1_POOLED_FAST_HPP
