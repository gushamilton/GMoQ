// Stage 1 fast panel kernel.
//
// ---------------------------------------------------------------------------
// Why this exists
// ---------------------------------------------------------------------------
// The exact panel kernel (panel_block_stats in stage1_pybind.cpp) is a pure
// re-nesting of the original per-SNP loops, so it is bit-for-bit identical to
// them - and it inherits their access pattern.  On a biobank-scale chromosome
// (N around 400,000, T = 45, K = 13, 1000 blocks) it dominates the Stage 1 run
// time and reaches only a small fraction of the SSE2 peak.  Three things hold
// it there:
//
//   * the genotype gather ``G_ptr[i + cols[p]*N]`` runs down the long axis and is
//     performed twice per (sample, SNP);
//   * the K = 13 covariate projection ``proj = sum_k qrow[k]*qt[k]`` is a serial
//     reduction recomputed for every (sample, SNP) pair;
//   * every sample writes into the block tile its block id selects, and with
//     i.i.d. uniform block ids that is an unrelated 5.8 kB tile per sample, so
//     the write-hot set is n_blocks * P * (1+T) * 8 bytes per thread and blows
//     past L3 as soon as more than a couple of threads run.
//
// ---------------------------------------------------------------------------
// What this kernel does instead
// ---------------------------------------------------------------------------
// Write r_ip = g_ip - Q_i . q_p for the residualised genotype, where
// q_p = Q^T g_p.  The block statistics the jackknife needs are
//
//   Out[b][p][0]   = sum_{i in b} r_ip^2
//   Out[b][p][1+t] = sum_{i in b} r_ip * RIF[i,t]
//
// Expanding r and collecting the terms that do not involve the genotype gives
//
//   Out[b][p][1+t] = W_b[p][t]  -  sum_k q_p[k] * S_b[k][t]
//   Out[b][p][0]   = s2_b[p] - 2 * sum_k q_p[k] * u_b[k][p]
//                            + sum_{k,l} q_p[k] A_b[k][l] q_p[l]
//
// with the per-block primitives
//
//   W_b[p][t] = sum_{i in b} g_ip RIF[i,t]     (the only O(N P T) term)
//   u_b[k][p] = sum_{i in b} Q[i,k] g_ip       (and q_p = sum_b u_b[:,p])
//   s2_b[p]   = sum_{i in b} g_ip^2
//   S_b[k][t] = sum_{i in b} Q[i,k] RIF[i,t]   <- SNP-independent, precomputed
//   A_b[k][l] = sum_{i in b} Q[i,k] Q[i,l]     <- SNP-independent, precomputed
//
// S and A depend only on Q, RIF and the block assignment, so they are computed
// once per Stage 1 run (BlockPlan) rather than once per (sample, SNP).  That is
// the hoist of the covariate projection: the serial K-length reduction that used
// to run N*P times per panel now runs n_blocks*P times, 400x less often at
// production N.
//
// The samples are counting-sorted by block id once, so each block owns a
// contiguous range of rows in RIF, Q and the packed genotype tile.  The three
// per-block primitives are then one fused streaming sweep over that range:
// sequential in every operand, one genotype read per (sample, SNP), and the
// accumulators are a small tile that stays in L1 - the scattered 5.8 kB-tile
// writes are gone entirely, which is why this kernel scales with thread count
// where the old one did not.
//
// ---------------------------------------------------------------------------
// Numerics
// ---------------------------------------------------------------------------
// This is NOT bit-identical to the exact kernel: the same quantity is formed
// from a different (algebraically equal) expression, so the rounding differs.
// The reassociation is benign because the cancellation is mild - at production
// scale |W_b| and |sum_k q S_b| are both O(10) while their difference is O(14),
// so roughly one decimal digit is lost, not ten.  Summation order is otherwise
// kept as close as the restructuring allows: within a block the samples are
// visited in ascending original row index (the counting sort is stable), which
// is exactly the order the exact kernel adds them in.
//
// Every result is independent of the thread count, of the panel width and of
// the batch size: a block's statistics depend only on that block's samples, and
// a SNP's statistics depend only on that SNP's column.
//
// The exact kernel remains available and untouched (--exact-kernel).

#ifndef GMOQ_STAGE1_FAST_HPP
#define GMOQ_STAGE1_FAST_HPP

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <new>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <numeric>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(__GNUC__) && defined(__x86_64__) && !defined(GMOQ_NO_AVX2_CLONE)
#define GMOQ_HAVE_AVX2_CLONE 1
#include <immintrin.h>
#endif

namespace gmoq_fast {

namespace py = pybind11;

// Register-blocking of the streaming sweep: PBLK genotype columns by CBLK
// (RIF | Q) columns.  PBLK*CBLK/4 = 8 ymm accumulators, one PBLK-wide genotype
// load and CBLK broadcasts per sample.
enum { PBLK = 8, CBLK = 6 };

// ---------------------------------------------------------------------------
// The hot leaves, compiled twice: once at the translation unit's baseline ISA
// and once for AVX2 + FMA, dispatched on __builtin_cpu_supports at first use.
// Raising the ISA for the whole TU would change the exact kernel's arithmetic
// (GCC contracts a*b+c into an FMA), so it is raised only here.
// ---------------------------------------------------------------------------

#define GMOQ_SUFFIX _base
#include "stage1_fast_hot.inc"
#undef GMOQ_SUFFIX

#ifdef GMOQ_HAVE_AVX2_CLONE
static_assert(PBLK == 8 && CBLK == 6, "the AVX2 sweep is written for PBLK=8, CBLK=6");
#pragma GCC push_options
#pragma GCC target("avx,avx2,fma")
#define GMOQ_SUFFIX _avx2
#define GMOQ_INTRIN 1
#include "stage1_fast_hot.inc"
#undef GMOQ_INTRIN
#undef GMOQ_SUFFIX
#pragma GCC pop_options
#endif

inline bool use_avx2() {
#ifdef GMOQ_HAVE_AVX2_CLONE
    static const bool ok = __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
    return ok;
#else
    return false;
#endif
}

// Uninitialised, 64-byte aligned scratch that grows monotonically and is reused
// across calls.  The buffers are fully overwritten before they are read, so
// zero-filling them (as std::vector would) is 300+ MB of pure waste per batch at
// production shapes - and reallocating them per batch would pay 300+ MB of
// first-touch page faults every time.  Alignment lets the tile packer use
// non-temporal stores, which halves its write traffic.
struct Scratch {
    double* p = nullptr;
    size_t n = 0;
    ~Scratch() { if (p) std::free(p); }
    double* get(size_t need) {
        if (need > n) {
            if (p) std::free(p);
            const size_t bytes = ((need * sizeof(double)) + 63) & ~static_cast<size_t>(63);
            void* q = nullptr;
            if (posix_memalign(&q, 64, bytes) != 0) throw std::bad_alloc();
            p = static_cast<double*>(q);
            n = bytes / sizeof(double);
        }
        return p;
    }
};

// ---------------------------------------------------------------------------
// BlockPlan: everything that depends on Q, RIF and the block assignment but not
// on the genotypes.  Built once per Stage 1 run and reused for every batch.
// ---------------------------------------------------------------------------

struct BlockPlan {
    size_t N = 0, T = 0, K = 0, Crow = 0;
    int n_blocks = 0;
    std::vector<int32_t> order;    // sorted position -> original sample row
    std::vector<int32_t> pos;      // original sample row -> sorted position
    std::vector<int64_t> offsets;  // n_blocks + 2 entries; [n_blocks] .. [n_blocks+1]
                                   // is the bucket of samples whose block id is out
                                   // of range (they contribute to q_p and to the
                                   // min_var gate, but to no block statistic)
    std::vector<double> RQ;        // (N, Crow) row-major in sorted order: columns
                                   // [0, T) are RIF, [T, T+K) are Q, column T+K is
                                   // e = 1 - Q Q^T 1, the rest is zero padding to a
                                   // multiple of CBLK so the sweep needs no
                                   // remainder branch
    std::vector<double> S;         // (n_blocks + 1, K, T)   Q_b^T RIF_b
    std::vector<double> A;         // (n_blocks + 1, K, K)   Q_b^T Q_b
    std::vector<double> E;         // (n_blocks + 1, T)      e_b^T RIF_b
    std::vector<double> F;         // (n_blocks + 1, K)      Q_b^T e_b
    std::vector<double> H;         // (n_blocks + 1)         e_b^T e_b

    // Per-batch working buffers.  They depend only on the panel shape, which is
    // constant across a run, so they are carried on the plan and allocated once.
    // Consequence: one plan serves one kernel call at a time.  The CLI's batch
    // pipeline runs a single batch at a time; concurrent callers need a plan each.
    mutable Scratch gt_buf, w_buf, u_buf, s2_buf, v_buf;

    // Optional bin-prefix data (set_bins).  RIF[i,t] = c_t - 1{y_i <= q_t}/f_t,
    // so each sample falls in one quantile bin, and the O(N P T) cross-product
    // W_b[p][t] = sum_{i in b} gc_ip RIF_resid[i,t] is recovered exactly from
    // per-(block, bin) genotype sums plus the U_b the sweep already forms:
    //
    //   W_b[p][t] = c_t * sum_{i in b} gc_ip
    //             - (1/f_t) * sum_{i in b, bin_i <= t} gc_ip
    //             - sum_k U_b[p][k] * (Q^T RIF)[k][t]
    //
    // The sweep then runs over [Q | e] only (K+1 columns instead of T+K+1).
    bool has_bins = false;
    size_t Crow2 = 0;
    std::vector<double> RQ2;           // (N, Crow2) [Q | e | 0 pad], sorted order
    std::vector<uint16_t> bin_sorted;  // (N,) bin of each sorted sample, in [0, T]
    std::vector<size_t> tau_ord;       // sorted-q position -> tau index
    std::vector<double> c_sorted;      // q + tau/f, sorted-q order
    std::vector<double> inv_f_sorted;  // 1/f, sorted-q order
    std::vector<double> QtR_tk;        // (T, K): (Q^T RIF)[k][t] stored [t*K + k]

    void set_bins(py::array_t<double, py::array::c_style | py::array::forcecast> y_arr,
                  py::array_t<double, py::array::c_style | py::array::forcecast> taus_arr,
                  py::array_t<double, py::array::c_style | py::array::forcecast> q_arr,
                  py::array_t<double, py::array::c_style | py::array::forcecast> f_arr);

    size_t n_buckets() const { return static_cast<size_t>(n_blocks) + 1; }

    BlockPlan(py::array_t<double, py::array::c_style | py::array::forcecast> RIF_arr,
              py::array_t<double, py::array::c_style | py::array::forcecast> Q_arr,
              py::array_t<int32_t> block_ids_arr,
              int n_blocks_,
              int n_threads);
};

inline BlockPlan::BlockPlan(
    py::array_t<double, py::array::c_style | py::array::forcecast> RIF_arr,
    py::array_t<double, py::array::c_style | py::array::forcecast> Q_arr,
    py::array_t<int32_t> block_ids_arr,
    int n_blocks_,
    int n_threads
) {
    auto RIF_buf = RIF_arr.request();
    auto Q_buf = Q_arr.request();
    auto block_buf = block_ids_arr.request();
    if (RIF_buf.ndim != 2 || Q_buf.ndim != 2 || block_buf.ndim != 1) {
        throw std::runtime_error("BlockPlan expects RIF (N, T), Q (N, K) and block_ids (N,).");
    }
    N = static_cast<size_t>(RIF_buf.shape[0]);
    T = static_cast<size_t>(RIF_buf.shape[1]);
    K = static_cast<size_t>(Q_buf.shape[1]);
    n_blocks = n_blocks_;
    if (static_cast<size_t>(Q_buf.shape[0]) != N || static_cast<size_t>(block_buf.shape[0]) != N) {
        throw std::runtime_error("RIF, Q and block_ids must all have N rows.");
    }
    if (n_blocks < 1) throw std::runtime_error("n_blocks must be >= 1.");

    const double* RIF_ptr = static_cast<const double*>(RIF_buf.ptr);
    const double* Q_ptr = static_cast<const double*>(Q_buf.ptr);
    const int32_t* block_ptr = static_cast<const int32_t*>(block_buf.ptr);

    const size_t NB = n_buckets();

    // Stable counting sort of the samples by block id.
    std::vector<int64_t> counts(NB, 0);
    for (size_t i = 0; i < N; i++) {
        const int b = block_ptr[i];
        const size_t bucket = (b >= 0 && b < n_blocks) ? static_cast<size_t>(b)
                                                       : static_cast<size_t>(n_blocks);
        counts[bucket]++;
    }
    offsets.assign(NB + 1, 0);
    int64_t running = 0;
    for (size_t b = 0; b < NB; b++) { offsets[b] = running; running += counts[b]; }
    offsets[NB] = running;

    std::vector<int64_t> cursor(offsets.begin(), offsets.begin() + NB);
    order.assign(N, 0);
    pos.assign(N, 0);
    for (size_t i = 0; i < N; i++) {
        const int b = block_ptr[i];
        const size_t bucket = (b >= 0 && b < n_blocks) ? static_cast<size_t>(b)
                                                       : static_cast<size_t>(n_blocks);
        const int64_t s = cursor[bucket]++;
        order[static_cast<size_t>(s)] = static_cast<int32_t>(i);
        pos[i] = static_cast<int32_t>(s);
    }

    Crow = ((T + K + 1 + CBLK - 1) / CBLK) * CBLK;
    RQ.assign(N * Crow, 0.0);
    S.assign(NB * K * T, 0.0);
    A.assign(NB * K * K, 0.0);
    E.assign(NB * T, 0.0);
    F.assign(NB * K, 0.0);
    H.assign(NB, 0.0);

#ifdef _OPENMP
    if (n_threads > 0) omp_set_num_threads(n_threads);
#else
    (void)n_threads;
#endif

    {
    py::gil_scoped_release gil_release;

    #pragma omp parallel for schedule(static)
    for (int64_t s = 0; s < static_cast<int64_t>(N); s++) {
        const size_t src = static_cast<size_t>(order[static_cast<size_t>(s)]);
        double* dst = &RQ[static_cast<size_t>(s) * Crow];
        std::memcpy(dst, RIF_ptr + src * T, T * sizeof(double));
        std::memcpy(dst + T, Q_ptr + src * K, K * sizeof(double));
    }

    // e = 1 - Q (Q^T 1) is the part of the intercept the basis does not span.
    // The kernel works on the mean-centred genotype, which is exact only when the
    // intercept lies in col(Q); e carries the correction for when it does not, and
    // is identically zero (to rounding) in every configuration the CLI builds,
    // where Q comes from a QR of [1 | covariates].
    {
        // w = Q^T 1, by compensated (Neumaier) summation.  A plain sum is not
        // good enough here: the intercept column of Q is N copies of 1/sqrt(N),
        // so a naive accumulation of it carries an O(N*eps) relative error -
        // 4e-12 at production N - and that error lands directly in e, whose true
        // value is zero.  The correction term would then inject more noise than
        // the centring removes.
        std::vector<double> w(K, 0.0), wc(K, 0.0);
        for (size_t s = 0; s < N; s++) {
            const double* q = &RQ[s * Crow + T];
            for (size_t k = 0; k < K; k++) {
                const double x = q[k], y = w[k], t = y + x;
                wc[k] += (std::abs(y) >= std::abs(x)) ? ((y - t) + x) : ((x - t) + y);
                w[k] = t;
            }
        }
        for (size_t k = 0; k < K; k++) w[k] += wc[k];
        #pragma omp parallel for schedule(static)
        for (int64_t s = 0; s < static_cast<int64_t>(N); s++) {
            const double* q = &RQ[static_cast<size_t>(s) * Crow + T];
            double proj = 0.0;
            for (size_t k = 0; k < K; k++) proj += q[k] * w[k];
            RQ[static_cast<size_t>(s) * Crow + T + K] = 1.0 - proj;
        }
    }

    // S_b = Q_b^T RIF_b and A_b = Q_b^T Q_b, accumulated over the block's samples
    // in ascending original row index.
    #pragma omp parallel for schedule(static)
    for (int64_t b = 0; b < static_cast<int64_t>(NB); b++) {
        const size_t bb = static_cast<size_t>(b);
        const size_t s0 = static_cast<size_t>(offsets[bb]);
        const size_t s1 = static_cast<size_t>(offsets[bb + 1]);
        double* Sb = &S[bb * K * T];
        double* Ab = &A[bb * K * K];
        double* Eb = &E[bb * T];
        double* Fb = &F[bb * K];
        double Hb = 0.0;
        for (size_t s = s0; s < s1; s++) {
            const double* r = &RQ[s * Crow];
            const double* q = r + T;
            const double ev = q[K];
            for (size_t k = 0; k < K; k++) {
                const double qk = q[k];
                double* Sk = Sb + k * T;
                for (size_t t = 0; t < T; t++) Sk[t] += qk * r[t];
                double* Ak = Ab + k * K;
                for (size_t l = 0; l < K; l++) Ak[l] += qk * q[l];
                Fb[k] += qk * ev;
            }
            for (size_t t = 0; t < T; t++) Eb[t] += ev * r[t];
            Hb += ev * ev;
        }
        H[bb] = Hb;
    }

    }  // gil_release
}

inline void BlockPlan::set_bins(
    py::array_t<double, py::array::c_style | py::array::forcecast> y_arr,
    py::array_t<double, py::array::c_style | py::array::forcecast> taus_arr,
    py::array_t<double, py::array::c_style | py::array::forcecast> q_arr,
    py::array_t<double, py::array::c_style | py::array::forcecast> f_arr
) {
    if (y_arr.ndim() != 1 || static_cast<size_t>(y_arr.shape(0)) != N) {
        throw std::runtime_error("set_bins: y must have N entries.");
    }
    if (taus_arr.shape(0) != static_cast<py::ssize_t>(T) || q_arr.shape(0) != static_cast<py::ssize_t>(T) ||
        f_arr.shape(0) != static_cast<py::ssize_t>(T)) {
        throw std::runtime_error("set_bins: taus, q_tau and fhat must have T entries.");
    }
    if (T >= 65535) throw std::runtime_error("set_bins requires T < 65535.");
    const double* y = y_arr.data();
    const double* taus = taus_arr.data();
    const double* q = q_arr.data();
    const double* f = f_arr.data();

    tau_ord.resize(T);
    std::iota(tau_ord.begin(), tau_ord.end(), size_t{0});
    std::stable_sort(tau_ord.begin(), tau_ord.end(),
                     [&](size_t a, size_t b) { return q[a] < q[b]; });
    std::vector<double> q_sorted(T);
    c_sorted.resize(T);
    inv_f_sorted.resize(T);
    for (size_t st = 0; st < T; st++) {
        const size_t t = tau_ord[st];
        if (!(f[t] > 0.0) || !std::isfinite(f[t])) throw std::runtime_error("set_bins: fhat must be positive.");
        q_sorted[st] = q[t];
        c_sorted[st] = q[t] + taus[t] / f[t];
        inv_f_sorted[st] = 1.0 / f[t];
    }

    // Bins in sorted-sample order: y_i <= q_sorted[st]  <=>  bin_i <= st.
    bin_sorted.assign(N, 0);
    for (size_t s = 0; s < N; s++) {
        const double yi = y[static_cast<size_t>(order[s])];
        bin_sorted[s] = static_cast<uint16_t>(
            std::lower_bound(q_sorted.begin(), q_sorted.end(), yi) - q_sorted.begin());
    }

    // Q^T RIF without the N x T RIF: per-bin sums of Q rows, then prefix sums.
    std::vector<double> qbin((T + 1) * K, 0.0), w(K, 0.0), wc(K, 0.0);
    for (size_t s = 0; s < N; s++) {
        const double* qr = &RQ[s * Crow + T];
        double* qb = &qbin[static_cast<size_t>(bin_sorted[s]) * K];
        for (size_t k = 0; k < K; k++) {
            qb[k] += qr[k];
            const double x = qr[k], yv = w[k], tt = yv + x;   // Neumaier, as for e
            wc[k] += (std::abs(yv) >= std::abs(x)) ? ((yv - tt) + x) : ((x - tt) + yv);
            w[k] = tt;
        }
    }
    for (size_t k = 0; k < K; k++) w[k] += wc[k];
    QtR_tk.assign(T * K, 0.0);
    std::vector<double> run(K, 0.0);
    for (size_t st = 0; st < T; st++) {
        for (size_t k = 0; k < K; k++) run[k] += qbin[st * K + k];
        const size_t t = tau_ord[st];
        for (size_t k = 0; k < K; k++) {
            QtR_tk[t * K + k] = c_sorted[st] * w[k] - run[k] * inv_f_sorted[st];
        }
    }

    // Slim sweep operand [Q | e | 0 pad].
    Crow2 = ((K + 1 + CBLK - 1) / CBLK) * CBLK;
    RQ2.assign(N * Crow2, 0.0);
    for (size_t s = 0; s < N; s++) {
        std::memcpy(&RQ2[s * Crow2], &RQ[s * Crow + T], (K + 1) * sizeof(double));
    }
    has_bins = true;
}

// ---------------------------------------------------------------------------
// Per-pass wall clock of the last run_fast call, for profiling only.  Not part
// of the numerical contract.
// ---------------------------------------------------------------------------
inline double* fast_timings() {
    static double t[6] = {0, 0, 0, 0, 0, 0};
    return t;
}

inline double wall_now() {
#ifdef _OPENMP
    return omp_get_wtime();
#else
    return 0.0;
#endif
}

// Jackknife beta and se only (no T x T covariance): the same definition as
// core.jackknife_beta_se_batch -- leave-one-block-out betas, 0 where a block's
// complement has a non-positive denominator, se^2 = (B-1)/B * sum (jk - mean)^2.
inline void jackknife_beta_se(const double* __restrict__ sp, size_t B, size_t T,
                              double* __restrict__ beta, double* __restrict__ se,
                              double* __restrict__ d_out, double* __restrict__ ntot,
                              double* __restrict__ mean, double* __restrict__ ss) {
    const size_t oc = 1 + T;
    double D = 0.0;
    for (size_t t = 0; t < T; t++) ntot[t] = 0.0;
    for (size_t b = 0; b < B; b++) {
        const double* r = sp + b * oc;
        D += r[0];
        for (size_t t = 0; t < T; t++) ntot[t] += r[1 + t];
    }
    *d_out = D;
    if (!(D > 0.0)) return;   // caller pre-filled NaN
    for (size_t t = 0; t < T; t++) { beta[t] = ntot[t] / D; mean[t] = 0.0; ss[t] = 0.0; }
    for (size_t b = 0; b < B; b++) {
        const double* r = sp + b * oc;
        const double den = D - r[0];
        if (!(den > 0.0)) continue;
        const double inv = 1.0 / den;
        for (size_t t = 0; t < T; t++) mean[t] += (ntot[t] - r[1 + t]) * inv;
    }
    const double invB = 1.0 / static_cast<double>(B);
    for (size_t t = 0; t < T; t++) mean[t] *= invB;
    for (size_t b = 0; b < B; b++) {
        const double* r = sp + b * oc;
        const double den = D - r[0];
        if (den > 0.0) {
            const double inv = 1.0 / den;
            for (size_t t = 0; t < T; t++) {
                const double v = (ntot[t] - r[1 + t]) * inv - mean[t];
                ss[t] += v * v;
            }
        } else {
            for (size_t t = 0; t < T; t++) ss[t] += mean[t] * mean[t];
        }
    }
    const double f = static_cast<double>(B - 1) / static_cast<double>(B);
    for (size_t t = 0; t < T; t++) se[t] = std::sqrt(f * ss[t]);
}

// ---------------------------------------------------------------------------
// Driver.  ``raw_out`` (M, n_blocks, 1+T) and the (betas, se, cov) triple are
// alternatives: exactly one set is non-null.
// ---------------------------------------------------------------------------

inline void run_fast(
    const double* __restrict__ G_ptr,
    size_t N, size_t M,
    const BlockPlan& plan,
    double min_mac,
    double min_var,
    int n_threads,
    int panel_width,
    double* raw_out,
    double* betas_ptr,
    double* se_ptr,
    double* cov_ptr,
    double* jk_beta = nullptr,   // (M, T) jackknife beta only, or null
    double* jk_se = nullptr,     // (M, T)
    double* jk_d = nullptr       // (M,)  sum of block denominators (d_resid)
) {
    const size_t T = plan.T, K = plan.K;
    const int n_blocks = plan.n_blocks;
    const size_t n_blocks_sz = static_cast<size_t>(n_blocks);
    const size_t NB = plan.n_buckets();
    const size_t out_cols = 1 + T;
    const size_t n_cov = T * (T + 1) / 2;
    const bool avx2 = use_avx2();

    size_t Pmax = (panel_width > 0) ? static_cast<size_t>(panel_width) : M;
    if (Pmax > M) Pmax = M;
    if (Pmax < 1) Pmax = 1;

#ifdef _OPENMP
    if (n_threads > 0) omp_set_num_threads(n_threads);
#endif

    const bool bins = plan.has_bins;
    const double* __restrict__ RQ = bins ? plan.RQ2.data() : plan.RQ.data();
    const size_t Crow = bins ? plan.Crow2 : plan.Crow;
    const size_t T_sweep = bins ? 0 : T;   // the sweep forms W only without bins
    const double* __restrict__ S_all = plan.S.data();
    const double* __restrict__ A_all = plan.A.data();
    const double* __restrict__ E_all = plan.E.data();
    const double* __restrict__ F_all = plan.F.data();
    const double* __restrict__ H_all = plan.H.data();
    const int32_t* __restrict__ posv = plan.pos.data();
    const double nan_val = std::numeric_limits<double>::quiet_NaN();

    std::vector<char> keep(M, 0);
    std::vector<double> gmean_all(M, 0.0);

    double* tm = fast_timings();
    for (int z = 0; z < 6; z++) tm[z] = 0.0;
    double t_mark = wall_now();

    // Pass 1: the MAC filter.  Reads one contiguous G column per SNP and forms the
    // finite count and sum in ascending sample order, exactly as the exact kernel
    // does, so the accept/reject decision and the imputation mean are bit-identical.
    #pragma omp parallel for schedule(static)
    for (int64_t j = 0; j < static_cast<int64_t>(M); j++) {
        const double* __restrict__ col = G_ptr + static_cast<size_t>(j) * N;
        size_t n_valid = 0;
        double g_sum = 0.0;
        for (size_t i = 0; i < N; i++) {
            const double v = col[i];
            if (std::isfinite(v)) { g_sum += v; n_valid++; }
        }
        if (n_valid < 100) continue;
        const double mac = std::min(g_sum, 2.0 * static_cast<double>(n_valid) - g_sum);
        if (mac < min_mac) continue;
        keep[static_cast<size_t>(j)] = 1;
        gmean_all[static_cast<size_t>(j)] = g_sum / static_cast<double>(n_valid);
    }

    std::vector<size_t> kept;
    kept.reserve(M);
    for (size_t j = 0; j < M; j++) if (keep[j]) kept.push_back(j);
    tm[0] += wall_now() - t_mark; t_mark = wall_now();

    std::vector<size_t> cols;
    std::vector<double> gmean, QTG;
    Scratch& Gt_buf = plan.gt_buf;
    Scratch& W_buf = plan.w_buf;
    Scratch& U_buf = plan.u_buf;
    Scratch& S2_buf = plan.s2_buf;
    Scratch& V_buf = plan.v_buf;

    for (size_t panel0 = 0; panel0 < kept.size(); panel0 += Pmax) {
        const size_t panel1 = std::min(kept.size(), panel0 + Pmax);
        const size_t P = panel1 - panel0;
        // 64-byte aligned tile rows, and a multiple of PBLK so the sweep runs
        // whole register blocks with no remainder (the padding columns are zero).
        const size_t LD = ((P + PBLK - 1) / PBLK) * PBLK;
        const size_t Ppad = LD;

        cols.assign(kept.begin() + panel0, kept.begin() + panel1);
        gmean.resize(P);
        for (size_t p = 0; p < P; p++) gmean[p] = gmean_all[cols[p]];

        // Every element of these is written before it is read (see pack_chunk and
        // sweep_block), so they are deliberately left uninitialised.
        double* __restrict__ Gt = Gt_buf.get(N * LD);
        double* __restrict__ W = W_buf.get(n_blocks_sz * P * T);
        double* __restrict__ U = U_buf.get(NB * P * K);
        double* __restrict__ S2 = S2_buf.get(NB * P);
        double* __restrict__ V = V_buf.get(NB * P);
        QTG.assign(P * K, 0.0);
        tm[1] += wall_now() - t_mark; t_mark = wall_now();

        // Pass 2: pack the panel's genotype columns into one (N, LD) row-major
        // tile in block-sorted sample order, mean-imputing on the way.  Reads run
        // down G's contiguous axis; the only scattered access is the row scatter,
        // which moves whole 64-byte-aligned rows.
        const size_t CHUNK = std::max<size_t>(64, 32768 / LD);
        const size_t n_chunks = (N + CHUNK - 1) / CHUNK;
        #pragma omp parallel
        {
            std::vector<double> tmp(CHUNK * LD, 0.0);  // padding columns stay zero
            #pragma omp for schedule(static)
            for (int64_t c = 0; c < static_cast<int64_t>(n_chunks); c++) {
                const size_t i0 = static_cast<size_t>(c) * CHUNK;
                const size_t i1 = std::min(N, i0 + CHUNK);
                if (avx2) {
#ifdef GMOQ_HAVE_AVX2_CLONE
                    pack_chunk_avx2(G_ptr, N, cols.data(), gmean.data(), P, LD,
                                    i0, i1, posv, tmp.data(), Gt);
#endif
                } else {
                    pack_chunk_base(G_ptr, N, cols.data(), gmean.data(), P, LD,
                                    i0, i1, posv, tmp.data(), Gt);
                }
            }
        }
        tm[2] += wall_now() - t_mark; t_mark = wall_now();

        // Pass 3: the streaming sweep, parallel over blocks.  Each block owns a
        // disjoint slice of W, U and S2, so there is no shared write set at all.
        // Samples per L1-resident sweep slice.  The slice holds PBLK genotype
        // columns and Crow (RIF | Q) columns per sample; the default target keeps
        // that plus the running accumulator inside a 32 kB 8-way L1 with room to
        // spare for conflict misses.  Overridable for tuning experiments.
        size_t stile_bytes = 49152;
        if (const char* e = std::getenv("GMOQ_STAGE1_STILE_BYTES")) {
            const long v = std::strtol(e, nullptr, 10);
            if (v >= 512 && v <= (1L << 20)) stile_bytes = static_cast<size_t>(v);
        }
        size_t STILE = stile_bytes / ((Crow + static_cast<size_t>(PBLK)) * sizeof(double));
        if (STILE < 8) STILE = 8;

        #pragma omp parallel
        {
            std::vector<double> accm(Crow * static_cast<size_t>(PBLK), 0.0);
            #pragma omp for schedule(static)
            for (int64_t b = 0; b < static_cast<int64_t>(NB); b++) {
                const size_t bb = static_cast<size_t>(b);
                const size_t s0 = static_cast<size_t>(plan.offsets[bb]);
                const size_t s1 = static_cast<size_t>(plan.offsets[bb + 1]);
                double* ubp = U + bb * P * K;
                double* vbp = V + bb * P;
                double* s2p = S2 + bb * P;
                double* Wbp = (!bins && bb < n_blocks_sz) ? W + bb * P * T : nullptr;
                if (avx2) {
#ifdef GMOQ_HAVE_AVX2_CLONE
                    sweep_block_avx2(Gt, LD, RQ, Crow, T_sweep, K, s0, s1, P, Ppad,
                                     Wbp, ubp, vbp, s2p, accm.data(), STILE);
#endif
                } else {
                    sweep_block_base(Gt, LD, RQ, Crow, T_sweep, K, s0, s1, P, Ppad,
                                     Wbp, ubp, vbp, s2p, accm.data(), STILE);
                }
            }
        }

        if (bins) {
            // Pass 3b: W from per-(block, bin) sums of the packed (centred) genotypes.
            const uint16_t* __restrict__ binv = plan.bin_sorted.data();
            const size_t* tord = plan.tau_ord.data();
            const double* __restrict__ cs = plan.c_sorted.data();
            const double* __restrict__ ifs = plan.inv_f_sorted.data();
            const double* __restrict__ QtR = plan.QtR_tk.data();
            #pragma omp parallel
            {
                std::vector<double> hacc((T + 1) * LD), tot(LD), cum(LD);
                #pragma omp for schedule(static)
                for (int64_t b = 0; b < static_cast<int64_t>(n_blocks_sz); b++) {
                    const size_t bb = static_cast<size_t>(b);
                    const size_t s0 = static_cast<size_t>(plan.offsets[bb]);
                    const size_t s1 = static_cast<size_t>(plan.offsets[bb + 1]);
                    std::fill(hacc.begin(), hacc.end(), 0.0);
                    for (size_t s = s0; s < s1; s++) {
                        const double* __restrict__ g = Gt + s * LD;
                        double* __restrict__ h = hacc.data() + static_cast<size_t>(binv[s]) * LD;
                        for (size_t p = 0; p < LD; p++) h[p] += g[p];
                    }
                    std::fill(tot.begin(), tot.end(), 0.0);
                    for (size_t st = 0; st <= T; st++) {
                        const double* __restrict__ h = hacc.data() + st * LD;
                        for (size_t p = 0; p < LD; p++) tot[p] += h[p];
                    }
                    std::fill(cum.begin(), cum.end(), 0.0);
                    double* __restrict__ Wb = W + bb * P * T;
                    const double* __restrict__ Ub = U + bb * P * K;
                    for (size_t st = 0; st < T; st++) {
                        const double* __restrict__ h = hacc.data() + st * LD;
                        for (size_t p = 0; p < LD; p++) cum[p] += h[p];
                        const size_t t = tord[st];
                        const double* __restrict__ qr = QtR + t * K;
                        for (size_t p = 0; p < P; p++) {
                            const double* __restrict__ up = Ub + p * K;
                            double proj = 0.0;
                            for (size_t k = 0; k < K; k++) proj += up[k] * qr[k];
                            Wb[p * T + t] = cs[st] * tot[p] - cum[p] * ifs[st] - proj;
                        }
                    }
                }
            }
        }
        tm[3] += wall_now() - t_mark; t_mark = wall_now();

        // Pass 4: q_p = sum_b u_b[:, p], blocks in ascending order.
        #pragma omp parallel for schedule(static)
        for (int64_t p = 0; p < static_cast<int64_t>(P); p++) {
            double* __restrict__ qp = QTG.data() + static_cast<size_t>(p) * K;
            for (size_t bb = 0; bb < NB; bb++) {
                const double* __restrict__ ubp = U + (bb * P + static_cast<size_t>(p)) * K;
                for (size_t k = 0; k < K; k++) qp[k] += ubp[k];
            }
        }
        tm[4] += wall_now() - t_mark; t_mark = wall_now();

        // Pass 5: apply the covariate correction and finalise, in groups of PG
        // SNPs so S_b and A_b are read once per group rather than once per SNP.
        const size_t PG = 8;
        const size_t n_groups = (P + PG - 1) / PG;
        #pragma omp parallel
        {
            std::vector<double> buf(PG * n_blocks_sz * out_cols, 0.0);
            std::vector<double> beta(T, nan_val), se(T, nan_val), covu(n_cov, nan_val);
            std::vector<double> g2(PG, 0.0);
            std::vector<double> jk(n_blocks_sz * T, 0.0), ntot(T, 0.0), mean(T, 0.0),
                                Cbuf(T * T, 0.0);

            #pragma omp for schedule(dynamic, 1)
            for (int64_t pg = 0; pg < static_cast<int64_t>(n_groups); pg++) {
                const size_t p0 = static_cast<size_t>(pg) * PG;
                const size_t ng = std::min(P, p0 + PG) - p0;

                if (avx2) {
#ifdef GMOQ_HAVE_AVX2_CLONE
                    correct_group_avx2(S_all, A_all, E_all, F_all, H_all,
                                       W, U, V, S2, QTG.data(), gmean.data(),
                                       NB, n_blocks_sz, P, T, K, p0, ng,
                                       buf.data(), g2.data());
#endif
                } else {
                    correct_group_base(S_all, A_all, E_all, F_all, H_all,
                                       W, U, V, S2, QTG.data(), gmean.data(),
                                       NB, n_blocks_sz, P, T, K, p0, ng,
                                       buf.data(), g2.data());
                }

                for (size_t jj = 0; jj < ng; jj++) {
                    if (!(g2[jj] >= min_var)) continue;
                    const size_t j = cols[p0 + jj];
                    const double* sp = buf.data() + jj * n_blocks_sz * out_cols;
                    if (jk_beta != nullptr) {
                        jackknife_beta_se(sp, n_blocks_sz, T, jk_beta + j * T, jk_se + j * T,
                                          jk_d + j, ntot.data(), mean.data(), jk.data());
                        if (raw_out == nullptr) continue;
                    }
                    if (raw_out != nullptr) {
                        std::memcpy(raw_out + j * n_blocks_sz * out_cols, sp,
                                    n_blocks_sz * out_cols * sizeof(double));
                        continue;
                    }
                    bool ok;
                    if (avx2) {
#ifdef GMOQ_HAVE_AVX2_CLONE
                        ok = finalize_stats_fast_avx2(sp, n_blocks, T, beta.data(), se.data(),
                                                      covu.data(), jk.data(), ntot.data(),
                                                      mean.data(), Cbuf.data());
#else
                        ok = false;
#endif
                    } else {
                        ok = finalize_stats_fast_base(sp, n_blocks, T, beta.data(), se.data(),
                                                      covu.data(), jk.data(), ntot.data(),
                                                      mean.data(), Cbuf.data());
                    }
                    if (!ok) continue;
                    for (size_t t = 0; t < T; t++) {
                        betas_ptr[j + M * t] = beta[t];
                        se_ptr[j + M * t] = se[t];
                    }
                    for (size_t c = 0; c < n_cov; c++) cov_ptr[j + M * c] = covu[c];
                }
            }
        }
        tm[5] += wall_now() - t_mark; t_mark = wall_now();
    }
}

inline py::tuple finalize_block_scores_fast(
    py::array_t<double, py::array::f_style | py::array::forcecast> G_arr,
    const BlockPlan& plan,
    double min_mac,
    double min_var,
    int n_threads,
    int panel_width
) {
    auto G_buf = G_arr.request();
    if (G_buf.ndim != 2) throw std::runtime_error("G must be 2-D.");
    const size_t N = static_cast<size_t>(G_buf.shape[0]);
    const size_t M = static_cast<size_t>(G_buf.shape[1]);
    if (N != plan.N) {
        throw std::runtime_error("G has a different number of rows than the BlockPlan was built for.");
    }
    const size_t T = plan.T;
    const size_t n_cov = T * (T + 1) / 2;

    py::array_t<double, py::array::f_style> betas_arr({M, T});
    py::array_t<double, py::array::f_style> se_arr({M, T});
    py::array_t<double, py::array::f_style> cov_arr({M, n_cov});
    double* betas_ptr = static_cast<double*>(betas_arr.request().ptr);
    double* se_ptr = static_cast<double*>(se_arr.request().ptr);
    double* cov_ptr = static_cast<double*>(cov_arr.request().ptr);
    const double nan_val = std::numeric_limits<double>::quiet_NaN();
    std::fill(betas_ptr, betas_ptr + M * T, nan_val);
    std::fill(se_ptr, se_ptr + M * T, nan_val);
    std::fill(cov_ptr, cov_ptr + M * n_cov, nan_val);
    if (M == 0) return py::make_tuple(betas_arr, se_arr, cov_arr);

    const double* G_ptr = static_cast<const double*>(G_buf.ptr);
    {
        py::gil_scoped_release gil_release;
        run_fast(G_ptr, N, M, plan, min_mac, min_var, n_threads, panel_width,
                 nullptr, betas_ptr, se_ptr, cov_ptr);
    }
    return py::make_tuple(betas_arr, se_arr, cov_arr);
}

inline py::array_t<double> block_scores_fast(
    py::array_t<double, py::array::f_style | py::array::forcecast> G_arr,
    const BlockPlan& plan,
    double min_mac,
    double min_var,
    int n_threads,
    int panel_width
) {
    auto G_buf = G_arr.request();
    if (G_buf.ndim != 2) throw std::runtime_error("G must be 2-D.");
    const size_t N = static_cast<size_t>(G_buf.shape[0]);
    const size_t M = static_cast<size_t>(G_buf.shape[1]);
    if (N != plan.N) {
        throw std::runtime_error("G has a different number of rows than the BlockPlan was built for.");
    }
    const size_t out_cols = 1 + plan.T;

    std::vector<size_t> shape = {M, static_cast<size_t>(plan.n_blocks), out_cols};
    py::array_t<double> result(shape);
    double* result_ptr = static_cast<double*>(result.request().ptr);
    std::fill(result_ptr,
              result_ptr + M * static_cast<size_t>(plan.n_blocks) * out_cols,
              std::numeric_limits<double>::quiet_NaN());
    if (M == 0) return result;

    const double* G_ptr = static_cast<const double*>(G_buf.ptr);
    {
        py::gil_scoped_release gil_release;
        run_fast(G_ptr, N, M, plan, min_mac, min_var, n_threads, panel_width,
                 result_ptr, nullptr, nullptr, nullptr);
    }
    return result;
}

// Block kernel + jackknife beta/se in C++: returns (raw or None, beta (M,T),
// se (M,T), d_resid (M,)), NaN for SNPs that fail the MAC/variance filters.
inline py::tuple block_jackknife_fast(
    py::array_t<double, py::array::f_style | py::array::forcecast> G_arr,
    const BlockPlan& plan,
    double min_mac,
    double min_var,
    int n_threads,
    int panel_width,
    bool want_raw
) {
    auto G_buf = G_arr.request();
    if (G_buf.ndim != 2) throw std::runtime_error("G must be 2-D.");
    const size_t N = static_cast<size_t>(G_buf.shape[0]);
    const size_t M = static_cast<size_t>(G_buf.shape[1]);
    if (N != plan.N) {
        throw std::runtime_error("G has a different number of rows than the BlockPlan was built for.");
    }
    const size_t T = plan.T;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    py::array_t<double> beta({M, T}), se({M, T}), d(std::vector<size_t>{M});
    std::fill(beta.mutable_data(), beta.mutable_data() + M * T, nan);
    std::fill(se.mutable_data(), se.mutable_data() + M * T, nan);
    std::fill(d.mutable_data(), d.mutable_data() + M, nan);
    py::object raw_obj = py::none();
    double* raw_ptr = nullptr;
    if (want_raw) {
        std::vector<size_t> shape = {M, static_cast<size_t>(plan.n_blocks), 1 + T};
        py::array_t<double> raw(shape);
        raw_ptr = raw.mutable_data();
        std::fill(raw_ptr, raw_ptr + M * static_cast<size_t>(plan.n_blocks) * (1 + T), nan);
        raw_obj = raw;
    }
    if (M > 0) {
        const double* G_ptr = static_cast<const double*>(G_buf.ptr);
        double* bp = beta.mutable_data();
        double* sp = se.mutable_data();
        double* dp = d.mutable_data();
        py::gil_scoped_release gil_release;
        run_fast(G_ptr, N, M, plan, min_mac, min_var, n_threads, panel_width,
                 raw_ptr, nullptr, nullptr, nullptr, bp, sp, dp);
    }
    return py::make_tuple(raw_obj, beta, se, d);
}

}  // namespace gmoq_fast

#endif  // GMOQ_STAGE1_FAST_HPP
