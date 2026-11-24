// Hot leaves of the fast Stage 1 panel kernel.
//
// Included TWICE by stage1_fast.hpp - once at the x86-64 baseline and once under
// `#pragma GCC target("arch=haswell")` - and dispatched at first use on
// __builtin_cpu_supports.  Doing it this way means the fast kernel gets AVX2+FMA
// without raising the ISA of the whole translation unit, so the exact kernel in
// stage1_pybind.cpp keeps its byte-for-byte arithmetic with the default (empty)
// GMOQ_STAGE1_ARCH.
//
// GMOQ_SUFFIX must be defined to the name suffix before including.

#define GMOQ_CAT2(a, b) a##b
#define GMOQ_CAT(a, b) GMOQ_CAT2(a, b)
#define GMOQ_FN(name) GMOQ_CAT(name, GMOQ_SUFFIX)

// One sweep over a block's samples producing, for PBLK genotype columns at a
// time, the cross-products against the combined (RIF | Q) row block and the
// genotype sum of squares.
//
// Two levels of blocking:
//
//   * register - a CBLK x PBLK tile of accumulators, i.e. CBLK*PBLK/4 = 8 ymm
//     registers held across the whole sample loop.  Each sample then costs two
//     vector loads, CBLK broadcasts and CBLK*PBLK/4 FMAs and touches nothing
//     else.  GCC's auto-vectoriser does not find this shape on its own - it
//     sinks the accumulator to the stack - which is worth ~4x here.
//
//   * sample (STILE) - the (pb, cb) tiles are swept over one STILE-row slice of
//     the block at a time, sized so the slice's genotype columns and (RIF | Q)
//     rows together fit in L1.  Without it the Crow/CBLK * Ppad/PBLK passes each
//     re-read the block's whole working set from L2, and the sweep runs at about
//     a third of its issue-limited rate.  The Crow x PBLK running accumulator
//     that carries across slices is a few kB and stays in L1 too.
//
// RQ is padded to Crow = round_up(T + K, CBLK) with zeros, so the c loop never
// needs a remainder branch; the padding columns accumulate zero and are dropped.
//
// Summation order: for every (SNP, column) pair the addends arrive in ascending
// sorted-sample order, which - the counting sort being stable - is ascending
// original row index, the same order the exact kernel uses.
static void GMOQ_FN(sweep_block)(
    const double* __restrict__ Gt, size_t LD,
    const double* __restrict__ RQ, size_t Crow, size_t T, size_t K,
    size_t s0, size_t s1, size_t P, size_t Ppad,
    double* __restrict__ Wb,      // P*T, or null to skip the RIF products
    double* __restrict__ ub,      // P*K
    double* __restrict__ vb,      // P   (cross product against the e column)
    double* __restrict__ s2b,     // P
    double* __restrict__ accm,    // Crow * PBLK running accumulator
    size_t STILE                  // samples per L1-resident slice
) {
    const size_t C = T + K + 1;   // RIF | Q | e

    for (size_t pb = 0; pb < Ppad; pb += PBLK) {
        for (size_t z = 0; z < Crow * static_cast<size_t>(PBLK); z++) accm[z] = 0.0;
        double st[PBLK];
        for (size_t j = 0; j < static_cast<size_t>(PBLK); j++) st[j] = 0.0;

        for (size_t t0 = s0; t0 < s1; t0 += STILE) {
            const size_t t1 = std::min(s1, t0 + STILE);

            for (size_t s = t0; s < t1; s++) {
                const double* __restrict__ g = Gt + s * LD + pb;
                for (size_t j = 0; j < static_cast<size_t>(PBLK); j++) st[j] += g[j] * g[j];
            }

            for (size_t cb = 0; cb < C; cb += CBLK) {
                double* __restrict__ ap = accm + cb * static_cast<size_t>(PBLK);
#ifdef GMOQ_INTRIN
                __m256d a0 = _mm256_loadu_pd(ap +  0), a1 = _mm256_loadu_pd(ap +  4);
                __m256d a2 = _mm256_loadu_pd(ap +  8), a3 = _mm256_loadu_pd(ap + 12);
                __m256d a4 = _mm256_loadu_pd(ap + 16), a5 = _mm256_loadu_pd(ap + 20);
                __m256d a6 = _mm256_loadu_pd(ap + 24), a7 = _mm256_loadu_pd(ap + 28);
                __m256d a8 = _mm256_loadu_pd(ap + 32), a9 = _mm256_loadu_pd(ap + 36);
                __m256d aa = _mm256_loadu_pd(ap + 40), ab = _mm256_loadu_pd(ap + 44);
                for (size_t s = t0; s < t1; s++) {
                    const double* g = Gt + s * LD + pb;
                    const double* r = RQ + s * Crow + cb;
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
#else
                for (size_t s = t0; s < t1; s++) {
                    const double* __restrict__ g = Gt + s * LD + pb;
                    const double* __restrict__ r = RQ + s * Crow + cb;
                    for (size_t c = 0; c < static_cast<size_t>(CBLK); c++) {
                        const double rv = r[c];
                        double* __restrict__ ac = ap + c * static_cast<size_t>(PBLK);
                        for (size_t j = 0; j < static_cast<size_t>(PBLK); j++) ac[j] += g[j] * rv;
                    }
                }
#endif
            }
        }

        for (size_t j = 0; j < static_cast<size_t>(PBLK) && pb + j < P; j++) s2b[pb + j] = st[j];
        for (size_t cc = 0; cc < C; cc++) {
            const double* __restrict__ ac = accm + cc * static_cast<size_t>(PBLK);
            if (cc < T) {
                if (Wb == nullptr) continue;
                for (size_t j = 0; j < static_cast<size_t>(PBLK) && pb + j < P; j++)
                    Wb[(pb + j) * T + cc] = ac[j];
            } else if (cc < T + K) {
                for (size_t j = 0; j < static_cast<size_t>(PBLK) && pb + j < P; j++)
                    ub[(pb + j) * K + (cc - T)] = ac[j];
            } else {
                for (size_t j = 0; j < static_cast<size_t>(PBLK) && pb + j < P; j++)
                    vb[pb + j] = ac[j];
            }
        }
    }
}

// Apply the covariate correction to one SNP group.  Kept out of line so it, too,
// gets the AVX2 clone: the inner loops run along the contiguous tau axis.
static void GMOQ_FN(correct_group)(
    const double* __restrict__ S_all, const double* __restrict__ A_all,
    const double* __restrict__ E_all, const double* __restrict__ F_all,
    const double* __restrict__ H_all,
    const double* __restrict__ W, const double* __restrict__ U,
    const double* __restrict__ V, const double* __restrict__ S2,
    const double* __restrict__ QTG, const double* __restrict__ gmean,
    size_t NB, size_t n_blocks_sz, size_t P, size_t T, size_t K,
    size_t p0, size_t ng,
    double* __restrict__ buf,    // ng * n_blocks * (1+T)
    double* __restrict__ g2      // ng
) {
    const size_t out_cols = 1 + T;
    for (size_t jj = 0; jj < ng; jj++) g2[jj] = 0.0;
    for (size_t bb = 0; bb < NB; bb++) {
        const double* __restrict__ Sb = S_all + bb * K * T;
        const double* __restrict__ Ab = A_all + bb * K * K;
        const double* __restrict__ Eb = E_all + bb * T;
        const double* __restrict__ Fb = F_all + bb * K;
        const double Hb = H_all[bb];
        const bool real_block = bb < n_blocks_sz;
        for (size_t jj = 0; jj < ng; jj++) {
            const size_t p = p0 + jj;
            const double* __restrict__ qp = QTG + p * K;
            const double* __restrict__ ubp = U + (bb * P + p) * K;
            const double gm = gmean[p];
            double t1 = 0.0, tf = 0.0;
            for (size_t k = 0; k < K; k++) { t1 += qp[k] * ubp[k]; tf += qp[k] * Fb[k]; }
            double t2 = 0.0;
            for (size_t k = 0; k < K; k++) {
                const double* __restrict__ Ak = Ab + k * K;
                double acc = 0.0;
                for (size_t l = 0; l < K; l++) acc += Ak[l] * qp[l];
                t2 += qp[k] * acc;
            }
            const double d0 = S2[bb * P + p] - 2.0 * t1 + t2
                            + 2.0 * gm * (V[bb * P + p] - tf) + gm * gm * Hb;
            g2[jj] += d0;
            if (!real_block) continue;
            double* __restrict__ o = buf + (jj * n_blocks_sz + bb) * out_cols;
            o[0] = d0;
            const double* __restrict__ Wbp = W + (bb * P + p) * T;
            for (size_t t = 0; t < T; t++) o[1 + t] = Wbp[t] + gm * Eb[t];
            for (size_t k = 0; k < K; k++) {
                const double qk = qp[k];
                const double* __restrict__ Sk = Sb + k * T;
                for (size_t t = 0; t < T; t++) o[1 + t] -= qk * Sk[t];
            }
        }
    }
}

// Jackknife finalisation.
//
// Identical arithmetic to finalize_stats_flat in stage1_pybind.cpp - every
// accumulator receives the same addends in the same order - but the covariance
// triple loop is re-nested with the block index outermost so the inner loop runs
// along the contiguous tau axis and vectorises.  The (a, b) accumulator tile is
// T*T doubles (16 kB at T = 45) and stays in L1/L2.
static bool GMOQ_FN(finalize_stats_fast)(
    const double* __restrict__ stats,
    int n_blocks,
    size_t T,
    double* __restrict__ beta,
    double* __restrict__ se,
    double* __restrict__ cov_upper,
    double* __restrict__ jk,       // n_blocks * T
    double* __restrict__ N_total,  // T
    double* __restrict__ mean,     // T
    double* __restrict__ C         // T * T
) {
    const size_t out_cols = 1 + T;
    double D_total = 0.0;
    for (size_t t = 0; t < T; t++) N_total[t] = 0.0;

    for (int b = 0; b < n_blocks; b++) {
        const double* __restrict__ s = stats + static_cast<size_t>(b) * out_cols;
        D_total += s[0];
        for (size_t t = 0; t < T; t++) N_total[t] += s[1 + t];
    }
    if (D_total <= 0.0) return false;
    for (size_t t = 0; t < T; t++) beta[t] = N_total[t] / D_total;

    for (size_t z = 0; z < static_cast<size_t>(n_blocks) * T; z++) jk[z] = 0.0;
    for (size_t t = 0; t < T; t++) mean[t] = 0.0;
    for (int b = 0; b < n_blocks; b++) {
        const double* __restrict__ s = stats + static_cast<size_t>(b) * out_cols;
        const double denom = D_total - s[0];
        if (denom <= 0.0) continue;
        double* __restrict__ jb = jk + static_cast<size_t>(b) * T;
        for (size_t t = 0; t < T; t++) {
            const double val = (N_total[t] - s[1 + t]) / denom;
            jb[t] = val;
            mean[t] += val;
        }
    }
    for (size_t t = 0; t < T; t++) mean[t] /= static_cast<double>(n_blocks);

    for (size_t z = 0; z < T * T; z++) C[z] = 0.0;
    for (int k = 0; k < n_blocks; k++) {
        const double* __restrict__ jb = jk + static_cast<size_t>(k) * T;
        for (size_t a = 0; a < T; a++) {
            const double da = jb[a] - mean[a];
            double* __restrict__ Ca = C + a * T;
            for (size_t b2 = a; b2 < T; b2++) Ca[b2] += da * (jb[b2] - mean[b2]);
        }
    }

    const double scale = (static_cast<double>(n_blocks) - 1.0) / static_cast<double>(n_blocks);
    size_t idx = 0;
    for (size_t a = 0; a < T; a++) {
        for (size_t b2 = a; b2 < T; b2++) {
            const double cov = scale * C[a * T + b2];
            cov_upper[idx++] = cov;
            if (a == b2) se[a] = std::sqrt(std::max(cov, 0.0));
        }
    }
    return true;
}

// Pack one chunk of samples into the block-sorted genotype tile.  The gather
// from G runs down G's contiguous axis; the scatter writes whole 64-byte-aligned
// tile rows.
static void GMOQ_FN(pack_chunk)(
    const double* __restrict__ G_ptr, size_t N,
    const size_t* __restrict__ cols, const double* __restrict__ gmean, size_t P,
    size_t LD, size_t i0, size_t i1,
    const int32_t* __restrict__ posv,
    double* __restrict__ tmp, double* __restrict__ Gt,
    const int64_t* __restrict__ rowmap, size_t N_src   // rowmap: sample i -> row of G (or null)
) {
    for (size_t p = 0; p < P; p++) {
        const double gm = gmean[p];
        double* __restrict__ dst = tmp + p;
        if (rowmap == nullptr) {
            const double* __restrict__ col = G_ptr + cols[p] * N;
            for (size_t i = i0; i < i1; i++) {
                const double v = col[i];
                // Mean-centred: a missing dosage is imputed to the mean, so it
                // contributes exactly zero.  Centring is what keeps the kernel's
                // cancellation ratio near one (see stage1_fast.hpp, Numerics).
                dst[(i - i0) * LD] = std::isfinite(v) ? (v - gm) : 0.0;
            }
        } else {
            const double* __restrict__ col = G_ptr + cols[p] * N_src;
            for (size_t i = i0; i < i1; i++) {
                const double v = col[rowmap[i]];
                dst[(i - i0) * LD] = std::isfinite(v) ? (v - gm) : 0.0;
            }
        }
    }
    for (size_t i = i0; i < i1; i++) {
        double* __restrict__ dst = Gt + static_cast<size_t>(posv[i]) * LD;
        const double* __restrict__ src = tmp + (i - i0) * LD;
#ifdef GMOQ_INTRIN
        // The tile is written once here and read once by the sweep, and it is far
        // larger than L3, so pulling each destination line into cache first (the
        // read-for-ownership a normal store pays) is wasted traffic.  Gt is
        // 64-byte aligned and LD is a multiple of 8, so whole rows stream.
        for (size_t z = 0; z < LD; z += 4) _mm256_stream_pd(dst + z, _mm256_loadu_pd(src + z));
#else
        std::memcpy(dst, src, LD * sizeof(double));
#endif
    }
#ifdef GMOQ_INTRIN
    _mm_sfence();
#endif
}

#undef GMOQ_FN
#undef GMOQ_CAT
#undef GMOQ_CAT2
