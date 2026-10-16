"""The panel kernel must reproduce the per-SNP kernel (bit for bit on x86 GCC; to ~1 ulp where the compiler contracts FMAs, e.g. clang on arm64).

The panel form re-nests the loops of ``prepare_genotype_raw`` +
``accumulate_block_stats_flat`` so that ``panel_size`` SNPs share one sweep over the
samples, and takes RIF/Q in C order instead of F order.  It is not a reassociation:
every accumulator still receives the same addends in the same order, so the output must
be *identical*, not merely close.  These tests assert agreement to 1e-12 relative.
"""

import numpy as np
import pytest

from fungwas_stage1 import core

pytestmark = pytest.mark.skipif(not core.HAVE_CPP, reason="C++ extension not built")


def _problem(seed=7, N=2003, M=23, T=9, K=5, B=17):
    rng = np.random.default_rng(seed)
    G = np.asfortranarray(rng.normal(1.0, 0.5, (N, M)))
    G[rng.random((N, M)) < 0.02] = np.nan          # missing dosages -> mean imputation
    G[:, 2] = 0.0
    G[:3, 2] = 1.0                                  # fails the MAC filter
    G[:, 4] = np.nan
    G[:50, 4] = 1.0                                 # fails n_valid >= 100
    X = np.column_stack([np.ones(N), rng.normal(0, 1, (N, K - 1))])
    Q, _ = np.linalg.qr(X)
    rif = rng.normal(0, 1, (N, T))
    rif_resid = rif - Q @ (Q.T @ rif)
    block_ids = rng.integers(0, B, N).astype(np.int32)
    block_ids[11] = -1                              # out-of-range block id is skipped
    return G, rif_resid, Q, block_ids, B


@pytest.mark.parametrize("panel_size", [1, 3, 8, 16, 64])
def test_panel_matches_per_snp_bitwise(panel_size):
    G, rif_resid, Q, block_ids, B = _problem()
    ref = core.compute_block_scores(G, rif_resid, Q, block_ids, B, n_threads=2)
    got = core.compute_block_scores_panel(
        G, rif_resid, Q, block_ids, B, n_threads=2, panel_size=panel_size
    )
    assert got.shape == ref.shape
    np.testing.assert_allclose(got, ref, rtol=1e-12, atol=1e-14, equal_nan=True)


@pytest.mark.parametrize("n_threads", [1, 2, 4])
def test_panel_is_thread_count_independent(n_threads):
    G, rif_resid, Q, block_ids, B = _problem()
    ref = core.compute_block_scores_panel(
        G, rif_resid, Q, block_ids, B, n_threads=1, panel_size=8
    )
    got = core.compute_block_scores_panel(
        G, rif_resid, Q, block_ids, B, n_threads=n_threads, panel_size=8
    )
    assert np.array_equal(ref, got, equal_nan=True)


def test_panel_rejects_the_same_snps():
    G, rif_resid, Q, block_ids, B = _problem()
    ref = core.compute_block_scores(G, rif_resid, Q, block_ids, B, n_threads=2)
    got = core.compute_block_scores_panel(
        G, rif_resid, Q, block_ids, B, n_threads=2, panel_size=8
    )
    ref_bad = np.isnan(ref[:, 0, 0])
    got_bad = np.isnan(got[:, 0, 0])
    assert ref_bad.sum() == 2          # the MAC and n_valid failures above
    assert np.array_equal(ref_bad, got_bad)


def test_fused_panel_matches_fused_per_snp_bitwise():
    G, rif_resid, Q, block_ids, B = _problem()
    b0, s0, c0 = core.compute_block_finalize(G, rif_resid, Q, block_ids, B, n_threads=2)
    b1, s1, c1 = core.compute_block_finalize_panel(
        G, rif_resid, Q, block_ids, B, n_threads=2, panel_size=8
    )
    np.testing.assert_allclose(b1, b0, rtol=1e-12, atol=1e-14, equal_nan=True)
    np.testing.assert_allclose(s1, s0, rtol=1e-12, atol=1e-14, equal_nan=True)
    np.testing.assert_allclose(c1, c0, rtol=1e-12, atol=1e-14, equal_nan=True)


def test_panel_does_not_mutate_its_inputs():
    G, rif_resid, Q, block_ids, B = _problem()
    G_copy, rif_copy, Q_copy = G.copy(), rif_resid.copy(), Q.copy()
    core.compute_block_scores_panel(G, rif_resid, Q, block_ids, B, n_threads=2, panel_size=8)
    assert np.array_equal(G, G_copy, equal_nan=True)
    assert np.array_equal(rif_resid, rif_copy)
    assert np.array_equal(Q, Q_copy)
