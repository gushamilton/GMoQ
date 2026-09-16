"""The fast panel kernel: what it must reproduce, and what it is allowed to change.

The fast kernel forms the same block statistics as the exact panel kernel from an
algebraically equivalent expansion (see ``src/stage1_fast.hpp``), so it reassociates
the summation and is *not* bit-identical.  What it must still guarantee:

  * agreement with the exact kernel to reassociation level;
  * the same accept/reject decisions - a SNP is NaN in one iff it is NaN in the other;
  * results that do not depend on the thread count, the panel width, or how the SNPs
    were split into batches.

The last three are what the restructuring rests on: a block's statistics depend only
on that block's samples, and a SNP's only on that SNP's column.  They are asserted
bit-exactly, because anything less would mean the decomposition leaks across
boundaries it is not supposed to.
"""

import numpy as np
import pytest

from fungwas_stage1 import core

pytestmark = pytest.mark.skipif(
    not core.HAVE_CPP or not hasattr(core.cpp_ext, "BlockPlan"),
    reason="C++ extension with the fast kernel not built",
)

# Reassociation-level tolerances, measured on the problems below and on full chr22.
# ``beta`` and the block statistics are compared on an absolute scale relative to the
# size of the quantity, because individual entries pass through zero.
ATOL_REL_SCALE = 1e-11


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


def _close(ref, got, scale=ATOL_REL_SCALE, name=""):
    """Agreement to reassociation level, on a scale set by the array itself."""
    assert ref.shape == got.shape
    assert np.array_equal(np.isnan(ref), np.isnan(got)), f"{name}: NaN pattern differs"
    m = ~np.isnan(ref)
    if not m.any():
        return
    tol = scale * max(float(np.max(np.abs(ref[m]))), 1e-300)
    err = float(np.max(np.abs(ref[m] - got[m])))
    assert err <= tol, f"{name}: max abs diff {err:.3e} exceeds {tol:.3e}"


def test_block_stats_match_the_exact_kernel():
    G, rif, Q, bid, B = _problem()
    ref = core.compute_block_scores_panel(G, rif, Q, bid, B, n_threads=2, panel_size=8)
    got = core.compute_block_scores_panel_fast(G, rif, Q, bid, B, n_threads=2)
    _close(ref, got, name="block stats")


def test_beta_se_cov_match_the_exact_kernel():
    G, rif, Q, bid, B = _problem()
    b0, s0, c0 = core.compute_block_finalize_panel(G, rif, Q, bid, B, n_threads=2, panel_size=8)
    b1, s1, c1 = core.compute_block_finalize_panel_fast(G, rif, Q, bid, B, n_threads=2)
    _close(b0, b1, name="beta")
    _close(s0, s1, name="se")
    _close(c0, c1, name="cov")
    # The SEs are strictly positive, so they also support a per-element relative check.
    m = ~np.isnan(s0)
    assert np.max(np.abs(s0[m] - s1[m]) / s0[m]) < 1e-12


def test_rejects_exactly_the_same_snps():
    G, rif, Q, bid, B = _problem()
    ref = core.compute_block_scores_panel(G, rif, Q, bid, B, n_threads=2, panel_size=8)
    got = core.compute_block_scores_panel_fast(G, rif, Q, bid, B, n_threads=2)
    ref_bad = np.isnan(ref[:, 0, 0])
    got_bad = np.isnan(got[:, 0, 0])
    assert ref_bad.sum() == 2          # the MAC and n_valid failures built above
    assert np.array_equal(ref_bad, got_bad)


@pytest.mark.parametrize("n_threads", [1, 2, 3, 4])
def test_thread_count_independence_is_bit_exact(n_threads):
    G, rif, Q, bid, B = _problem()
    plan = core.make_block_plan(rif, Q, bid, B, n_threads=1)
    ref = core.compute_block_scores_panel_fast(G, rif, Q, bid, B, n_threads=1, plan=plan)
    got = core.compute_block_scores_panel_fast(G, rif, Q, bid, B, n_threads=n_threads, plan=plan)
    assert np.array_equal(ref, got, equal_nan=True)


@pytest.mark.parametrize("panel_width", [0, 1, 2, 7, 8, 16, 64])
def test_panel_width_independence_is_bit_exact(panel_width):
    G, rif, Q, bid, B = _problem()
    plan = core.make_block_plan(rif, Q, bid, B, n_threads=2)
    ref = core.compute_block_scores_panel_fast(G, rif, Q, bid, B, n_threads=2,
                                               panel_width=0, plan=plan)
    got = core.compute_block_scores_panel_fast(G, rif, Q, bid, B, n_threads=2,
                                               panel_width=panel_width, plan=plan)
    assert np.array_equal(ref, got, equal_nan=True)


@pytest.mark.parametrize("batch", [1, 5, 11, 23])
def test_batch_split_independence_is_bit_exact(batch):
    """A SNP's answer must not depend on which SNPs it was batched with."""
    G, rif, Q, bid, B = _problem()
    plan = core.make_block_plan(rif, Q, bid, B, n_threads=2)
    ref = core.compute_block_scores_panel_fast(G, rif, Q, bid, B, n_threads=2, plan=plan)
    pieces = [
        core.compute_block_scores_panel_fast(
            np.asfortranarray(G[:, j0:j0 + batch]), rif, Q, bid, B, n_threads=2, plan=plan
        )
        for j0 in range(0, G.shape[1], batch)
    ]
    got = np.concatenate(pieces, axis=0)
    assert np.array_equal(ref, got, equal_nan=True)


def test_prebuilt_plan_matches_an_implicit_one():
    G, rif, Q, bid, B = _problem()
    plan = core.make_block_plan(rif, Q, bid, B, n_threads=2)
    a = core.compute_block_scores_panel_fast(G, rif, Q, bid, B, n_threads=2, plan=plan)
    b = core.compute_block_scores_panel_fast(G, rif, Q, bid, B, n_threads=2)
    assert np.array_equal(a, b, equal_nan=True)


def test_fused_output_matches_finalising_the_raw_stats():
    """The fused path must agree with running the jackknife on the raw statistics."""
    G, rif, Q, bid, B = _problem()
    plan = core.make_block_plan(rif, Q, bid, B, n_threads=2)
    stats = core.compute_block_scores_panel_fast(G, rif, Q, bid, B, n_threads=2, plan=plan)
    b1, s1, _ = core.compute_block_finalize_panel_fast(G, rif, Q, bid, B, n_threads=2, plan=plan)
    for j in range(G.shape[1]):
        if np.isnan(stats[j, 0, 0]):
            assert np.isnan(b1[j]).all()
            continue
        D = stats[j, :, 0].sum()
        Nt = stats[j, :, 1:].sum(axis=0)
        np.testing.assert_allclose(b1[j], Nt / D, rtol=1e-12, atol=0)


def test_out_of_range_block_ids_are_excluded_from_the_statistics():
    """Samples with an out-of-range block id feed the min_var gate but no block."""
    G, rif, Q, bid, B = _problem()
    bid2 = bid.copy()
    bid2[:200] = -1
    ref = core.compute_block_scores_panel(G, rif, Q, bid2, B, n_threads=2, panel_size=8)
    got = core.compute_block_scores_panel_fast(G, rif, Q, bid2, B, n_threads=2)
    _close(ref, got, name="block stats with out-of-range ids")


def test_does_not_mutate_its_inputs():
    G, rif, Q, bid, B = _problem()
    G_copy, rif_copy, Q_copy, bid_copy = G.copy(), rif.copy(), Q.copy(), bid.copy()
    core.compute_block_scores_panel_fast(G, rif, Q, bid, B, n_threads=2)
    assert np.array_equal(G, G_copy, equal_nan=True)
    assert np.array_equal(rif, rif_copy)
    assert np.array_equal(Q, Q_copy)
    assert np.array_equal(bid, bid_copy)


def test_single_tau_and_single_covariate_shapes():
    """T and K below the kernel's register-block widths must still be handled."""
    G, rif, Q, bid, B = _problem(seed=3, N=901, M=5, T=1, K=1, B=7)
    ref = core.compute_block_scores_panel(G, rif, Q, bid, B, n_threads=2, panel_size=4)
    got = core.compute_block_scores_panel_fast(G, rif, Q, bid, B, n_threads=2)
    _close(ref, got, name="block stats (T=1, K=1)")


def test_wide_panel_beyond_one_register_block():
    """More SNPs than the tile padding width, and a non-multiple of it."""
    G, rif, Q, bid, B = _problem(seed=11, N=1511, M=37, T=45, K=13, B=23)
    ref = core.compute_block_finalize_panel(G, rif, Q, bid, B, n_threads=2, panel_size=16)
    got = core.compute_block_finalize_panel_fast(G, rif, Q, bid, B, n_threads=2)
    for name, a, b in zip(("beta", "se", "cov"), ref, got):
        _close(a, b, name=name)
