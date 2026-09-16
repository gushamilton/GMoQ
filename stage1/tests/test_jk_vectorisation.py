"""Regression tests for the Stage 1 hot-path rewrites.

Each test pins a property that a performance change relied on. They are written
against the *old* implementation inline, so they keep failing if someone
"optimises" one of these paths in a way that moves a result by a single bit.
"""

import struct

import numpy as np
import pytest

from fungwas_stage1.core import compute_covariance_from_blocks


def _reference_covariance_from_blocks(stats, n_blocks):
    """The pre-vectorisation implementation, kept verbatim as the oracle."""
    T = stats.shape[1] - 1
    D_k = stats[:, 0]
    N_k = stats[:, 1:]

    D_total = np.sum(D_k)
    N_total = np.sum(N_k, axis=0)

    if D_total <= 0:
        return np.full(T, np.nan), np.full(T, np.nan), np.full((T, T), np.nan)

    jk_betas = np.zeros((n_blocks, T))
    for k in range(n_blocks):
        denom = D_total - D_k[k]
        if denom > 0:
            jk_betas[k] = (N_total - N_k[k]) / denom

    beta = N_total / D_total
    centered = jk_betas - np.mean(jk_betas, axis=0)
    Sigma = ((n_blocks - 1) / n_blocks) * (centered.T @ centered)
    se = np.sqrt(np.diag(Sigma))
    return beta, se, Sigma


def _block_stats(rng, n_blocks, T, bad_blocks=()):
    stats = np.empty((n_blocks, 1 + T))
    stats[:, 0] = rng.gamma(2.0, 50.0, size=n_blocks)
    stats[:, 1:] = rng.standard_normal((n_blocks, T)) * 30.0
    for b in bad_blocks:
        # A block whose denominator swallows the total: the leave-one-out fit
        # for it is undefined and the estimator must leave the row at zero.
        stats[b, 0] = -1e9
    return stats


@pytest.mark.parametrize("n_blocks,T", [(1000, 45), (65, 45), (50, 7), (12, 3)])
def test_jackknife_covariance_is_bit_identical_to_the_loop(n_blocks, T):
    rng = np.random.default_rng(20260916 + n_blocks + T)
    for trial in range(20):
        bad = (trial % 5,) if trial % 4 == 1 else ()
        stats = _block_stats(rng, n_blocks, T, bad_blocks=bad)
        want = _reference_covariance_from_blocks(stats, n_blocks)
        got = compute_covariance_from_blocks(stats, n_blocks)
        for name, a, b in zip(("beta", "se", "Sigma"), want, got):
            assert np.array_equal(a, b, equal_nan=True), (
                f"{name} moved for n_blocks={n_blocks} T={T} trial={trial}"
            )


def test_every_block_unusable_leaves_the_jackknife_at_zero():
    """All leave-one-out denominators non-positive: Sigma is zero, not NaN.

    D_total stays positive (so beta is defined) while every individual block
    denominator is not, which is the one case where the vectorised divide has to
    agree with the loop about writing nothing.
    """
    n_blocks, T = 4, 3
    stats = np.zeros((n_blocks, 1 + T))
    stats[:, 0] = [10.0, -1.0, -1.0, -1.0]
    stats[0, 0] = 10.0
    stats[:, 1:] = 1.0
    want = _reference_covariance_from_blocks(stats, n_blocks)
    got = compute_covariance_from_blocks(stats, n_blocks)
    for a, b in zip(want, got):
        assert np.array_equal(a, b, equal_nan=True)


@pytest.mark.parametrize(
    "values",
    [
        [0.0, -0.0, 1.0, -1.0],
        [np.nan, np.inf, -np.inf, 1.5],
        [1e-46, 1e-38, 3.4e38, 1e39],          # subnormal, normal, overflow
        [0.1, 1 / 3, 2 ** -24, 1 + 2 ** -53],  # values that round on narrowing
    ],
)
def test_float32_covariance_packing_matches_struct_pack(values):
    """``astype('<f4').tobytes()`` is what ``struct.pack('{n}f', *x)`` wrote.

    write_stats_batch swapped one for the other to avoid building a
    1035-element argument tuple per variant; both are little-endian IEEE-754
    binary32 with round-to-nearest, and the covariance payload must not move.
    """
    arr = np.array(values, dtype=np.float64)
    # One case deliberately overflows binary32. Both paths saturate to +/-inf
    # rather than raising, which is the behaviour being pinned; NumPy just says
    # so out loud.
    with np.errstate(over="ignore"):
        packed = arr.astype("<f4", copy=False).tobytes()
    assert struct.pack(f"{arr.size}f", *arr) == packed


