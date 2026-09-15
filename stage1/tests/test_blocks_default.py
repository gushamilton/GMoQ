"""Tests for the jackknife block-count default and its guard rails.

The default moved from ``max(n_taus + 20, 64)`` (65 on the production 45-tau
grid) to 1000.  Under the fixed-RIF, regular-block approximation, what changes
is the noise in ``SE_hat``, whose design-variance approximation is
``1/sqrt(2*(B-1))`` -- 8.8% at B = 65 against 2.2% at B = 1000.  This is not a
universal calibration guarantee for the full plug-in RIF estimator.
"""

import logging

import numpy as np
import pytest

from fungwas_stage1 import cov_meta
from fungwas_stage1.core import (
    DEFAULT_N_BLOCKS,
    check_blocks_vs_samples,
    count_nonempty_blocks,
    resolve_n_blocks,
)


def test_default_block_count_is_1000():
    assert DEFAULT_N_BLOCKS == 1000
    assert resolve_n_blocks(None, n_taus=45) == 1000


def test_default_is_used_for_the_production_tau_grid():
    """The old default gave 65 here, which is the configuration that is 17x
    inflated at p = 5e-8 under a normal reference."""
    assert resolve_n_blocks(None, n_taus=45) != max(45 + 20, 64)


def test_blocks_remains_overridable():
    assert resolve_n_blocks(200, n_taus=45) == 200
    assert resolve_n_blocks(65, n_taus=45) == 65
    assert resolve_n_blocks(2500, n_taus=45) == 2500


def test_blocks_must_exceed_taus():
    with pytest.raises(ValueError, match="greater than the number of taus"):
        resolve_n_blocks(45, n_taus=45)
    with pytest.raises(ValueError, match="greater than the number of taus"):
        resolve_n_blocks(10, n_taus=45)
    assert resolve_n_blocks(46, n_taus=45) == 46


def test_below_default_warns_about_se_noise(caplog):
    with caplog.at_level(logging.WARNING, logger="fungwas_stage1.core"):
        resolve_n_blocks(65, n_taus=45)
    text = caplog.text
    assert "below the default" in text
    # 1/sqrt(2*64) = 8.84%
    assert "8.8%" in text


def test_at_or_above_default_does_not_warn(caplog):
    with caplog.at_level(logging.WARNING, logger="fungwas_stage1.core"):
        resolve_n_blocks(None, n_taus=45)
        resolve_n_blocks(2500, n_taus=45)
    assert "below the default" not in caplog.text


def test_more_blocks_than_samples_is_refused():
    """Reachable only now that the default is 1000: at B = 65 no realistic
    sample was smaller than the block count."""
    with pytest.raises(ValueError, match="exceeds the number of analysed samples"):
        check_blocks_vs_samples(DEFAULT_N_BLOCKS, n_samples=500)
    check_blocks_vs_samples(DEFAULT_N_BLOCKS, n_samples=50000)


def test_empty_blocks_are_counted_and_reduce_nu(caplog):
    """Deleting an empty block reproduces the full estimate, so empty blocks
    contribute no replicate and must not be counted towards nu."""
    rng = np.random.default_rng(1)
    n_blocks = 200
    block_ids = rng.integers(0, n_blocks, 300)
    with caplog.at_level(logging.WARNING, logger="fungwas_stage1.core"):
        n_nonempty = count_nonempty_blocks(block_ids, n_blocks)
    assert n_nonempty == len(np.unique(block_ids))
    assert n_nonempty < n_blocks
    assert "empty" in caplog.text

    meta = cov_meta.build_cov_meta(
        cov_dtype="int8", n_taus=45, n_blocks=n_blocks, n_blocks_nonempty=n_nonempty
    )
    assert meta["jackknife"]["n_blocks"] == n_blocks
    assert meta["jackknife"]["n_blocks_nonempty"] == n_nonempty
    assert cov_meta.resolve_block_count(meta) == n_nonempty


def test_no_empty_blocks_at_production_scale(caplog):
    rng = np.random.default_rng(2)
    block_ids = rng.integers(0, DEFAULT_N_BLOCKS, 50000)
    with caplog.at_level(logging.WARNING, logger="fungwas_stage1.core"):
        n_nonempty = count_nonempty_blocks(block_ids, DEFAULT_N_BLOCKS)
    assert n_nonempty == DEFAULT_N_BLOCKS
    assert "empty" not in caplog.text


def test_nonempty_count_cannot_exceed_requested():
    with pytest.raises(ValueError, match="cannot exceed"):
        cov_meta.build_cov_meta(
            cov_dtype="int8", n_taus=45, n_blocks=100, n_blocks_nonempty=101
        )


def test_default_block_count_reaches_the_metadata():
    meta = cov_meta.build_cov_meta(
        cov_dtype="int8", n_taus=45, n_blocks=resolve_n_blocks(None, 45), seed=42
    )
    assert cov_meta.resolve_block_count(meta) == 1000
    # 1/(B-1) is only the first term; the realised nu is lower and SNP-specific.
    assert cov_meta.effective_nu(1000, 9134, cov_meta.C_LOCATION) < 999
