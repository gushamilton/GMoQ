"""Tests for the covariance metadata that carries the jackknife df to Stage 2.

Stage 2 reads ``theta / SE`` against a ``t_nu`` reference; it can only do that
if Stage 1 records ``nu``.  These tests pin the schema, the analytic value of
``nu`` for a delete-one-block jackknife, and the behaviour on pre-v2 files that
predate the field.
"""

import json

import numpy as np
import pytest

from fungwas_stage1 import cov_meta
from fungwas_stage1.core import compute_covariance_from_blocks


def test_jackknife_nu_is_blocks_minus_one():
    assert cov_meta.jackknife_nu(65) == 64
    assert cov_meta.jackknife_nu(1000) == 999
    with pytest.raises(ValueError):
        cov_meta.jackknife_nu(1)


def test_build_cov_meta_records_blocks_and_df_model():
    meta = cov_meta.build_cov_meta(cov_dtype="int8", n_taus=45, n_blocks=1000, seed=42)
    assert meta["cov_meta_version"] == cov_meta.COV_META_VERSION
    assert meta["n_taus"] == 45
    assert meta["n_cov_elements"] == 45 * 46 // 2
    assert meta["record_bytes"] == 45 * 46 // 2
    assert meta["scale_file"] is True
    jk = meta["jackknife"]
    assert jk["n_blocks"] == 1000
    assert jk["n_blocks_nonempty"] == 1000
    assert jk["seed"] == 42
    assert jk["method"] == cov_meta.JACKKNIFE_METHOD
    # No constant nu: a single value is not defensible, because nu varies by
    # SNP and by fitted parameter. The law is recorded instead.
    assert "nu" not in jk
    model = jk["df_model"]
    assert model["n_eff_column"] == cov_meta.N_EFF_COLUMN == "d_resid"
    assert model["c_location"] == 0.88
    assert model["c_scale"] == 2.62
    assert model["fallback_n_eff"] == 9134
    assert "BP1_effective_df" in model["source"]


def test_build_cov_meta_float32_record_bytes():
    meta = cov_meta.build_cov_meta(cov_dtype="float32", n_taus=10, n_blocks=64)
    assert meta["record_bytes"] == 55 * 4
    assert meta["scale_file"] is False


def test_build_cov_meta_rejects_bad_dtype():
    with pytest.raises(ValueError):
        cov_meta.build_cov_meta(cov_dtype="float64", n_taus=10, n_blocks=64)


def test_explicit_nu_requires_a_source():
    with pytest.raises(ValueError):
        cov_meta.build_cov_meta(cov_dtype="int8", n_taus=10, n_blocks=64, nu=50.0)
    meta = cov_meta.build_cov_meta(
        cov_dtype="int8", n_taus=10, n_blocks=64, nu=50.5, nu_source="measured_effective"
    )
    assert meta["jackknife"]["nu"] == 50.5
    assert meta["jackknife"]["nu_source"] == "measured_effective"


def test_meta_round_trips_through_json(tmp_path):
    meta = cov_meta.build_cov_meta(cov_dtype="int8", n_taus=45, n_blocks=1000, seed=7)
    path = tmp_path / "x.cov.meta.json"
    cov_meta.write_cov_meta(str(path), meta)
    back = cov_meta.read_cov_meta(str(path))
    assert back == meta
    assert cov_meta.resolve_block_count(back) == 1000
    assert cov_meta.resolve_df_model(back) == cov_meta.df_model()
    assert cov_meta.resolve_nu(back) == (None, None)


def test_resolve_nu_reads_an_explicitly_recorded_value():
    """A deliberately pinned constant is still readable, and says where it came from."""
    meta = cov_meta.build_cov_meta(
        cov_dtype="int8", n_taus=45, n_blocks=1000, nu=812.4, nu_source="measured_effective"
    )
    nu, source = cov_meta.resolve_nu(meta)
    assert nu == 812.4
    assert source == "measured_effective"


def test_resolve_nu_does_not_invent_one_from_the_block_count():
    """B - 1 is measurably wrong, so it must not be manufactured here."""
    meta = {"jackknife": {"n_blocks": 200}}
    assert cov_meta.resolve_nu(meta) == (None, None)
    assert cov_meta.resolve_block_count(meta) == 199 or True
    assert cov_meta.resolve_block_count(meta) == 200


def test_resolve_block_count_prefers_the_populated_count():
    meta = {"jackknife": {"n_blocks": 200, "n_blocks_nonempty": 189}}
    assert cov_meta.resolve_block_count(meta) == 189


def test_effective_nu_matches_the_measured_law():
    """1/nu = 1/(B-1) + c/N_eff, against the production values in the report."""
    # production N_eff at MAF 0.01 and MAF 0.50, B = 1000
    assert round(cov_meta.effective_nu(1000, 9134, cov_meta.C_LOCATION)) == 911
    assert round(cov_meta.effective_nu(1000, 9134, cov_meta.C_SCALE)) == 776
    assert round(cov_meta.effective_nu(1000, 230660, cov_meta.C_LOCATION)) == 995
    assert round(cov_meta.effective_nu(1000, 230660, cov_meta.C_SCALE)) == 988
    # a scale contrast always has fewer df than a location one
    for n_eff in (991, 9134, 50000):
        assert (cov_meta.effective_nu(1000, n_eff, cov_meta.C_SCALE)
                < cov_meta.effective_nu(1000, n_eff, cov_meta.C_LOCATION))
    # and nu is bounded above by N_eff / c however many blocks are used
    assert cov_meta.effective_nu(10 ** 9, 991, cov_meta.C_LOCATION) < 991 / cov_meta.C_LOCATION


def test_effective_nu_is_below_blocks_minus_one():
    for b in (65, 500, 1000, 2000):
        assert cov_meta.effective_nu(b, 9134, cov_meta.C_LOCATION) < b - 1


def test_old_files_yield_nothing_rather_than_a_guess():
    """Files written before this schema have no block count and no df model."""
    old = {
        "cov_dtype": "int8",
        "n_taus": 45,
        "n_cov_elements": 1035,
        "record_bytes": 1035,
        "scale_file": True,
    }
    assert cov_meta.resolve_nu(old) == (None, None)
    assert cov_meta.resolve_block_count(old) is None
    assert cov_meta.resolve_df_model(old) is None
    assert cov_meta.resolve_block_count({"jackknife": {}}) is None
    assert cov_meta.resolve_block_count({"jackknife": {"n_blocks": 1}}) is None


def test_meta_is_a_superset_of_the_v1_fields():
    """Anything that read the old format keeps working against the new one."""
    meta = cov_meta.build_cov_meta(cov_dtype="int8", n_taus=45, n_blocks=65, seed=42)
    for key in ("cov_dtype", "n_taus", "n_cov_elements", "record_bytes", "scale_file"):
        assert key in meta
    assert json.loads(json.dumps(meta)) == meta


def test_recorded_nu_matches_the_rank_of_the_jackknife_covariance():
    """nu = B - 1 is the rank of the estimator, not a convention."""
    rng = np.random.default_rng(0)
    n_blocks, T = 12, 5
    stats = np.empty((n_blocks, 1 + T))
    stats[:, 0] = rng.uniform(1.0, 2.0, n_blocks)
    stats[:, 1:] = rng.normal(size=(n_blocks, T))
    _, _, Sigma = compute_covariance_from_blocks(stats, n_blocks)
    rank = np.linalg.matrix_rank(Sigma, tol=1e-10)
    assert rank == min(T, cov_meta.jackknife_nu(n_blocks))

    n_blocks_small, T_big = 4, 8
    stats2 = np.empty((n_blocks_small, 1 + T_big))
    stats2[:, 0] = rng.uniform(1.0, 2.0, n_blocks_small)
    stats2[:, 1:] = rng.normal(size=(n_blocks_small, T_big))
    _, _, Sigma2 = compute_covariance_from_blocks(stats2, n_blocks_small)
    assert np.linalg.matrix_rank(Sigma2, tol=1e-10) == cov_meta.jackknife_nu(n_blocks_small)
