"""Tests for the explicit Stage 1 density-kernel choice."""

import numpy as np
import pytest

from fungwas_stage1.core import (
    DENSITY_METHOD_FOURTH_ORDER_GAUSSIAN,
    DENSITY_METHOD_SILVERMAN_GAUSSIAN,
    compute_rif,
)
from fungwas_stage1.pooled_cov import analysis_fingerprint, build_sigma0_meta


def test_default_density_method_is_the_historical_gaussian_path():
    rng = np.random.default_rng(20260915)
    y = rng.normal(size=401)
    taus = np.array([0.02, 0.5, 0.98])
    rif_default, q_default, f_default = compute_rif(y, taus)
    rif_gaussian, q_gaussian, f_gaussian = compute_rif(
        y, taus, density_method=DENSITY_METHOD_SILVERMAN_GAUSSIAN
    )
    assert np.array_equal(q_default, q_gaussian)
    assert np.array_equal(f_default, f_gaussian)
    assert np.array_equal(rif_default, rif_gaussian)


def test_fourth_order_kernel_changes_density_but_keeps_quantiles():
    rng = np.random.default_rng(20260915)
    y = rng.normal(size=2001)
    taus = np.linspace(0.02, 0.98, 9)
    _, q_gaussian, f_gaussian = compute_rif(y, taus)
    rif_fourth, q_fourth, f_fourth = compute_rif(
        y, taus, density_method=DENSITY_METHOD_FOURTH_ORDER_GAUSSIAN
    )
    assert np.array_equal(q_gaussian, q_fourth)
    assert np.all(np.isfinite(rif_fourth))
    assert np.all(np.isfinite(f_fourth))
    assert np.all(f_fourth > 0)
    assert np.any(np.abs(f_fourth - f_gaussian) > 1e-12)


def test_invalid_density_method_is_rejected():
    with pytest.raises(ValueError, match="density_method"):
        compute_rif(np.arange(10.0), np.array([0.5]), density_method="unknown")


def test_density_method_is_part_of_pooled_provenance_and_fingerprint():
    common = dict(
        phenotype="testosterone",
        covariates=["age"],
        taus=np.linspace(0.02, 0.98, 45),
        n_samples=1000,
        chromosomes=["ALL"],
    )
    gaussian = build_sigma0_meta(
        **common, density_method=DENSITY_METHOD_SILVERMAN_GAUSSIAN
    )
    fourth = build_sigma0_meta(
        **common, density_method=DENSITY_METHOD_FOURTH_ORDER_GAUSSIAN
    )
    assert gaussian["density_method"] == DENSITY_METHOD_SILVERMAN_GAUSSIAN
    assert fourth["density_method"] == DENSITY_METHOD_FOURTH_ORDER_GAUSSIAN
    assert gaussian["fingerprint"] != fourth["fingerprint"]
    assert gaussian["fingerprint"] == analysis_fingerprint(
        "testosterone", ["age"], common["taus"], 1000, None,
        density_method=DENSITY_METHOD_SILVERMAN_GAUSSIAN
    )
