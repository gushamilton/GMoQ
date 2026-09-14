"""The pooled cross-quantile covariance, and the refusal that protects it.

``Sigma_v = Sigma_0 / d_resid`` is only the sampling covariance of ``beta_tau``
when ``Sigma_0`` was built from the same phenotype, covariate set, tau grid,
LOCO predictor and sample as the scan it is applied to.  Getting that wrong is
silent: the point estimates are untouched and the standard errors stay
plausible.  So half of these tests are about the arithmetic and half are about
the refusal.

The statistical claim itself -- that the pooled estimator reproduces the exact
per-individual HC0 sandwich -- is tested here on simulated data at n = 20,000
and was measured on the production chr22 inputs at n = 461,320 (median diagonal
ratio 1.00000, median correlation deviation 0.0053, both traits).
"""

import gzip
import json

import numpy as np
import pytest

from fungwas_stage1 import cli, core, cov_meta, pooled_cov


# ---------------------------------------------------------------------------
# the arithmetic
# ---------------------------------------------------------------------------

def test_compute_sigma0_is_the_gram_matrix_over_n():
    rng = np.random.default_rng(1)
    Psi = rng.normal(size=(500, 7))
    S = pooled_cov.compute_sigma0(Psi)
    assert S.shape == (7, 7)
    assert np.allclose(S, Psi.T @ Psi / 500.0, rtol=0, atol=1e-12)
    # Exactly symmetric, not just symmetric to rounding: Stage 2 reads the
    # upper triangle and must not depend on which one.
    assert np.array_equal(S, S.T)
    # PSD by construction -- one of the reasons this replaces the jackknife,
    # which needs a PSD repair and an eigenvalue sign rule.
    assert np.min(np.linalg.eigvalsh(S)) >= -1e-12


def test_compute_sigma0_rejects_degenerate_input():
    with pytest.raises(ValueError):
        pooled_cov.compute_sigma0(np.zeros(5))
    with pytest.raises(ValueError):
        pooled_cov.compute_sigma0(np.zeros((1, 4)))


def test_sigma_for_variant_rejects_a_non_positive_scale():
    S = np.eye(3)
    assert np.allclose(pooled_cov.sigma_for_variant(S, 4.0), np.eye(3) / 4.0)
    for bad in (0.0, -1.0, np.nan, np.inf):
        with pytest.raises(ValueError):
            pooled_cov.sigma_for_variant(S, bad)


def test_pooled_covariance_reproduces_the_exact_hc0_sandwich():
    """The estimator claim, on data where the exact answer can be computed.

    ``S_hc0 = sum_i g_i^2 e_i e_i' / (g'g)^2`` with ``e_i`` the per-individual
    residual is the textbook robust variance of the partialled-out regression.
    The pooled estimator is its working-independence form, exact when ``g_i^2``
    is uncorrelated with ``e_i e_i'``, which is what an i.i.d. sample gives.
    """
    rng = np.random.default_rng(20260914)
    n, T, K = 20000, 9, 4
    taus = np.linspace(0.1, 0.9, T)

    X = np.column_stack([np.ones(n), rng.normal(size=(n, K - 1))])
    y = X @ np.array([1.0, 0.4, -0.2, 0.1]) + rng.gamma(3.0, 1.5, n)
    rif, _, _ = core.compute_rif(y, taus)
    Q, _ = np.linalg.qr(X)
    Psi_r = rif - Q @ (Q.T @ rif)
    Sigma0 = pooled_cov.compute_sigma0(Psi_r)

    ratios = []
    corr_dev = []
    for _ in range(12):
        g = rng.binomial(2, rng.uniform(0.05, 0.5), size=n).astype(float)
        gr = g - Q @ (Q.T @ g)
        D = float(gr @ gr)
        b = (gr @ Psi_r) / D
        e = Psi_r - np.outer(gr, b)
        S_hc0 = (e * (gr * gr)[:, None]).T @ e / D ** 2
        S_pool = pooled_cov.sigma_for_variant(Sigma0, D)
        ratios.append(np.diag(S_pool) / np.diag(S_hc0))

        def corr(S):
            d = np.sqrt(np.diag(S))
            return S / np.outer(d, d)

        corr_dev.append(np.max(np.abs(corr(S_pool) - corr(S_hc0))))

    ratios = np.concatenate(ratios)
    assert abs(np.median(ratios) - 1.0) < 0.01, np.median(ratios)
    assert np.max(np.abs(ratios - 1.0)) < 0.10, np.max(np.abs(ratios - 1.0))
    assert np.median(corr_dev) < 0.05, np.median(corr_dev)


# ---------------------------------------------------------------------------
# the on-disk format
# ---------------------------------------------------------------------------

def _meta(**kw):
    base = dict(phenotype="trait", covariates=["age", "sex"],
                taus=np.linspace(0.1, 0.9, 5), n_samples=1234,
                chromosomes=["22"], loco_predictor=None)
    base.update(kw)
    return pooled_cov.build_sigma0_meta(**base)


def test_round_trip_is_bit_exact(tmp_path):
    rng = np.random.default_rng(7)
    A = rng.normal(size=(400, 5))
    S = pooled_cov.compute_sigma0(A)
    meta = _meta()
    prefix = str(tmp_path / "out")
    pooled_cov.write_sigma0(prefix, {"22": S}, meta)

    back_meta, back = pooled_cov.read_sigma0(prefix)
    # %.17g round-trips float64 exactly. A lossy sidecar would reintroduce the
    # quantisation error the int8 covariance already costs (df_eff 0.85-0.93 of
    # truth on the published files).
    assert np.array_equal(back["22"], S)
    assert back_meta["fingerprint"] == meta["fingerprint"]
    assert back_meta["matrix_file"].endswith(".sigma0.tsv.gz")


def test_write_rejects_a_missing_or_misshapen_matrix(tmp_path):
    meta = _meta(chromosomes=["21", "22"])
    S = np.eye(5)
    with pytest.raises(ValueError, match="no matrix"):
        pooled_cov.write_sigma0(str(tmp_path / "a"), {"22": S}, meta)
    with pytest.raises(ValueError, match="shape"):
        pooled_cov.write_sigma0(str(tmp_path / "b"), {"21": np.eye(4), "22": S}, meta)
    with pytest.raises(ValueError, match="finite"):
        bad = np.eye(5) * np.nan
        pooled_cov.write_sigma0(str(tmp_path / "c"), {"21": bad, "22": S}, meta)


# ---------------------------------------------------------------------------
# the refusal
# ---------------------------------------------------------------------------

def test_fingerprint_moves_with_everything_sigma0_depends_on():
    base = _meta()["fingerprint"]
    assert _meta()["fingerprint"] == base  # stable
    assert _meta(phenotype="other")["fingerprint"] != base
    assert _meta(covariates=["age", "sex", "PC1"])["fingerprint"] != base
    assert _meta(covariates=["sex", "age"])["fingerprint"] != base  # order matters
    assert _meta(taus=np.linspace(0.1, 0.9, 6))["fingerprint"] != base
    assert _meta(n_samples=1235)["fingerprint"] != base
    assert _meta(loco_predictor="/p/bmi.loco")["fingerprint"] != base
    # ...and not with anything it does not depend on.
    assert _meta(chromosomes=["1", "22"])["fingerprint"] == base


def test_cov_meta_carries_the_analysis_block_and_check_passes():
    meta = _meta()
    cm = cov_meta.build_cov_meta(
        cov_dtype="int8", n_taus=5, n_blocks=1000, seed=42,
        analysis=pooled_cov.sigma0_analysis_block(meta),
    )
    assert cm["cov_meta_version"] == 4
    assert cov_meta.resolve_analysis(cm)["fingerprint"] == meta["fingerprint"]
    pooled_cov.check_fingerprint(meta, cm)  # does not raise


@pytest.mark.parametrize("field,value", [
    ("phenotype", "other"),
    ("covariates", ["age"]),
    ("n_samples", 9999),
    ("loco_predictor", "/p/bmi.loco"),
])
def test_check_fingerprint_refuses_a_mismatched_scan(field, value):
    scan = cov_meta.build_cov_meta(
        cov_dtype="int8", n_taus=5, n_blocks=1000,
        analysis=pooled_cov.sigma0_analysis_block(_meta()),
    )
    other = _meta(**{field: value})
    with pytest.raises(ValueError, match="does not match the scan"):
        pooled_cov.check_fingerprint(other, scan)


def test_check_fingerprint_refuses_metadata_with_no_analysis_block():
    old = cov_meta.build_cov_meta(cov_dtype="int8", n_taus=5, n_blocks=1000)
    old.pop("analysis", None)
    with pytest.raises(ValueError, match="no 'analysis' block"):
        pooled_cov.check_fingerprint(_meta(), old)


def test_check_taus_refuses_a_different_grid():
    meta = _meta()
    pooled_cov.check_taus(meta, np.linspace(0.1, 0.9, 5))
    with pytest.raises(ValueError, match="tau grid"):
        pooled_cov.check_taus(meta, np.linspace(0.1, 0.9, 6))
    with pytest.raises(ValueError, match="tau grid"):
        pooled_cov.check_taus(meta, np.linspace(0.05, 0.95, 5))


# ---------------------------------------------------------------------------
# end to end through main()
# ---------------------------------------------------------------------------

N_KEEP = 900
N_GENO = 700


class _FakeVariant:
    def __init__(self, rsid, pos, alt_dosage, chrom="22"):
        self.rsid = rsid
        self.chrom = chrom
        self.pos = pos
        self.alleles = ("A", "G")
        self.alt_dosage = alt_dosage


def _inputs(tmp_path, seed=11):
    rng = np.random.default_rng(seed)
    ids = [f"ID{i:05d}" for i in range(N_KEEP)]
    age = rng.normal(57.0, 8.0, N_KEEP)
    sex = rng.integers(0, 2, N_KEEP).astype(float)
    y = 0.02 * age + 0.5 * sex + rng.gamma(4.0, 1.2, N_KEEP)
    with open(tmp_path / "pheno.txt", "w") as fh:
        fh.write("FID IID trait\n")
        for sid, v in zip(ids, y):
            fh.write(f"{sid} {sid} {v:.10g}\n")
    with open(tmp_path / "covar.txt", "w") as fh:
        fh.write("FID IID age sex\n")
        for sid, a, s in zip(ids, age, sex):
            fh.write(f"{sid} {sid} {a:.10g} {s:.10g}\n")
    with open(tmp_path / "keep.txt", "w") as fh:
        fh.write("\n".join(ids) + "\n")
    return ids, y


def _run_main(tmp_path, monkeypatch, ids, extra_argv=(), chroms=("22",)):
    bgen_samples = ids[:N_GENO]
    rng = np.random.default_rng(5)
    variants = []
    for chrom in chroms:
        variants.extend(
            _FakeVariant(f"rs{chrom}_{j}", 1000 + j,
                         rng.binomial(2, 0.3, size=N_GENO).astype(float), chrom)
            for j in range(4)
        )
    monkeypatch.setattr(cli, "get_variants_bgen",
                        lambda *a, **k: (variants, bgen_samples, None))
    out = tmp_path / "o"
    monkeypatch.setattr("sys.argv", [
        "fungwas-stage1", "--bgen", str(tmp_path / "fake.bgen"),
        "--pheno", str(tmp_path / "pheno.txt"), "--pheno-col", "trait",
        "--covar", str(tmp_path / "covar.txt"), "--keep", str(tmp_path / "keep.txt"),
        "--taus", "0.25,0.5,0.75", "--blocks", "20", "--batch-size", "3",
        "--out", str(out), *extra_argv,
    ])
    cli.main()
    return str(out)


def test_main_writes_a_sigma0_sidecar_matching_the_analysed_rows(tmp_path, monkeypatch):
    ids, _ = _inputs(tmp_path)
    out = _run_main(tmp_path, monkeypatch, ids)

    meta, mats = pooled_cov.read_sigma0(out)
    assert list(mats) == [pooled_cov.SIGMA0_ALL_CHROM]
    assert mats[pooled_cov.SIGMA0_ALL_CHROM].shape == (3, 3)
    # n is the genotype-matched count, not the keep count. Sigma_0 has to be
    # built on the rows that are actually swept -- the same requirement, and the
    # same failure mode, as the covariate basis (fixed at 4754794).
    assert meta["n_samples"] == N_GENO
    assert meta["covariates"] == ["age", "sex"]
    assert meta["taus"] == [0.25, 0.5, 0.75]
    assert meta["loco_predictor"] is None
    assert meta["estimator"] == pooled_cov.SIGMA0_ESTIMATOR

    with open(f"{out}.cov.meta.json") as fh:
        cm = json.load(fh)
    pooled_cov.check_fingerprint(meta, cm)


def test_sidecar_reproduces_the_tsv_standard_errors(tmp_path, monkeypatch):
    """The whole point, checked end to end: Sigma_0 / d_resid is the covariance.

    Against the jackknife written in the same run, so the comparison is on real
    Stage 1 output rather than on a reimplementation.  The tolerance is loose
    because the jackknife is the noisy one: at 20 blocks its own relative error
    on an SE is ~16%.
    """
    ids, _ = _inputs(tmp_path)
    out = _run_main(tmp_path, monkeypatch, ids)

    _, mats = pooled_cov.read_sigma0(out)
    S0 = mats[pooled_cov.SIGMA0_ALL_CHROM]
    with gzip.open(f"{out}.stage1.tsv.gz", "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [ln.rstrip("\n").split("\t") for ln in fh]
    se_idx = [i for i, c in enumerate(header) if c.startswith("se_tau_")]
    d_idx = header.index("d_resid")

    ratios = []
    for row in rows:
        d = float(row[d_idx])
        se_pool = np.sqrt(np.diag(S0) / d)
        se_jk = np.array([float(row[i]) for i in se_idx])
        ratios.append(se_pool / se_jk)
    ratios = np.concatenate(ratios)
    assert 0.7 < np.median(ratios) < 1.4, np.median(ratios)


def test_no_sigma0_suppresses_the_sidecar(tmp_path, monkeypatch):
    ids, _ = _inputs(tmp_path)
    out = _run_main(tmp_path, monkeypatch, ids, extra_argv=("--no-sigma0",))
    matrix_path, meta_path = pooled_cov.sigma0_paths(out)
    import os
    assert not os.path.exists(matrix_path)
    assert not os.path.exists(meta_path)
    # The jackknife path is untouched: this must stay a two-way door.
    assert os.path.exists(f"{out}.cov.gz")
    with open(f"{out}.cov.meta.json") as fh:
        cm = json.load(fh)
    assert cm["jackknife"]["n_blocks"] == 20


def test_loco_run_writes_one_matrix_per_chromosome(tmp_path, monkeypatch):
    ids, y = _inputs(tmp_path)
    loco = tmp_path / "pred.loco"
    rng = np.random.default_rng(3)
    with open(loco, "w") as fh:
        fh.write("FID " + " ".join(ids) + "\n")
        for chrom in ("21", "22"):
            vals = rng.normal(0.0, 0.3, len(ids))
            fh.write(f"{chrom} " + " ".join(f"{v:.8g}" for v in vals) + "\n")
    out = _run_main(tmp_path, monkeypatch, ids, chroms=("21", "22"),
                    extra_argv=("--loco-preds", str(loco)))

    meta, mats = pooled_cov.read_sigma0(out)
    # Per chromosome because the LOCO predictor changes the phenotype, and the
    # influence functions with it. One matrix for the whole genome would be the
    # wrong object for every chromosome but one.
    assert sorted(mats) == ["21", "22"]
    assert meta["chromosomes"] == ["21", "22"]
    assert meta["loco_predictor"] == str(loco)
    assert not np.allclose(mats["21"], mats["22"])


# --- the LOCO scale must be recorded and must enter the fingerprint ---------

def test_loco_scale_is_recorded_in_the_sidecar():
    """
    A scan made under the pre-v0.4.0 unscaled behaviour must be distinguishable
    from a corrected one. Effects differ by a factor of 2.2 on log_bmi, so a
    sidecar that does not say which mode was used cannot be audited later.
    """
    from fungwas_stage1 import pooled_cov
    meta = pooled_cov.build_sigma0_meta(
        phenotype="log_bmi", covariates=["age", "sex"],
        taus=[0.25, 0.5, 0.75], n_samples=1000, chromosomes=["22"],
        loco_predictor="/path/to/log_bmi.loco",
        loco_scale_mode="regression", loco_scale_factor=0.160696,
    )
    assert meta["loco_scale_mode"] == "regression"
    assert meta["loco_scale_factor"] == pytest.approx(0.160696)
    # and it must travel in the block Stage 2 compares
    blk = pooled_cov.sigma0_analysis_block(meta)
    assert blk["loco_scale_mode"] == "regression"
    assert blk["loco_scale_factor"] == pytest.approx(0.160696)


def test_loco_scale_changes_the_fingerprint():
    """
    Two runs identical but for --loco-scale analyse different variables, so a
    Sigma_0 from one is invalid for the other. Without this they hash alike and
    the mismatch guard would pass them.
    """
    from fungwas_stage1 import pooled_cov
    kw = dict(phenotype="log_bmi", covariates=["age"], taus=[0.25, 0.5, 0.75],
              n_samples=1000, loco_predictor="/p/log_bmi.loco")
    fp_reg = pooled_cov.analysis_fingerprint(**kw, loco_scale_mode="regression",
                                             loco_scale_factor=0.1607)
    fp_none = pooled_cov.analysis_fingerprint(**kw, loco_scale_mode="none",
                                              loco_scale_factor=1.0)
    fp_sd = pooled_cov.analysis_fingerprint(**kw, loco_scale_mode="sd",
                                            loco_scale_factor=0.1665)
    assert len({fp_reg, fp_none, fp_sd}) == 3


def test_fingerprint_unchanged_when_no_loco_adjustment():
    """Backward compatibility: a run with no LOCO keeps the hash it always had."""
    from fungwas_stage1 import pooled_cov
    kw = dict(phenotype="bmi", covariates=["age"], taus=[0.25, 0.5, 0.75],
              n_samples=1000, loco_predictor=None)
    assert pooled_cov.analysis_fingerprint(**kw) == \
           pooled_cov.analysis_fingerprint(**kw, loco_scale_mode=None,
                                           loco_scale_factor=None)
