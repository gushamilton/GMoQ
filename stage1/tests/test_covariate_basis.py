"""The covariate basis must be built on the rows that are actually analysed.

Stage 1 residualises every genotype as ``g_resid = g - Q (Q' g)``. That is the
orthogonal projection onto the complement of the covariate span *only if*
``Q'Q = I`` over the rows being swept. Historically the QR was computed on the
keep sample and then row-subset to the samples that matched the genotype file,
which does not preserve orthonormality: on the UK Biobank BMI inputs (486,254
keep rows, 461,320 matched) the subset basis reached ``max |Q'Q - I| = 0.97`` and
shifted ``beta_tau`` by a median of 0.12 SE, up to 2.3 SE, in a consistent
direction at all 45 taus, while leaving ``se_tau`` unchanged - a silent bias.

These tests pin the ordering (match first, then build) and the guard that makes
a mis-built basis fail loudly instead of quietly biasing the output.
"""

import gzip
import numpy as np
import pytest

from fungwas_stage1 import cli, core


# ---------------------------------------------------------------------------
# fixtures: a keep sample strictly larger than the genotype sample
# ---------------------------------------------------------------------------

N_KEEP = 4000
N_GENO = 2500
N_PC = 6


def _make_sample(seed=20260914):
    """Keep sample of N_KEEP; only the first N_GENO are in the genotype file.

    The unmatched rows are given a shifted, inflated PC distribution. This is not
    decoration: in the real case the genotype release is ancestry-restricted while
    the keep file is built on phenotype availability alone, so the unmatched
    samples carry most of the ancestry signal and the PC columns of a row-subset
    basis collapse to a few percent of their norm. A uniformly random split would
    make the defect look like a benign ~n/N rescaling.
    """
    rng = np.random.default_rng(seed)
    ids = [f"ID{i:06d}" for i in range(N_KEEP)]
    matched = np.zeros(N_KEEP, dtype=bool)
    matched[:N_GENO] = True

    age = rng.normal(57.0, 8.0, N_KEEP)
    sex = rng.integers(0, 2, N_KEEP).astype(float)
    pcs = rng.normal(0.0, 1.0, size=(N_KEEP, N_PC))
    # Ancestry-like structure: the unmatched rows dominate the leading PCs.
    pcs[~matched, :3] += 40.0 * rng.normal(1.0, 0.3, size=((~matched).sum(), 3))
    pcs[~matched, :3] *= 25.0

    y = (0.02 * age + 0.5 * sex + 0.01 * pcs[:, 0]
         + rng.gamma(shape=4.0, scale=1.2, size=N_KEEP))

    X = np.column_stack([np.ones(N_KEEP), age, sex, pcs])
    return ids, matched, X, y


def _write_inputs(tmp_path, ids, X, y):
    covar_names = ["age", "sex"] + [f"PC{i + 1}" for i in range(N_PC)]
    pheno = tmp_path / "pheno.txt"
    covar = tmp_path / "covar.txt"
    keep = tmp_path / "keep.txt"
    with open(pheno, "w") as fh:
        fh.write("FID IID trait\n")
        for sid, val in zip(ids, y):
            fh.write(f"{sid} {sid} {val:.10g}\n")
    with open(covar, "w") as fh:
        fh.write("FID IID " + " ".join(covar_names) + "\n")
        for i, sid in enumerate(ids):
            vals = " ".join(f"{v:.10g}" for v in X[i, 1:])
            fh.write(f"{sid} {sid} {vals}\n")
    with open(keep, "w") as fh:
        fh.write("\n".join(ids) + "\n")
    return str(pheno), str(covar), str(keep)


def _max_gram_deviation(Q):
    G = np.asarray(Q, dtype=np.float64).T @ np.asarray(Q, dtype=np.float64)
    return float(np.max(np.abs(G - np.eye(G.shape[0]))))


# ---------------------------------------------------------------------------
# the defect itself, and the corrected construction
# ---------------------------------------------------------------------------

def test_row_subset_of_a_full_sample_basis_is_not_orthonormal():
    _, matched, X, _ = _make_sample()

    Q_full, _ = np.linalg.qr(X)
    assert _max_gram_deviation(Q_full) < 1e-12  # orthonormal over all N_KEEP rows

    Q_subset = Q_full[matched]
    # Not a small perturbation: the leading columns lose most of their norm,
    # exactly as measured on the production covariate matrix.
    assert _max_gram_deviation(Q_subset) > 0.1

    with pytest.raises(ValueError, match="not orthonormal"):
        core.check_orthonormal_basis(Q_subset)


def test_guard_tolerance_separates_qr_noise_from_a_single_dropped_row():
    """ORTHONORMAL_BASIS_TOL must sit between float64 QR noise and the smallest
    realisable instance of the defect - a keep sample one row larger than the
    genotype sample."""
    rng = np.random.default_rng(11)
    n, k = 50_000, 13
    X = np.column_stack([np.ones(n), rng.normal(size=(n, k - 1)) * rng.uniform(1, 500, k - 1)])
    Q, _ = np.linalg.qr(X)

    noise = core.check_orthonormal_basis(Q)
    assert noise < core.ORTHONORMAL_BASIS_TOL / 1e4, noise

    drop_one = np.ones(n, dtype=bool)
    drop_one[int(np.argmax(np.abs(Q[:, 1])))] = False
    assert _max_gram_deviation(Q[drop_one]) > core.ORTHONORMAL_BASIS_TOL

    # ... and even dropping an unremarkable row is caught.
    drop_first = np.ones(n, dtype=bool)
    drop_first[0] = False
    assert _max_gram_deviation(Q[drop_first]) > core.ORTHONORMAL_BASIS_TOL


@pytest.mark.parametrize("call", ["scores", "multi", "rss"])
def test_kernel_entry_points_reject_a_row_subset_basis(call):
    _, matched, X, _ = _make_sample()
    rng = np.random.default_rng(3)
    n = int(matched.sum())
    G = rng.binomial(2, 0.3, size=(n, 4)).astype(float)
    rif = rng.normal(size=(n, 3))
    block_ids = rng.integers(0, 8, n).astype(np.int32)
    Q_bad = np.linalg.qr(X)[0][matched]

    with pytest.raises(ValueError, match="not orthonormal"):
        if call == "scores":
            core.compute_block_scores(G, rif, Q_bad, block_ids, 8)
        elif call == "multi":
            core.compute_block_scores_multi(G, [rif], Q_bad, block_ids, 8)
        else:
            core.residual_sum_squares(G, Q_bad)


# ---------------------------------------------------------------------------
# the pipeline: match first, then build
# ---------------------------------------------------------------------------

def test_compute_chr_rif_residuals_builds_on_the_matched_rows():
    _, matched, X, y = _make_sample()
    taus = np.linspace(0.1, 0.9, 5)

    rif_list, Q, _ = cli.compute_chr_rif_residuals(
        ["trait"], {"trait": y}, X, taus, None, None, matched
    )

    n = int(matched.sum())
    assert Q.shape == (n, X.shape[1])
    assert rif_list[0].shape == (n, len(taus))
    assert core.check_orthonormal_basis(Q) < 1e-12


class _FakeVariant:
    def __init__(self, rsid, pos, alt_dosage):
        self.rsid = rsid
        self.chrom = "22"
        self.pos = pos
        self.alleles = ("A", "G")
        self.alt_dosage = alt_dosage
        self.bit_depth = 8


def test_cli_builds_the_projector_on_the_genotype_matched_rows(tmp_path, monkeypatch):
    """End-to-end through main(): keep sample of 4000, genotype sample of 2500.

    This is the ordering test. The unit tests above would still pass if main()
    went back to building the basis on the keep sample and row-subsetting it.
    """
    ids, matched, X, y = _make_sample()
    pheno, covar, keep = _write_inputs(tmp_path, ids, X, y)

    bgen_samples = [sid for sid, m in zip(ids, matched) if m]
    rng = np.random.default_rng(99)
    variants = [
        _FakeVariant(f"rs{j}", 1000 + j, rng.binomial(2, 0.3, size=len(bgen_samples)).astype(float))
        for j in range(6)
    ]
    monkeypatch.setattr(cli, "get_variants_bgen",
                        lambda *a, **k: (variants, bgen_samples, None))

    seen = []
    real_process = cli._process_batch

    def spy(G_view, RIF_resid_list, Q, *a, **k):
        seen.append((G_view.shape, Q.shape, _max_gram_deviation(Q), RIF_resid_list[0].shape))
        return real_process(G_view, RIF_resid_list, Q, *a, **k)

    monkeypatch.setattr(cli, "_process_batch", spy)

    out = tmp_path / "o"
    monkeypatch.setattr("sys.argv", [
        "fungwas-stage1", "--bgen", str(tmp_path / "fake.bgen"),
        "--pheno", pheno, "--pheno-col", "trait", "--covar", covar, "--keep", keep,
        "--taus", "0.25,0.5,0.75", "--blocks", "20", "--batch-size", "4",
        "--out", str(out),
    ])
    cli.main()

    assert seen, "no batch was processed"
    for g_shape, q_shape, deviation, rif_shape in seen:
        assert g_shape[0] == N_GENO
        assert q_shape == (N_GENO, X.shape[1])
        assert rif_shape[0] == N_GENO
        # The projector is orthonormal over exactly the rows it is applied to.
        assert deviation < 1e-12, deviation

    with gzip.open(f"{out}.stage1.tsv.gz", "rt") as fh:
        lines = fh.read().strip().split("\n")
    assert len(lines) == 1 + len(variants)
