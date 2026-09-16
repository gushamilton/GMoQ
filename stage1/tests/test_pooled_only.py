"""``--cov-mode pooled-only`` must agree with the full path everywhere it claims to.

The pooled-only path drops the per-block accumulator entirely, so it cannot be
bit-identical to the block path and does not claim to be: the totals are summed
over samples rather than over blocks, which is the same addends in a different
association.  What it does claim is pinned here.

1. The same variants are kept and the same variants are rejected -- the MAC and
   variance filters are the same two tests on inputs that are bit-for-bit equal
   up to the filters themselves.
2. ``beta_tau`` and ``d_resid`` agree with the block path to rounding, and agree
   *exactly* once written at the ``%.6g`` / ``%.9g`` precision of the TSV.
3. The kernel is independent of the panel size and of the thread count, bit for
   bit -- nothing is reassociated by parallelism.
4. End to end through the CLI: the ``.sigma0.tsv.gz`` sidecar is byte-identical
   to a full-mode run's, the shared TSV columns are identical, the jackknife
   payload is absent, and ``se_tau_*`` are the pooled standard errors
   ``sqrt(Sigma_0[t,t] / d_resid)`` rather than jackknife estimates.
"""

import gzip
import sys
from pathlib import Path

import numpy as np
import pytest

from fungwas_stage1 import cli, core
from fungwas_stage1.pooled_cov import compute_sigma0, read_sigma0

pytestmark = pytest.mark.skipif(not core.HAVE_CPP, reason="C++ extension not built")


def _problem(seed=11, N=2003, M=23, T=9, K=5, B=17):
    rng = np.random.default_rng(seed)
    G = np.asfortranarray(rng.normal(1.0, 0.5, (N, M)))
    G[rng.random((N, M)) < 0.02] = np.nan
    G[:, 2] = 0.0
    G[:3, 2] = 1.0                                  # fails the MAC filter
    G[:, 4] = np.nan
    G[:50, 4] = 1.0                                 # fails n_valid >= 100
    X = np.column_stack([np.ones(N), rng.normal(0, 1, (N, K - 1))])
    Q, _ = np.linalg.qr(X)
    rif = rng.normal(0, 1, (N, T))
    rif_resid = rif - Q @ (Q.T @ rif)
    block_ids = rng.integers(0, B, N).astype(np.int32)
    return G, rif_resid, Q, block_ids, B


def _block_reference(G, rif_resid, Q, block_ids, B):
    """beta and d_resid as the production block path derives them."""
    stats = core.compute_block_scores(G, rif_resid, Q, block_ids, B, n_threads=2)
    D = stats[:, :, 0].sum(axis=1)
    num = stats[:, :, 1:].sum(axis=1)
    with np.errstate(invalid="ignore"):
        beta = num / D[:, None]
    return beta, D


def test_pooled_keeps_and_rejects_exactly_the_same_variants():
    G, rif_resid, Q, block_ids, B = _problem()
    _, D_ref = _block_reference(G, rif_resid, Q, block_ids, B)
    _, d_pooled = core.compute_pooled_panel(G, rif_resid, Q, n_threads=2, panel_size=8)
    rejected_ref = np.isnan(D_ref)
    rejected_pooled = np.isnan(np.asarray(d_pooled))
    assert rejected_ref.sum() == 2          # the MAC and n_valid failures above
    assert np.array_equal(rejected_ref, rejected_pooled)


def test_pooled_agrees_with_the_block_path_to_rounding():
    G, rif_resid, Q, block_ids, B = _problem()
    beta_ref, D_ref = _block_reference(G, rif_resid, Q, block_ids, B)
    beta, d = core.compute_pooled_panel(G, rif_resid, Q, n_threads=2, panel_size=8)
    ok = np.isfinite(D_ref)

    # Same association would be bit-equality; different association is rounding.
    assert np.max(np.abs(d[ok] - D_ref[ok]) / D_ref[ok]) < 1e-12
    scale = np.maximum(np.abs(beta_ref[ok]), np.sqrt(np.diag(rif_resid.T @ rif_resid) / D_ref[ok, None]))
    assert np.max(np.abs(beta[ok] - beta_ref[ok]) / scale) < 1e-9


def test_pooled_matches_the_block_path_at_tsv_precision():
    """What actually reaches disk is identical, which is the claim that matters."""
    G, rif_resid, Q, block_ids, B = _problem()
    beta_ref, D_ref = _block_reference(G, rif_resid, Q, block_ids, B)
    beta, d = core.compute_pooled_panel(G, rif_resid, Q, n_threads=2, panel_size=8)
    ok = np.isfinite(D_ref)
    assert [f"{v:.9g}" for v in d[ok]] == [f"{v:.9g}" for v in D_ref[ok]]
    assert ([f"{v:.6g}" for v in beta[ok].ravel()]
            == [f"{v:.6g}" for v in beta_ref[ok].ravel()])


@pytest.mark.parametrize("panel_size", [1, 3, 8, 16, 64])
def test_pooled_is_panel_size_independent_bitwise(panel_size):
    G, rif_resid, Q, _, _ = _problem()
    b0, d0 = core.compute_pooled_panel(G, rif_resid, Q, n_threads=2, panel_size=8)
    b1, d1 = core.compute_pooled_panel(G, rif_resid, Q, n_threads=2, panel_size=panel_size)
    assert np.array_equal(b0, b1, equal_nan=True)
    assert np.array_equal(d0, d1, equal_nan=True)


@pytest.mark.parametrize("n_threads", [1, 2, 4])
def test_pooled_is_thread_count_independent_bitwise(n_threads):
    G, rif_resid, Q, _, _ = _problem()
    b0, d0 = core.compute_pooled_panel(G, rif_resid, Q, n_threads=1, panel_size=8)
    b1, d1 = core.compute_pooled_panel(G, rif_resid, Q, n_threads=n_threads, panel_size=8)
    assert np.array_equal(b0, b1, equal_nan=True)
    assert np.array_equal(d0, d1, equal_nan=True)


# ---------------------------------------------------------------------------
# End to end through the CLI, on the same real-BGEN fixture the input-boundary
# test uses.
# ---------------------------------------------------------------------------

def _cli_run(bgen, pheno, snps, out, cov_mode, monkeypatch, taus):
    monkeypatch.setattr(sys, "argv", [
        "fungwas-stage1",
        "--bgen", str(bgen),
        "--pheno", str(pheno),
        "--pheno-col", "trait",
        "--snps", str(snps),
        "--taus", taus,
        "--blocks", "16",
        "--batch-size", "1",
        "--cov-dtype", "float32",
        "--threads", "1",
        "--cov-mode", cov_mode,
        "--out", str(out),
    ])
    cli.main()


def _rows(prefix):
    with gzip.open(f"{prefix}.stage1.tsv.gz", "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return head, [line.rstrip("\n").split("\t") for line in fh]


def test_cli_pooled_only_against_full_mode(tmp_path, monkeypatch):
    fixture = pytest.importorskip("tests.test_bgen_fixture", reason="fixture helpers")
    dosages = fixture._known_dosages()
    bgen = tmp_path / "fixture.bgen"
    fixture._write_bgen(bgen, fixture.SAMPLES, dosages)
    pheno = tmp_path / "pheno.txt"
    fixture._write_phenotype(pheno, dosages)
    snps = tmp_path / "snps.txt"
    fixture._write_snp_list(snps)

    full = tmp_path / "full"
    pooled = tmp_path / "pooled_only"
    _cli_run(bgen, pheno, snps, full, "full", monkeypatch, fixture.TAUS)
    _cli_run(bgen, pheno, snps, pooled, "pooled-only", monkeypatch, fixture.TAUS)

    # Sigma_0 has no genotype in it and no block structure, so it must be
    # byte-identical once decompressed.
    with gzip.open(f"{full}.sigma0.tsv.gz", "rb") as a, \
            gzip.open(f"{pooled}.sigma0.tsv.gz", "rb") as b:
        assert a.read() == b.read()

    # The jackknife payload is not written at all, rather than written empty.
    for suffix in (".cov.gz", ".cov.ids.tsv.gz"):
        assert not Path(f"{pooled}{suffix}").exists()
        assert Path(f"{full}{suffix}").exists()

    head_f, rows_f = _rows(full)
    head_p, rows_p = _rows(pooled)
    assert head_f == head_p
    assert len(rows_f) == len(rows_p) > 0

    beta_cols = [i for i, c in enumerate(head_f) if c.startswith("beta_tau_")]
    se_cols = [i for i, c in enumerate(head_f) if c.startswith("se_tau_")]
    d_col = head_f.index("d_resid")
    id_cols = [i for i, c in enumerate(head_f)
               if not c.startswith(("beta_tau_", "se_tau_")) and c != "d_resid"]

    for rf, rp in zip(rows_f, rows_p):
        assert [rf[i] for i in id_cols] == [rp[i] for i in id_cols]
        assert [rf[i] for i in beta_cols] == [rp[i] for i in beta_cols]
        assert rf[d_col] == rp[d_col]

    # se_tau_* is a different estimator, and specifically this one.
    _, sigma0 = read_sigma0(str(pooled))
    S = next(iter(sigma0.values()))
    sd = np.sqrt(np.diag(S))
    for rp in rows_p:
        d = float(rp[d_col])
        expected = [f"{v:.6g}" for v in sd / np.sqrt(d)]
        assert [rp[i] for i in se_cols] == expected

    # And it is not the jackknife: the two must actually differ, or this mode
    # would be claiming a saving it did not make.
    assert any(rf[i] != rp[i] for rf, rp in zip(rows_f, rows_p) for i in se_cols)


def test_cli_pooled_only_records_the_se_source(tmp_path, monkeypatch):
    import json

    fixture = pytest.importorskip("tests.test_bgen_fixture", reason="fixture helpers")
    dosages = fixture._known_dosages()
    bgen = tmp_path / "fixture.bgen"
    fixture._write_bgen(bgen, fixture.SAMPLES, dosages)
    pheno = tmp_path / "pheno.txt"
    fixture._write_phenotype(pheno, dosages)
    snps = tmp_path / "snps.txt"
    fixture._write_snp_list(snps)

    full = tmp_path / "full"
    pooled = tmp_path / "pooled_only"
    _cli_run(bgen, pheno, snps, full, "full", monkeypatch, fixture.TAUS)
    _cli_run(bgen, pheno, snps, pooled, "pooled-only", monkeypatch, fixture.TAUS)

    with open(f"{pooled}.cov.meta.json") as fh:
        meta_p = json.load(fh)
    with open(f"{full}.cov.meta.json") as fh:
        meta_f = json.load(fh)
    assert meta_p["se_source"] == "pooled_sigma0"
    # Absent, not null, on the jackknife path: files written before the field
    # existed must stay readable as "jackknife" without a special case.
    assert "se_source" not in meta_f


def test_sigma0_is_a_by_product_of_the_rif_only():
    """The claim the mode rests on: Sigma_0 needs neither genotype nor blocks."""
    _, rif_resid, _, _, _ = _problem()
    S = compute_sigma0(rif_resid)
    assert np.allclose(S, S.T, rtol=0, atol=0)
    assert np.all(np.linalg.eigvalsh(S) > -1e-12)


def test_pooled_only_is_opt_in_and_the_default_is_unchanged():
    """Adding a third cov mode must not move the default off ``full``.

    The whole point of making this a separate mode rather than a change to
    ``pooled`` is that ``se_tau_*`` means something different under it.  If the
    default ever moved, every existing Stage 2 run reading ``se_tau_*`` would
    silently change its Q without anyone asking for it.
    """
    parser = cli.build_arg_parser() if hasattr(cli, "build_arg_parser") else None
    if parser is None:
        # The CLI builds its parser inline in main(); fall back to the constants,
        # which is what main()'s `default=` is written from.
        assert cli.COV_MODE_FULL == "full"
        assert cli.COV_MODES[0] == cli.COV_MODE_FULL
    else:
        assert parser.get_default("cov_mode") == cli.COV_MODE_FULL
    assert cli.COV_MODE_POOLED_ONLY == "pooled-only"
    assert cli.COV_MODE_POOLED_ONLY in cli.COV_MODES
    # ``pooled`` and ``pooled-only`` both drop the payload; only ``full`` keeps it.
    assert cli.COV_MODE_FULL not in cli.COV_MODES_NO_PAYLOAD
    assert cli.COV_MODE_POOLED in cli.COV_MODES_NO_PAYLOAD
    assert cli.COV_MODE_POOLED_ONLY in cli.COV_MODES_NO_PAYLOAD


def test_pooled_only_moves_se_tau_and_nothing_else(tmp_path, monkeypatch):
    """Measured on the full chr22 testosterone scan, pinned here on the fixture.

    On chr22 (17,550 variants, N=399,710, T=45, --blocks 1000) the full-mode and
    pooled-only TSVs agreed on every ``beta_tau_*`` and every ``d_resid`` value as
    written text -- 0 of 789,750 and 0 of 17,550 differing -- and on the whole
    ``.sigma0.tsv.gz``, while 789,704 of 789,750 ``se_tau_*`` values differed.
    Downstream that made Stage 2's gamma/beta_left/beta_right, their SEs, their
    p-values and their dfs byte-identical, and moved only Q, Q_pval, Q_df_eff and
    Q_scale.  This test pins the Stage 1 half of that claim per column class.
    """
    fixture = pytest.importorskip("tests.test_bgen_fixture", reason="fixture helpers")
    dosages = fixture._known_dosages()
    bgen = tmp_path / "fixture.bgen"
    fixture._write_bgen(bgen, fixture.SAMPLES, dosages)
    pheno = tmp_path / "pheno.txt"
    fixture._write_phenotype(pheno, dosages)
    snps = tmp_path / "snps.txt"
    fixture._write_snp_list(snps)

    full = tmp_path / "full"
    pooled = tmp_path / "pooled_only"
    _cli_run(bgen, pheno, snps, full, "full", monkeypatch, fixture.TAUS)
    _cli_run(bgen, pheno, snps, pooled, "pooled-only", monkeypatch, fixture.TAUS)

    head_f, rows_f = _rows(full)
    head_p, rows_p = _rows(pooled)
    assert head_f == head_p

    beta_cols = [i for i, c in enumerate(head_f) if c.startswith("beta_tau_")]
    se_cols = [i for i, c in enumerate(head_f) if c.startswith("se_tau_")]
    d_col = head_f.index("d_resid")
    id_cols = [i for i, c in enumerate(head_f)
               if not c.startswith(("beta_tau_", "se_tau_")) and c != "d_resid"]

    n_beta_diff = sum(rf[i] != rp[i] for rf, rp in zip(rows_f, rows_p) for i in beta_cols)
    n_id_diff = sum(rf[i] != rp[i] for rf, rp in zip(rows_f, rows_p) for i in id_cols)
    n_d_diff = sum(rf[d_col] != rp[d_col] for rf, rp in zip(rows_f, rows_p))
    n_se_diff = sum(rf[i] != rp[i] for rf, rp in zip(rows_f, rows_p) for i in se_cols)

    assert n_id_diff == 0, "identity columns must not move"
    assert n_beta_diff == 0, "beta_tau_* must agree at TSV precision"
    assert n_d_diff == 0, "d_resid must agree at TSV precision"
    assert n_se_diff > 0, "se_tau_* is the one column class that is a different estimator"

    # There is no frequency column in the Stage 1 schema, so the correctness
    # comparison has no allele-frequency class to check -- the allele columns are
    # the only allele-level output, and they are covered by id_cols above.
    assert not [c for c in head_f
                if any(k in c.lower() for k in ("freq", "eaf", "maf"))]
