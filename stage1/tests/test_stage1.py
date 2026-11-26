"""Core checks for G-MoQ Stage 1.

Every test builds its own small simulated BGEN, so the suite needs no data files.
Expected values are computed directly from the documented formulas, not by
calling the code under test:

    RIF_i(tau) = q_tau + (tau - 1{y_i <= q_tau}) / f(q_tau)

with q_tau the type-8 (median-unbiased) sample quantile and f a Gaussian kernel
density with Silverman bandwidth 1.06 * sd(y) * N^(-1/5), and the Stage 1 slope
at each tau is the OLS coefficient of the genotype in RIF(tau) ~ G + covariates.
"""

import csv
import gzip
import logging
import sqlite3
import sys

import numpy as np
import pytest

from gmoq import cli, core, pooled_cov

N = 600
TAUS = np.array([0.1, 0.25, 0.5, 0.75, 0.9])
TAU_ARG = ",".join(f"{t:g}" for t in TAUS)
RSIDS = [f"rs{j}" for j in range(1, 7)]


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

def _probs(dosage_allele1):
    """[P(11), P(12), P(22)] for an integer allele-1 dosage 0, 1 or 2."""
    p = np.zeros((len(dosage_allele1), 3))
    p[np.arange(len(dosage_allele1)), (2 - dosage_allele1).astype(int)] = 1.0
    return p


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    """A BGEN with six SNPs plus one rsID split over two records, and inputs."""
    from bgen import BgenWriter

    d = tmp_path_factory.mktemp("stage1")
    rng = np.random.default_rng(2026)
    samples = [f"S{i:04d}" for i in range(N)]
    G = np.column_stack([rng.binomial(2, p, N) for p in (0.1, 0.2, 0.3, 0.4, 0.45, 0.5)])
    C = rng.normal(size=(N, 2))
    loco = 0.5 * C[:, 0] + rng.normal(scale=0.3, size=N)
    # a mean effect at rs1 and a variance effect at rs2
    y = (1.0 + 0.4 * G[:, 0] + 0.5 * C[:, 0] - 0.3 * C[:, 1] + 0.4 * loco
         + rng.standard_normal(N) * np.sqrt(1.0 + 0.6 * G[:, 1]))

    bgen = d / "geno.bgen"
    with BgenWriter(str(bgen), n_samples=N, samples=samples, compression="zlib") as w:
        for j, rsid in enumerate(RSIDS):
            w.add_variant(varid=f"v{j}", rsid=rsid, chrom="22", pos=1000 * (j + 1),
                          alleles=["A", "G"], genotypes=_probs(G[:, j]), bit_depth=8)
        # one rsID on two records (a split multiallelic site): excluded by default
        for k, alt in enumerate(("C", "T")):
            w.add_variant(varid=f"m{k}", rsid="rs_multi", chrom="22", pos=9000,
                          alleles=["A", alt], genotypes=_probs(rng.binomial(2, 0.3, N)),
                          bit_depth=8)

    pheno = d / "pheno.txt"
    with pheno.open("w") as fh:
        fh.write("FID IID y\n")
        for s, v in zip(samples, y):
            fh.write(f"{s} {s} {v:.17g}\n")
    covar = d / "covar.txt"
    with covar.open("w") as fh:
        fh.write("FID IID c1 c2\n")
        for s, (a, b) in zip(samples, C):
            fh.write(f"{s} {s} {a:.17g} {b:.17g}\n")
    loco_file = d / "y.loco"
    with loco_file.open("w") as fh:
        fh.write("FID_IID " + " ".join(f"{s}_{s}" for s in samples) + "\n")
        fh.write("22 " + " ".join(f"{v:.17g}" for v in loco) + "\n")
    return dict(dir=d, bgen=bgen, pheno=pheno, covar=covar, loco=loco_file,
                G=G.astype(float), C=C, y=y, L=loco)


def _run(monkeypatch, data, out, *extra):
    argv = ["gmoq-stage1", "--bgen", str(data["bgen"]), "--pheno", str(data["pheno"]),
            "--pheno-col", "y", "--covar", str(data["covar"]), "--taus", TAU_ARG,
            "--blocks", "20", "--threads", "1", "--seed", "1", "--out", str(out), *extra]
    monkeypatch.setattr(sys, "argv", argv)
    cli.main()
    with gzip.open(f"{out}.stage1.tsv.gz", "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _cols(rows, prefix):
    return np.array([[float(r[f"{prefix}_{t:.3f}"]) for t in TAUS] for r in rows])


def _effects_from_z(rows, out):
    """beta and se from a pooled-only file: se = sqrt(diag(Sigma_0) / d_resid), beta = z * se."""
    _, mats = pooled_cov.read_sigma0(str(out))
    sd0 = np.sqrt(np.diag(next(iter(mats.values()))))
    d = np.array([float(r["d_resid"]) for r in rows])
    se = sd0[None, :] / np.sqrt(d)[:, None]
    return _cols(rows, "z_tau") * se, se


def _assert_beta_at_output_precision(rows, beta):
    """beta_tau_* is written to 0.001 of its standard error."""
    got, se = _cols(rows, "beta_tau"), _cols(rows, "se_tau")
    assert np.all(np.abs(got - beta) <= 1e-3 * se), np.max(np.abs(got - beta) / se)


# ---------------------------------------------------------------------------
# direct computation
# ---------------------------------------------------------------------------

def direct_rif(y, taus):
    q = np.quantile(y, taus, method="median_unbiased")
    h = 1.06 * np.std(y) * len(y) ** -0.2
    f = np.array([np.mean(np.exp(-0.5 * ((y - qt) / h) ** 2)) / (h * np.sqrt(2 * np.pi))
                  for qt in q])
    return q + (taus - (y[:, None] <= q).astype(float)) / f


def direct_fit(y, G, X, taus):
    """Per-tau OLS slope of each SNP, residualised RIF, and D = sum g_resid^2."""
    rif = direct_rif(y, taus)
    X = np.column_stack([np.ones(len(y)), X])
    resid = lambda A: A - X @ np.linalg.lstsq(X, A, rcond=None)[0]
    psi = resid(rif)
    beta = np.array([np.linalg.lstsq(np.column_stack([G[:, j], X]), rif, rcond=None)[0][0]
                     for j in range(G.shape[1])])
    Gr = resid(G)
    return beta, psi, (Gr ** 2).sum(0), Gr


# ---------------------------------------------------------------------------
# estimates and covariance
# ---------------------------------------------------------------------------

def test_reproduces_direct_rif_ols_at_every_quantile(data, monkeypatch, tmp_path):
    rows = _run(monkeypatch, data, tmp_path / "full")
    assert [r["snp_id"] for r in rows] == RSIDS            # rs_multi excluded
    assert {(r["effect_allele"], r["other_allele"]) for r in rows} == {("A", "G")}

    beta, psi, D, Gr = direct_fit(data["y"], data["G"], data["C"], TAUS)
    _assert_beta_at_output_precision(rows, beta)
    np.testing.assert_allclose([float(r["d_resid"]) for r in rows], D, rtol=1e-8)
    assert np.all(beta[0] > 0.2)                          # the mean effect, every tau

    # jackknife SEs are a noisy estimate of the HC0 standard error
    hc0 = np.array([np.sqrt(np.sum((psi - np.outer(g, b)) ** 2 * g[:, None] ** 2, 0)) / Dj
                    for g, b, Dj in zip(Gr.T, beta, D)])
    ratio = _cols(rows, "se_tau") / hc0
    assert 0.7 < np.median(ratio) < 1.3


def test_pooled_covariance_matches_direct_computation(data, monkeypatch, tmp_path):
    out = tmp_path / "pooled"
    rows = _run(monkeypatch, data, out, "--cov-mode", "pooled-only")
    beta, psi, D, _ = direct_fit(data["y"], data["G"], data["C"], TAUS)

    meta, mats = pooled_cov.read_sigma0(str(out))
    sigma0 = mats["ALL"]
    np.testing.assert_allclose(sigma0, psi.T @ psi / N, rtol=1e-9)
    assert meta["n_samples"] == N
    np.testing.assert_allclose(meta["taus"], TAUS)

    # pooled-only writes z_tau; beta and se_tau = sqrt(diag(Sigma0) / D) are rebuilt
    got_beta, se = _effects_from_z(rows, out)
    np.testing.assert_allclose(se, np.sqrt(np.diag(sigma0)[None, :] / D[:, None]), rtol=1e-8)
    assert np.all(np.abs(got_beta - beta) <= 1e-3 * se)


def test_pooled_covariance_reproduces_the_hc0_sandwich():
    """Sigma0 / D equals the exact per-individual HC0 covariance in an i.i.d. sample."""
    rng = np.random.default_rng(20260914)
    n = 20000
    taus = np.linspace(0.1, 0.9, 9)
    X = np.column_stack([np.ones(n), rng.normal(size=(n, 3))])
    y = X @ np.array([1.0, 0.4, -0.2, 0.1]) + rng.gamma(3.0, 1.5, n)
    Q, _ = np.linalg.qr(X)
    psi = direct_rif(y, taus)
    psi -= Q @ (Q.T @ psi)
    sigma0 = pooled_cov.compute_sigma0(psi)
    ratios = []
    for _ in range(10):
        g = rng.binomial(2, rng.uniform(0.05, 0.5), n).astype(float)
        g -= Q @ (Q.T @ g)
        D = g @ g
        e = psi - np.outer(g, (g @ psi) / D)
        hc0 = (e * (g * g)[:, None]).T @ e / D ** 2
        ratios.append(np.diag(pooled_cov.sigma_for_variant(sigma0, D)) / np.diag(hc0))
    ratios = np.concatenate(ratios)
    assert abs(np.median(ratios) - 1) < 0.01
    assert np.max(np.abs(ratios - 1)) < 0.10


def test_loco_prediction_enters_as_a_covariate(data, monkeypatch, tmp_path):
    rows = _run(monkeypatch, data, tmp_path / "loco", "--loco-preds", str(data["loco"]),
                "--loco-mode", "covariate", "--cov-mode", "pooled-only")
    X = np.column_stack([data["C"], data["L"]])
    beta, psi, D, _ = direct_fit(data["y"], data["G"], X, TAUS)
    np.testing.assert_allclose([float(r["d_resid"]) for r in rows], D, rtol=1e-8)
    got_beta, se = _effects_from_z(rows, tmp_path / "loco")
    assert np.all(np.abs(got_beta - beta) <= 1e-3 * se)
    _, mats = pooled_cov.read_sigma0(str(tmp_path / "loco"))
    np.testing.assert_allclose(mats["22"], psi.T @ psi / N, rtol=1e-9)


# ---------------------------------------------------------------------------
# BGEN reading
# ---------------------------------------------------------------------------

def _write_bgi(bgen):
    """A bgenix-style .bgi index for a BGEN, from the record offsets."""
    from bgen import BgenReader

    size = bgen.stat().st_size
    with BgenReader(str(bgen), delay_parsing=True) as r:
        recs = [(v.chrom, v.pos, v.rsid, v.alleles, v.fileoffset) for v in r]
    con = sqlite3.connect(f"{bgen}.bgi")
    con.execute("CREATE TABLE Variant (chromosome TEXT NOT NULL, position INT NOT NULL, "
                "rsid TEXT NOT NULL, number_of_alleles INT NOT NULL, allele1 TEXT NOT NULL, "
                "allele2 TEXT NULL, file_start_position INT NOT NULL, "
                "size_in_bytes INT NOT NULL, PRIMARY KEY (chromosome, position, rsid, "
                "allele1, allele2, file_start_position))")
    ends = [rec[4] for rec in recs[1:]] + [size]
    con.executemany("INSERT INTO Variant VALUES (?,?,?,?,?,?,?,?)",
                    [(c, p, s, len(a), a[0], a[1], off, end - off)
                     for (c, p, s, a, off), end in zip(recs, ends)])
    con.commit()
    con.close()


def test_bgen_with_and_without_index_give_identical_output(data, monkeypatch, tmp_path):
    import shutil

    bgen = tmp_path / "copy.bgen"
    shutil.copy(data["bgen"], bgen)
    snps = tmp_path / "snps.txt"
    snps.write_text("rs5\nrs2\nrs_multi\nrs_absent\n")    # unordered, with a missing ID
    local = dict(data, bgen=bgen)

    no_index = _run(monkeypatch, local, tmp_path / "noidx", "--snps", str(snps))
    _write_bgi(bgen)
    indexed = _run(monkeypatch, local, tmp_path / "idx", "--snps", str(snps))
    whole = _run(monkeypatch, local, tmp_path / "whole")

    assert [r["snp_id"] for r in indexed] == ["rs2", "rs5"]   # file order, rs_multi dropped
    assert indexed == no_index
    assert indexed == [r for r in whole if r["snp_id"] in ("rs2", "rs5")]


def test_dosage_is_the_expected_count_of_allele_1(data):
    variants, samples, handle = cli.get_variants_bgen(str(data["bgen"]), None, ["rs3"])
    (v,) = list(variants)
    assert handle is not None and samples[0] == "S0000"
    np.testing.assert_allclose(cli.extract_dosage(v), data["G"][:, 2], atol=1e-6)


# ---------------------------------------------------------------------------
# input validation
# ---------------------------------------------------------------------------

def _write(path, header, rows):
    path.write_text("\n".join([header] + rows) + "\n")
    return str(path)


def test_dropped_rows_and_missing_codes_are_reported(tmp_path, caplog):
    pheno = _write(tmp_path / "p.txt", "FID IID y",
                   ["a a 1.0", "b b NA", "c c .", "d d -9", "e e 2.0", "f f 3.0"])
    covar = _write(tmp_path / "c.txt", "FID IID age",
                   ["a a 40", "b b 41", "c c 42", "d d 43", "e e NA"])  # f has no row
    with caplog.at_level(logging.WARNING):
        ids, *_ = cli.load_analysis_inputs(pheno, ["y"], covar)
    assert ids == ["a", "d"]
    for msg in ("Dropped 4 phenotype rows", "1 with a missing (NA) phenotype",
                "1 with a non-numeric or empty phenotype", "1 with a missing (NA) covariate",
                "1 with no row in the covariate file", "1 samples have the value -9 in column 'y'"):
        assert msg in caplog.text


def test_ambiguous_sample_ids_are_refused():
    with pytest.raises(ValueError, match="BGEN sample IDs are not unique"):
        cli.match_samples(["S1", "S1"], ["S1"])
    with pytest.raises(ValueError, match="Phenotype sample IDs are not unique"):
        cli.match_samples(["S1", "S2"], ["S1", "S1"])


def test_invalid_rif_inputs_are_refused():
    with pytest.raises(ValueError):
        core.compute_rif(np.array([1.0, np.nan, 2.0]), TAUS)
    with pytest.raises(ValueError):
        core.compute_rif(np.arange(10.0), np.array([0.0, 0.5]))
    with pytest.raises(ValueError):
        core.resolve_n_blocks(5, n_taus=len(TAUS) + 10)


# ---------------------------------------------------------------------------
# fast paths: same answers as the reference paths
# ---------------------------------------------------------------------------

def test_pooled_only_cli_fast_matches_exact(data, monkeypatch, tmp_path):
    loco_mode = "covariate"
    loco = [] if loco_mode is None else ["--loco-preds", str(data["loco"]), "--loco-mode", loco_mode]
    rows_fast = _run(monkeypatch, data, tmp_path / "p_fast", "--cov-mode", "pooled-only", *loco)
    rows_exact = _run(monkeypatch, data, tmp_path / "p_exact", "--cov-mode", "pooled-only",
                      "--exact-kernel", *loco)
    assert [r["snp_id"] for r in rows_fast] == [r["snp_id"] for r in rows_exact]
    assert np.max(np.abs(_cols(rows_fast, "z_tau") - _cols(rows_exact, "z_tau"))) <= 1e-3
    np.testing.assert_allclose([float(r["d_resid"]) for r in rows_fast],
                               [float(r["d_resid"]) for r in rows_exact], rtol=1e-8)


def test_block_kernel_bins_match_dense():
    case = "tied"
    rng = np.random.default_rng({"cov": 1, "int_only": 2, "tied": 3, "unsorted": 4, "no_int": 5}[case])
    n, m, B = 2500, 13, 40
    taus = np.linspace(0.05, 0.95, 15)
    if case == "unsorted":
        taus = taus[rng.permutation(len(taus))]
    y = rng.normal(size=n)
    if case == "tied":
        y = np.round(2 * y)
    if case == "int_only":
        X = np.ones((n, 1))
    elif case == "no_int":
        X = rng.normal(size=(n, 3))
    else:
        X = np.column_stack([np.ones(n), rng.normal(size=(n, 4))])
    Q, _ = np.linalg.qr(X)
    G = rng.binomial(2, 0.3, size=(n, m)).astype(float)
    G[rng.random(G.shape) < 0.02] = np.nan
    G[:, 0] = 0.0
    block_ids = rng.integers(0, B, size=n).astype(np.int32)
    rif, q_tau, fhat = core.compute_rif(y, taus)
    rif_resid = rif - Q @ (Q.T @ rif)

    dense = core.make_block_plan(rif_resid, Q, block_ids, B, n_threads=2)
    binned = core.make_block_plan(rif_resid, Q, block_ids, B, n_threads=2)
    binned.set_bins(y, taus, q_tau, fhat)
    assert binned.has_bins and not dense.has_bins
    ref = core.compute_block_scores_panel_fast(G, rif_resid, Q, block_ids, B, n_threads=2, plan=dense)
    got = core.compute_block_scores_panel_fast(G, rif_resid, Q, block_ids, B, n_threads=2, plan=binned)
    ref, got = np.asarray(ref), np.asarray(got)
    assert np.array_equal(np.isnan(ref), np.isnan(got))
    ok = ~np.isnan(ref)
    scale = np.nanmax(np.abs(ref), axis=(0, 1), keepdims=True)
    np.testing.assert_allclose((got / scale)[ok], (ref / scale)[ok], rtol=0, atol=1e-11)


def _write_bgen(path, probs_list, bit_depths, compression):
    from bgen import BgenWriter
    n = probs_list[0].shape[0]
    samples = [f"S{i:04d}" for i in range(n)]
    with BgenWriter(str(path), n_samples=n, samples=samples, compression=compression) as w:
        for j, (p, b) in enumerate(zip(probs_list, bit_depths)):
            w.add_variant(varid=f"v{j}", rsid=f"rs{j}", chrom="7", pos=100 * (j + 1),
                          alleles=["A", "T"], genotypes=p, bit_depth=b)
    return samples


def _indexed_copy(data, tmp_path):
    import shutil
    bgen = tmp_path / "indexed.bgen"
    shutil.copy(data["bgen"], bgen)
    _write_bgi(bgen)
    return dict(data, bgen=bgen)


def _snps_file(tmp_path):
    p = tmp_path / "snps.txt"
    p.write_text("\n".join(RSIDS + ["rs_multi"]) + "\n")
    return p


def _write_keep(path, ids):
    path.write_text("\n".join(f"{s} {s}" for s in ids) + "\n")
    return str(path)


def test_native_decoder_matches_bgen_package(tmp_path):
    if not hasattr(core.cpp_ext, "bgen_decode_block"):
        pytest.skip("extension built without the native BGEN decoder")
    rng = np.random.default_rng(7)
    n = 257
    bit_depths = [1, 3, 8, 8, 12, 16, 23, 32]
    probs_list = []
    for b in bit_depths:
        p = rng.dirichlet([0.6, 0.8, 1.0], size=n)
        p[rng.random(n) < 0.05] = np.nan           # missing samples
        probs_list.append(p)
    path = tmp_path / "g.bgen"
    _write_bgen(path, probs_list, bit_depths, "zlib")      # BgenWriter writes the .bgi too

    variants, _, handle = cli.get_variants_bgen(str(path), None, [f"rs{j}" for j in range(len(bit_depths))])
    assert isinstance(variants, cli.SeekingVariantGenerator)
    assert cli.native_bgen_decoder_reason(str(path), variants) is None

    idx = rng.permutation(n)[: n - 20]               # subset, out of order
    ref = {}
    for v in variants:
        bits = cli._variant_bit_depth(v, bgen_path=str(path))
        ref[v.rsid] = (bits, cli.extract_dosage(v, probability_bit_depth=bits)[idx])

    skipped = cli._SkipCount()
    # Yielded columns are views into a reused block buffer: copy before the next item.
    got = [(i, c, d.copy()) for i, c, d in
           cli._decode_variants_native(variants.records, str(path), idx, skipped,
                                       n_threads=2, block_size=3)]
    assert [g[0][0] for g in got] == sorted(ref, key=lambda r: int(r[2:]))
    for info, chrom, dos in got:
        bits, want = ref[info[0]]
        assert info[1:] == ("7", 100 * (int(info[0][2:]) + 1), "A", "T") and chrom == "7"
        assert np.array_equal(np.isnan(dos), np.isnan(want))
        ok = ~np.isnan(want)
        # bgen returns float32; native is the exact float64 of the stored integers
        np.testing.assert_allclose(dos[ok], want[ok], rtol=0, atol=1e-6 if bits <= 23 else 2e-6)
        assert np.all((dos[ok] >= 0) & (dos[ok] <= 2))


def test_cli_native_and_bgen_decoders_agree(data, monkeypatch, tmp_path):
    cov_mode, extra = "full", ["--batch-size", "4", "LOCO"]
    data = _indexed_copy(data, tmp_path)
    extra = [x for x in extra if x != "LOCO"] + (
        ["--loco-preds", str(data["loco"]), "--loco-mode", "covariate"] if "LOCO" in extra else [])
    rows_nat = _run(monkeypatch, data, tmp_path / "nat", "--cov-mode", cov_mode, *extra,
                    "--snps", str(_snps_file(tmp_path)), "--bgen-decoder", "native")
    rows_bgen = _run(monkeypatch, data, tmp_path / "pkg", "--cov-mode", cov_mode, *extra,
                     "--snps", str(_snps_file(tmp_path)), "--bgen-decoder", "bgen")
    if cov_mode == "full":
        import gzip as _gz
        for suffix in (".cov.gz", ".cov.ids.tsv.gz"):
            with _gz.open(f"{tmp_path / 'nat'}{suffix}", "rb") as a, _gz.open(f"{tmp_path / 'pkg'}{suffix}", "rb") as b:
                assert a.read() == b.read()
    assert rows_nat == rows_bgen        # hard-call genotypes: identical text output


def test_jobs_outputs_identical_to_separate_runs(data, monkeypatch, tmp_path):
    threads = "2"
    local = _indexed_copy(data, tmp_path)
    rng = np.random.default_rng(5)
    samples = [f"S{i:04d}" for i in range(N)]
    keep_a = _write_keep(tmp_path / "keep_a.txt", [s for s in samples if rng.random() < 0.8])
    keep_b = _write_keep(tmp_path / "keep_b.txt", [s for s in samples if rng.random() < 0.7])
    pheno2 = tmp_path / "pheno2.txt"
    with pheno2.open("w") as fh:
        fh.write("FID IID y z\n")
        z = rng.gamma(2.0, size=N)
        for s, v, w in zip(samples, data["y"], z):
            fh.write(f"{s} {s} {v:.17g} {'NA' if rng.random() < 0.1 else f'{w:.17g}'}\n")
    snps = _snps_file(tmp_path)
    shared = ["--bgen", str(local["bgen"]), "--snps", str(snps), "--taus", TAU_ARG,
              "--blocks", "20", "--seed", "1", "--batch-size", "2", "--threads", threads]
    jobs = [
        f"--pheno {pheno2} --pheno-col y --covar {local['covar']} --keep {keep_a} "
        f"--loco-preds {local['loco']} --loco-mode covariate --cov-mode pooled --out {tmp_path / 'j1'}",
        f"--pheno {pheno2} --pheno-col z --covar {local['covar']} --covar-col c1 "
        f"--keep {keep_b} --cov-mode full --out {tmp_path / 'j2'}",
        f"--pheno {pheno2} --pheno-col z --cov-mode pooled-only --out {tmp_path / 'j3'}",
    ]
    jf = tmp_path / "jobs.txt"
    jf.write_text("# three analyses\n" + "\n".join(jobs) + "\n")
    monkeypatch.setattr(sys, "argv", ["gmoq-stage1", *shared, "--jobs", str(jf)])
    cli.main()
    import shlex, gzip as _gz, glob, os
    for j, line in enumerate(jobs, 1):
        out_sep = tmp_path / f"s{j}"
        argv = shared + shlex.split(line.replace(f"--out {tmp_path / f'j{j}'}", f"--out {out_sep}"))
        monkeypatch.setattr(sys, "argv", ["gmoq-stage1", *argv])
        cli.main()
        produced = sorted(glob.glob(f"{tmp_path / f'j{j}'}.*"))
        assert produced
        for path in produced:
            suffix = path[len(str(tmp_path / f"j{j}")):]
            other = f"{out_sep}{suffix}"
            assert os.path.exists(other), suffix
            if path.endswith(".gz"):
                with _gz.open(path, "rb") as a, _gz.open(other, "rb") as b:
                    assert a.read() == b.read(), suffix
            elif path.endswith(".json"):
                import json
                da, db = json.load(open(path)), json.load(open(other))
                for d in (da, db):
                    for key in list(d):
                        if any(t in key for t in ("time", "date", "path", "out", "command", "host", "seconds", "file")):
                            d.pop(key)
                assert da == db, suffix
