"""Independent, real-BGEN checks for the Stage 1 input boundary.

The expected dosages in this file are derived from the genotype probabilities supplied
to ``BgenWriter``.  They are deliberately not obtained by reading the fixture back with
the same package that the CLI uses.  That leaves a reader allele or missingness bug
visible to the test.
"""

import csv
import gzip
import sys
from pathlib import Path

import numpy as np
import pytest

from fungwas_stage1 import cli


N_GENO = 160
SAMPLES = np.asarray([f"S{i:03d}" for i in range(N_GENO)])
TAUS = "0.2,0.5,0.8"
ALLELES = {
    "rs_frac": ("A", "C"),
    "rs_missing": ("G", "T"),
}


def _probabilities_for_first_allele(dosage):
    """Return [AA, Aa, aa] probabilities with the requested allele-0 dosage."""
    alt = 2.0 - np.asarray(dosage, dtype=float)
    probs = np.zeros((len(alt), 3), dtype=float)
    lower = alt <= 1.0
    probs[lower, 0] = 1.0 - alt[lower]
    probs[lower, 1] = alt[lower]
    probs[~lower, 1] = 2.0 - alt[~lower]
    probs[~lower, 2] = alt[~lower] - 1.0
    return probs


def _known_dosages():
    """Known allele-0 dosages in phenotype/sample order, before BGEN encoding."""
    # All probabilities below are exact multiples of 1/255, so the standard
    # 8-bit BGEN representation retains the stated fractional truth exactly.
    cycle = np.asarray([0.2, 0.8, 1.2, 1.8], dtype=float)
    fractional = np.resize(cycle, N_GENO)
    missing = np.resize(np.asarray([0.0, 0.4, 1.2, 1.8]), N_GENO)
    missing[17] = np.nan
    return {"rs_frac": fractional, "rs_missing": missing}


def _write_bgen(path: Path, sample_order, dosages):
    from bgen import BgenWriter

    order = [int(s[1:]) for s in sample_order]
    with BgenWriter(str(path), n_samples=N_GENO, samples=list(sample_order), compression="zlib") as writer:
        for j, (rsid, (effect, other)) in enumerate(ALLELES.items()):
            encoded = np.asarray(dosages[rsid])[order]
            probs = _probabilities_for_first_allele(np.nan_to_num(encoded, nan=0.0))
            missing = ~np.isfinite(encoded)
            # BgenWriter recognises an all-NaN row and sets the BGEN missing flag.
            # Zeroes are valid probabilities for a real genotype and would decode
            # as a called zero dosage.
            probs[missing, :] = np.nan
            writer.add_variant(
                varid=f"v{j + 1}",
                rsid=rsid,
                chrom="22",
                pos=1000 + j,
                alleles=[effect, other],
                genotypes=probs,
                bit_depth=8,
            )


def _write_phenotype(path: Path, dosages):
    # The phenotype is tied to the known allele-0 dosage.  The extra row exercises
    # the supported case where the phenotype file covers a wider cohort than BGEN.
    rng = np.random.default_rng(20260915)
    y = 2.0 + 1.4 * dosages["rs_frac"] + rng.normal(0.0, 0.12, N_GENO)
    with path.open("w") as fh:
        fh.write("FID IID trait\n")
        for sid, value in zip(SAMPLES, y):
            fh.write(f"{sid} {sid} {value:.12g}\n")
        fh.write("S_absent S_absent 3.5\n")
    return y


def _write_snp_list(path: Path):
    # The index path returns records in BGEN file order even when the request is
    # listed in the reverse order; checking that is part of this boundary test.
    path.write_text("rs_missing\nrs_frac\n")


def _run_cli(bgen, pheno, snps, out, monkeypatch):
    monkeypatch.setattr(sys, "argv", [
        "fungwas-stage1",
        "--bgen", str(bgen),
        "--pheno", str(pheno),
        "--pheno-col", "trait",
        "--snps", str(snps),
        "--taus", TAUS,
        "--blocks", "16",
        "--batch-size", "1",
        "--cov-dtype", "float32",
        "--threads", "1",
        "--out", str(out),
    ])
    cli.main()


def _read_stage1(path):
    with gzip.open(f"{path}.stage1.tsv.gz", "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_real_bgen_cli_matches_independent_dosage_truth_and_sample_order(tmp_path, monkeypatch):
    """The real CLI preserves alleles, fractional values, missingness and IDs."""
    dosages = _known_dosages()
    permutation = np.roll(np.arange(N_GENO), 17)
    permuted_samples = SAMPLES[permutation]

    bgen_natural = tmp_path / "natural.bgen"
    bgen_permuted = tmp_path / "permuted.bgen"
    _write_bgen(bgen_natural, SAMPLES, dosages)
    _write_bgen(bgen_permuted, permuted_samples, dosages)
    pheno = tmp_path / "pheno.txt"
    _write_phenotype(pheno, dosages)
    snps = tmp_path / "snps.txt"
    _write_snp_list(snps)

    # This is an independent truth check for the reader's dosage and missingness
    # convention.  The expected arrays come from the probabilities above, not from
    # BgenReader or extract_dosage.
    variants, samples, handle = cli.get_variants_bgen(
        str(bgen_permuted), None, ["rs_frac", "rs_missing"], seek_variants=True
    )
    assert samples == list(permuted_samples)
    observed = {variant.rsid: cli.extract_dosage(variant) for variant in variants}
    expected = {
        rsid: np.asarray(dosages[rsid])[permutation]
        for rsid in dosages
    }
    for rsid in expected:
        np.testing.assert_allclose(observed[rsid], expected[rsid], atol=2e-6, equal_nan=True)
    assert not np.isnan(observed["rs_missing"][17])
    assert np.isnan(observed["rs_missing"][int(np.flatnonzero(permutation == 17)[0])])
    assert handle is not None

    # Decode the same bytes with the independent bgen-reader implementation.  This
    # confirms that the writer fixture really contains the probabilities specified
    # above; a test that only reads with bgen would make writer and reader errors
    # cancel each other.
    from bgen_reader import read_bgen

    decoded = read_bgen(str(bgen_permuted), verbose=False)
    for j, rsid in enumerate(("rs_frac", "rs_missing")):
        probabilities = decoded["genotype"][j].compute()["probs"]
        encoded = np.asarray(dosages[rsid])[permutation]
        expected_probabilities = _probabilities_for_first_allele(
            np.nan_to_num(encoded, nan=0.0)
        )
        expected_probabilities[~np.isfinite(encoded), :] = np.nan
        np.testing.assert_allclose(
            probabilities, expected_probabilities, atol=2e-6, equal_nan=True
        )

    out_natural = tmp_path / "natural"
    out_permuted = tmp_path / "permuted"
    _run_cli(bgen_natural, pheno, snps, out_natural, monkeypatch)
    _run_cli(bgen_permuted, pheno, snps, out_permuted, monkeypatch)
    rows_natural = _read_stage1(out_natural)
    rows_permuted = _read_stage1(out_permuted)

    assert [row["snp_id"] for row in rows_natural] == ["rs_frac", "rs_missing"]
    assert [row["snp_id"] for row in rows_permuted] == ["rs_frac", "rs_missing"]
    assert [(row["effect_allele"], row["other_allele"]) for row in rows_natural] == [
        ("A", "C"), ("G", "T")
    ]
    assert rows_permuted == rows_natural

    # The phenotype was generated from allele-0 dosage.  A flipped allele would
    # reverse this sign while leaving many downstream diagnostics plausible.
    assert float(rows_natural[0]["beta_tau_0.500"]) > 0.0

    # Missing genotypes are mean-imputed by the Stage 1 kernel.  With an intercept
    # only, d_resid is the centered sum of squares after that imputation.
    d = np.nan_to_num(dosages["rs_missing"], nan=np.nanmean(dosages["rs_missing"]))
    expected_d_resid = float(np.sum((d - np.mean(d)) ** 2))
    assert float(rows_natural[1]["d_resid"]) == pytest.approx(expected_d_resid, rel=2e-6)


def test_match_samples_rejects_ambiguous_ids():
    with pytest.raises(ValueError, match="BGEN sample IDs are not unique"):
        cli.match_samples(["S1", "S1"], ["S1"])
    with pytest.raises(ValueError, match="Phenotype sample IDs are not unique"):
        cli.match_samples(["S1", "S2"], ["S1", "S1"])
