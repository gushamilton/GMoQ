"""Dropped input rows and likely missing codes are reported, not silent."""

import logging

from fungwas_stage1.cli import load_analysis_inputs


def _write(path, header, rows):
    path.write_text("\n".join([header] + rows) + "\n")
    return str(path)


def test_dropped_rows_and_minus9_are_warned(tmp_path, caplog):
    pheno = _write(tmp_path / "pheno.txt", "FID IID y",
                   ["a a 1.0", "b b NA", "c c .", "d d -9", "e e 2.0", "f f 3.0"])
    covar = _write(tmp_path / "covar.txt", "FID IID age",
                   ["a a 40", "b b 41", "c c 42", "d d 43", "e e NA"])  # f has no row
    with caplog.at_level(logging.WARNING):
        ids, *_ = load_analysis_inputs(pheno, ["y"], covar)
    assert ids == ["a", "d"]                      # rows kept exactly as before
    text = caplog.text
    assert "Dropped 4 phenotype rows" in text
    assert "1 with a missing (NA) phenotype" in text
    assert "1 with a non-numeric or empty phenotype" in text
    assert "1 with a missing (NA) covariate" in text
    assert "1 with no row in the covariate file" in text
    assert "1 samples have the value -9 in column 'y'" in text


def test_clean_input_is_quiet(tmp_path, caplog):
    pheno = _write(tmp_path / "pheno.txt", "FID IID y", ["a a 1.0", "b b 2.0"])
    with caplog.at_level(logging.WARNING):
        ids, *_ = load_analysis_inputs(pheno, ["y"], None)
    assert ids == ["a", "b"]
    assert "Dropped" not in caplog.text and "value -9" not in caplog.text
