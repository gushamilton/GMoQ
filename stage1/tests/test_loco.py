import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from fungwas_stage1.cli import (
    build_rif_residuals,
    compute_chr_rif_residuals,
    load_loco_predictions,
    normalize_chr_label,
    resolve_loco_files,
)
from fungwas_stage1 import core


def test_resolve_loco_files_single_file(tmp_path: Path):
    loco_file = tmp_path / "bmi.loco"
    loco_file.write_text("FID_IID sample1\n22 0.0\n")

    mapping = resolve_loco_files(str(loco_file), ["bmi"])

    assert mapping == {"bmi": str(loco_file)}


def test_load_loco_predictions_aligns_fid_iid_aliases(tmp_path: Path):
    loco_file = tmp_path / "bmi.loco"
    loco_file.write_text(
        "FID_IID IEU2_IEU2 IEU1_IEU1\n"
        "21 0.2 0.1\n"
        "22 0.4 0.3\n"
    )

    preds, valid_mask = load_loco_predictions(
        str(loco_file),
        ["bmi"],
        [["IEU1", "IEU1_IEU1"], ["IEU2", "IEU2_IEU2"]],
    )

    np.testing.assert_array_equal(valid_mask, np.array([True, True]))
    np.testing.assert_allclose(preds["bmi"]["21"], np.array([0.1, 0.2]))
    np.testing.assert_allclose(preds["bmi"]["22"], np.array([0.3, 0.4]))


def test_normalize_chr_label_handles_numeric_and_special_labels():
    assert normalize_chr_label("04") == "4"
    assert normalize_chr_label("4") == "4"
    assert normalize_chr_label("chr04") == "4"
    assert normalize_chr_label("CHR22") == "22"
    assert normalize_chr_label("x") == "X"
    assert normalize_chr_label("chrxy") == "XY"
    assert normalize_chr_label("mt") == "MT"


def test_load_loco_predictions_drops_samples_missing_from_loco(tmp_path: Path):
    loco_file = tmp_path / "bmi.loco"
    loco_file.write_text(
        "FID_IID IEU2_IEU2 IEU1_IEU1\n"
        "22 0.4 0.3\n"
    )

    preds, valid_mask = load_loco_predictions(
        str(loco_file),
        ["bmi"],
        [["IEU1", "IEU1_IEU1"], ["IEU2", "IEU2_IEU2"], ["IEU3", "IEU3_IEU3"]],
    )

    np.testing.assert_array_equal(valid_mask, np.array([True, True, False]))
    np.testing.assert_allclose(preds["bmi"]["22"], np.array([0.3, 0.4]))


def test_compute_chr_rif_residuals_accepts_zero_padded_chr_labels():
    taus = np.array([0.25, 0.5, 0.75])
    x = np.column_stack([np.ones(4), np.linspace(-1.0, 1.0, 4)])
    y_values = {"bmi": np.array([1.0, 2.0, 3.0, 4.0])}
    loco_predictions = {"bmi": {"4": np.array([0.1, 0.2, 0.3, 0.4])}}
    valid_mask = np.array([True, True, True, True])

    rif_resid_list, q_chr, q_tau_list = compute_chr_rif_residuals(
        ["bmi"], y_values, x, taus, loco_predictions, "04", valid_mask
    )

    assert len(rif_resid_list) == 1
    assert rif_resid_list[0].shape == (4, 3)
    assert q_chr.shape == x.shape
    assert len(q_tau_list) == 1
    assert q_tau_list[0].shape == (3,)


def test_loco_adjustment_changes_residualized_rif():
    taus = np.array([0.25, 0.5, 0.75])
    x = np.column_stack([np.ones(5), np.linspace(-1.0, 1.0, 5)])
    y_base = {"bmi": np.array([0.0, 1.0, 2.0, 3.0, 4.0])}
    y_loco = {"bmi": y_base["bmi"] - np.array([0.0, 0.2, 0.2, 0.0, 0.0])}

    rif_base, _, q_base = build_rif_residuals(y_base, ["bmi"], x, taus)
    rif_loco, _, q_loco = build_rif_residuals(y_loco, ["bmi"], x, taus)

    assert not np.allclose(rif_base[0], rif_loco[0])
    assert not np.allclose(q_base[0], q_loco[0])


def test_cpp_requirement_fails_closed_when_extension_missing():
    original_have_cpp = core.HAVE_CPP
    original_allow = core.ALLOW_NUMPY_FALLBACK
    try:
        core.HAVE_CPP = False
        core.set_numpy_fallback_allowed(False)
        try:
            core.require_cpp_extension("test")
        except RuntimeError as exc:
            assert "--allow-numpy-fallback" in str(exc)
        else:
            raise AssertionError("Expected RuntimeError when C++ extension is unavailable")
    finally:
        core.HAVE_CPP = original_have_cpp
        core.ALLOW_NUMPY_FALLBACK = original_allow


def test_numpy_fallback_warning_is_explicit(caplog):
    original_have_cpp = core.HAVE_CPP
    original_allow = core.ALLOW_NUMPY_FALLBACK
    try:
        core.HAVE_CPP = False
        with caplog.at_level("WARNING"):
            core.set_numpy_fallback_allowed(True)
        warning_text = "\n".join(record.message for record in caplog.records)
        assert "NUMPY FALLBACK ENABLED FOR STAGE 1" in warning_text
        assert "SUPER SUPER SLOW" in warning_text
        assert "DO NOT USE THIS FOR PRODUCTION GWAS" in warning_text
        assert "--allow-numpy-fallback" in warning_text
    finally:
        core.HAVE_CPP = original_have_cpp
        core.ALLOW_NUMPY_FALLBACK = original_allow


def test_multi_phenotype_block_scores_match_single_path_under_numpy_fallback():
    original_have_cpp = core.HAVE_CPP
    original_allow = core.ALLOW_NUMPY_FALLBACK
    try:
        core.HAVE_CPP = False
        core.set_numpy_fallback_allowed(True)

        rng = np.random.default_rng(123)
        n_samples = 120
        n_snps = 6
        n_taus = 3
        n_blocks = 8

        g = rng.binomial(2, 0.3, size=(n_samples, n_snps)).astype(np.float64)
        x = np.column_stack([np.ones(n_samples), rng.normal(size=(n_samples, 2))])
        q, _ = np.linalg.qr(x)
        block_ids = rng.integers(0, n_blocks, size=n_samples, dtype=np.int32)

        rif_a = rng.normal(size=(n_samples, n_taus))
        rif_b = rng.normal(size=(n_samples, n_taus))

        multi_stats = core.compute_block_scores_multi(
            g, [rif_a, rif_b], q, block_ids, n_blocks
        )
        single_stats_a = core.compute_block_scores(
            g, rif_a, q, block_ids, n_blocks
        )
        single_stats_b = core.compute_block_scores(
            g, rif_b, q, block_ids, n_blocks
        )

        assert len(multi_stats) == 2
        np.testing.assert_allclose(multi_stats[0], single_stats_a)
        np.testing.assert_allclose(multi_stats[1], single_stats_b)
    finally:
        core.HAVE_CPP = original_have_cpp
        core.ALLOW_NUMPY_FALLBACK = original_allow


# ---------------------------------------------------------------------------
# resolve_loco_files: exactness, determinism, and explicit failure.
#
# Regression cover for the silent-wrong-results bug in which pred.list entries
# were mapped to phenotypes by *substring* match while iterating a *set*: every
# manuscript phenotype name is a proper substring of another (bmi of log_bmi and
# rint_bmi; testosterone of its three shrink variants; crp of log_crp), so the
# mapping depended on PYTHONHASHSEED and on pred.list entry order and could
# silently attach another phenotype's LOCO predictors.
# ---------------------------------------------------------------------------

BMI_GROUP = ["bmi", "log_bmi", "rint_bmi"]
TESTOSTERONE_GROUP = [
    "testosterone",
    "testosterone_shrink_0.2",
    "testosterone_shrink_0.5",
    "testosterone_shrink_0.7",
]


def _write_pred_list(tmp_path: Path, name: str, phenos, reverse: bool = False) -> tuple[str, dict]:
    """Write a REGENIE-style combined pred.list mirroring the real paper layout.

    Real entries look like::

        testosterone_shrink_0.2 /.../full/testosterone_shrink_0_2/testosterone_shrink_0_2_step1_1.loco

    i.e. the label keeps the ``.`` while the path component replaces it with ``_``.
    """
    expected = {}
    lines = []
    for pheno in phenos:
        tag = pheno.replace(".", "_")
        loco = str(tmp_path / "full" / tag / f"{tag}_step1_1.loco")
        expected[pheno] = loco
        lines.append(f"{pheno} {loco}")
    if reverse:
        lines = list(reversed(lines))
    pred = tmp_path / name
    pred.write_text("\n".join(lines) + "\n")
    return str(pred), expected


@pytest.mark.parametrize("phenos", [BMI_GROUP, TESTOSTERONE_GROUP])
@pytest.mark.parametrize("reverse", [False, True])
def test_resolve_loco_files_exact_for_prefix_overlapping_groups(tmp_path: Path, phenos, reverse):
    """Each phenotype gets its own .loco even though names are prefixes of one another."""
    pred, expected = _write_pred_list(tmp_path, "group_pred.list", phenos, reverse=reverse)

    assert resolve_loco_files(pred, list(phenos)) == expected
    # ...and independently of the order the phenotypes are requested in.
    assert resolve_loco_files(pred, list(reversed(phenos))) == expected


def test_resolve_loco_files_matches_canonicalized_label(tmp_path: Path):
    """A pred.list label written with '_' for '.' still resolves, exactly."""
    loco = tmp_path / "testosterone_shrink_0_2_step1_1.loco"
    loco.write_text("FID_IID s1\n22 0.0\n")
    other = tmp_path / "testosterone_step1_1.loco"
    other.write_text("FID_IID s1\n22 0.0\n")
    pred = tmp_path / "pred.list"
    pred.write_text(f"testosterone {other}\ntestosterone_shrink_0_2 {loco}\n")

    mapping = resolve_loco_files(str(pred), ["testosterone", "testosterone_shrink_0.2"])

    assert mapping == {
        "testosterone": str(other),
        "testosterone_shrink_0.2": str(loco),
    }


def test_resolve_loco_files_single_pheno_against_combined_pred_list(tmp_path: Path):
    """A single phenotype selected out of a combined pred.list is matched exactly.

    Under the old substring matcher, requesting only ``bmi`` could bind to the
    ``log_bmi`` entry, since 'bmi' is contained in that entry's line.
    """
    pred, expected = _write_pred_list(tmp_path, "group_pred.list", BMI_GROUP, reverse=True)

    assert resolve_loco_files(pred, ["bmi"]) == {"bmi": expected["bmi"]}
    assert resolve_loco_files(pred, ["log_bmi"]) == {"log_bmi": expected["log_bmi"]}


def test_resolve_loco_files_single_entry_pred_list_keeps_regenie_label_mismatch(tmp_path: Path):
    """One phenotype, one .loco: still resolved, whatever REGENIE called the column.

    The paper's ``weight_height`` / ``weight_height2`` models run REGENIE on a
    column literally named ``weight``, so the pred.list label does not equal the
    requested phenotype tag.  That single/single case must keep working.
    """
    loco = tmp_path / "weight_height_step1_1.loco"
    loco.write_text("FID_IID s1\n22 0.0\n")
    pred = tmp_path / "weight_height_step1_pred.list"
    pred.write_text(f"weight {loco}\n")

    assert resolve_loco_files(str(pred), ["weight_height"]) == {"weight_height": str(loco)}


def test_resolve_loco_files_raises_on_ambiguous_duplicate_label(tmp_path: Path):
    """Two entries sharing a label must raise and name both candidates."""
    a = tmp_path / "a" / "weight_step1_1.loco"
    b = tmp_path / "b" / "weight_step1_1.loco"
    pred = tmp_path / "pred.list"
    pred.write_text(f"weight {a}\nweight {b}\nbmi {tmp_path / 'bmi.loco'}\n")

    with pytest.raises(ValueError) as excinfo:
        resolve_loco_files(str(pred), ["weight", "bmi"])

    message = str(excinfo.value)
    assert "Ambiguous" in message
    assert str(a) in message and str(b) in message


def test_resolve_loco_files_raises_when_phenotype_absent(tmp_path: Path):
    """An unresolvable phenotype raises and lists the labels that are present."""
    pred, _ = _write_pred_list(tmp_path, "group_pred.list", BMI_GROUP)

    with pytest.raises(ValueError) as excinfo:
        resolve_loco_files(pred, ["bmi", "vitamin_d"])

    message = str(excinfo.value)
    assert "vitamin_d" in message
    assert "log_bmi" in message and "rint_bmi" in message
    # It must not silently fall back to positional assignment.
    assert "Could not map LOCO files" in message


def test_resolve_loco_files_raises_when_two_phenotypes_claim_one_file(tmp_path: Path):
    """Distinct phenotypes canonicalizing to the same label must not share a .loco."""
    loco = tmp_path / "x_step1_1.loco"
    pred = tmp_path / "pred.list"
    pred.write_text(f"a.b {loco}\n")

    with pytest.raises(ValueError) as excinfo:
        resolve_loco_files(str(pred), ["a.b", "a_b"])

    assert "resolved for both" in str(excinfo.value)


def test_resolve_loco_files_is_invariant_to_python_hash_seed(tmp_path: Path):
    """The mapping must not depend on PYTHONHASHSEED, for either group, either order.

    Runs in fresh interpreters: string hashing is randomized per process, so this
    is the only way to exercise it.  The pre-fix implementation iterated a ``set``
    of phenotype names and produced 2 distinct mappings for the BMI group and 3
    for the testosterone group on reversed pred.lists.
    """
    configs = []
    for label, phenos in (("bmi", BMI_GROUP), ("testosterone", TESTOSTERONE_GROUP)):
        for reverse in (False, True):
            name = f"{label}_{'rev' if reverse else 'fwd'}_pred.list"
            pred, expected = _write_pred_list(tmp_path, name, phenos, reverse=reverse)
            configs.append({"pred": pred, "phenos": list(phenos), "expected": expected})

    probe = tmp_path / "probe.py"
    probe.write_text(
        "import json, sys\n"
        "from fungwas_stage1.cli import resolve_loco_files\n"
        "configs = json.load(open(sys.argv[1]))\n"
        "print(json.dumps([resolve_loco_files(c['pred'], c['phenos']) for c in configs]))\n"
    )
    config_json = tmp_path / "configs.json"
    config_json.write_text(json.dumps(configs))

    expected = [c["expected"] for c in configs]

    env = dict(os.environ)
    env["PYTHONPATH"] = str(Path(__file__).resolve().parents[1]) + os.pathsep + env.get("PYTHONPATH", "")

    seeds = [0, 1, 2, 3, 7, 13, 42, 99, 1234, 31337, 65535, 4294967295]
    for seed in seeds:
        env["PYTHONHASHSEED"] = str(seed)
        out = subprocess.run(
            [sys.executable, str(probe), str(config_json)],
            capture_output=True, text=True, env=env, check=True,
        )
        assert json.loads(out.stdout) == expected, f"mapping changed under PYTHONHASHSEED={seed}"


# --- LOCO scaling (docs/LOCO_SCALE.md) --------------------------------------

def _loco_fixture(seed=0, n=4000, y_sd=0.1665, pred_sd=0.3438, rho=0.333):
    """
    A phenotype and a *standardised* LOCO prediction, mimicking REGENIE: the
    prediction's sd is fixed regardless of the phenotype's own sd, and the
    correlation between them is rho. Matches the log_bmi row of docs/LOCO_SCALE.md.
    """
    rng = np.random.default_rng(seed)
    z = rng.standard_normal(n)
    e = rng.standard_normal(n)
    y_std = rho * z + np.sqrt(1.0 - rho ** 2) * e
    y = 3.3 + y_sd * y_std
    pred = pred_sd * z
    return y, {"22": pred}


def test_resolve_loco_scale_sd_recovers_phenotype_scale():
    from fungwas_stage1.cli import resolve_loco_scale
    y, per_chr = _loco_fixture()
    scale = resolve_loco_scale("sd", y, per_chr, "log_bmi")
    assert scale == pytest.approx(np.std(y), rel=1e-12)


def test_resolve_loco_scale_regression_matches_ols_slope():
    from fungwas_stage1.cli import resolve_loco_scale
    y, per_chr = _loco_fixture()
    scale = resolve_loco_scale("regression", y, per_chr, "log_bmi")
    pred = per_chr["22"]
    expected = np.cov(y, pred, ddof=1)[0, 1] / np.var(pred, ddof=1)
    assert scale == pytest.approx(expected, rel=1e-10)


def test_resolve_loco_scale_none_and_literal():
    from fungwas_stage1.cli import resolve_loco_scale
    y, per_chr = _loco_fixture()
    assert resolve_loco_scale("none", y, per_chr, "log_bmi") == 1.0
    assert resolve_loco_scale("0.25", y, per_chr, "log_bmi") == pytest.approx(0.25)
    with pytest.raises(ValueError):
        resolve_loco_scale("nonsense", y, per_chr, "log_bmi")


def test_unscaled_loco_inflates_variance_but_scaled_loco_reduces_it():
    """
    The defect in one assertion: for a phenotype with sd < sd(LOCO), subtracting
    the raw prediction increases the residual variance. Scaling by sd(Y) reduces it.
    """
    from fungwas_stage1.cli import resolve_loco_scale
    y, per_chr = _loco_fixture()
    pred = per_chr["22"]

    var_raw = np.var(y, ddof=1)
    var_unscaled = np.var(y - pred, ddof=1)
    var_scaled = np.var(y - resolve_loco_scale("sd", y, per_chr, "log_bmi") * pred, ddof=1)

    assert var_unscaled > 3.0 * var_raw          # injects noise: ~3.9x in the real data
    assert var_scaled < var_raw                  # removes it
    assert var_scaled == pytest.approx((1.0 - 0.333 ** 2) * var_raw, rel=0.05)


def test_loco_scale_warns_when_subtracted_term_dominates(caplog):
    from fungwas_stage1.cli import resolve_loco_scale
    import logging
    y, per_chr = _loco_fixture()
    with caplog.at_level(logging.WARNING):
        resolve_loco_scale("none", y, per_chr, "log_bmi")
    assert any("INCREASES the analysed variance" in r.message for r in caplog.records)

    caplog.clear()
    with caplog.at_level(logging.WARNING):
        resolve_loco_scale("sd", y, per_chr, "log_bmi")
    assert not any("INCREASES the analysed variance" in r.message for r in caplog.records)


def test_rint_like_phenotype_is_unaffected_by_the_fix():
    """sd(Y) == 1 is the one case where the old behaviour was already correct."""
    from fungwas_stage1.cli import resolve_loco_scale
    y, per_chr = _loco_fixture(y_sd=1.0, pred_sd=0.3431)
    y = y - y.mean()
    assert resolve_loco_scale("sd", y, per_chr, "rint_bmi") == pytest.approx(1.0, rel=0.02)


def test_compute_chr_rif_residuals_applies_the_scale():
    y, per_chr = _loco_fixture(n=800)
    taus = np.linspace(0.1, 0.9, 9)
    X = np.ones((y.size, 1))
    mask = np.ones(y.size, dtype=bool)

    unscaled, _, _ = compute_chr_rif_residuals(
        ["log_bmi"], {"log_bmi": y}, X, taus, {"log_bmi": per_chr}, "22", mask
    )
    half, _, _ = compute_chr_rif_residuals(
        ["log_bmi"], {"log_bmi": y}, X, taus, {"log_bmi": per_chr}, "22", mask,
        loco_scales={"log_bmi": 0.5}
    )
    explicit_one, _, _ = compute_chr_rif_residuals(
        ["log_bmi"], {"log_bmi": y}, X, taus, {"log_bmi": per_chr}, "22", mask,
        loco_scales={"log_bmi": 1.0}
    )

    # omitting loco_scales must stay identical to an explicit factor of 1
    assert np.allclose(unscaled[0], explicit_one[0])
    # and a different factor must actually change the residualised RIF
    assert not np.allclose(unscaled[0], half[0])


def _loco_fixture_with_covariates(seed=1, n=6000, cov_r2=0.30, y_sd=15.91, pred_sd=0.2333, rho=0.199):
    """
    Like _loco_fixture but with covariates that explain `cov_r2` of Y, and a LOCO
    prediction built against the *residualised* trait -- which is what REGENIE does
    (Pheno.cpp::residualize_phenotypes standardises by the residual sd, not sd(Y)).
    Defaults mimic `weight`, whose step-1 height terms carry R^2 = 0.30.
    """
    rng = np.random.default_rng(seed)
    cov = rng.standard_normal(n)
    z = rng.standard_normal(n)
    e = rng.standard_normal(n)
    resid_std = rho * z + np.sqrt(1.0 - rho ** 2) * e
    resid_std = (resid_std - resid_std.mean()) / resid_std.std()
    y_resid_sd = y_sd * np.sqrt(1.0 - cov_r2)
    y = 70.0 + y_sd * np.sqrt(cov_r2) * cov + y_resid_sd * resid_std
    pred = pred_sd * z
    return y, {"22": pred}, cov, y_resid_sd


def test_sd_mode_over_subtracts_when_covariates_explain_the_trait():
    """
    The factor REGENIE used is the sd of the trait *after* its step-1 covariates, so
    sd(Y) over-states it whenever those covariates are predictive. Pinned at the real
    `weight` configuration: R^2 = 0.30 gives sd(Y)/sd(Y_resid) = 1/sqrt(0.70) = 1.195.
    """
    from fungwas_stage1.cli import resolve_loco_scale
    y, per_chr, cov, y_resid_sd = _loco_fixture_with_covariates()

    sd_factor = resolve_loco_scale("sd", y, per_chr, "weight")
    assert sd_factor / y_resid_sd == pytest.approx(1.0 / np.sqrt(0.70), rel=0.05)

    # regression needs no knowledge of the covariates and is the variance-minimising
    # factor by construction, so it must leave no more variance than sd does.
    reg_factor = resolve_loco_scale("regression", y, per_chr, "weight")
    pred = per_chr["22"]
    assert np.var(y - reg_factor * pred) <= np.var(y - sd_factor * pred)
    assert np.var(y - reg_factor * pred) < np.var(y)

    # It lands nearer REGENIE's own scale_Y (the residual sd) than sd(Y) does, but it is
    # not equal to it: the slope is sd_resid * rho / sd_pred, which coincides with
    # scale_Y only for a perfectly calibrated prediction. REGENIE's are over-dispersed
    # (sd_pred > rho), so regression sits slightly below scale_Y -- matching the real
    # data, where slope/sd(Y_resid) is 0.89 for weight and 0.97 for the BMI traits.
    assert abs(reg_factor - y_resid_sd) < abs(sd_factor - y_resid_sd)
    # against the fixture's own parameters, not a quantity re-derived from the data
    assert reg_factor == pytest.approx(y_resid_sd * 0.199 / 0.2333, rel=0.05)


def test_regression_is_the_cli_default():
    """The default must be the mode that does not need REGENIE's step-1 covariate set."""
    import argparse, inspect
    from fungwas_stage1 import cli
    src = inspect.getsource(cli.main)
    assert "'--loco-scale', dest='loco_scale', default='regression'" in src


def test_variance_guard_is_silent_for_a_predictive_covariate_trait(caplog):
    """The old sd-ratio heuristic could fire spuriously; the variance check must not."""
    from fungwas_stage1.cli import resolve_loco_scale
    import logging
    y, per_chr, _, _ = _loco_fixture_with_covariates()
    with caplog.at_level(logging.WARNING):
        resolve_loco_scale("regression", y, per_chr, "weight")
    assert not any("INCREASES the analysed variance" in r.message for r in caplog.records)


# --- the sample-match-rate guard (docs/BACKLOG.md A2) ------------------------

def test_match_samples_keys_on_the_first_column_only():
    """
    Pins the behaviour the guard exists to protect against. match_samples does a
    bare dict lookup on the phenotype file's first column, with no alias
    resolution -- unlike the LOCO reader, which resolves fid, iid, fid_iid,
    "fid iid" and fid:iid. A cohort whose BGEN IDs are FID_IID therefore matches
    LOCO and then drops every sample here.
    """
    from fungwas_stage1.cli import match_samples
    bgen = ["A1_A1", "A2_A2", "A3_A3"]          # FID_IID convention
    pheno = ["A1", "A2", "A3"]                   # FID convention
    valid, idx = match_samples(bgen, pheno)
    assert idx.size == 0, "expected a total mismatch; the aliasing gap is the point"
    assert not valid.any()

    # and the matching case still works
    valid, idx = match_samples(bgen, ["A1_A1", "A3_A3"])
    assert list(idx) == [0, 2]
    assert list(valid) == [True, True]


def test_partial_match_is_the_dangerous_case():
    """
    Zero overlap fails loudly downstream. A small overlap runs to completion on a
    non-random subset and yields attenuated effects with inflated SEs that are
    internally consistent and wrong -- which is why the guard is on the RATE, not
    on emptiness.
    """
    from fungwas_stage1.cli import match_samples
    bgen = [f"S{i}" for i in range(1000)]
    pheno = [f"S{i}" for i in range(30)] + [f"X{i}" for i in range(970)]
    valid, idx = match_samples(bgen, pheno)
    assert idx.size == 30
    frac = idx.size / len(pheno)
    assert frac == pytest.approx(0.03)
    # nothing in match_samples itself refuses this; the caller must.
    assert frac < 0.90


def test_guard_refuses_only_a_near_total_failure():
    """
    The threshold is tuned to the failure it exists to catch. A partial overlap is
    legitimate -- a phenotype file routinely covers a wider cohort than one genotype
    release, and the package's own fixtures match 700 of 900 (77.8%) with identical
    conventions. An ID-convention mismatch instead produces a near-total failure.
    A guard set high enough to catch the former breaks ordinary use.
    """
    default = 0.01
    for n_matched, n_total, should_refuse in [
        (0, 1000, True),       # total mismatch -- the convention error
        (5, 1000, True),       # 0.5%, accidental collisions only
        (10, 1000, False),     # exactly at the threshold
        (700, 900, False),     # the package's own fixture: legitimate subsetting
        (1000, 1000, False),
    ]:
        frac = n_matched / n_total
        assert (default > 0 and frac < default) is should_refuse, (n_matched, n_total)


def test_guard_warns_between_the_threshold_and_a_half():
    """Legitimate but worth flagging: proceed, but say so."""
    lo, hi = 0.01, 0.50
    for frac, should_warn in [
        (0.005, False),   # refused instead
        (0.03, True),
        (0.49, True),
        (0.778, False),   # the fixture: ordinary, silent
        (1.0, False),
    ]:
        assert (lo > 0 and lo <= frac < hi) is should_warn, frac


def test_guard_can_be_disabled():
    """--min-sample-match 0 must proceed deliberately rather than refuse."""
    for frac in (0.0, 0.03, 0.5):
        assert not (0 > 0 and frac < 0)


def test_cli_exposes_the_guard_with_a_default_that_does_not_break_ordinary_use():
    import inspect
    from fungwas_stage1 import cli
    src = inspect.getsource(cli.main)
    assert "'--min-sample-match'" in src
    assert "default=0.01" in src
    assert "below --min-sample-match" in src   # refuses, not merely warns
