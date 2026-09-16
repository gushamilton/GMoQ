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


def _loco_mode_fixture(seed: int = 7, n: int = 4000):
    rng = np.random.default_rng(seed)
    X = np.column_stack([np.ones(n), rng.normal(size=n), rng.normal(size=n)])
    L = 0.2 * rng.normal(size=n)
    y = 2.0 + 0.5 * X[:, 1] + 1.5 * L + rng.standard_exponential(n)
    taus = np.linspace(0.1, 0.9, 5)
    return X, L, y, taus, np.ones(n, dtype=bool)


def test_loco_mode_covariate_leaves_the_unconditional_quantiles_alone():
    X, L, y, taus, mask = _loco_mode_fixture()
    preds = {"y": {"22": L}}

    _, _, q_none = compute_chr_rif_residuals(
        ["y"], {"y": y}, X, taus, None, None, mask
    )
    _, _, q_cov = compute_chr_rif_residuals(
        ["y"], {"y": y}, X, taus, preds, "22", mask, None, loco_mode="covariate"
    )
    _, _, q_sub = compute_chr_rif_residuals(
        ["y"], {"y": y}, X, taus, preds, "22", mask, {"y": 1.0}, loco_mode="subtract"
    )

    np.testing.assert_allclose(q_cov[0], q_none[0])
    # The historical mode estimates something else: the quantiles themselves move.
    assert not np.allclose(q_sub[0], q_none[0])


def test_loco_mode_covariate_residualises_on_x_and_l():
    X, L, y, taus, mask = _loco_mode_fixture()

    rif_resid, Q, _ = compute_chr_rif_residuals(
        ["y"], {"y": y}, X, taus, {"y": {"22": L}}, "22", mask,
        None, loco_mode="covariate"
    )

    rif, _, _ = core.compute_rif(y, taus)
    Q_ref, _ = np.linalg.qr(np.column_stack([X, L]))
    np.testing.assert_allclose(rif_resid[0], rif - Q_ref @ (Q_ref.T @ rif), atol=1e-10)
    # Q also projects the genotype in the scoring kernel, so the extra column has
    # to reach it: L is controlled for on both sides of the regression.
    assert Q.shape[1] == X.shape[1] + 1


def test_loco_mode_covariate_ignores_the_loco_scale():
    X, L, y, taus, mask = _loco_mode_fixture()

    a, _, _ = compute_chr_rif_residuals(
        ["y"], {"y": y}, X, taus, {"y": {"22": L}}, "22", mask,
        None, loco_mode="covariate"
    )
    b, _, _ = compute_chr_rif_residuals(
        ["y"], {"y": y}, X, taus, {"y": {"22": 7.3 * L}}, "22", mask,
        {"y": 1.088}, loco_mode="covariate"
    )

    np.testing.assert_allclose(a[0], b[0], atol=1e-8)


def test_loco_mode_subtract_is_the_untouched_default():
    X, L, y, taus, mask = _loco_mode_fixture()

    default, _, _ = compute_chr_rif_residuals(
        ["y"], {"y": y}, X, taus, {"y": {"22": L}}, "22", mask, {"y": 1.088}
    )
    expected, _, _ = build_rif_residuals({"y": y - 1.088 * L}, ["y"], X, taus)

    np.testing.assert_array_equal(default[0], expected[0])


def test_covariate_mode_changes_the_sigma0_fingerprint_but_subtract_keeps_it():
    from fungwas_stage1 import pooled_cov

    kw = dict(phenotype="y", covariates=["age"], taus=[0.1, 0.5, 0.9],
              n_samples=1000, loco_predictor="/p/y.loco")
    historical = pooled_cov.analysis_fingerprint(
        **kw, loco_scale_mode="regression", loco_scale_factor=1.088,
        density_method="silverman_gaussian")

    assert historical == pooled_cov.analysis_fingerprint(
        **kw, loco_scale_mode="regression", loco_scale_factor=1.088,
        density_method="silverman_gaussian", loco_mode=None)
    assert historical != pooled_cov.analysis_fingerprint(
        **kw, loco_scale_mode=None, loco_scale_factor=None,
        density_method="silverman_gaussian", loco_mode="covariate")


def test_sample_match_guard_refuses_warns_and_passes_at_the_right_fractions():
    """The guard itself, not a copy of its arithmetic.

    The three tests this replaces asserted `not (0 > 0 and frac < 0)`, grepped
    the source of `cli.main` for a literal argparse string, and recomputed
    `n_matched / n_total` inline -- none of them called the package, so all
    three passed whatever the guard did.
    """
    from fungwas_stage1.cli import sample_match_verdict, SAMPLE_MATCH_WARN_BELOW

    default = 0.01

    # A near-total failure is the ID-convention error the guard exists to catch.
    assert sample_match_verdict(0.0, default) == "refuse"
    assert sample_match_verdict(0.005, default) == "refuse"

    # Exactly at the threshold is not below it.
    assert sample_match_verdict(default, default) == "warn"

    # Legitimate partial subsetting proceeds, loudly. 700/900 is this package's
    # own fixture, with identical ID conventions on both sides.
    assert sample_match_verdict(0.3, default) == "warn"
    assert sample_match_verdict(700 / 900, default) == "ok"
    assert sample_match_verdict(1.0, default) == "ok"

    # The warning band is half-open at the top.
    assert sample_match_verdict(SAMPLE_MATCH_WARN_BELOW, default) == "ok"
    assert sample_match_verdict(SAMPLE_MATCH_WARN_BELOW - 1e-9, default) == "warn"

    # Disabling it proceeds deliberately, even on a total failure.
    for frac in (0.0, 0.03, 0.5, 1.0):
        assert sample_match_verdict(frac, 0) == "ok"


def test_loco_mode_default_follows_whether_loco_preds_were_given():
    """The default is resolved after argparse, so the parser default is None.

    The test this replaces asserted that the string
    "'--loco-scale', dest='loco_scale', default='regression'" appears in the
    source of `cli.main`, which says nothing about what the CLI resolves to and
    breaks on reformatting.
    """
    import argparse
    from fungwas_stage1.cli import (
        LOCO_MODE_COVARIATE, LOCO_MODE_SUBTRACT, resolve_loco_mode,
    )

    # With LOCO predictions, the covariate mode is the 0.4.0 default.
    args = argparse.Namespace(loco_mode=None, loco_preds="preds.list")
    assert resolve_loco_mode(args) == LOCO_MODE_COVARIATE

    # Without them there is nothing to add to the basis, so subtract is the only
    # mode that means anything.
    args = argparse.Namespace(loco_mode=None, loco_preds=None)
    assert resolve_loco_mode(args) == LOCO_MODE_SUBTRACT

    # An explicit choice is never overridden.
    args = argparse.Namespace(loco_mode=LOCO_MODE_SUBTRACT, loco_preds="preds.list")
    assert resolve_loco_mode(args) == LOCO_MODE_SUBTRACT

    # Covariate mode without predictions is a contradiction, not a silent fallback.
    args = argparse.Namespace(loco_mode=LOCO_MODE_COVARIATE, loco_preds=None)
    try:
        resolve_loco_mode(args)
    except ValueError:
        pass
    else:
        raise AssertionError("covariate mode without --loco-preds must raise")
