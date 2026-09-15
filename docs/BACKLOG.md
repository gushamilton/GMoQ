# G-MoQ backlog

One ordered list of everything open in **the method and the package**.
Updated 2026-09-15 at `ce9fb92`.

Nothing here is manuscript work. Regenerating any published result is out of scope by
definition: the method defect that would require it is already fixed and committed, and
re-running an analysis is not a method task.

Status key: **DONE** (committed, with evidence) / **OPEN** (known, not done) /
**UNVERIFIED** (claimed, evidence not confirmed) / **UNKNOWN** (never looked).

---

## A. Blocking - a user would get a wrong answer

| # | Item | Status |
|---|---|---|
| A1 | **The input pipeline has never been checked against ground truth.** `extract_dosage`, `match_samples`, `load_analysis_inputs`, `read_keep_ids`, `write_stats_batch` have zero test references; both tests touching `get_variants_bgen` monkeypatch it away; all six simulation campaigns call `core.run_stage1()` directly. Everything passes today and nothing would notice if it were wrong. Tier 4 of `test/suite` is built and unrun. | **OPEN** |
| A2 | **No sample-match-rate guard.** `cli.py:1262` logs `Matched N/M` and proceeds. `match_samples` keys on FID alone while LOCO matching uses five aliases, so an ID-convention mismatch silently analyses a small non-random subset. Three lines to fix; protects every future user rather than every future test. | **OPEN** |

## B. Decisions the evidence does not yet settle

| # | Item | Status |
|---|---|---|
| B1 | **Should `cov_source` default to pooled?** Currently jackknife. The pooled-OLS vs jackknife-OLS comparison - the one that decides it - has never been measured. `FINDINGS_pooled_dispersion.md`'s jackknife claim rests on a table with no pooled-OLS column. Simulation running on `sim/covsource`. | **OPEN** |
| B2 | **Does empirical `W` help anywhere?** Settled as no overall (`EMPIRICAL_W.md`), but it helps chi2(4) under OLS by 83x while still leaving 108x nominal. Whether that is worth a mode is undecided. | **OPEN** |
| B3 | **Q calibration on the jackknife path.** `Q_CALIBRATION` verified with an exact Wishart Sigma independent of `b`. A block jackknife is neither exactly Wishart nor independent of `b`, and its `nu` comes from a law fitted for a scalar contrast. The default path rests on derivation, not measurement. | **UNVERIFIED** |
| B4 | **Q on int8-quantised Sigma.** GLS needs an inverse the stored covariance is documented as unfit for. Bounds where the guard can be trusted. | **UNKNOWN** |

## C. Claims in the docs that need checking

Audit running. Each of these has already produced one wrong claim of the same shape.

| # | Item | Status |
|---|---|---|
| C1 | `FINDINGS_pooled_dispersion.md` - "the jackknife is a genuine and complete fix": which estimator produced 5.00e-08? | **UNVERIFIED** |
| C2 | `WHY_GLS_DEGRADES.md` - the restated `Q/(T-K)` thresholds. The first version was wrong; is the restatement right? | **UNVERIFIED** |
| C3 | `BP2_lumpy_phenotypes.md` - "KDE error cancels exactly in z" (5e-14); "ceiling censoring 270x at 1e-6" carries ~14% MC error and the extrapolation to 5e-8 is unmeasured. | **UNVERIFIED** |
| C4 | `RELATEDNESS.md` - `lambda_sigma2 - 1` proportional to `rho^2`; the claim that Stage 1's covariate step is stronger than REGENIE's. | **UNVERIFIED** |
| C5 | `N_SCALING.md` - the genome-wide extrapolation revised from 0.99-1.03 to 0.87-1.19. Is the revision sound? | **UNVERIFIED** |
| C6 | `LOCO_SCALE.md` - cor 0.500 -> 0.993 against REGENIE: confirm the allele harmonisation makes it like-for-like. | **UNVERIFIED** |

## D. Code that does not do what it says

| # | Item | Status |
|---|---|---|
| D1 | `make_weights_vqtl` validates `q_tau` and `mu`, then discards both and hardcodes `qnorm(taus)`. Both are mandatory with no defaults and the `@details` presents them as load-bearing. | **OPEN** |
| D2 | `make_weights_generic` declares `dist_cdf` as mandatory and never references it. The documented example passes it positionally. | **OPEN** |
| D3 | `dispersion_report` computes p-values with `2*pnorm()` - the normal reference `tail-reference.R` argues at length must not be used for Stage 2 parameters. | **OPEN** |
| D4 | `estimate_rav(se_sigma2)` used only in a finiteness filter, never in the estimator. `estimate_rav(trim)`'s comment states the opposite of what the line beneath it does. | **OPEN** |
| D5 | `quantile_gwas` uses `bw.nrd0` where Python Stage 1 uses `1.06*sd*N^-0.2`, so the two paths differ by construction. Exported, used by nine demos, asserted by none. | **OPEN** |
| D6 | `quantile_gwas(covar_mode="include")` is documented, accepted by `match.arg`, and `stop()`s on the first line of the branch. | **OPEN** |
| D7 | `dispersion_skew_guard(N)` documented as "required unless `maf` is supplied"; the body requires it unconditionally. | **OPEN** |

## E. Testing

| # | Item | Status |
|---|---|---|
| E1 | **Nine of 23 Stage 2 scripts assert nothing** - including `test-qstat.R`, the only file named after the Q statistic, and `proof_stage2_se.R`, cited as "step 3 of the validation chain". They consume 89% of suite runtime (591s of 667s) while the three files carrying 248 of ~300 assertions run in under two seconds. They already compute the right quantities and `cat()` them; turning prints into `stopifnot()` is a small edit. | **OPEN** |
| E2 | **Run tiers 0 and 1 of `test/suite`** - ~2 minutes, additive-only, closes the largest Stage 2 gaps: every weight builder against a finite difference of the real quantile function, and the Q tail against `CompQuadForm::davies`. | **OPEN** |
| E3 | Nine test scripts call `devtools::document()` at line 5 - the exact mechanism that silently dropped `useDynLib` for three commits. | **OPEN** |
| E4 | Jackknife SE tolerances of +/-20-25% pass a systematic 19% error. Replace with the exact HC0 sandwich comparison already implemented at `test_pooled_cov.py:58`. | **OPEN** |
| E5 | `test_multi_phenotype_perf.py` collects zero tests - a `__main__` benchmark in the test directory. | **OPEN** |
| E6 | CI runs 3 of 12 Python files and 2 of 23 R scripts, `R CMD INSTALL` but never `R CMD check`. | **OPEN** |

## F. Packaging

| # | Item | Status |
|---|---|---|
| F1 | No top-level `LICENSE`; Stage 1 ships no license text though `pyproject.toml` declares MIT. | **OPEN** |
| F2 | Version numbers mean nothing: both components say `0.1.0`, tags say `v0.2.0`/`v0.3.0`, code advertises "pre-v0.4.0". No `NEWS`/`CHANGELOG`. | **OPEN** |
| F3 | READMEs never mention `--loco-scale`, `tau_cov_estimator`, `cov_source` or the dispersion guards - the substance of the last 31 commits. | **OPEN** |
| F4 | `param_gwas_from_file` is exported with no man page; its roxygen is attached to the internal `.validate_cov_ids`. `q_method` undocumented. | **OPEN** |
| F5 | Phantom functions in docs: `make_weights_lognormal`, `make_weights_mixture_vqtl`. | **OPEN** |
| F6 | `docs/slurm_templates/submit_stage1.sh` fails three ways against the current CLI. | **OPEN** |
| F7 | No end-to-end synthetic example that runs as written; `docs/examples/` predates the pooled path. | **OPEN** |

## G. Unmerged branches

`sim/dispersion` (1), `sim/nscale` (2), `test/suite` (3), `wt/density` (3) carry evidence and
scripts not on master. Decide: merge the evidence, or archive the branch.

---

## Done this session

| Item | Commit |
|---|---|
| LOCO scaling bug - found, fixed, validated against REGENIE | `8ca9e0b`..`b648947` |
| LOCO scale recorded and hashed into the fingerprint | `24eef02` |
| Pooled covariance Sigma_0 | `dc19377`..`4658a17` |
| False vQTLs from genotype skewness - corrected | `c5b2395` |
| Three guardrails that refuse rather than mislead | `0ec64c5` |
| OLS default settled on measured coverage | `497d074`, `b5024d8`, `44681df` |
| `Q_pval` calibrated on the GLS path (Hotelling's T2) | `ce9fb92` |
| NAMESPACE root fix - roxygen regenerates the C++ block | `db7bff4` + package doc |
| Package installs and builds outside this cluster | `e96f3ef` |
| Allele convention documented | `bace033` |
| Relatedness, density, tau grid - closed | `6807d6e`, `56f1059` |
