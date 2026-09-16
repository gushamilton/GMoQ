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
| ~~A2~~ | ~~**No sample-match-rate guard.** `cli.py:1262` logs `Matched N/M` and proceeds. `match_samples` keys on FID alone while LOCO matching uses five aliases, so an ID-convention mismatch silently analyses a small non-random subset. Three lines to fix; protects every future user rather than every future test. | **OPEN** |

## B. Decisions the evidence does not yet settle

| # | Item | Status |
|---|---|---|
| B1 | **REOPENED 2026-09-16 by the B5 measurement. `cov_source` default.** Currently `pooled`, set on the view that the 780,000-scan campaign used a W we do not ship. That is now measured on the W we DO ship and the result is worse, not better: `beta_left`/`beta_right` inflated 3510x at MAF 0.005/Z30 and 236x at MAF 0.01/Z30, against `jackknife` nominal in every regime. Estimates are bit-identical between sources so switching costs nothing in effect sizes and under 1% in SE on real chr22 data. The counter-argument is that common variants are safe at any effect size reachable in practice (onset \|Z\| = 23.2 at MAF 0.05, 87.9 at MAF 0.30) and pooled is much cheaper. **Decision pending.** | **REOPENED - evidence favours jackknife** |
| B2 | **Does empirical `W` help anywhere?** Settled as no overall (`EMPIRICAL_W.md`), but it helps chi2(4) under OLS by 83x while still leaving 108x nominal. Whether that is worth a mode is undecided. | **OPEN** |
| B3 | **Q calibration on the jackknife path.** `Q_CALIBRATION` verified with an exact Wishart Sigma independent of `b`. A block jackknife is neither exactly Wishart nor independent of `b`, and its `nu` comes from a law fitted for a scalar contrast. The explicit jackknife path rests on derivation, not measurement. | **UNVERIFIED** |
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

## H. Opened 2026-09-16 - LOCO and efficiency

| id | item | status |
|---|---|---|
| H1 | **REOPENED then CLOSED differently, 2026-09-16. Genomic inflation without LOCO.** The earlier resolution claimed the 5-9% lambda rise on removing LOCO was real lost relatedness control. That was wrong. Masking the signal (1 Mb around every REGENIE |z|>=5 hit, 5,114 variants including SHBG) leaves lambda unchanged or higher in every arm, so at N=400,000 this lambda is polygenicity, not confounding, and it never had the power to adjudicate. The rise was the arithmetic consequence of the subtraction inflating SEs -- confirmed by gamma, whose SE ratio between arms was 1.0037 and whose lambda was likewise unchanged. See LOCO_AS_COVARIATE.md. | **CLOSED - lambda is the wrong instrument** |
| H2 | **Run the blind LOCO diagnostic on the antibody traits.** `jcv_log_all` and `mcv_log_all` are mixture-family, bimodal, and went through the same step-1 workflow with no membership covariate. Regress their LOCO predictors on the Stage 2 posterior membership weights. | **OPEN** |
| H3 | **SUPERSEDED 2026-09-16. Two-pass blind LOCO.** Built and it works (SE 1.0113 / 0.9134), but `--loco-mode covariate` reaches the same place with a one-line Stage 1 change instead of a second REGENIE step 1, and additionally retires the LOCO_SCALE hazard since a projection is invariant to rescaling. Two-pass retains a marginal SE edge and remains on branch `fix/two-pass-loco` if that 1-2% is ever wanted. | **SUPERSEDED by --loco-mode covariate** |
| H4 | **RESOLVED 2026-09-16. GLS is the only arm with an efficiency guarantee.** Zhao & Xiao Thms 8-9 cover the optimally-weighted projection and the OLS arm has no such theorem, but the gap is now quantified rather than open: with the membership column present `col(W)` is Sigma-invariant, so OLS is already efficient and the measured penalty is 0-6% for every shipped W. The guarantee is missing but the quantity it would bound is negligible. GLS remains the wrong trade because it is anti-conservative under Sigma misspecification (0.402-0.818 of its true SE at testosterone's 16% density error) while OLS's point estimate never touches Sigma. | **RESOLVED** |
| H5 | **RESOLVED 2026-09-16. The OLS projection penalty: 57% in simulation, 0.8% on real testosterone.** The cause is neither the distribution family nor the covariance path: it is the **membership column**. The simulation's W had two columns (`beta_1`,`beta_2`); testosterone's shipped W has three (`gamma`,`beta_left`,`beta_right`). At matched separation and matched K, skew and normal W give the same penalty (S=4: 1.593 normal vs 1.568/1.732 skew); deleting `gamma` from the real testosterone W moves it from 1.013 to 2.13/2.51. Mechanism: `Sigma = diag(1/f) B diag(1/f)` collapses onto a single spike at the bimodal trough (98.5% of the spectrum at S=6), and `W_gamma = -(p1(1-p1)/f)(F1-F2)` *is* that spike - it sits 0.011 deg from the leading eigenvector, so including it makes `col(W)` Sigma-invariant and OLS already efficient. Nothing is unclaimed for jcv/mcv: they use the same 3-column builder, predicted penalty 1.0-9.2% (median 1.0%), confirmed by simulation. Diagnostic and evidence: `sim/projeff/src` with aggregates in `sim/projeff/tables`, register `projection-efficiency-20260916-01`. | **RESOLVED - no default change** |
| H6 | **`jackknife + GLS` is unsafe.** 1.4-5x more liberal than the other three arms on the Type I proxy in six of six component/alpha combinations over 85,000+ null variants, and the only arm producing a false positive at 5e-8 (2 of 69). Corroborates the simulated null inflation. Should it be refused rather than merely not-default? | **OPEN** |
| H7 | **`flexmix` at default settings silently returns a degenerate fit at biobank N.** 200/200 replicates in all four simulation cells collapsed to one component; on real data it stops at iteration 2 reporting `converged = TRUE`, with both component SDs equal to the whole-phenotype SD and NaN SEs on a third of variants. Cause: default `tolerance = 1e-6` is relative to a log-likelihood of ~-5e5. Informed initialisation fixes it. Worth a warning in any doc that mentions mixture-model comparisons, since a reader who tries flexmix will get a confident wrong answer. | **OPEN - documentation** |
| H8 | **Stage 1 throughput: 43 SNPs/s on 8 threads, against ~200 SNPs/s previously observed.** The measured run carried `--blocks 1000` and wrote a 66 MB per-chromosome covariance payload that the production `pooled` + `ols` path never reads. Profile the split between BGEN decode, RIF construction (once per chromosome, should not scale with variant count), the per-SNP regression, the block jackknife, and covariance serialisation. Check specifically whether the jackknife uses per-block sufficient statistics (leave-one-out by subtraction, O(B) per tau) or recomputes each block-deleted fit (O(B*T) per SNP) -- the naive form alone would explain a 4-5x gap. Design question that falls out: if pooled is the default, Stage 1 should be able to skip the jackknife on request, which is both the optimisation and a much smaller output. | **OPEN** |
| H9 | **Add a Stage 2 guard warning when `rho_max > 0.6`.** One 45x45 solve. `rho_max` is the largest canonical correlation between in-subspace and complement noise; it predicts the OLS projection penalty at Spearman 0.993 across 38 configurations. Every shipped W sits at 0.03-0.29, but a user who builds a means-only W (no membership column) lands at 0.67-0.99 and silently loses 2-5x efficiency. The guard would catch that. | **OPEN - recommended** |
| H10 | **The membership column is load-bearing and undocumented as such.** Deleting `gamma` from the real testosterone W raises the OLS penalty from 1.013 to 2.13/2.51. Nothing in the weight-builder documentation says that omitting it costs a factor of two. Document it, and consider refusing a mixture W without one. | **OPEN** |
| H11 | **Should `--loco-mode covariate` become the default?** It reproduces no-LOCO point estimates to r = 0.9999, removes the 9.5-14.7% SE inflation the subtraction imposed, raises REGENIE agreement from 0.898/0.817 to 0.970/0.919, and makes `--loco-scale` inert. Blocking question: every REGENIE reference used so far was run `--ignore-pred`, so the ratios measure efficiency against an unadjusted comparator and cannot certify structure control. Decide once the stratified-LOCO comparator lands. | **OPEN - blocked on the fair comparator** |
| H12 | **Is the gamma deconfounding signal real?** Under `--loco-mode covariate` gamma is the only parameter that moves: lambda 1.2650 -> 1.2239, estimates shrinking toward null (0.985) while SEs also shrink (0.987), and the high-z slope against REGENIE moving 1.124 -> 1.090. That signature is the opposite of SE inflation so it cannot be the same artefact, but it is ~3% on one chromosome and one trait. Needs a second chromosome or trait before it is claimed. | **OPEN** |
