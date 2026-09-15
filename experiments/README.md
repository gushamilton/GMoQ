# Experiment register

Check this register **before launching an experiment**. A completed experiment
may have found a failure; “completed” does not mean “passed”. Reuse its evidence
unless a new source change, scientific question or independent confirmation
justifies a new run.

Small aggregate results, scripts, settings and conclusions belong in Git.
Participant-level inputs and large binary payloads stay on BluePebble. The
existing private repository is `gushamilton/GMoQ`; development work is on
`work/complete-gmoq`.

## Runs

| Stable ID | State | Question and settings | Source and preserved evidence | BP output |
|---|---|---|---|---|
| `testosterone-gap-20260915-01` | Registered; diagnostic pending | Bounded reproduction of the testosterone two-density interface: 60 paired datasets, N=10,000, separation 8, upper-membership gamma .08, MAF .30, 45 taus, seed 20260915. Compare production Stage 1 KDE density with historical fitted-mixture density and a same-density repair. | `tests/validity_pilot/testosterone_gap_diagnostic.py`; this README. The diagnostic is oracle-mixture evidence and does not claim to replay UK Biobank. | Reserved: `gmoq_finish_20260915/validity/testosterone-gap-20260915-01` |
| `final-package-check-20260915-01` | Completed | Final archive built with vignettes, checked and installed into a fresh library; focused multi-file alignment regression and existing pooled/ID checks passed. | Source `ccc2d97`; [commands, logs and outcome](final-package-check-20260915-01/). Zero errors/notes, one existing GNU Makevars warning. Allocation 18966618 on bp1-compute019 released. | `gmoq_finish_20260915/root/final-package-check-20260915-01` |
| `nuisance-bootstrap-20260915-01` | Deferred; not run | Proposed N=10,000, 100 independent datasets, MAF .30, variance effects 0/.24, seed 20260916; Gaussian versus fourth-order KDE with fixed oracle Normal OLS weights. Compare 199 full-recomputation bootstrap resamples against fixed-RIF JK and pooled SEs. | User prioritised the fast pooled workflow and proportionate validation. Preserve this design for a future non-null uncertainty investigation; it is not a release gate. Freeze script/configuration before any future submission. | Reserved: `gmoq_finish_20260915/inference/nuisance-bootstrap-20260915-01` |
| `density-repair-20260915-01` | Completed; bias improved, uncertainty unresolved | Paired production Gaussian KDE versus fourth-order Gaussian kernel at the same bandwidth; fresh seed 20260915, N=10,000, 100 datasets, variance effects 0/.06/.24. | Test-only source `28c5316`; [configuration](../tests/validity_pilot/density_repair_config.json), [results](../tests/validity_pilot/evidence/DENSITY_REPAIR.md). Allocation 18966618, bp1-compute019. Strong-cell paired shift +.003761 ± .000101 removes measured projection bias; pooled coverage .91 ± .029 and SE/empirical SD .845. All 2,400 Stage 2 fits succeeded; no repaired density floors. No production repair accepted. | `gmoq_finish_20260915/validity/density-repair-20260915-01` |
| `pooled-memory-repair-20260915-01` | Completed | Shared chromosome covariance with per-variant scales. Dense/compact agreement passed for estimates, SEs, p-values, Q and returned covariances, including multiple models/chromosomes and invalid inputs. Ten installed-package regression scripts also passed. | Source `e1bc65b`; [fixture, manifest and memory profile](../stage2/tests/pooled-memory-repair-20260915-01/). At 45 taus/10,000 variants, covariance storage decreased from 82,880,096 to 168,424 bytes (~492-fold). Construction timing is a single bounded measurement, not genome-wide throughput. | `gmoq_finish_20260915/inputs/pooled-memory-repair-20260915-01` |
| `allele-alignment-20260915-01` | Completed | Do saved chr22 log-BMI/REGENIE effect alleles agree? 17,550 matched variants. | [Result, inputs and command](../docs/ALLELE_CONVENTION.md); script introduced at `39a3ebc`. Same pairs throughout; correlation .99285. | `gmoq_finish_20260915/root/out/allele_alignment.json` |
| `location-scale-pilot-20260915-01` | Completed; superseded | Initial 100-dataset, N=10,000 location/scale/covariate pilot. | Source `f1b4a2c`; superseded by the expanded run below. Do not rerun for missing summaries. | `gmoq_finish_20260915/validity/pilot` |
| `location-scale-pilot-20260915-02` | Completed; limitations found | Expanded 15-cell pilot, 100 datasets, N=10,000, 45 taus, 1,000 blocks; pooled/JK, oracle/estimated W. | Source `48171e3`; [metadata, aggregate results, per-cell provenance and commands](../tests/validity_pilot/evidence/VALIDITY_EVIDENCE.md). Strong variance projection bias detected. | `gmoq_finish_20260915/validity/pilot_linear` |
| `covariate-basis-20260915-01` | Completed; limitations found | Paired linear C versus C,C²,C³; 100 datasets, N=10,000, common/low-frequency correlated null SNP. | Source `be70f0a`; [summary and provenance](../tests/validity_pilot/evidence/VALIDITY_EVIDENCE.md). Rich basis removes most bias; pooled low-frequency inference still inflated. | `gmoq_finish_20260915/validity/cov_probe` |
| `biobank-covariance-20260915-01` | Completed; limitations found | Reanalyse six completed cells, 15,000 datasets/cell, N=461,320, 45 taus, 1,000 blocks; exact OLS projection target. | Current Stage 2 `b1be4d7`; [frozen truth, selected chunks, source compatibility, commands and aggregates](../stage2/tests/validity_biobank_check/README.md). Null/location calibrated; low-frequency dispersion pooled coverage .9285 versus JK .9447. | `gmoq_finish_20260915/validity/biobank_check/{out,out_compare}`; original raw data `covsource_run/raw` |
| `stage2-gzip-reader-20260915-01` | Completed | Ensure Stage 2 reads gzipped Stage 1 and covariance-ID sidecars on BP nodes without R.utils. | Package commit `92fef3a` streams `.gz` text through `gzip -dc` for Stage 1 and covariance-ID reads; matched full chr22 pooled/jackknife benchmark passed on BP compute053 (Slurm 18970068), with no estimator change. | Paper evidence: `Gmoq_paper/outputs/stage2-benchmark-20260915-05` |
| `real-chr22-smoke-20260915-01` | Completed | Actual CLI, 24 real variants, N=463,005, known mean effect .25, 45 taus, 65 blocks; independent rerun and decode timing. | Input source `d231e27`; [scripts and result](../scripts/validation/README.md), [timings and interpretation](../docs/VALIDATION_2026-09-15.md). Mean recovered .250635; repeated decompressed files identical. This is not replicated coverage evidence. | `gmoq_finish_20260915/inputs/real_chr22_smoke` |
| `empirical-logbmi-20260915-01` | Completed | Fresh 24-variant empirical run, N=461,320, age/sex/10 PCs, regression LOCO scaling; matched Normal mean/variance weights. | Source `ba82da8`; [exact preprocessing and aggregate result](../scripts/validation/real_chr22_empirical_provenance.md). Earlier cached raw-BMI weight fit excluded. | `gmoq_finish_20260915/inputs/real_chr22_smoke/empirical/log_bmi_matched_both_pooled.tsv` |
| `pooled-memory-20260915-01` | Completed | Measure existing covariance expansion at 45 taus. | Installed candidate `root/lib`, R 4.5.1; [script and output](pooled-memory-20260915-01/). 1,000 variants = 8,280,048 covariance bytes/.037 s; 10,000 = 82,800,048 bytes/.443 s. This measures allocation, not genome-wide throughput. | `gmoq_finish_20260915/root/profile_pooled_memory.R` |
| `hc0-covariate-20260915-01` | Completed | Independent HC0 covariance pilot, 100 datasets, N=10,000; OLS estimates unchanged. | [Exact executed script snapshot, manifest and results](../stage2/tests/robust_probe/README.md), commit `490857e`. Fixed Normal cutoff covariance diagnostic, not the full shipped t-reference calculation. | `gmoq_finish_20260915/inference/robust_probe/reps100` |
| `hc0-covariate-20260915-02` | Completed; residual limitation | Fresh 1,000-dataset confirmation, seed 19152026, N=10,000. | [Manifest and aggregates](../stage2/tests/robust_probe/evidence/hc0-covariate-20260915-02/). Correlated low-frequency variance-null rejection: pooled .268, HC0 .070, JK .067 with fixed Normal cutoff. Substantial improvement, not proof of nominal calibration. | `gmoq_finish_20260915/inference/robust_probe/hc0-covariate-20260915-02` |
| `density-bias-20260915-01` | Completed; cause identified | Paired empirical quantiles/KDE versus exact mixture density/quantiles; N=10,000, 100 datasets, v=0/.06/.24. | Source `47ec83d`; [commands, hashes and results](../tests/validity_pilot/evidence/DENSITY_DIAGNOSIS.md). At v=.24, KDE versus exact density changes the variance estimate by -.003989 ± .000324. | `gmoq_finish_20260915/validity/density-bias-20260915-01` |
| `density-bias-20260915-02` | Completed; cause confirmed | Larger-N paired diagnostic, N=100,000, 50 datasets, v=.24. | Source `eadf3fa`; [evidence](../tests/validity_pilot/evidence/DENSITY_DIAGNOSIS.md). KDE contrast -.001834 ± .000146; exact-density/population-quantile arm has negligible measured projection bias. No production correction implemented. | `gmoq_finish_20260915/validity/density-bias-20260915-02` |

BP paths above are relative to `/user/work/fh6520/`. Legacy directories retain
their existing names; the stable IDs index them without moving or overwriting
them. `validity/pilot_smoke` was a one-replicate execution check, not validity
evidence. `validity/pilot_final` was an earlier Stage-1-only product; use the
completed `pilot_linear` result instead.

The unexecuted nuisance-bootstrap prototype is preserved at source `8a2cdfd`
and tag `experiments/2026-09-15/nuisance-bootstrap-deferred-source`; it is not
merged into the package branch. Review its output guard against the exclusive
directory rule below before any future execution. Preparation is not evidence
that the prototype has been tested.

The independent `covsource_run`/`cs1`/`cs2` campaign has existing jobs and outputs.
Inspect its live job handles and completed chunks before scheduling anything
similar; this development work does not own or restart those jobs.

## Earlier investigations to consult before proposing new runs

These are inherited reports, not newly reproduced results. Preserve that
distinction when using their conclusions or deciding whether a confirmation
is needed.

- [Empirical model weights](../docs/EMPIRICAL_W.md): exact-quadrature comparison
  of shipped versus empirical dispersion columns; working files at
  `/user/work/fh6520/empw2/`. It reports poor general calibration from simply
  substituting empirical quantiles into W. It explicitly lacks a Monte Carlo
  reconciliation through the shipped Stage 2 pipeline.
- [Density bandwidth investigation](../docs/KNOWN_ISSUES.md#8-kernel-density-bandwidth-as-of-2026-09-14)
  and [lumpy phenotypes](../docs/BP2_lumpy_phenotypes.md): prior bandwidth and
  discrete/mixture phenotype investigations. Their Q/model-selection findings
  do not establish unbiased non-null parameter estimation or interval coverage.
  Consult these before repeating a bandwidth sweep or assuming a smaller
  bandwidth fixes tied outcomes.

## Rules for new runs

1. Assign a new stable ID and register the scientific question before submission.
   A changed seed, configuration, source, or confirmation dataset gets a new ID.
2. Commit the script/configuration first. Record source commit (and any dirty
   patch/hash), seed, dataset replication unit, input identifiers, environment,
   tau grid, sample size, blocks, threads, host and Slurm job ID.
3. Create the output directory exclusively; refuse an existing directory. Never
   overwrite a completed run. Resume only a verified incomplete run with the same
   configuration and a checked manifest; record the resumed job handle.
4. Commit aggregate results and the conclusion, including failures and Monte
   Carlo uncertainty. Record exact commands and large-output locations. Preserve
   source commits referenced by the evidence with experiment-source Git tags.
5. Mark the run completed only after checking its outputs. Reanalysis of existing
   raw data should get its own report/run ID and retain the original raw data.

All BluePebble execution, including small validation and manifest-check scripts,
belongs on a compute node. Reuse at most one held allocation and release it when
finished. Git and file inspection can use the login node.
