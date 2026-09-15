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
| `density-repair-20260915-01` | Dispatched | Paired production Gaussian KDE versus fourth-order Gaussian kernel at the same bandwidth; fresh seed 20260915, N=10,000, 100 datasets, variance effects 0/.06/.24. Check bias, shipped-reference coverage, SE calibration, density failures and runtime. | Test-only source `28c5316`; [configuration](../tests/validity_pilot/density_repair_config.json). Allocation 18966618 on bp1-compute019. Density nuisance uncertainty remains a separate question; improved bias alone is not a passing result. | `gmoq_finish_20260915/validity/density-repair-20260915-01` |
| `pooled-memory-repair-20260915-01` | In preparation; not run | Replace repeated pooled covariance storage with shared chromosome covariance and per-variant scales. Require dense/compact agreement for estimates, SEs, p-values, Q and returned covariances, including multiple models/chromosomes and invalid inputs; measure memory. | Production mathematics must remain unchanged. Freeze implementation and validation source before submission. | Reserved: `gmoq_finish_20260915/root/pooled-memory-repair-20260915-01` |
| `allele-alignment-20260915-01` | Completed | Do saved chr22 log-BMI/REGENIE effect alleles agree? 17,550 matched variants. | [Result, inputs and command](../docs/ALLELE_CONVENTION.md); script introduced at `39a3ebc`. Same pairs throughout; correlation .99285. | `gmoq_finish_20260915/root/out/allele_alignment.json` |
| `location-scale-pilot-20260915-01` | Completed; superseded | Initial 100-dataset, N=10,000 location/scale/covariate pilot. | Source `f1b4a2c`; superseded by the expanded run below. Do not rerun for missing summaries. | `gmoq_finish_20260915/validity/pilot` |
| `location-scale-pilot-20260915-02` | Completed; limitations found | Expanded 15-cell pilot, 100 datasets, N=10,000, 45 taus, 1,000 blocks; pooled/JK, oracle/estimated W. | Source `48171e3`; [metadata, aggregate results, per-cell provenance and commands](../tests/validity_pilot/evidence/VALIDITY_EVIDENCE.md). Strong variance projection bias detected. | `gmoq_finish_20260915/validity/pilot_linear` |
| `covariate-basis-20260915-01` | Completed; limitations found | Paired linear C versus C,C²,C³; 100 datasets, N=10,000, common/low-frequency correlated null SNP. | Source `be70f0a`; [summary and provenance](../tests/validity_pilot/evidence/VALIDITY_EVIDENCE.md). Rich basis removes most bias; pooled low-frequency inference still inflated. | `gmoq_finish_20260915/validity/cov_probe` |
| `biobank-covariance-20260915-01` | Completed; limitations found | Reanalyse six completed cells, 15,000 datasets/cell, N=461,320, 45 taus, 1,000 blocks; exact OLS projection target. | Current Stage 2 `b1be4d7`; [frozen truth, selected chunks, source compatibility, commands and aggregates](../stage2/tests/validity_biobank_check/README.md). Null/location calibrated; low-frequency dispersion pooled coverage .9285 versus JK .9447. | `gmoq_finish_20260915/validity/biobank_check/{out,out_compare}`; original raw data `covsource_run/raw` |
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
