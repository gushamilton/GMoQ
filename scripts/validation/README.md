# Bounded input and real-chr22 validation

These scripts exercise the installed Stage 1 CLI and the Stage 2 OLS + pooled
covariance path on compute nodes. They do not copy participant-level inputs into
the repository. The cluster outputs are kept under
`/user/work/fh6520/gmoq_finish_20260915/inputs/real_chr22_smoke/`.

`real_chr22_smoke.py` selects 24 common, called HM3 chr22 SNPs, generates a
known-truth phenotype from the first SNP's allele-0 dosage, and runs the actual
CLI. `real_chr22_smoke.R` fits the vQTL weights with pooled and jackknife
covariance. The known-truth run recovered beta_mu = 0.250635 for an expected
0.25; the pooled and jackknife runs both had the correct sign. Non-causal SNPs
are descriptive because LD prevents a null assertion.

`run_real_chr22_empirical.sh` uses the paper's exact empirical definition:
`log_bmi = log(bmi)`, covariates `age, sex, PC1..PC10`, complete-case
`keep_log_bmi`, and the phenotype-specific REGENIE `pred.list` with
`--loco-scale regression`. `real_chr22_empirical_stage2.R` fits the paper's
additive or mean-plus-variance weights with pooled covariance. The empirical
outputs are descriptive and carry no known-effect truth claim.

`time_bgen_decode.py` compares the native 8-bit dosage path with the checked
per-record preamble path. `compare_repeat.py` compares decompressed output
content across repeated runs, avoiding gzip timestamp differences.
