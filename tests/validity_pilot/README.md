# Compact validity pilot

This pilot is a bounded, paired check of the shipped G-MoQ engines. It is not
the confirmatory tier-2 grid and its numbers are not release criteria.

The pilot includes null, mean-only, scale-only, combined mean+scale, and a
covariate stress case at MAF 0.30 and 0.01. Each dataset has six SNPs. It also
contains a common-MAF Normal variance-linear arm with

```
E[Y | g] = a * (g - 2p),    Var[Y | g] = 1 + v * (g - 2p)
```

at `(a,v) = (0.03,0)`, `(0,0.06)`, `(0.03,0.06)`, and `(0,0.24)`. The
conditional variances are positive for every genotype. The exponential-log-SD
cells remain a separate nonlinear stress. The first SNP is causal in the four
direct-effect cells; the remaining five are null. In `cov_both`, all six are
conditionally null. The covariate is

```
C ~ Normal(0,1),   G_null | C ~ Binomial(2, expit(logit(p) + alpha*C))
Y = 0.50*C + exp(0.45*C)*E
```

so the null SNP is correlated with a covariate that changes both the mean and
the log standard deviation of the phenotype. The production arm uses
`fungwas_stage1.core.run_stage1`, including its RIF-before-covariate order. The
`y_first` comparator uses the shipped `compute_rif`, `residualize`, block-score,
and covariance kernels after residualising the outcome first. Stage 2 fits the
installed package with pooled+OLS and jackknife+OLS, each with an oracle W and
an estimated W.

The `a=b=0.12` direct-effect cells are strong-effect stress contrasts; one
additional common-MAF cell uses `a=b=0.03` as a routine-strength check.

The scientific target for the exponential stress is the generative
conditional mean coefficient `a` and the log-SD coefficient `b` in
`s(g)=exp(b*g)` (the baseline is `s0=1`). The shipped `beta_sigma2` basis
estimates the local variance coefficient `2*b`, while the corresponding
one-dosage variance ratio is `exp(2*b)`. The Normal variance-linear arm's
scientific target is `(a,v)` directly, matching `make_weights_vqtl`'s
variance column. Both parameterisations are stored explicitly in each truth
file so an SD effect is not silently used as a variance effect. Both are zero
for the covariate SNP. The projection target is
the population slope of the Stage-1 RIF after the named preprocessing,
projected onto the fitted W. Direct-effect projection targets come from the
existing conditional-CDF truth machinery or the exact three-component Normal
mixture calculation; covariate targets are computed once from a separate
reference draw and labelled reference-MC. The independent
dataset is the Monte Carlo unit. SNP rows are retained for diagnostics but are
not treated as independent phenotype replicates.

`pilot_covariate_snp1.csv` reports the correlated conditional-null SNP1
separately for `beta_mu` and `beta_sigma2`, for current and `y_first`, each
covariance source and W arm. It includes mean scientific Z, rejection rate,
effect bias, and the package-reported Q statistic with dataset-level MCSEs.

The historical tier report in `tests/report.R` is labelled a regression check;
its old SOUND wording and permissive checks cannot certify scientific validity.
Use this pilot's paired known-truth outputs for that question.

Run from a BP compute node after setting `GMOQ_PYTHON` and the Stage 2 library:

```bash
bash tests/validity_pilot/run_pilot.sh /user/work/fh6520/gmoq_finish_20260915/validity/pilot 10 5000 1000 8
```

The output includes `pilot_meta.json`, one truth JSON per cell, binary Stage 1
artefacts, `pilot_rows.csv`, and `pilot_rows_datasets.csv`. The metadata stores
the source commit, working-tree status, environment, host, seeds and the
target definitions. `summarise_pilot.R` writes `pilot_summary.csv`,
`pilot_bias.csv`, SNP1 covariate summaries and paired covariate deltas;
`report_pilot.R` renders these to `pilot_report.txt` with dataset-level MCSEs.
