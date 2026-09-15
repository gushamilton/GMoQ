# Compact validity pilot

This pilot is a bounded, paired check of the shipped G-MoQ engines. It is not
the confirmatory tier-2 grid and its numbers are not release criteria.

The eleven cells are null, mean-only, scale-only, combined mean+scale, and a
covariate stress case at MAF 0.30 and 0.01. Each dataset has six SNPs. The
first SNP is causal in the four direct-effect cells; the remaining five are
null. In `cov_both`, all six are conditionally null. The covariate is

```
C = 0.60 * standardized(G_null) + sqrt(1 - 0.60^2) * U
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

The scientific target is the generative conditional mean coefficient `a` and
the log-SD coefficient `b` in `s(g)=exp(b*g)` (the baseline is `s0=1`). The
shipped `beta_sigma2` basis estimates the local variance coefficient `2*b`,
while the corresponding one-dosage variance ratio is `exp(2*b)`; both are
stored explicitly in each truth file so an SD effect is not silently used as a
variance effect. Both are zero
for the covariate SNP. The projection target is
the population slope of the Stage-1 RIF after the named preprocessing,
projected onto the fitted W. Direct-effect projection targets come from the
existing conditional-CDF truth machinery; covariate targets are computed once
from a separate reference draw and labelled reference-MC. The independent
dataset is the Monte Carlo unit. SNP rows are retained for diagnostics but are
not treated as independent phenotype replicates.

Run from a BP compute node after setting `GMOQ_PYTHON` and the Stage 2 library:

```bash
bash tests/validity_pilot/run_pilot.sh /user/work/fh6520/gmoq_finish_20260915/validity/pilot 10 5000 1000 8
```

The output includes `pilot_meta.json`, one truth JSON per cell, binary Stage 1
artefacts, `pilot_rows.csv`, and `pilot_rows_datasets.csv`. The metadata stores
the source commit, working-tree status, environment, host, seeds and the
target definitions. `summarise_pilot.R` writes `pilot_summary.csv`,
`pilot_bias.csv`, and paired covariate deltas; `report_pilot.R` renders these to
`pilot_report.txt`.
