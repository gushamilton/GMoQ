# Independent HC0 covariance probe

This directory records a bounded diagnostic of the rich covariate-basis
conditional-null failure. It does not change Stage 1 or Stage 2 production
code. `hc0_covariate_probe.py` regenerates the exact `cov_both_common` and
`cov_both_rare` cells from the paired validity pilot, using the pilot's
dataset seeds, rich basis `C, C^2, C^3`, `N=10,000`, 45 taus and 1,000 blocks.
The fixed two-parameter OLS design is

```text
W(tau) = [1, qnorm(tau) / 2]
```

For each SNP it forms the independent sandwich requested for the diagnostic:

```text
score_i = g_resid_i * (RIF_resid_i - g_resid_i * beta)
V_HC0  = score' score / d_resid^2
```

The first SNP is correlated with `C` and remains conditionally null. SNPs 2--6
are independent marginal-null controls. Summaries treat the dataset as the
Monte Carlo unit; the SNP-group summary for SNPs 2--6 averages within dataset
before calculating MCSE, while `empirical_sd_row` reports the individual-SNP
spread for SE-ratio interpretation.

The immutable run records are:

* `hc0-covariate-20260915-01`: 100 datasets, seed `9142026`. The completed
  remote output was deliberately retained at
  `/user/work/fh6520/gmoq_finish_20260915/inference/robust_probe/reps100`;
  its exact pre-guard script snapshot is
  `hc0_covariate_probe_pilot_exec.py`.
* `hc0-covariate-20260915-02`: 1,000 fresh datasets, seed `19152026`, output
  `/user/work/fh6520/gmoq_finish_20260915/inference/robust_probe/hc0-covariate-20260915-02`.

Both runs used bp1-compute053, two OpenMP/BLAS threads, the Stage 1 snapshot
at source commit `be70f0a`, and the covariate-probe source hash recorded in
each `meta.json`. The aggregate `summary.csv` and `timing.csv` files are
checked in; raw per-dataset rows remain on BP and are not needed to reproduce
the reported aggregate evidence.

The 1,000-dataset production-order (`current`) result for the correlated rare
SNP is the key diagnostic. For `beta_sigma2`, HC0 has rejection `0.070` and
coverage `0.930`, JK has `0.067` and `0.933`, while pooled has `0.268` and
`0.732`; mean SE / individual-SNP empirical SD is `0.956`, `0.956`, and
`0.581`, respectively. For `beta_mu`, HC0/JK rejection is `0.049`/`0.050`
versus pooled `0.169`. The common-MAF correlated SNP is near nominal for HC0
and JK, with pooled modestly low on beta_mu and beta_sigma2. Independent
SNPs 2--6 remain near nominal for all sources, so the defect is tied to the
genotype-correlated heteroscedastic null rather than a universal variance
inflation.

On this small probe, HC0 formation averaged `0.0145` seconds per arm, versus
`0.0256` seconds for the current Stage 1 call and `0.0145` versus `0.0212`
seconds for the y-first comparator. This is about 57--69% extra CPU time in
the six-SNP N=10,000 diagnostic; full-scan storage and compute would require a
separate design review.
