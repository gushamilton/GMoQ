# Stage 2: map quantile effects onto model parameters (`fungwasStage2`)

Stage 1 gives each SNP's effect on 45 quantiles of the phenotype. Stage 2 asks which change in
the parameters of a model you choose — a mean and a variance, a two-component mixture, or your
own — best explains those 45 effects. It fits that by least squares through a weight matrix `W`
(one row per quantile, one column per parameter: `W[t, k]` is how quantile `t` moves when
parameter `k` moves), with standard errors from the Stage 1 covariance.

## Install

```bash
R CMD INSTALL stage2
```

## Use

1. Read the phenotype values Stage 1 analysed (the samples it kept) and their quantiles on the
   Stage 1 grid.
2. Build `W` for your model.
3. Map.

```r
library(fungwasStage2)
y     <- data.table::fread("pheno.txt")$bmi              # restrict to the samples Stage 1 kept
taus  <- round(seq(0.02, 0.98, length.out = 45), 3)      # Stage 1's default grid
q_tau <- as.numeric(quantile(y, taus, type = 8))

W <- make_weights_vqtl(taus, q_tau, mu = mean(y), sd = sd(y))

res <- param_gwas_from_file(
  stage1_file = "results/chr22.stage1.tsv.gz",
  W           = W,
  sigma0_file = "results/chr22.sigma0.tsv.gz"   # pooled covariance (the default source)
)
```

For the per-SNP jackknife covariance instead, pass `cov_source = "jackknife", cov_file =
"results/chr22.cov.gz"`. The estimates are identical; the standard errors and the Q test differ.
[../docs/BACKGROUND.md](../docs/BACKGROUND.md) says when to prefer it.

### Weight builders

| function | parameters (columns of `W`) |
|---|---|
| `make_weights_vqtl(taus, q_tau, mu, sd)` | `beta_mu` (mean), `beta_sigma2` (variance, per allele, in squared phenotype units) |
| `make_weights_normal_mixture(taus, q_tau, p1, mu1, sd1, mu2, sd2, include_membership = TRUE, stage1_density = ...)` | `gamma` (log-odds of belonging to component 1; only with `include_membership = TRUE`), `beta_1`, `beta_2` (component means) |
| `make_weights_normal_mixture_vqtl(...)` | `gamma`, `beta_mu1`, `beta_mu2`, `beta_sigma1`, `beta_sigma2` |
| `make_weights_generic(taus, q_tau, dist_cdf, dist_pdf, params, grad_funcs)` | any family: one column per gradient you supply |

The mixture parameters (`p1`, `mu1`, `sd1`, `mu2`, `sd2`) describe the baseline distribution of
`y`, which you fit first. A W with a membership column needs the Stage 1 density:

```r
fit <- mixtools::normalmixEM(y, k = 2)     # or any two-component fit you trust
W <- make_weights_normal_mixture(
  taus, q_tau, p1 = fit$lambda[1],
  mu1 = fit$mu[1], sd1 = fit$sigma[1], mu2 = fit$mu[2], sd2 = fit$sigma[2],
  include_membership = TRUE,
  stage1_density = stage1_kde_density(y, q_tau))
```

`mixtools` is not a dependency; install it separately. More:
[../docs/MIXTURE_WEIGHT_DENSITY.md](../docs/MIXTURE_WEIGHT_DENSITY.md).

### Output

One row per SNP: `snp_id chr bp effect_allele other_allele`; for each parameter `<param>`,
`<param>_se`, `<param>_p` and `<param>_df`; then `Q` and `Q_pval`, a test of whether the model's
parameters fully explain the SNP's quantile profile; then diagnostic columns.

### Other functions

`param_gwas()` (in-memory data), `param_gwas_multi()` (compare several `W` by AIC/BIC),
`functional_grad()` (SEs for functions of the parameters) and `prepare_regenie()`. Each has a help
page, e.g. `?param_gwas_from_file`.

## Tests

```bash
for t in stage2/tests/test-*.R; do Rscript --vanilla "$t" || exit 1; done   # from the repo root
Rscript -e 'testthat::test_dir("stage2/tests/testthat", package = "fungwasStage2", load_package = "installed")'
```
