# FungWas Stage 2: Parametric Mapping (R Package)

This R package implements "Stage 2" of the FungWas framework: mapping quantile GWAS results (from Stage 1) to parameters of a user-chosen structural model.

**Core idea:** FungWAS is a general framework for Quantile GWAS. You fit **your** null phenotype model externally (normal, mixture, skewed, etc.), then pass the resulting parameters (e.g., posterior weights, residuals, or model-derived weights `W`) into Stage 2. Mixture models are a common use case (e.g., Alfred's serology data), but they are **not** required.

## Key Features
-   **Flexible Mapping:** Map quantile effects to ANY parametric model where you can derive $dQ(\tau) / d\theta$.
-   **Accurate SEs:** Uses the validated pooled Stage 1 covariance sidecar by default (`cov_source = "pooled"`) or the per-SNP delete-one-block jackknife when selected explicitly.
-   **Pre-built Models:** Includes helpers for:
    -   vQTLs (`make_weights_vqtl`)
    -   Normal Mixtures (`make_weights_normal_mixture`,
        `make_weights_normal_mixture_vqtl`)
    -   Anything else, from a user-supplied quantile function
        (`make_weights_generic`)
-   **Mixture weight density (0.5.0):** `density = "hybrid"` divides only the
    membership column by the Stage 1 density (`stage1_kde_density`,
    `read_stage1_density`); the default `"model"` is unchanged. See
    `../docs/MIXTURE_WEIGHT_DENSITY.md`.

## Installation

```r
devtools::install(".")
library(fungwasStage2)
```

## Workflow: From Raw Data to GWAS

### Step 0: The Research Question
Examples of **why** you would use Quantile GWAS:

- **Example A (Standard Model):** "Does a SNP affect the variance of BMI?"
- **Example B (Mixture Model):** "Does a SNP affect the risk of being seropositive in a continuous antibody distribution?"

### Step 1: Fit Your Model (External to FungWAS)
You must first fit your **null phenotype model** in R/Python. The goal is to extract the inputs your model implies, such as **posterior weights** (for mixture models) and **residuals** (for normal or generalized models).

**Example A: Normal model (single-component)**
```r
fit <- lm(y ~ age + sex + batch, data = pheno)
resid <- resid(fit)
# Single-component model has no mixing: weights are all 1.
weights <- rep(1, length(resid))
```

**Example B: Two-component mixture model (mixtools)**

`mixtools` is not a dependency of the package; install it separately for this example.

```r
library(mixtools)

mix <- normalmixEM(y, k = 2)
posterior <- mix$posterior
# Posterior probability of component 1 for each sample
weights <- posterior[, 1]
# Model-implied mean for each sample
mu_hat <- posterior[, 1] * mix$mu[1] + posterior[, 2] * mix$mu[2]
resid <- y - mu_hat
```

You then use these outputs (or other model-derived quantities) to build your **weight matrix** `W` for Stage 2. See `make_weights_*` helpers for common cases, or use `make_weights_generic(...)` for custom models.

**Important:** **Mean-center non-intercept columns of `W` (e.g., multiplicative/variance columns)** to reduce collinearity with the intercept and avoid power loss. This is a re-parameterization (it does not change the model fit), but it can dramatically improve identifiability.

### Step 2: Understanding the “Cutoff”
In a mixture model, there is **no hard cutoff**. The relevant quantity is the **posterior probability** (the weight) computed in Step 1. FungWAS is designed to use these **continuous weights** rather than a hard dichotomy.

If a binary cutoff is required for interpretation, the natural threshold is the point where:

$$P(\text{Class A}) = P(\text{Class B})$$

But for analysis, fungwas prefers the continuous posterior probabilities.

### Step 3: Run Stage 1 & Stage 2
**Stage 1 (Quantile GWAS):**
```bash
fungwas-stage1 \
  --bgen data.bgen \
  --pheno pheno.txt \
  --pheno-col phenotype \
  --covar covar.txt \
  --covar-col age sex batch \
  --out results/chr22
```

**Stage 2 (Parametric mapping):**
```r
W <- make_weights_vqtl(taus, q_tau, mu = mean(y), sd = sd(y))

results <- param_gwas_from_file(
  stage1_file = "results/chr22.stage1.tsv.gz",
  W = W,
  se_mode = "tau_cov",
  sigma0_file = "results/chr22.sigma0.tsv.gz"
)
```

The pooled covariance is the default for `se_mode = "tau_cov"`; it is checked
against the Stage 1 phenotype, covariates, tau grid and sample metadata before
use. When `sigma0_file` is omitted, Stage 2 looks for the conventional sibling
`*.sigma0.tsv.gz` sidecar. To reproduce an analysis that relies on the per-SNP delete-one-block
covariance, select it explicitly:

```r
results <- param_gwas_from_file(
  stage1_file = "results/chr22.stage1.tsv.gz",
  W = W,
  se_mode = "tau_cov",
  cov_source = "jackknife",
  cov_file = "results/chr22.cov.gz"
)
```

## Core Functions

### `param_gwas_from_file()`
Streamlined wrapper to run Stage 2 directly from Stage 1 output files.

```r
results <- param_gwas_from_file(
  stage1_file = "path/to/stage1.tsv.gz",
  W = my_weight_matrix,
  se_mode = "tau_cov",
  sigma0_file = "path/to/stage1.sigma0.tsv.gz"
)
```

If your covariance file predates `.cov.ids.tsv.gz`, you must opt in to unsafe mode:

```r
results <- param_gwas_from_file(
  stage1_file = "path/to/stage1.tsv.gz",
  W = my_weight_matrix,
  se_mode = "tau_cov",
  cov_source = "jackknife",
  cov_file = "path/to/stage1.cov.gz",
  unsafe_skip_cov_ids = TRUE
)
```

Covariance storage compatibility:
- Stage 2 auto-detects and reads both formats:
- Legacy float32 `.cov.gz`
- Compact int8 `.cov.gz` with per-SNP scales in `.cov.scale.gz`

### Parameter covariance (generic, any model)
If you run with `se_mode = "tau_cov"`, Stage 2 can emit per-SNP covariance of the
mapped parameters (for mashr or other downstream tools).

```r
out <- param_gwas_from_file(
  stage1_file = "path/to/stage1.tsv.gz",
  W = my_weight_matrix,
  se_mode = "tau_cov",
  sigma0_file = "path/to/stage1.sigma0.tsv.gz",
  return_cov = "upper",
  param_cov_file = "stage2.paramcov.gz"
)

# Upper-triangle covariance per SNP (P x K*(K+1)/2)
param_cov_upper <- out$param_cov_upper
```

For mashr-style inputs:

```r
fit <- param_gwas(
  stage1_object,
  transform = "custom_W",
  transform_args = list(W = my_weight_matrix),
  se_mode = "tau_cov",
  cov_data = my_cov_data,
  return_cov = "mashr"
)

mashr_inputs <- fit$mashr
# mashr_inputs$Bhat  (P x K), mashr_inputs$Shat (P x K), mashr_inputs$V (K x K x P)
```

### `param_gwas()`
Low-level function working on in-memory objects (useful for simulations).

```r
stage2_fit <- param_gwas(
  stage1_object,
  transform = "custom_W",
  transform_args = list(W = W),
  se_mode = "tau_cov",
  cov_data = cov_data
)
```

Use `se_mode = "diagonal"` only when no tau-covariance object is available and the run is explicitly a fallback or diagnostic check.

### `param_gwas_multi()` (model tournament)
Fit multiple models in one pass and return per-SNP AIC/BIC for model selection.

```r
W_list <- list(
  add = matrix(1, nrow = length(taus), ncol = 1),
  add_mult = cbind(1, seq(-1, 1, length.out = length(taus)))
)

model_scores <- param_gwas_multi(stage1_object, W_list, N = 20000)
head(model_scores)
```

### `param_gwas_multi_from_file()`
Convenience wrapper that reads Stage 1 and (optionally) covariance files.

```r
model_scores <- param_gwas_multi_from_file(
  stage1_file = "path/to/stage1.tsv.gz",
  sigma0_file = "path/to/stage1.sigma0.tsv.gz",
  W_list = W_list,
  N = 20000
)
```

If your covariance file predates `.cov.ids.tsv.gz`, you must opt in to unsafe mode:

```r
model_scores <- param_gwas_multi_from_file(
  stage1_file = "path/to/stage1.tsv.gz",
  cov_file = "path/to/stage1.cov.gz",
  cov_source = "jackknife",
  W_list = W_list,
  N = 20000,
  unsafe_skip_cov_ids = TRUE
)
```

### Weight Builders
Functions to generate the `W` matrix ($T \times K$: one row per quantile, one column per parameter) mapping $T$ quantiles to $K$ parameters.

-   `make_weights_vqtl(taus, q_tau, mu, sd)`
-   `make_weights_normal_mixture(taus, q_tau, p1, mu1, sd1, ...)`
-   `make_weights_generic(...)`

## Outputs and Definitions

The Stage 2 output includes per-SNP parameter estimates and diagnostics. For score/LRT-style inference, the key statistics are:

- **beta_mu:** Estimated effect of the SNP on the mean of the target distribution parameterization.
- **beta_mu_se:** Standard error of `beta_mu`.
- **beta_sigma2:** Estimated effect of the SNP on the variance $\sigma^2$ of the target distribution, per allele, in squared phenotype units (not the log-variance).
- **beta_sigma2_se:** Standard error of `beta_sigma2`.
- **Q:** Heterogeneity (misfit) statistic, computed as a quadratic form of residual quantile slopes, testing whether the SNP affects *any* parameters of the chosen model. Its reference depends on the estimator and covariance source.
- **Q_pval:** Calibrated lack-of-fit p-value when the selected `q_method` has a supported reference. The default OLS path uses a covariance-aware weighted chi-square mixture; the GLS path uses the Hotelling F reference when covariance degrees of freedom are known. An explicitly naive or unsupported reference omits `Q_pval` rather than silently applying $\chi^2_{T-K}$.

## Testing

```bash
R CMD INSTALL stage2
for t in stage2/tests/test-*.R; do Rscript --vanilla "$t" || exit 1; done   # from the repo root
Rscript -e 'testthat::test_dir("stage2/tests/testthat", package = "fungwasStage2", load_package = "installed")'
```

For recovery of known effects end to end, run `docs/examples/1kg_chr22/`. The older
simulation drivers and weight diagnostics (`stage2/scripts/`, `diagnose_*.R`) are at the
tag `archive/research-records-20261016`.

## REGENIE inputs (fast Stage 1)

If you want to run Stage 1 with REGENIE, use `prepare_regenie()` to generate
RIF phenotypes and covariates.

```r
library(fungwasStage2)
library(data.table)

pheno <- fread("inputs/regenie_pheno_shrink_0.7.txt")
covar <- fread("inputs/regenie_covar.txt")

prepare_regenie(
  pheno_df = pheno,
  covar_df = covar,
  phenotypes = c("testosterone_shrink_0.7_male"),
  taus = seq(0.1, 0.9, 0.05),
  covar_cols = setdiff(names(covar), c("FID", "IID")),
  out_inputs = "regenie_inputs",
  out_rtau = "rtau"
)
```
