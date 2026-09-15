# G-MoQ: Genomic Method of Quantiles

*(the method and this codebase were formerly called "FungWAS"; the Python package is still
`fungwas_stage1` and the R package `fungwasStage2`, and those names are used throughout below.)*

A framework for performing Genome-Wide Association Studies (GWAS) on distributional parameters,
enabling discovery beyond standard mean effects.

**G-MoQ is a two-stage pipeline.** Stage 1 (Python/C++) runs a quantile GWAS and emits per-SNP
effects at a grid of quantiles together with the full covariance of those effects across quantiles.
Stage 2 (R/C++) maps that quantile profile onto the parameters of a distributional model you
choose, with standard errors that respect the Stage 1 covariance. The two stages are separate
packages in one repository and are installed separately.

The primary workflow is **OLS with pooled covariance**. Stage 1 writes the pooled
`Sigma_0` sidecar by default; Stage 2 scales it by each SNP's residual genotype
variation. Per-SNP jackknife covariance remains available as an explicit comparison.
Pooled covariance is a working approximation: validation must include true variance
effects and covariate structure, not just null phenotypes.
See [what the primary analysis estimates](docs/ESTIMAND.md) for effect units,
the local model interpretation and the scope of validation.

**FungWas** allows you to test SNP associations with:
- **Variance** (vQTLs)
- **Mixture components** (e.g., latent subtypes in a disease)
- **Skewness / Kurtosis**
- **Any parametric system** definable by a map from quantiles.

## How It Works

FungWas operates in two stages:

See [Stage 2 Documentation](stage2/README.md) for full details.

1.  **Stage 1 (Python/C++):**
    -   Performs **Recentered Influence Function (RIF) regression** on the phenotype at multiple quantile levels (taus).
    -   This is a highly optimized, parallelized quantile GWAS.
    -   **Output:** Beta estimates and standard errors for each tau, plus the full covariance matrix of estimates across taus for each SNP.
    -   **Key Tech:** Python CLI, `bgen-reader` for fast I/O, C++ kernel with OpenMP.

2.  **Stage 2 (R/C++):**
    -   Maps the quantile effects from Stage 1 to parameters of interest (e.g., mean/variance of a mixture model).
    -   Uses the **Delta Method** or efficient **GMM-like** estimators.
    -   Input: Stage 1 results + a "Weight Matrix" describing the derivative of quantiles with respect to parameters.
    -   **Output:** GWAS summary statistics (Beta, SE, P-value) for your custom parameters.

## Installation

### 1. Create Environment (Recommended)

We provide a Conda environment that handles Python dependencies and system libraries (Armadillo, OpenMP).

```bash
# Create and activate the environment
conda env create -f stage1/environment.yml
conda activate fungwas-stage1
```

### 2. Install Stage 1 (Python)

```bash
cd stage1
pip install .
```
*Note: This builds the optimized C++ extension.*

### 3. Install Stage 2 (R)

```r
# From R (ensure devtools is installed)
devtools::install("stage2")
```

## Detailed Build Instructions

### System Requirements
- C/C++ compiler with OpenMP (e.g. GCC >= 9).
- R >= 4.2 with `Rcpp` and `RcppArmadillo` available.
- Python >= 3.10 with `pip` and `numpy`.

### Build Stage 1 (Python/C++)

```bash
conda activate fungwas-stage1
cd stage1
pip install -e .
python -c "from fungwas_stage1 import core; print(core.HAVE_CPP)"
```

If `core.HAVE_CPP` is `False`, ensure your compiler and OpenMP runtime are available,
then reinstall with `pip install -e .`.

### Build Stage 2 (R/C++)

```bash
R CMD INSTALL stage2
```

If you install into a user library, set:

```bash
export R_LIBS_USER=/path/to/your/R/library
```

### Example simulation helper

The Stage 1 simulation helper used by the Stage 2 proof script lives in
`docs/examples/generate_stage1_sim.py`. It writes `.stage1.tsv.gz` and `.cov.gz`
files compatible with `param_gwas_from_file`.

## Testing

Both stages ship their own checks. Tests should use the installed package and must not
regenerate package documentation or modify `NAMESPACE`.

### Stage 1 (pytest)

```bash
conda activate fungwas-stage1
cd stage1 && python setup.py build_ext --inplace && cd ..
PYTHONPATH=stage1 python -m pytest stage1/tests/ -q
```

### Stage 2 (R proof scripts)

`stage2/tests/` holds standalone R proof scripts (not `testthat` files); each is run with
`Rscript` and exits non-zero on failure.

```bash
R CMD INSTALL -l "$R_LIBS_USER" stage2/
cd stage2/tests
for t in test-*.R; do Rscript --vanilla "$t" || exit 1; done
```

Some scripts are large simulation demonstrations or require optional dependencies/data;
the CI workflow lists the fast assertion-based regression checks separately.
`test-plink-files.R` requires `GMOQ_PLINK_PREFIX`; `test-qstat.R` needs `fitdistrplus`.

### End-to-end check (synthetic data only)

`docs/examples/generate_stage1_sim.py` and `stage2/scripts/simulate_full_pipeline_stage1.py`
simulate data, run Stage 1, and write Stage 1 outputs that Stage 2 can consume. For a check that
also exercises the real BGEN path, write a small synthetic BGEN with `bgen.BgenWriter`, run the
`fungwas-stage1` CLI over it, and feed the resulting `.stage1.tsv.gz` plus the covariance quartet
(`.cov.gz`, `.cov.scale.gz`, `.cov.ids.tsv.gz`, `.cov.meta.json`) into `param_gwas_from_file()`
with `se_mode = "tau_cov"`.

**Never run these inside a directory holding, or writing to, individual-level data.** Write all
run outputs outside this repository.

## Workflow Guide

### Step 0: Data Preparation
-   **Genotypes:** BGEN format (indexed with `bgenix`).
-   **Phenotypes:** Standard space-delimited text file (FID IID Pheno).
-   **Covariates:** Standard space-delimited text file (FID IID Age Sex ...).

See [Phenotype Preparation Guide](docs/phenotype_preparation.md) for details.

### Allele Handling
- Stage 1 outputs `effect_allele`/`other_allele` and uses BGEN allele1 dosage as the effect allele.
- Non-biallelic variants are dropped during Stage 1 and the CLI reports the count.

### Step 1: Run Quantile GWAS

Use the `fungwas-stage1` CLI to run the genome-wide scan.

```bash
# Example: Chromosome 22
fungwas-stage1 \
    --bgen data/chr22.bgen \
    --pheno inputs/phenotypes.txt \
    --pheno-col testosterone \
    --covar inputs/covariates.txt \
    --out results/chr22 \
    --threads 16
```
*Produces `results/chr22.stage1.tsv.gz` (estimates), per-SNP covariance files, and
`results/chr22.sigma0.tsv.gz` with metadata (pooled covariance).*
By default, Stage 1 writes covariance in compact int8 format with sidecar
`results/chr22.cov.scale.gz` (dequantization scales).

The default grid is 45 quantiles from 0.02 to 0.98, with 1,000 jackknife blocks.
When supplying REGENIE LOCO predictions, Stage 1 uses `--loco-scale regression`
by default to convert predictions to the phenotype's scale before subtraction.

See [Stage 1 Documentation](stage1/README.md) for full CLI options.

### Optional: Use REGENIE for ultra-fast Stage 1

If you want speed over detailed jackknife SEs, you can run Stage 1 via
REGENIE by transforming phenotypes into RIFs and running a standard GWAS.
The `prepare_regenie()` helper in the Stage 2 R package writes the required
phenotype and covariate files.

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

Then run REGENIE step 2 with each `RIF_tau*.txt` phenotype column. The
resulting per-tau betas/SEs can be combined for Stage 2 mapping using
`param_gwas_from_file()`.

### Step 2: Parametric Mapping

In R, load the results and map them to your model of interest.

See [Stage 2 Documentation](stage2/README.md) for full details.

```r
library(fungwasStage2)

# 1. Define a Normal location/variance model on the Stage 1 grid.
#    This example assumes a baseline SD of 1; use the analysed phenotype's
#    baseline SD for your model (after the same preprocessing as Stage 1).
taus <- round(seq(0.02, 0.98, length.out = 45), 3)
# This helper builds the weight matrix W
W <- make_weights_vqtl(taus, qnorm(taus), mu = 0, sd = 1)

# 2. Run Stage 2
results <- param_gwas_from_file(
  stage1_file = "results/chr22.stage1.tsv.gz",
  W = W,
  se_mode = "tau_cov",
  tau_cov_estimator = "ols",
  cov_source = "pooled",
  sigma0_file = "results/chr22.sigma0.tsv.gz"
)

# 3. Save or analyze
head(results) # Contains betas and SEs for 'beta_mu' and 'beta_sigma2'
```

### HPC Submission
See [SLURM Templates](docs/slurm_templates/) for example submission scripts.

## Advanced Topics

### Standard Error Modes
Stage 2 supports different methods for SE calculation (`param_gwas` argument `se_mode`):
-   `tau_cov`: **Recommended**. Uses full cross-quantile covariance. With `cov_source="pooled"`, this is `Sigma_0 / d_resid`; with `cov_source="jackknife"`, it is the estimated per-SNP covariance. The latter supports float32 and int8 storage. These are covariance estimates, not exact sampling covariances.
-   `plugin_cor`: Uses a single correlation matrix for all SNPs (from null distribution). Faster, good for initial scans.
-   `diagonal`: Assumes independence between quantiles. Least accurate, not recommended for final results.

### Proof & Verification
The repository includes proof scripts demonstrating the statistical correctness of the methods:
-   `stage1/tests/test_jk_equivalence.py`: Proves Stage 1 jackknife SEs match analytic theory.
-   `stage2/tests/test-vQTL.R`: Proves vQTL recovery on simulated data.

## License
MIT License
