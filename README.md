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

**What a G-MoQ effect means.** Stage 1 is a *recentered influence function* (RIF)
regression, so a coefficient at quantile `tau` is an **unconditional** quantile
effect: how the `tau`-th quantile of Y in the population would move under a
marginal shift in the genotype distribution, holding the covariate basis fixed.
It is *not* the conditional quantile effect that classical quantile regression
estimates, and the two do not coincide. Stage 2 then reports the parameters of
your chosen distributional model that best explain the whole quantile profile.
Two consequences bite in practice: anything you subtract from Y before Stage 1
(a LOCO prediction, a transformation) **redefines the estimand** rather than
removing a nuisance, which is why `--loco-mode covariate` is the default; and
the estimand is defined against the unconditional distribution of the analysed
phenotype, so the analysed phenotype must be recorded with the result. See
[`docs/ESTIMAND.md`](docs/ESTIMAND.md).

The primary workflow is **OLS with pooled covariance** (`cov_source = "pooled"`,
`tau_cov_estimator = "ols"`). Both defaults are measured, not assumed:

- Under OLS the point estimate does not touch the covariance, so `pooled` and
  `jackknife` give **bit-identical** estimates and differ only in the reported
  SE, by under 1% on real chr22 testosterone. The choice is therefore purely
  about which SE is better calibrated.
- `pooled` is the default: it reads only the `Sigma_0` sidecar rather than the
  per-SNP payload, and on real data the two are indistinguishable -- under 1%
  apart on chr22 testosterone, with correlations identical to four decimals in
  all 24 z-strata.
- **Where they differ, and it is narrower than it is small.** Across 780,000
  simulated scans, pooled is anti-conservative under a TRUE DISPERSION effect at
  LOW frequency. On the 3-column mixture `W` shipped here the defect is
  CONCENTRATED rather than avoided: `gamma` is protected (SE ratio 1.0136 at
  MAF 0.005 / Z=30) because it takes 99.96% of its variance on the tau-covariance's
  leading eigendirection, but `beta_left`/`beta_right` sit at exactly 90 degrees
  to that direction with zero share and are inflated **3510x at MAF 0.005 / Z=30
  and 236x at MAF 0.01 / Z=30**. Under the null the two sources agree EXACTLY
  (provably: the discrepancy term vanishes), which is why null simulations never
  detect this.

  The onset of a 5% SE shortfall on `beta_left` is `|Z|` = 6.9 at MAF 0.005, 9.8
  at 0.01, 23.2 at 0.05 and 87.9 at 0.30 -- so common variants are safe at any
  effect size reachable in practice, and rare variants with strong effects are
  not. **Use `cov_source = "jackknife"` for rare variants, for component-specific
  parameters, or when mixture separation is below about 4**; it is nominal in
  every regime measured and costs nothing in the estimates. See
  [FINDINGS_covsource_default.md](docs/FINDINGS_covsource_default.md).

Stage 1 writes the pooled `Sigma_0` sidecar by default and Stage 2 scales it by
each SNP's residual genotype variation. The per-SNP jackknife remains available
as an explicit comparison; note that `se_tau_*` in the Stage 1 output are the
jackknife SEs and the pooled path DOES read them, for the `Q` lack-of-fit
statistic.

### LOCO predictions go in the basis, not in the phenotype

**Changed in 0.4.0.** `--loco-mode` now defaults to `covariate`: Stage 1 leaves
the phenotype unadjusted and appends `L(chr)` to the basis that residualises the
RIF columns and the genotype. It used to subtract `c * L(chr)` from the
phenotype before building the RIF.

This matters because RIF/UQR estimands are defined against the **unconditional
distribution of Y**. Subtracting a predictor changes `q_tau`, `f(q_tau)` and
every indicator -- it redefines the target rather than removing a nuisance,
unlike a mean GWAS where any `g`-independent subtraction leaves the slope alone.
Stage 1 already residualised covariates out of the RIF columns, which is
correct; LOCO was the inconsistency.

On testosterone chr22, against a sex-stratified REGENIE comparator carrying its
own within-stratum LOCO, subtracting cost 9.5-14.7% in standard error and moved
the median effect by a third of a standard error away from the reference.
`--loco-scale` is inert under `covariate`, since a projection is invariant to
rescaling the column. Pass `--loco-mode subtract` to reproduce pre-0.4.0 runs
byte-identically.

A step-1 ridge fitted without the mixture variable as a covariate will learn the
mixture rather than polygenic signal, which is what happened here. Check any
LOCO predictor with `scripts/loco_blind_diagnostic.R` before trusting it; it
needs no membership label. See [LOCO_MIXTURE_HAZARD.md](docs/LOCO_MIXTURE_HAZARD.md).

See [what the primary analysis estimates](docs/ESTIMAND.md) for effect units,
the local model interpretation and the scope of validation.

Current validation is intentionally narrow: testosterone on chromosome 22 first.
The primary result is Stage 2 OLS with pooled covariance; jackknife covariance is
retained as a sensitivity analysis. Commit `afb59ab` adds an optional
`fourth_order_gaussian` RIF density kernel and records its method in the
Stage 1/Stage 2 provenance, while preserving `silverman_gaussian` as the default.
That kernel was measured and **not taken forward**: it is unreliable for
multimodal, rounded and count phenotypes (issue #1), so do not use it for
results. Do not treat genome-wide,
other-phenotype or whole-paper runs as part of this validation phase.

**FungWas** allows you to test SNP associations with:
- **Variance** (vQTLs)
- **Mixture components** (e.g., latent subtypes in a disease)
- **Skewness / Kurtosis**
- **Any parametric system** definable by a map from quantiles.

## What is in this repository

| path | what it is |
|---|---|
| `stage1/` | the Python/C++ quantile GWAS. Installed on its own. |
| `stage2/` | the R/C++ projection onto a distributional model. Installed on its own. |
| `docs/` | user guidance; see [`docs/README.md`](docs/README.md). `docs/examples/1kg_chr22/` is a one-minute end-to-end check. |
| `scripts/` | `loco_blind_diagnostic.R`, a check on any LOCO predictor before trusting it. |
| `tests/` | tier 0, the fast cross-package unit suite (`bash tests/run_tier0.sh`), separate from the unit tests in `stage1/tests` and `stage2/tests`. |

The research record behind the method -- simulation campaigns (`sim/`), dated run
records (`experiments/`), the BluePebble/UK Biobank runners, the tier 1-4 simulation
suites and the audit and validation reports -- is not needed to use the packages and
lives at the git tag
[`archive/research-records-20261016`](https://github.com/gushamilton/GMoQ/tree/archive/research-records-20261016)
(`git checkout archive/research-records-20261016`). A file named in these docs or in a
code comment that is not on `main` is there.

## How It Works

FungWas operates in two stages:

See [Stage 2 Documentation](stage2/README.md) for full details.

1.  **Stage 1 (Python/C++):**
    -   Performs **Recentered Influence Function (RIF) regression** on the phenotype at multiple quantile levels (taus).
    -   This is a highly optimized, parallelized quantile GWAS.
    -   **Output:** Beta estimates and standard errors for each tau, plus the full covariance matrix of estimates across taus for each SNP.
    -   **Key Tech:** Python CLI, the `bgen` package for BGEN I/O, C++ kernel with OpenMP.

2.  **Stage 2 (R/C++):**
    -   Maps the quantile effects from Stage 1 to parameters of interest (e.g., mean/variance of a mixture model).
    -   Uses the **Delta Method** or efficient **GMM-like** estimators.
    -   Input: Stage 1 results + a "Weight Matrix" describing the derivative of quantiles with respect to parameters.
    -   **Output:** GWAS summary statistics (Beta, SE, P-value) for your custom parameters.

## Installation

### 1. Create Environment (Recommended)

We provide a Conda environment that handles Python dependencies and system libraries
(Armadillo, a C++ compiler and OpenMP). It works on Linux and on macOS, including Apple Silicon.
Both packages share one version (currently 0.6.1); see [stage2/NEWS.md](stage2/NEWS.md).

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

```bash
R CMD INSTALL stage2
```

### 4. Check it works

`docs/examples/1kg_chr22/` runs the whole pipeline on public 1000 Genomes genotypes with
simulated phenotypes in about a minute and prints PASS/FAIL against the simulated truth.

## Detailed Build Instructions

### System Requirements
- A C++ compiler with OpenMP: GCC >= 9, or LLVM clang with libomp. The conda
  environment provides one. Apple's own clang has no OpenMP; the build then
  looks for a separately installed libomp (conda-forge `llvm-openmp` or
  `brew install libomp`). If it finds none, it still builds but runs
  single-threaded, and says so.
- R >= 4.2 with `Rcpp`, `RcppArmadillo` and `data.table`. On macOS, Homebrew R
  can build `data.table` without zlib, which breaks writing `.gz` files; if
  `data.table::fwrite(data.frame(a = 1), tempfile(fileext = ".gz"))` fails,
  reinstall it with `install.packages("data.table", type = "source")`.
- Python >= 3.10 with `pip` and `numpy`.

### Build Stage 1 (Python/C++)

**You probably do not need a build flag.** Since the restructured kernel landed,
Stage 1's hot loops are compiled twice -- once at the portable SSE2 baseline and
once under AVX2+FMA -- and the right one is chosen at run time from
`__builtin_cpu_supports`. A single build is therefore both portable and fast,
and it is the default.

`FUNGWAS_STAGE1_ARCH` still exists, but it now buys 10-14% on the **exact**
kernel only (the `--exact-kernel` path), which is not run-time dispatched:

```bash
FUNGWAS_STAGE1_ARCH="-march=haswell -mtune=skylake-avx512" pip install -e stage1
```

Measured on full chr22 testosterone with LOCO, 4 physical cores, `--cov-mode
pooled-only`: 172.4s to 142.2s, i.e. 101.8 to 123.4 SNPs/s, with the outputs
byte-identical to the SSE2 build. It is not the default because
`-march=haswell` requires AVX2 (Intel 2013+ / AMD Zen+) and a binary built with
it will SIGILL on older hardware, with no graceful diagnostic -- and the binary
outlives the build host when the package sits on a shared filesystem. Setting it
also changes the exact kernel's arithmetic in the last bits, so a run that must
reproduce historical output byte for byte should leave it alone.

### Which Stage 1 kernel runs

The **restructured (fast) kernel is the default** and is what the timings in
this README refer to. On full chr22 at 4 threads it takes the kernel from 138.7s
to 33.0s (4.2x) and the wall clock from 190s to 104s; the gap widens with cores,
to 8.8x at 16, because the old kernel stops scaling past 8 and this one does
not.

It is **not bit-identical** to the old path -- the summation is reassociated --
but it is not meaningfully different either: across all 176 real chr22 batches
the largest disagreement is 5.4e-14 in |dbeta|/SE and 8.7e-13 relative in the
SE, no variant changes accept/reject, `stage1.tsv.gz` is byte-identical at the
`%.6g` resolution Stage 1 writes, and Stage 2 output is byte-identical to the
accepted chr22 gate.

Pass `--exact-kernel` when a run must reproduce historical output byte for byte
including `cov.gz`. See [docs/STAGE1_FAST_KERNEL.md](docs/STAGE1_FAST_KERNEL.md)
for the derivation, the numerics and the per-change timings.

Note also that **more threads is not always faster**: on a node whose "8 CPUs"
are 4 physical cores plus SMT siblings, 8 threads is slower than 4 in every
mode, because filling the siblings starves the decode and write threads the
pipeline overlaps. Match `--threads` to PHYSICAL cores.


```bash
conda activate fungwas-stage1
cd stage1
pip install -e .
python -c "from fungwas_stage1 import core; print(core.HAVE_CPP)"
```

If `core.HAVE_CPP` is `False`, the extension did not build: check the `pip install`
output for compiler errors and reinstall. A build without OpenMP prints a warning
saying so; set `FUNGWAS_STAGE1_OPENMP=0` to ask for one deliberately.

### Build Stage 2 (R/C++)

```bash
R CMD INSTALL stage2
```

If you install into a user library, set:

```bash
export R_LIBS_USER=/path/to/your/R/library
```

## Testing

Both stages ship their own checks. Tests should use the installed package and must not
regenerate package documentation or modify `NAMESPACE`.

### Stage 1 (pytest)

```bash
conda activate fungwas-stage1
cd stage1 && python setup.py build_ext --inplace && cd ..
PYTHONPATH=stage1 python -m pytest stage1/tests/ -q
```

164 tests, a few seconds.

**Rebuild the extension after pulling anything that touches `stage1/src/`.**
A stale `_stage1_cpp` produces a wall of failures that look like broken code but
are only a missing symbol -- pulling the `--cov-mode pooled-only` work without
rebuilding fails 13 tests, because it added `panel_pooled_stats` to the pybind
layer, and pulling the restructured kernel without rebuilding fails on a missing
`BlockPlan`. The `build_ext --inplace` line above is not optional after such a
pull.

### Stage 2 (R)

`stage2/tests/test-*.R` are standalone assertion scripts against the installed package;
each exits non-zero on failure. `stage2/tests/testthat/` holds the testthat suite, which
`R CMD check` runs.

```bash
R CMD INSTALL stage2
for t in stage2/tests/test-*.R; do Rscript --vanilla "$t" || exit 1; done
Rscript -e 'testthat::test_dir("stage2/tests/testthat", package = "fungwasStage2", load_package = "installed")'
bash tests/run_tier0.sh     # tier 0: both of the above plus the cross-package unit checks
```

### End-to-end check

[`docs/examples/1kg_chr22/`](docs/examples/1kg_chr22/README.md) is the quickest
check: real 1000 Genomes chr22 genotypes, phenotypes simulated under a null, a
mean shift, a variance effect and a two-component mixture, then Stage 1 and
Stage 2, then a PASS/FAIL comparison with the truth.

**Never run these inside a directory holding, or writing to, individual-level data.** Write all
run outputs outside this repository.

## Workflow Guide

### Step 0: Data Preparation
-   **Genotypes:** BGEN format. A `bgenix` index (`.bgen.bgi`) is optional: it
    speeds up `--snps` subsets, and without one Stage 1 reads the variant
    headers instead.
-   **Phenotypes:** Standard space-delimited text file (FID IID Pheno).
-   **Covariates:** Standard space-delimited text file (FID IID Age Sex ...).

See [Phenotype Preparation Guide](docs/phenotype_preparation.md) for details.

**Sample identifiers -- get this right before anything else.** Stage 1 matches
the phenotype file's **first column (FID)** against the BGEN sample identifiers,
with **no alias resolution**. It does not try IID, `FID_IID`, or `FID:IID`. If
your BGEN sample IDs are not what sits in column 1 of the phenotype file,
rewrite that column; there is no flag that will do it for you.

Pass `--sample <file>` when the BGEN does not embed sample identifiers. Without
it, and without embedded IDs, the run cannot align anything.

This asymmetry is the one worth knowing: the LOCO reader **does** resolve five
aliases (`fid`, `iid`, `fid_iid`, `"fid iid"`, `fid:iid`), so a cohort whose
BGEN IDs are `FID_IID` will match the LOCO predictions happily and then drop
every sample here. Stage 1 refuses below `--min-sample-match` (1% by default)
and warns between there and 50%, because a *partial* match is the dangerous
case -- it runs to completion on a small non-random subset and produces
attenuated effects with inflated SEs that are internally consistent and wrong.
Check the matched-sample count in the log before trusting a run.

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
    --sample data/chr22.sample \
    --out results/chr22 \
    --threads 16
```
*Drop `--sample` only if the BGEN embeds its sample identifiers. Match
`--threads` to PHYSICAL cores, not to the CPU count the scheduler reports.*
*Produces `results/chr22.stage1.tsv.gz` (estimates), per-SNP covariance files, and
`results/chr22.sigma0.tsv.gz` with metadata (pooled covariance).*
By default, Stage 1 writes covariance in compact int8 format with sidecar
`results/chr22.cov.scale.gz` (dequantization scales).

The default grid is 45 quantiles from 0.02 to 0.98, with 1,000 jackknife blocks.
When supplying REGENIE LOCO predictions, Stage 1 uses `--loco-scale regression`
by default to convert predictions to the phenotype's scale before subtraction.

`--loco-mode covariate` instead leaves the phenotype unadjusted and appends the
per-chromosome prediction to the basis the RIF columns and the genotype are
residualised on. A RIF/UQR estimand is defined against the *unconditional*
distribution of `Y`, so subtracting a prediction from `Y` redefines the target
rather than removing a nuisance; residualising the influence-function columns
does not. On testosterone chr22 this recovers the precision and the REGENIE
agreement of a no-LOCO scan. See [`docs/LOCO_AS_COVARIATE.md`](docs/LOCO_AS_COVARIATE.md).
Since 0.4.0 the default is `covariate` whenever `--loco-preds` is given; with
no LOCO predictions the mode is moot and resolves to `subtract`. Pass
`--loco-mode subtract` explicitly to reproduce pre-0.4.0 output.

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

Keep pooled covariance for the primary scan. To compare a result with the
retained per-SNP jackknife, use the same Stage 1 file and weights:

```r
jackknife_results <- param_gwas_from_file(
  stage1_file = "results/chr22.stage1.tsv.gz",
  W = W,
  se_mode = "tau_cov",
  tau_cov_estimator = "ols",
  cov_source = "jackknife",
  cov_file = "results/chr22.cov.gz"
)
```

The OLS point estimates are unchanged; covariance and reference calibration
can change the standard errors and p-values. Keep the covariance sidecars
beside the `.cov.gz` file so the reader can check variant alignment and metadata.

**Mixture weights: which density.** The mixture and generic weight
builders take `density = "model"`, `"hybrid"` or `"stage1"`. **Since 0.6.0 the
default is `"hybrid"` whenever W has a membership (`gamma`) column**, and
`"model"` (every column divided by the fitted density) otherwise. Pass
`density = "model"` to reproduce 0.5.0. For a
two-component mixture with a membership parameter, `"hybrid"` changes only the
`gamma` column. It divides that column by the Stage 1 kernel density at
`q_tau`, the density the Stage 1 coefficients carry, and the component shares
still sum to 1. In separated mixtures this removes gamma's attenuation and its
leakage into the component locations.

```r
fhat <- stage1_kde_density(y_stage1, q_tau)   # or read_stage1_density("rif_nuisance.json")
W <- make_weights_normal_mixture(taus, q_tau, p1, mu1, sd1, mu2, sd2,
                                 include_membership = TRUE,
                                 stage1_density = fhat)   # density = "hybrid" by default
res <- param_gwas_from_file(..., W = W, require_density_mode = "hybrid")
```

See [docs/MIXTURE_WEIGHT_DENSITY.md](docs/MIXTURE_WEIGHT_DENSITY.md) for when
to use each mode.

## Advanced Topics

### Standard Error Modes
Stage 2 supports different methods for SE calculation (`param_gwas` argument `se_mode`):
-   `tau_cov`: **Recommended**. Uses full cross-quantile covariance. With `cov_source="pooled"`, this is `Sigma_0 / d_resid`; with `cov_source="jackknife"`, it is the estimated per-SNP covariance. The latter supports float32 and int8 storage. These are covariance estimates, not exact sampling covariances.
-   `plugin_cor`: Uses a single correlation matrix for all SNPs (from null distribution). Faster, good for initial scans.
-   `diagonal`: Assumes independence between quantiles. Least accurate, not recommended for final results.

### Proof & Verification
-   `stage1/tests/test_jk_equivalence.py`: Stage 1 jackknife SEs match analytic theory.
-   `docs/examples/1kg_chr22/`: recovery of known shift, variance and mixture effects, and
    null calibration, on real genotypes.

## License
MIT; see [LICENSE](LICENSE).
