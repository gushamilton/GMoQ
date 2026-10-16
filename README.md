# G-MoQ: Genomic Method of Quantiles

GWAS on the whole distribution of a phenotype, not just its mean: variance effects (vQTLs),
mixture membership, and any parametric model you can write as a map from quantiles.

The Python package is `fungwas_stage1` (command `fungwas-stage1`) and the R package is
`fungwasStage2`.

## How it works

Two stages, installed separately:

1. **Stage 1** (`fungwas-stage1`, Python/C++). A quantile GWAS: for every SNP, the effect on 45
   quantiles of the phenotype (unconditional quantile regression via the recentered influence
   function, RIF), plus the covariance of those effects across quantiles.
2. **Stage 2** (`fungwasStage2`, R). Maps each SNP's 45 quantile effects onto the parameters of
   a model you choose — a mean and a variance, or a two-component mixture, or your own — by
   least squares through a weight matrix `W`, with standard errors from the Stage 1 covariance.

A G-MoQ effect is an *unconditional* quantile effect, defined on the distribution of the
phenotype you analysed; see [docs/ESTIMAND.md](docs/ESTIMAND.md). More background is in
[docs/BACKGROUND.md](docs/BACKGROUND.md).

## Install

With conda (Linux or macOS, including Apple Silicon):

```bash
conda env create -f stage1/environment.yml
conda activate fungwas-stage1
pip install ./stage1
R CMD INSTALL stage2
```

Without conda, install Armadillo (and libomp on macOS) first — `brew install armadillo libomp`
or `apt install libarmadillo-dev` — then `pip install ./stage1` in a virtual environment.
Requirements, troubleshooting and running the tests: [docs/BUILD.md](docs/BUILD.md).

## Check it works

```bash
bash docs/examples/1kg_chr22/00_get_data.sh   # public 1000 Genomes chr22, ~200 MB
bash docs/examples/1kg_chr22/run_all.sh       # ~1 minute; prints RECOVERY / CALIBRATION PASS
```

It simulates a null, a mean shift, a variance effect and a two-component mixture on real
genotypes, runs both stages and compares the estimates with the truth.
See [docs/examples/1kg_chr22/](docs/examples/1kg_chr22/README.md).

## Quick start

**Stage 1** — one BGEN, a phenotype file and a covariate file (space-delimited, `FID IID ...`):

```bash
fungwas-stage1 \
    --bgen data/chr22.bgen \
    --pheno inputs/pheno.txt --pheno-col bmi \
    --covar inputs/covar.txt \
    --out results/chr22 \
    --threads 8
```

This writes `results/chr22.stage1.tsv.gz` (quantile effects), `results/chr22.sigma0.tsv.gz`
(pooled covariance) and the per-SNP covariance files `results/chr22.cov.*`. Run it once per
BGEN (e.g. per chromosome); Stage 2 then reads each output separately.
All options: [stage1/README.md](stage1/README.md).

**Stage 2** — build `W` for your model from the phenotype values Stage 1 analysed (the samples
it kept: non-missing, with covariates, present in the BGEN), then map:

```r
library(fungwasStage2)
y     <- data.table::fread("inputs/pheno.txt")$bmi      # restrict to the samples Stage 1 kept
taus  <- round(seq(0.02, 0.98, length.out = 45), 3)   # Stage 1's default grid
q_tau <- as.numeric(quantile(y, taus, type = 8))

W <- make_weights_vqtl(taus, q_tau, mu = mean(y), sd = sd(y))   # mean and variance

res <- param_gwas_from_file(
  stage1_file = "results/chr22.stage1.tsv.gz", W = W,
  sigma0_file = "results/chr22.sigma0.tsv.gz"
)
head(res)   # beta_mu, beta_sigma2 (per allele, on the variance scale), SEs, p-values, Q
```

Mixtures and other models: [stage2/README.md](stage2/README.md).

## Before you run

- **Sample IDs.** Stage 1 matches the phenotype file's **first column (FID)** to the BGEN sample
  IDs exactly. If your BGEN IDs look like `FID_IID`, rewrite the FID column to match. Then check
  the "Matched N samples" line in the log: a partial match only warns.
- **Covariates.** Every non-ID column of `--covar` is used unless you pass `--covar-col`.
  Do not pre-adjust the phenotype for covariates; Stage 1 handles them.
- **Missing values.** Only `NA` means missing. Stage 1 warns about dropped rows and about `-9`,
  which it otherwise analyses as a number.
- **Effect allele.** Effects are per copy of BGEN allele 1, reported as `effect_allele`; see
  [docs/ALLELE_CONVENTION.md](docs/ALLELE_CONVENTION.md).

## More

- [docs/README.md](docs/README.md) — index of the user documentation.
- [docs/KNOWN_ISSUES.md](docs/KNOWN_ISSUES.md) — what the package does not do well.
- [stage2/NEWS.md](stage2/NEWS.md) — changes; both packages share one version.
- The research record (simulation campaigns, run records, validation reports) is at the tag
  [`archive/research-records-20261016`](https://github.com/gushamilton/GMoQ/tree/archive/research-records-20261016).

## License

MIT; see [LICENSE](LICENSE).
