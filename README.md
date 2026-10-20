<p align="center"><img src="docs/logo/gmoq_logo.png" alt="G-MoQ" width="420"></p>

# G-MoQ: Genomic Method of Quantiles

GWAS on the whole distribution of a continuous phenotype, not just its mean: variance effects (vQTLs),
GWAS of mixture distributions, the means of the mixture components, and their membership. G-MoQ enables GWAS of any statistics, biological, or evolutionary parametric model you can write in terms of its effects on the quantiles of the phenotype distribution.

## How it works

Two stages, installed separately:

1. **Stage 1** (Python/C++ package `gmoq`, command `gmoq-stage1`). An unconditional quantile GWAS: for every
   SNP, the effect on 45 quantiles of the phenotype (unconditional quantile regression via the
   recentered influence function, RIF), plus the covariance of those effects across quantiles.
2. **Stage 2** (R package `gmoq`). Maps each SNP's 45 quantile effects onto the parameters of a
   model you choose — a mean and a variance, or a two-component mixture, or your own — by least
   squares through a weight matrix `W`, with standard errors from the Stage 1 covariance.

A G-MoQ effect is an *unconditional* quantile effect, defined on the distribution of the
phenotype you analysed; see [docs/ESTIMAND.md](docs/ESTIMAND.md). The method, its defaults and
its known limitations are described in [docs/METHOD.md](docs/METHOD.md).

## Install

With conda (Linux or macOS, including Apple Silicon):

```bash
conda env create -f stage1/environment.yml
conda activate gmoq
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
gmoq-stage1 \
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
library(gmoq)
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

More in [docs/INPUT_PREPARATION.md](docs/INPUT_PREPARATION.md).

## Repository layout

| path | contents |
|---|---|
| `stage1/` | Stage 1 Python package with its C++ kernel (`gmoq/`, `src/`) and tests |
| `stage2/` | Stage 2 R package, with tests and vignettes |
| `docs/` | user documentation ([index](docs/README.md)) and the 1000 Genomes example |
| `stage2/NEWS.md` | changes; both packages share one version |

## License

MIT; see [LICENSE](LICENSE).
