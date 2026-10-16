# G-MoQ on 1000 Genomes with simulated phenotypes

A one-minute end-to-end check: real genotypes (1000 Genomes phase 3 chr22, all 2,504 samples,
2,000 random common SNPs), simulated phenotypes with known effects, Stage 1 then Stage 2, then a
PASS/FAIL check that the truth is recovered and the null SNPs are calibrated.

Only public 1000 Genomes data is used. Data and outputs go to `~/gmoq_1kg_example/` by default,
outside the repository; set `GMOQ_1KG_DIR` to put them elsewhere.

## Prerequisites

- Stage 1 installed (`fungwas-stage1` on PATH) and Stage 2 installed (`library(fungwasStage2)`);
  see the top-level README.
- [`plink2`](https://www.cog-genomics.org/plink/2.0/) on PATH (or `PLINK2=/path/to/plink2`), `curl`,
  and R with `data.table` and `jsonlite`.
- About 300 MB of disk.

## Run

```bash
bash docs/examples/1kg_chr22/00_get_data.sh   # once: chr22 VCF + index, ~200 MB
bash docs/examples/1kg_chr22/run_all.sh       # ~1 min: BGEN, simulate, Stage 1, Stage 2, checks
SEED=7 bash docs/examples/1kg_chr22/run_all.sh  # same genotypes, new causal SNPs and phenotypes
```

Settings (paths, threads) are in `config.sh` and can be overridden from the environment.

| step | what it does |
|---|---|
| `01_make_bgen.sh` | biallelic SNPs, MAF ≥ 0.05; 5 PCs from LD-pruned SNPs; 2,000 random SNPs as BGEN v1.2 (`ref-first`, so allele 1 = REF = Stage 1's effect allele) |
| `02_simulate.R` | four phenotypes plus `truth.tsv` (below) |
| `03_stage1.sh` | `fungwas-stage1`, 45 quantiles, PCs as covariates, `--blocks 100` |
| `04_stage2.R` | `param_gwas_from_file()` with pooled and jackknife covariance |
| `05_check.R` | `recovery.tsv`, `null_calibration.tsv`, `qq.png`, PASS/FAIL |

## The simulated models

`G` is the allele-1 dosage and `e ~ N(0,1)`. Every model except the mixture adds `0.3·PC1`.

| phenotype | model | Stage 2 weights | truth |
|---|---|---|---|
| `null` | no SNP effect | Normal mean/variance | everything 0 |
| `shift` | `Y = 0.25·G + e`, 5 SNPs | Normal mean/variance | `beta_mu = 0.25`, `beta_sigma2 = 0` |
| `variance` | `Y = sqrt(1 + 0.3·(M − mean M))·e`, 3 SNPs, where M counts the minor allele | Normal mean/variance | `beta_mu = 0`, `beta_sigma2 = ±0.3` per copy of allele 1 (+ if allele 1 is the minor allele), **on the variance scale** |
| `mixture` | N(0,1) with probability `p1`, else N(4,1); `logit p1 = logit 0.7 + 0.5·(G − mean G)`, 5 SNPs | Normal mixture with membership: `density = "hybrid"` (default) and `"model"` | `gamma = 0.5`, `beta_1 = beta_2 = 0` |

Notes:
- Stage 1 matches the phenotype file's first column (FID) to the BGEN sample IDs, so FID = IID here.
- With N = 2,504, the default `--blocks 1000` leaves about 2.5 samples per block, so the example
  uses 100. Scale `--blocks` with N.
- Ancestry is mixed across all five superpopulations, so the 5 PCs are in the model as covariates.
- A SNP counts as null for a model if its r² with every causal SNP of that model is below 0.01.
  Chr22 LD in a mixed-ancestry sample is long-range, so this keeps roughly 600 to 900 of the
  2,000 SNPs.

## What "pass" looks like

`05_check.R` prints `RECOVERY PASS` and `CALIBRATION PASS` when:
- every causal SNP's estimate is within 3 SEs of its truth, on every parameter (with a truth of 0
  off target), and its target parameter has p < 1e-4;
- for every parameter, the share of null SNPs with p < 0.05 is within 3 binomial SDs of 0.05.
  λ_GC is reported too, but not judged, because a median over a few hundred SNPs is too noisy to
  gate on.

The mixture is fitted twice. The default `density = "hybrid"` is judged. The `density = "model"`
fit (the default before 0.6.0) is reported but not judged: in this well-separated mixture it
underestimates `gamma` by about 30%, which is why the default changed.

Expected result for seed 2026, pooled covariance. Seeds 7 and 99 also pass.

| model | parameter | truth | estimates (5 or 3 causal SNPs) |
|---|---|---|---|
| shift | `beta_mu` | 0.25 | 0.21, 0.25, 0.19, 0.25, 0.28 (SE ≈ 0.033) |
| variance | `beta_sigma2` | ±0.30 | −0.39, 0.28, −0.36 (SE ≈ 0.05) |
| mixture, hybrid | `gamma` | 0.50 | 0.38, 0.42, 0.45, 0.43, 0.56 (SE ≈ 0.07) |
| mixture, model | `gamma` | 0.50 | 0.31, 0.34, 0.36, 0.35, 0.46 (attenuated) |

Across every model, parameter and covariance source, 4.5% to 6.5% of null SNPs have p < 0.05,
and the Q lack-of-fit test rejects 4.5% to 6.5% of them. Pooled and jackknife standard errors
agree, as expected at MAF ≥ 0.05.
