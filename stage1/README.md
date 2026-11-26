# Stage 1: RIF quantile GWAS (`gmoq-stage1`)

For every SNP in a BGEN file, Stage 1 estimates its effect on 45 quantiles of the phenotype
(unconditional quantile regression via the recentered influence function, RIF) and the
covariance of those effects across quantiles. Stage 2 reads the output.

## Install

```bash
conda env create -f environment.yml && conda activate gmoq
pip install .
```

Without conda, or if something goes wrong: [../docs/BUILD.md](../docs/BUILD.md).

## Run

```bash
gmoq-stage1 \
    --bgen data/chr22.bgen \
    --pheno pheno.txt --pheno-col bmi \
    --covar covar.txt \
    --out results/chr22 \
    --threads 8
```

### Inputs

- **`--bgen`**: BGEN v1.2 dosages. Effects are per copy of allele 1 (`effect_allele` in the
  output). Index it with `bgenix -index` and pass `--snps`: Stage 1 then reads only those
  variants and decodes them in parallel, which is much faster.
- **`--pheno`**: space- or tab-delimited, header `FID IID <columns>`, the raw phenotype (Stage 1
  builds the RIF itself). The first column (FID) must equal the BGEN sample IDs. Only `NA` means
  missing.
- **`--covar`** (optional): same layout. All non-ID columns are used unless you pass
  `--covar-col`. Covariates are projected out of the RIF and the genotypes, not out of Y.
- **`--sample`**: only if the BGEN does not embed sample IDs.

The log reports the samples loaded, the rows dropped and why, and the samples matched to the
BGEN. Read it.

### Options you are likely to use

| option | default | what it does |
|---|---|---|
| `--pheno-col` / `--pheno-cols a,b` | | one phenotype, or several in one genotype pass (samples complete for all of them) |
| `--covar-col` | all | which covariates to use |
| `--snps FILE` | all | restrict to these rsIDs (one per line) |
| `--keep FILE` | all | restrict to these FIDs (one per line) |
| `--threads` | 1 | set to the number of physical cores |
| `--blocks` | 1000 | jackknife blocks; use fewer (e.g. 100) below N of about 5,000 |
| `--taus` | 45 from 0.02 to 0.98 | comma-separated quantile grid |
| `--loco-preds FILE` | | REGENIE step 1 `.loco` or `_pred.list`; added to the covariate basis |

`gmoq-stage1 --help` lists every option.

### Several phenotypes in one pass

Each phenotype can have its own sample list, covariates and LOCO file. Put one analysis per
line in a text file, with the options that differ, and give the shared options once:

```text
--pheno pheno.txt --pheno-col bmi --keep keep_bmi.txt --covar covar.txt --out results/bmi_chr22
--pheno pheno.txt --pheno-col ldl --keep keep_ldl.txt --covar covar.txt --out results/ldl_chr22
```

```bash
gmoq-stage1 --jobs jobs.txt --bgen data/chr22.bgen --snps hm3_chr22.txt --threads 8
```

Stage 1 reads each variant once and stores dosages only for samples matched to at least one
analysis. Each line uses its own matched samples and produces the same output as running it alone.

### Outputs

| file | contents |
|---|---|
| `{out}.stage1.tsv.gz` | per SNP: `snp_id chr bp effect_allele other_allele`, `beta_tau_*` and `se_tau_*` for each quantile, and `d_resid`. With `--cov-mode pooled-only`, `z_tau_*` instead (`se_tau = sqrt(diag(Sigma_0) / d_resid)`, `beta_tau = z_tau * se_tau`; Stage 2 does this for you) |
| `{out}.sigma0.tsv.gz`, `{out}.sigma0.meta.json` | pooled cross-quantile covariance, used by Stage 2 by default |
| `{out}.cov.gz`, `{out}.cov.scale.gz`, `{out}.cov.ids.tsv.gz`, `{out}.cov.meta.json` | per-SNP jackknife covariance, for Stage 2's `cov_source = "jackknife"` |

Keep these files together; Stage 2 reads the sidecars next to the file you give it. Variants
with more than two alleles, rsIDs with several BGEN records, and SNPs with minor allele count
below 5 are skipped.

## Python API

For small in-memory data:

```python
import numpy as np
from gmoq import run_stage1

G = np.random.binomial(2, 0.3, (5000, 100)).astype(float)   # N x M dosages
Y = np.random.randn(5000)                                    # raw phenotype
res = run_stage1(G, Y, taus=np.round(np.linspace(0.02, 0.98, 45), 3), n_threads=4)
res["betas"].shape, res["se"].shape, res["covariances"][0].shape   # (M, T), (M, T), (T, T)
```

## Tests

From the repository root (so the installed package, not the source tree, is imported):

```bash
python -m pytest stage1/tests -q
```

How the method works: [../docs/METHOD.md](../docs/METHOD.md).
