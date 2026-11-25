# Preparing the inputs

Stage 1 takes the **raw phenotype** and computes the recentered influence functions itself. Do
not supply a residualised or pre-adjusted phenotype: the estimand is defined on the distribution
of the phenotype you give it (see [ESTIMAND.md](ESTIMAND.md)).

## Genotypes

BGEN v1.2 files with genotype probabilities (8-bit is enough). One file per chromosome is usual;
run Stage 1 once per file. Index each file with `bgenix -index` (`.bgen.bgi`): with an index and
`--snps`, Stage 1 reads only the requested variants and decodes them in parallel. Effects are per copy of the
first allele of each record ([ALLELE_CONVENTION.md](ALLELE_CONVENTION.md)).

If the BGEN does not embed sample IDs, pass the `.sample` file with `--sample`.

## Phenotype file

Space- or tab-delimited, with a header whose first two columns are `FID` and `IID`:

```
FID IID trait_a trait_b
1001 1001 15.3 24.5
1002 1002 12.1 26.2
1003 1003 NA 23.8
```

- Select the column with `--pheno-col`, or several with `--pheno-cols a,b` (analysed in one pass
  over the genotypes, on the samples complete for all of them).
- **Sample matching.** Stage 1 matches the **first column (FID)** to the BGEN sample IDs, exactly.
  If your BGEN IDs have the form `FID_IID`, rewrite the FID column to match. Check the
  `Matched N/M samples` line in the log: a near-total failure to match stops the run
  (`--min-sample-match`), and a partial match only warns.
- **Missing values.** Only `NA` means missing. Rows with a missing or non-numeric phenotype or
  covariate, or with no row in the covariate file, are dropped, and the log says how many and why.
  Codes such as `-9` are analysed as numbers; Stage 1 warns when it sees `-9`.
- `--keep FILE` restricts the analysis to the FIDs listed (one per line).

## Covariate file

Same layout as the phenotype file:

```
FID IID age sex PC1 PC2 PC3 PC4 PC5
1001 1001 45 1 0.012 -0.034 0.056 -0.023 0.089
```

Every non-ID column is used unless you choose with `--covar-col`. Typical covariates are age, sex,
genetic principal components and genotyping batch or assessment centre. Covariates are projected
out of the RIF and the genotypes, not out of the phenotype. Because each quantile gets a linear
covariate coefficient, include the terms (for example interactions or nonlinear terms) that the
phenotype's distribution depends on.

## Transformations

Choose the scale on which you want to interpret effects before Stage 1, and analyse that. For
right-skewed measurements, a log transform usually removes most of the mean-variance coupling
that otherwise shows up as spurious variance effects (see [METHOD.md](METHOD.md)). A rank-based
inverse normal transform removes it entirely, but also removes the distributional shape that
mixture models try to describe:

```r
int_transform <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
```

Remove clear measurement errors. Heavily tied phenotypes (few distinct values) are outside what
has been validated.

## LOCO predictions

To account for polygenic background and relatedness, pass REGENIE step 1 predictions with
`--loco-preds` (a `.loco` file, or the `_pred.list` for several phenotypes). The prediction is
added to the covariates per chromosome. For phenotypes that are mixtures of groups, read the LOCO
section of [METHOD.md](METHOD.md) first.

## Stage 2 inputs

Stage 2 needs the Stage 1 output files (keep them together) and a weight matrix `W` built from
the quantiles of the phenotype values Stage 1 analysed: the same transformation and the same
samples (non-missing, with complete covariates, present in the BGEN). Use the Stage 1 quantile
grid; it is recorded in `<out>.sigma0.meta.json`.
