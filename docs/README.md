# Documentation

| | |
|---|---|
| [BUILD.md](BUILD.md) | Installing and building both packages, and running the tests. |
| [METHOD.md](METHOD.md) | How the method works: Stage 1 RIF regression, the pooled and jackknife covariances, the Stage 2 mapping, p-values and the Q test, LOCO predictions, dispersion artefacts, and known limitations. |
| [ESTIMAND.md](ESTIMAND.md) | What a G-MoQ effect means. Read this before interpreting results. |
| [INPUT_PREPARATION.md](INPUT_PREPARATION.md) | Genotype, phenotype and covariate files, sample matching, missing values, transformations and LOCO predictions. |
| [ALLELE_CONVENTION.md](ALLELE_CONVENTION.md) | Which allele the effects are counted against, and how to align results with other software. |
| [MIXTURE_WEIGHT_DENSITY.md](MIXTURE_WEIGHT_DENSITY.md) | The `density` option of the mixture weight builders, and why the default for a membership parameter is `"hybrid"`. |
| [examples/1kg_chr22/](examples/1kg_chr22/README.md) | A one-minute worked example on public 1000 Genomes genotypes with simulated phenotypes, checking that known effects are recovered. Run it to check an installation end to end. |

Function-level documentation for Stage 2 is in the R help pages (`?param_gwas_from_file`,
`?make_weights_normal_mixture`, ...) and the package vignettes. For Stage 1, run
`gmoq-stage1 --help`.
