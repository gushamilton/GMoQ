# Documentation index

## If you are going to run G-MoQ

Read these four, in order:

| | |
|---|---|
| [ESTIMAND.md](ESTIMAND.md) | **What a G-MoQ effect means.** Unconditional (RIF) quantile effects, not conditional ones. Anything subtracted from Y redefines the estimand. Read this first or you will misinterpret the output. |
| [phenotype_preparation.md](phenotype_preparation.md) | How to prepare the phenotype. The estimand is defined against its unconditional distribution, so this is not boilerplate. |
| [ALLELE_CONVENTION.md](ALLELE_CONVENTION.md) | Which allele the effect is counted against. Get this wrong and every sign flips. |
| [KNOWN_ISSUES.md](KNOWN_ISSUES.md) | What the package does not do well, stated plainly. |

Then run [examples/1kg_chr22/](examples/1kg_chr22/README.md) to check your install end to end.

Also: [BUILD.md](BUILD.md) (requirements, non-conda install, troubleshooting, tests) and
[BACKGROUND.md](BACKGROUND.md) (short notes on the estimand and the defaults).

## If you hit a specific question

| Question | File |
|---|---|
| Which covariance source, and when does it matter? | [FINDINGS_covsource_default.md](FINDINGS_covsource_default.md) |
| Why is the estimator OLS and not GLS? | [OLS_VS_GLS.md](OLS_VS_GLS.md), [WHY_GLS_DEGRADES.md](WHY_GLS_DEGRADES.md) |
| How should I handle LOCO / REGENIE predictions? | [LOCO_AS_COVARIATE.md](LOCO_AS_COVARIATE.md), [LOCO_SCALE.md](LOCO_SCALE.md), [LOCO_MIXTURE_HAZARD.md](LOCO_MIXTURE_HAZARD.md) |
| Which density divides a mixture weight column? | [MIXTURE_WEIGHT_DENSITY.md](MIXTURE_WEIGHT_DENSITY.md) |
| How do I get SEs for a function of the fitted parameters? | [FUNCTIONAL_GRADIENTS.md](FUNCTIONAL_GRADIENTS.md) |
| Where does the RIF's own uncertainty go? | [DENSITY_RIF_NUISANCE.md](DENSITY_RIF_NUISANCE.md) |

## Research record

The research notes, simulation campaigns and validation reports behind the method are not
needed to use the packages. They are at the tag
[`archive/research-records-20261016`](https://github.com/gushamilton/GMoQ/tree/archive/research-records-20261016/docs).
