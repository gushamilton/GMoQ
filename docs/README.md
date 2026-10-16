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

## If you hit a specific question

| Question | File |
|---|---|
| Which covariance source, and when does it matter? | [FINDINGS_covsource_default.md](FINDINGS_covsource_default.md) |
| Why is the estimator OLS and not GLS? | [OLS_VS_GLS.md](OLS_VS_GLS.md), [WHY_GLS_DEGRADES.md](WHY_GLS_DEGRADES.md) |
| How should I handle LOCO / REGENIE predictions? | [LOCO_AS_COVARIATE.md](LOCO_AS_COVARIATE.md), [LOCO_SCALE.md](LOCO_SCALE.md), [LOCO_MIXTURE_HAZARD.md](LOCO_MIXTURE_HAZARD.md) |
| Which density divides a mixture weight column? | [MIXTURE_WEIGHT_DENSITY.md](MIXTURE_WEIGHT_DENSITY.md) |
| How do I get SEs for a function of the fitted parameters? | [FUNCTIONAL_GRADIENTS.md](FUNCTIONAL_GRADIENTS.md) |
| How fast is Stage 1, and which kernel runs? | [STAGE1_FAST_KERNEL.md](STAGE1_FAST_KERNEL.md) |
| Where does the RIF's own uncertainty go? | [DENSITY_RIF_NUISANCE.md](DENSITY_RIF_NUISANCE.md) |

## Archived research notes

These are not needed to run the packages. They live at the tag
[`archive/research-records-20261016`](https://github.com/gushamilton/GMoQ/tree/archive/research-records-20261016/docs),
together with the simulation campaigns, run records and audit reports they cite.

| Question | File |
|---|---|
| My phenotype is lumpy / has a point mass | [BP2_lumpy_phenotypes.md](https://github.com/gushamilton/GMoQ/blob/archive/research-records-20261016/docs/BP2_lumpy_phenotypes.md) |
| I found a variance effect -- is it real? | [MEAN_VARIANCE_ARTEFACT.md](https://github.com/gushamilton/GMoQ/blob/archive/research-records-20261016/docs/MEAN_VARIANCE_ARTEFACT.md), [DISPERSION_SKEW_CORRECTION.md](https://github.com/gushamilton/GMoQ/blob/archive/research-records-20261016/docs/DISPERSION_SKEW_CORRECTION.md) |
| How does this behave as N grows? | [N_SCALING.md](https://github.com/gushamilton/GMoQ/blob/archive/research-records-20261016/docs/N_SCALING.md), [EFFICIENCY_BOUND.md](https://github.com/gushamilton/GMoQ/blob/archive/research-records-20261016/docs/EFFICIENCY_BOUND.md) |
| Related samples? | [RELATEDNESS.md](https://github.com/gushamilton/GMoQ/blob/archive/research-records-20261016/docs/RELATEDNESS.md) |
| Two-pass LOCO? | [TWO_PASS_LOCO.md](https://github.com/gushamilton/GMoQ/blob/archive/research-records-20261016/docs/TWO_PASS_LOCO.md) |
| Can I run several phenotypes in one pass? | [MULTI_PHENOTYPE_FUSION.md](https://github.com/gushamilton/GMoQ/blob/archive/research-records-20261016/docs/MULTI_PHENOTYPE_FUSION.md) (answer: no, it is slower) |
| Where do the shipped weight matrices come from? | [EMPIRICAL_W.md](https://github.com/gushamilton/GMoQ/blob/archive/research-records-20261016/docs/EMPIRICAL_W.md) |
| Every G-MoQ-vs-REGENIE number | [COMPARATOR_LEDGER.md](https://github.com/gushamilton/GMoQ/blob/archive/research-records-20261016/docs/COMPARATOR_LEDGER.md), [POWER_AND_CONSERVATISM.md](https://github.com/gushamilton/GMoQ/blob/archive/research-records-20261016/docs/POWER_AND_CONSERVATISM.md) |
| Every simulation campaign and what it concluded | [SIMULATION_INVENTORY.md](https://github.com/gushamilton/GMoQ/blob/archive/research-records-20261016/docs/SIMULATION_INVENTORY.md), [TEST_SUITE.md](https://github.com/gushamilton/GMoQ/blob/archive/research-records-20261016/docs/TEST_SUITE.md) |
| Dated audits and validation reports | [AUDIT_2026-09-16.md](https://github.com/gushamilton/GMoQ/blob/archive/research-records-20261016/docs/AUDIT_2026-09-16.md), [VALIDATION_2026-09-15.md](https://github.com/gushamilton/GMoQ/blob/archive/research-records-20261016/docs/VALIDATION_2026-09-15.md), [VERIFICATION_v0.3.0.md](https://github.com/gushamilton/GMoQ/blob/archive/research-records-20261016/docs/VERIFICATION_v0.3.0.md), [CHR22_VALIDATION.md](https://github.com/gushamilton/GMoQ/blob/archive/research-records-20261016/docs/CHR22_VALIDATION.md), [FINDINGS_pooled_dispersion.md](https://github.com/gushamilton/GMoQ/blob/archive/research-records-20261016/docs/FINDINGS_pooled_dispersion.md), [BACKLOG.md](https://github.com/gushamilton/GMoQ/blob/archive/research-records-20261016/docs/BACKLOG.md) |
