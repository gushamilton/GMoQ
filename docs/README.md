# Documentation index

Thirty-odd files live here and they are not all for the same reader. Start in
the right column.

## If you are going to RUN G-MoQ

Read these four, in order:

| | |
|---|---|
| [ESTIMAND.md](ESTIMAND.md) | **What a G-MoQ effect means.** Unconditional (RIF) quantile effects, not conditional ones. Anything subtracted from Y redefines the estimand. Read this first or you will misinterpret the output. |
| [phenotype_preparation.md](phenotype_preparation.md) | How to prepare the phenotype. The estimand is defined against its unconditional distribution, so this is not boilerplate. |
| [ALLELE_CONVENTION.md](ALLELE_CONVENTION.md) | Which allele the effect is counted against. Get this wrong and every sign flips. |
| [KNOWN_ISSUES.md](KNOWN_ISSUES.md) | What the package does not do well, stated plainly. |

## If you hit a specific question

| Question | File |
|---|---|
| Which covariance source, and when does it matter? | [FINDINGS_covsource_default.md](FINDINGS_covsource_default.md) |
| Why is the estimator OLS and not GLS? | [OLS_VS_GLS.md](OLS_VS_GLS.md), [WHY_GLS_DEGRADES.md](WHY_GLS_DEGRADES.md) |
| How should I handle LOCO / REGENIE predictions? | [LOCO_AS_COVARIATE.md](LOCO_AS_COVARIATE.md), [LOCO_SCALE.md](LOCO_SCALE.md), [LOCO_MIXTURE_HAZARD.md](LOCO_MIXTURE_HAZARD.md), [TWO_PASS_LOCO.md](TWO_PASS_LOCO.md) |
| My phenotype is lumpy / has a point mass | [BP2_lumpy_phenotypes.md](BP2_lumpy_phenotypes.md) |
| I found a variance effect -- is it real? | [MEAN_VARIANCE_ARTEFACT.md](MEAN_VARIANCE_ARTEFACT.md), [DISPERSION_SKEW_CORRECTION.md](DISPERSION_SKEW_CORRECTION.md) |
| How does this behave as N grows? | [N_SCALING.md](N_SCALING.md), [EFFICIENCY_BOUND.md](EFFICIENCY_BOUND.md) |
| Related samples? | [RELATEDNESS.md](RELATEDNESS.md) |
| How fast is Stage 1, and which kernel runs? | [STAGE1_FAST_KERNEL.md](STAGE1_FAST_KERNEL.md) |
| Can I run several phenotypes in one pass? | [MULTI_PHENOTYPE_FUSION.md](MULTI_PHENOTYPE_FUSION.md) (answer: no, it is slower) |

## If you are developing or reviewing the package

| | |
|---|---|
| [BACKLOG.md](BACKLOG.md) | Open defects and decisions, with status. |
| [TEST_SUITE.md](TEST_SUITE.md) | What the tiered simulation suite covers. |
| [SIMULATION_INVENTORY.md](SIMULATION_INVENTORY.md) | Every simulation campaign: what it asked, what it concluded, whether that still holds under the 0.4.0 defaults, and where the archived ones went. Read this before reusing a number from `sim/`. |
| [COMPARATOR_LEDGER.md](COMPARATOR_LEDGER.md) | Every G-MoQ-vs-REGENIE number, and what a fair comparator is. |
| [POWER_AND_CONSERVATISM.md](POWER_AND_CONSERVATISM.md) | What the LOCO change bought; null calibration at 5e-8. |
| [EMPIRICAL_W.md](EMPIRICAL_W.md) | Where the shipped weight matrices come from. |

## Audit trail -- point-in-time records, not current guidance

These say what was true on a date. They are kept for provenance and are **not**
maintained against the current code: [AUDIT_2026-09-16.md](AUDIT_2026-09-16.md),
[VALIDATION_2026-09-15.md](VALIDATION_2026-09-15.md),
[VERIFICATION_v0.3.0.md](VERIFICATION_v0.3.0.md),
[CHR22_VALIDATION.md](CHR22_VALIDATION.md),
[FINDINGS_pooled_dispersion.md](FINDINGS_pooled_dispersion.md) (section 8 item 1
is withdrawn; see FINDINGS_covsource_default.md),
[PAPER_FIGURES_AND_TABLES.md](PAPER_FIGURES_AND_TABLES.md).

When one of these disagrees with the README or with the four "going to run it"
documents above, the README wins.
