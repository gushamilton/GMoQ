# OLS and pooled covariance: implementation and validation

This report records a bounded development pass starting from commit `46cdbc0`.
OLS with pooled covariance is the primary file-based workflow; jackknife is an
explicit comparison. Execution ran on BluePebble compute nodes, with at most
two threads per concurrent worker. The initial checks used `bp1-compute001`
(allocation 18964919); subsequent registered checks used `bp1-compute053`
(18966436) and `bp1-compute019` (18966618). See the
[experiment register](../experiments/README.md) for per-run provenance. The
original cluster checkout and its staged work were preserved.

The user confirmed that pooled speed and retaining both covariance options
are priorities. The stress-case limitations below remain documented; the
optional full-nuisance bootstrap investigation is deferred. No claim that
pooled is preferable in a quantified fraction of all applications is supported
by these bounded experiments.

## Repairs

- Correct the high-bit BGEN dosage path using decoded genotype probabilities.
  Independent fixtures cover allele labels, missingness, fractional dosages,
  mixed record precision, sample ordering and the actual CLI. Duplicate sample
  IDs now fail explicitly.
- Make pooled covariance the file-wrapper default and discover its conventional
  sibling sidecar. Missing pooled inputs fail rather than silently switching
  covariance source.
- Preserve parameter p-value references in dispersion reporting; use available
  degrees of freedom for corrected statistics and return unavailable inference
  explicitly. Corrected statistics still use an approximate reference.
- Trim per-variant correction ratios by the ratio itself, and remove the silent
  naive-Q fallback.
- Correct the saved REGENIE allele-alignment claim and withdraw unsupported
  universal claims about covariate adjustment and relatedness.

## Reproduced checks

| Check | Observed result | What it establishes |
|---|---|---|
| Full Stage 1 Python suite | 126 passed, 19.04 seconds | Tested input, kernel and CLI behavior |
| Mathematical regression checks | 131 assertions passed | RIF, weights, covariance propagation, Q and reference calculations |
| Stage 1 to Stage 2 integration | 18 end-to-end, 45 family and 22 refusal assertions passed | File handoff, fitting and explicit error paths |
| Targeted Stage 2 tests | Passed | Covariance IDs, pooled handling, dispersion guards, t references, GLS Q, negative eigenvalues and namespace |
| R package build and check | Zero errors; one GNU Makevars portability warning | Installation, documentation, examples and vignettes; tests run separately |
| Saved empirical allele comparison | 17,550 matching allele pairs; correlation 0.99285 | Alignment of those saved log-BMI results |

The historical family checks use permissive regression thresholds. Their passing
does not certify bias, coverage or model selection. The validation suite now
labels them accordingly. The public BGEN fixtures are independent of the reader
used to analyse them; a phenotype simulated from the production reader alone
cannot expose a consistently reversed dosage convention.

## Real-genotype and empirical runs

A current-CLI check used 24 real chromosome-22 variants and 463,005 people,
45 quantiles and 65 blocks. A simulated mean effect of +0.25 at rs9617528
produced an OLS estimate of +0.250635, with pooled SE 0.003961. Stage 1 took
5.2 seconds on two threads and approximately 958 MB peak memory. This is a
bounded smoke test, not replicated coverage evidence. The fitted raw variance
effect was -0.01533 (pooled SE 0.00577); a successful mean recovery must not be
presented as proof of unbiased variance recovery.
An independent rerun produced identical decompressed Stage 1, jackknife covariance
and pooled sidecar contents. Across these 24 8-bit records, decoding took 0.12016
seconds with the original native path and 0.19432 seconds with the per-record
precision check, with zero dosage difference. The measured added cost was 0.07416
seconds on this panel; file compression timestamps were excluded from the
reproducibility comparison.

A fresh empirical log-BMI run used the same 24 variants, 461,320 matched people,
age, sex and ten PCs, the paper's keep file and LOCO predictions, and explicit
regression LOCO scaling. Stage 1 took 22.1 seconds and approximately 2.674 GB
peak memory; additive pooled Stage 2 took about 1.8 seconds. These are timings
for the stated small variant panel, not genome-wide throughput estimates.
Empirical effect truth is unknown, so this establishes execution and traceable
preprocessing rather than unbiasedness.
An exploratory two-parameter fit using the paper's cached BMI weights is
excluded from validated empirical evidence: those weights were built on raw
BMI without this log-BMI LOCO preprocessing. The corrected fit has now run
using the exact matched log-BMI baseline and serialized tau grid. Its
[provenance and aggregate result](../scripts/validation/real_chr22_empirical_provenance.md)
are preserved separately from that excluded fit.

## Statistical findings

Current-package reanalysis of completed independent simulations at N=461,320
used 15,000 datasets per selected cell and the exact OLS projection target.
Pooled null false-positive rates were about 0.05, interval coverage about 0.95,
and mean SE close to empirical SD at MAF 0.30 and 0.05. Common-variant location
and dispersion checks were also near nominal.

A strong dispersion signal at MAF 0.05 exposed a specific limitation. With
`Y=(1+c*g)*E`, Normal E and `c=0.0604686271168916`, pooled variance-effect
coverage was 0.92847 (Monte Carlo SE 0.00210), versus 0.94473 (0.00187) for
jackknife on the same datasets. The mean-SE/empirical-SD ratios were 0.92058
and 0.98597. Both use the same OLS estimate and target. Pooled covariance can
therefore understate uncertainty under genotype-dependent dispersion. Use
jackknife confirmation for strong dispersion signals, particularly at lower
allele frequencies; it does not repair a misspecified effect interpretation.

The 100-dataset, N=10,000 pilot separates scientific coefficients from population
RIF projection targets and compares oracle with estimated weights. Its centered
linear-variance cells are the relevant check of the Normal variance parameter.
Strong exponential-scale contrasts are separate stress tests; treating their
log-SD coefficient as the package's variance coefficient would be incorrect.

The covariate stress test demonstrates that a correlated, conditionally null
SNP can remain associated after an inadequate linear covariate adjustment.
Residualising the phenotype first is not a universal repair. This is an
estimand/nuisance-model issue as well as an uncertainty issue. The probe in
`tests/validity_pilot` compares a richer supplied covariate basis on paired data.
Adding C-squared and C-cubed removed most of the effect bias in that probe.
However, in its low-frequency correlated-SNP cell, pooled variance-effect
false-positive frequency remained 0.30 (MCSE 0.046), versus 0.09 (0.029) with
jackknife. The common-frequency pooled value was 0.05. Thus an adequate basis
does not itself validate pooled uncertainty under strong covariate-driven
heteroscedasticity. This remains an unresolved correctness limitation for that
setting, not a passing result or merely a signal-confirmation concern.

A follow-up [independent covariance diagnostic](../stage2/tests/robust_probe/README.md)
used 1,000 fresh datasets. Genotype-weighted HC0 substantially reduced the
correlated low-frequency variance-null rejection from .268 (pooled) to .070,
similar to .067 for jackknife. These comparisons use a fixed Normal cutoff
for all methods to diagnose covariance scaling; they are not a replay of the
package's full t-reference inference. The residual .070 rejection is still
above .05, so this does not yet certify a production repair. The experiment
does not change the package covariance path.

The strong centered-linear variance pilot also has a measurable finite-sample
projection bias: -0.00912 with MCSE 0.00218, about 0.42 empirical SD. It must
not be described as negligible relative to Monte Carlo uncertainty or hidden
by a coverage average across parameters. A subsequent
[paired density diagnosis](../tests/validity_pilot/evidence/DENSITY_DIAGNOSIS.md)
identified KDE smoothing as a cause: at v=.24, empirical KDE versus exact
density changed the estimate by -.003989 ± .000324 at N=10,000 and
-.001834 ± .000146 at N=100,000. Using population quantiles and exact density
gave negligible measured projection bias. The diagnostic uses fresh single-SNP
datasets, so its absolute bias need not equal the earlier six-SNP pilot.
This identifies a repair target; no density correction has yet been validated
or introduced into production.

See [the estimand contract](ESTIMAND.md) for units and interpretation. These
checks do not establish calibration at genome-wide tails, unbiased effects in
every phenotype regime, or validity of all mixture/shape-model selections.

## Reproduction

The checked-in `tests/validity_pilot` scripts record seeds, source and target
definitions. `tests/run_tier0.sh` and `tests/run_tier1.sh` exercise the mathematical
and file-integration checks. Stage 1's complete suite is under `stage1/tests`;
the individual Stage 2 regression scripts are under `stage2/tests`.

On BluePebble, run all commands inside a compute allocation. Use an isolated R
library, `Rscript --vanilla`, an explicitly selected Python interpreter, and the
recorded thread limits. Local development evidence and exact wrapper commands
are retained in the project's `review/` directory; aggregate simulation evidence
and data-specific smoke scripts accompany this report as they complete.

Preserved evidence:

- [Pilot and covariate probe](../tests/validity_pilot/evidence/VALIDITY_EVIDENCE.md)
- [Biobank-size covariance audit](../stage2/tests/validity_biobank_check/README.md)
- [Real-genotype and empirical scripts](../scripts/validation/README.md)
