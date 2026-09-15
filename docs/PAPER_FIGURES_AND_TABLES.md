# G-MoQ paper: figures and supplementary evidence

Status: agreed scope translated into a production plan, 15 September 2026.
Figures and tables below are planned outputs, not newly generated results.

The paper should establish three things: the implementation computes the
intended estimator, the estimator behaves as claimed in defined scenarios,
and the package makes distributional analysis practical. OLS with pooled
covariance is the primary analysis; jackknife remains an uncertainty sensitivity.
Their point estimates are identical when all other settings match.

## Main figures

| ID and proposed title | Panels and question | Evidence and remaining work |
|---|---|---|
| Figure 1. Reconstruction of mean genetic effects | Matched G-MoQ versus direct OLS and REGENIE effect plots; identity lines; differences against direct OLS by MAF; SE agreement. Report slope, intercept and absolute error as well as correlation. | Saved 17,550-variant comparison has r=.99285 for mean quantile coefficients versus REGENIE. Audit its provenance and reconstruct current Stage 2 results before using that label. Add independent direct OLS. Match samples, phenotype, covariates, alleles and LOCO scale. |
| Figure 2. Recovery and calibration under known truth | Null rejection, standardized bias, mean SE/empirical SD, and 95% coverage across null, location and dispersion scenarios at MAF .30/.05. Reference lines and Monte Carlo intervals. | Reuse six verified biobank cells, N=461,320 and 15,000 datasets/cell. Use exact OLS projection targets; list generating coefficients separately. Include pooled and jackknife uncertainty, including the low-frequency dispersion limitation. |
| Figure 3. Distributional information beyond the mean | Genotype-specific distributions, quantile-effect profiles and recovered mean/variance effects for pure location and pure dispersion; power curves across effect size. | Reuse validated scenarios. Inventory completed weaker-effect outputs before any new runs. REGENIE tests the mean; it is not a competing variance test. A superiority claim over variance methods would require an additional matched comparator and is outside this core scope. |
| Figure 4. Computational performance | Stage-specific wall time and peak RSS against variant count; pooled versus jackknife Stage 2; compact covariance storage. | Short matched benchmark on identical artifacts, host and threads with repeated timings. Existing 492-fold reduction concerns a covariance object only. Stage 1 currently computes both covariance representations, so do not attribute Stage 1 savings to pooled Stage 2. |

G-MoQ's finite-grid mean reconstruction is not assumed to equal direct OLS
exactly. Figure 1 measures that approximation. Empirical agreement cannot
establish unbiasedness because empirical truth is unknown.

## Supplementary figures

| ID and proposed title | Content |
|---|---|
| Figure S1. Independent verification of Stage 1 RIF regression | Independently construct RIF outcomes and fit ordinary least squares at each tau. First hold quantiles, densities, sample mask and nuisance treatment fixed to isolate regression correctness; then compare independently estimated quantiles/densities under matched conventions. Plot coefficient agreement and absolute error versus tau, with missingness/covariate cases. Check pooled covariance and delete-block reference separately under the same definitions; do not conflate different standard-error estimators. |
| Figure S2. Separation of location and dispersion | MAF-by-location-effect display of spurious dispersion and correction behaviour, with genuine dispersion retained as a sensitivity panel. Distinguish analytic calculations from simulation and audit inherited evidence against the current source. |
| Figure S3. Statistical boundaries | Compact panels for rare-variant covariate-dependent variance and density-estimation bias. Show pooled and jackknife limitations. Label HC0 and fourth-order KDE as diagnostic prototypes, not shipped repairs. |

Stage 1 deserves this small independent check because it supplies every
downstream result. It does not need a broad UQR tournament. The direct comparator
is conventional RIF–OLS UQR with the same nuisance definitions. Conditional
quantile regression targets a different quantity and is not a numerical
equivalence reference. A named external UQR implementation can be an optional
confirmation once its documentation, defaults and estimand have been checked;
it should not become a dependency for finishing this figure set.

## Supplementary tables

1. **Table S1. Design and provenance.** DGP, sample size, independent dataset
   count, MAF, effect settings, tau grid, blocks, covariates, density method,
   exact projection target and generating coefficient, source commit, seed,
   input identifiers and experiment ID. Include the Stage 1 comparator design.
2. **Table S2. Complete statistical results.** Bias and Monte Carlo SE,
   empirical SD, mean reported SE, coverage, null rejection, power and failure
   counts, for every reported scenario and method. State the uncertainty
   reference used. Include numerical Stage 1 agreement metrics.
3. **Table S3. Empirical and computational reproducibility.** Comparator
   matching details, effect/SE agreement, representative empirical results,
   software versions, host/resources, timings, peak RSS, output hashes and
   artifact locations. Distinguish object storage from process memory.

## Production and presentation

- First collate the completed evidence in [the experiment register](../experiments/README.md).
  Configured cells are not completed cells. Do not repeat an existing campaign.
- New work is bounded to independent Stage 1/direct OLS checks, resolving the
  matched empirical comparison, fair timings, and necessary power-curve gaps.
  Do not restart the deferred nuisance bootstrap.
- Register and commit scripts/configuration before computation. Use unique
  output directories and preserve failed runs. All R/Python, plotting and
  analyses run on BP compute nodes, with one coordinated allocation.
- Use Luna xhigh workers for bounded preparation and implementation tasks.
- Produce vector PDF and SVG, high-resolution PNG previews, source CSVs,
  editable plotting scripts and standalone captions. Use consistent method
  colours and panel lettering; show reference lines and Monte Carlo uncertainty.
- Each caption specifies the estimand, independent replication unit, sample
  size, inference method and source experiment. Do not treat correlated SNPs
  as independent simulated datasets. Fifteen thousand null datasets cannot
  validate a 5e-8 false-positive rate empirically.
- A figure is ready when every plotted value maps to preserved evidence, its
  caption states relevant limitations, and the rendering has been inspected.
  Statistical limitations remain findings; plotting does not repair them.

The bounded target is **four main figures, three supplementary figures and
three supplementary tables**. The additional supplement relative to the initial
proposal is the independent Stage 1 check.
