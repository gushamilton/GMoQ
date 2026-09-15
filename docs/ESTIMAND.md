# What the primary analysis estimates

The primary file-based analysis uses OLS with pooled covariance. Jackknife
covariance remains an explicit comparison. This choice concerns uncertainty;
with the same quantile coefficients and model weights, the OLS point estimate
is the same under both covariance sources.

Stage 1 constructs quantile recentered influence functions (RIFs) from the
analysed phenotype, removes the supplied covariate basis from each RIF and
genotype, and estimates one genotype coefficient per quantile. Genotypes count
the first BGEN allele; see [the allele convention](ALLELE_CONVENTION.md).
LOCO adjustment and phenotype transformations change the analysed outcome and
must be recorded with the result.

Stage 2 projects that coefficient vector onto the selected model's weight
matrix. For OLS this is `(W'W)^(-1) W' beta`. Its covariance is obtained by
applying the same linear transformation to the supplied quantile covariance.
The pooled path uses the shared residual RIF covariance divided by the
variant's residual genotype sum of squares. It does not estimate arbitrary
genotype-dependent residual covariance separately at every variant.

Use jackknife covariance as a confirmation for strong dispersion signals,
particularly at lower allele frequencies. A paired 15,000-dataset Normal
simulation at N=461,320 and MAF=0.05 gave 92.85% coverage with pooled versus
94.47% with jackknife for a nominal 95% interval around the same OLS dispersion
projection target. Their point estimates were identical. This is evidence of
a pooled uncertainty limitation under that heteroscedastic signal; changing
covariance does not repair confounding or change the effect's interpretation.
Pooled inference also remains unvalidated under strong covariate-driven
heteroscedasticity correlated with genotype: the richer-basis probe removes
effect bias but still shows inflated pooled false positives at low frequency.
Use the jackknife comparison and the setting-specific validation report before
interpreting pooled p-values in that regime.

The current jackknife also holds the fitted quantiles and density fixed. It
therefore does not propagate all nuisance-estimation uncertainty for nonzero
effects. A lower-bias density estimate alone does not close that gap; see the
[RIF uncertainty derivation](../stage2/docs/DENSITY_RIF_NUISANCE.md). The
jackknife comparison above is useful evidence, not a universal calibration
guarantee.

For the Normal mean/variance model, `beta_mu` is on the outcome scale and
`beta_sigma2` is on the squared outcome scale per counted allele. Interpretation
as a conditional mean or variance coefficient is a local model interpretation,
not an identity for every finite genotype contrast. For example, generating
`SD(Y|g) = exp(b*g)` does not make the reported variance coefficient equal to
`b`: the local variance derivative at zero is `2*b`, and a one-allele variance
ratio is `exp(2*b)`. Strong effects can depart from the linear approximation.

Validation therefore reports both the population RIF projection target and
the scientific coefficient used to generate a phenotype. Null calibration,
agreement with another program, and a non-significant lack-of-fit Q statistic
do not by themselves establish unbiased recovery of the scientific coefficient.

Covariate adjustment requires an adequate basis for the relevant conditional
expectations. Separate linear coefficients at each quantile do not guarantee
control of nonlinear mean and variance confounding. The paired pilot tests
this explicitly. Tied phenotypes, strong effects, relatedness and alternative
model families need their own evidence; results from the smooth, small-effect
location/variance setting should not be generalized to these automatically.

The reproducible pilot in `tests/validity_pilot` is a diagnostic with independent
datasets as its Monte Carlo unit. Its report distinguishes scientific and
projection bias, coverage, covariance source and estimated versus oracle
weights. A small pilot is not evidence of calibration at genome-wide tails.
