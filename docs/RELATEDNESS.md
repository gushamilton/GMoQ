# Relatedness and structure: G-MoQ inherits REGENIE's control, and does better on dispersion

Status: derived and confirmed numerically 2026-09-14. **Not a concern at UK Biobank scale.**
The mean parameter carries exactly REGENIE's exposure; the dispersion parameter carries its
square. What does need attention is the mean-variance scale artefact, three to four orders
of magnitude larger - see the last section.

## The result

Let `rho` be the residual post-LOCO phenotypic correlation between relatives. Then

    lambda_mu - 1      is proportional to rho
    lambda_sigma2 - 1  is proportional to rho^2

Residual relatedness hits the mean at first order and dispersion only at second. Since
LOCO's whole job is to make `rho` small, squaring it settles the question.

## Why

Write the between-relative dependence as a shared additive component of variance `s2`.
To first order, for any such component (no normality needed):

    Cov(1{r_i <= q_t}, 1{r_j <= q_u})  ~  s2 * f(q_t) * f(q_u)

The RIF divides by `f(q_t)`, so

    Cov(RIF_i(t), RIF_j(u))  =  s2 * f(q_t) f(q_u) / (f(q_t) f(q_u))  =  s2

**The `1/f` in the RIF cancels the density exactly.** Additive relatedness therefore
contaminates `Cov(b_hat)` by a rank-one matrix along the constant vector, `s2 * 11'`,
identically at every tau. The nonlinearity of the RIF is not the problem here - it is the
reason the result holds, because `1/f` is exactly the Jacobian that cancels.

Stage 2's `W_mu(tau) = 1` exactly, so the projection sends that direction entirely to the
mean parameter and nothing to the others. This holds for every weight function the package
offers: in the mixture parameterisations `W_mu1 + W_mu2 = (p1 f1 + (1-p1) f2)/f = 1`, so the
constant vector is in the column space of `W` in all of them.

Exact computation (`/user/work/fh6520/scratch_rel/exact.R`):

| rho | A_mu | A_sigma2 | ratio |
|---|---|---|---|
| 0.30 | 0.30006 | 0.180 | 0.281 |
| 0.20 | 0.20003 | 0.0797 | 0.187 |
| 0.10 | 0.10002 | 0.0196 | 0.092 |
| 0.05 | 0.05001 | 0.00462 | 0.043 |
| 0.02 | 0.02001 | 0.00042 | 0.010 |

`A_mu = rho` to five digits - the mean parameter inherits exactly REGENIE's exposure, no
more and no less. `A_sigma2` scales as `rho^2` (0.180/0.0797 = 2.26 = 1.5^2), reproducing
the elementary identity `Corr(X^2, Y^2) = rho^2` for bivariate normals.

**UK Biobank**: ~22,000 full-sib pairs in ~460,000, generous post-LOCO `rho = 0.15`:

    lambda_mu - 1     = 0.0071     (REGENIE's own residual)
    lambda_sigma2 - 1 = 0.0010     (0.1%)

Per-tau statistics are inflated by `(lambda_mean - 1) * phi(z_t)^2 / [t(1-t)]`, maximised at
the median at exactly `2/pi = 0.637` and falling to 0.12 at tau = 0.02. Every per-tau
statistic is therefore inflated *less* than REGENIE's mean statistic.

## The covariate step is stronger than REGENIE's, and the order of operations is why

Stage 1 residualises the RIF matrix on covariates **column by column** (`cli.py:667-677`:
`Q.T @ rif` is K x T, a separate coefficient vector per tau). The PC adjustment therefore has
**40 free coefficients at each of 45 quantiles - 1,800 parameters, against REGENIE's 40.**
Those extra degrees of freedom are exactly what describes how the *shape* of the conditional
distribution varies with ancestry: a location shift across strata has a flat tau-profile, a
scale shift an odd `q_tau` profile, and each tau having its own coefficient absorbs both.

This depends on the order. Residualising `Y` on PCs and *then* building the RIF would give
one coefficient set, not 45, and would not absorb second-moment structure. **The pipeline's
order is the correct one.**

Stress test (`pc.R`; sd varying as `exp(0.15 * PC)`, SNP with `cor(g, PC) = 0.05`, no true effect):

| adjustment | spurious `theta_sigma2` bias |
|---|---|
| none | 1.85 SD |
| PC | 0.03 SD |
| PC + PC^2 | 0.03 SD |

A factor of 60, to zero within Monte Carlo error, with the SE essentially unchanged
(0.0068 -> 0.0065) - genuine bias removal, not SE inflation.

## A genuine second-moment shared component

The mechanism exists even though the magnitude does not. If relatives share a *variance*
level rather than a mean, `Cov(RIF_i(t), RIF_j(u)) = s2_d * q_t * q_u` - rank-one along `q`,
not along `1`. Since `q = mu * W_mu + 2 sigma^2 * W_sigma2` lies in the column space of `W`,
that does load onto the vQTL parameter. Simulation confirms it: the leading eigenvector
correlates 0.9988 with `q_tau`, and 100.3% of the induced variance lands on `beta_sigma2`,
with the diagonal contribution falling to 0.0009 at the median against a predicted 0.0000.

It does not bite in UK Biobank because only ~9.6% of the sample sits in related clusters and
the exposure goes as `rho^2`. **It would bite in a fully clustered sample** - ALSPAC
mother-child pairs, or an admixed cohort. Shriner et al. 2025 found admixture explained ~53%
of a chr2 vQTL signal despite PC1 plus a GRM. This analysis does not transfer to those.

## A free diagnostic

The two mechanisms have opposite tau-profiles, so the *shape* of `lambda_GC(tau)` identifies
which is present - no LD reference, no assumption about polygenicity:

| tau | 0.02 | 0.10 | 0.25 | 0.50 | 0.75 | 0.90 | 0.98 |
|---|---|---|---|---|---|---|---|
| additive relatedness | 0.120 | 0.356 | 0.547 | **0.637** | 0.547 | 0.356 | 0.120 |
| second moment | 0.831 | 0.904 | 0.373 | **0.000** | 0.373 | 0.904 | 0.831 |

Inverted-U peaking at the median means ordinary relatedness, confined to `beta_mu`, harmless.
U-shaped and zero at the median means genuine second-moment contamination. Compute
`lambda_GC(tau)` per quantile from Stage 1 output already on disk and regress it on these two
basis vectors. Being a shape test rather than a level test, it separates confounding from
real signal without LDSC - which matters, since nobody has justified applying mean-effect
LDSC to variance statistics.

## What to worry about instead: the mean-variance scale artefact

Three to four orders of magnitude larger than anything above. A SNP with a **pure mean effect
and no variance effect at all**, N = 200,000:

| scale | spurious `theta_sigma2` bias |
|---|---|
| raw (lognormal) | **28.1 SD** |
| log-transformed | 2.7 SD |
| RINT | 3.2 SD |

**RINT does not remove it**, independently reproducing Young et al. 2018: variance effects
arising from a mean-variance relationship are not in general removed by an inverse normal
transform. The residual after the correct transformation is `W`-misspecification: the
Gaussian `W` is built on a marginal that is a mixture over genotypes, so any departure from
Gaussianity couples mean into variance.

Stage 2's GLS orthogonalises `theta_sigma2` against `theta_mu` *within* a variant under an
assumed shape. That is not the same as Young's genome-wide residualisation
`d_l = alpha_vl - r_av * alpha_l` (they estimate `r_av = 0.135` for UK Biobank BMI). Since
this project analyses raw `bmi` and `weight` alongside `rint_bmi`, **that residualisation is
worth adding.**

## Context

QUAIL (Miao et al. 2022, PNAS) states this as an open problem: whether a mixed model
accounting for relatedness is compatible with a quantile-integral phenotype. The rank-one
argument above answers it - the constant direction *is* the quantile-integral direction.

Separately, Stage 1's sample-level block jackknife is a sandwich estimator, so it is robust
to exactly the heteroskedasticity a vQTL test looks for - more robust than REGENIE's
model-based SE. Iachine et al. 2010 use precisely a clustered Huber-White sandwich for this.

Working files: `/user/work/fh6520/scratch_rel/` (`exact.R`, `rel.R`, `pc.R`, `prof.R`).
