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

## What to worry about instead: a mean effect that manufactures a vQTL

Three to four orders of magnitude larger than anything above, and **two separate
mechanisms, not one**. An earlier version of this section attributed the whole effect
to the phenotype's mean-variance relationship and blamed the residual on
`W`-misspecification. Both attributions were wrong. The exact population calculation
is in `docs/MEAN_VARIANCE_ARTEFACT.md`; the numbers below supersede the ones that
stood here.

A SNP with a **pure mean effect and no variance effect at all**, N = 200,000,
MAF 0.30, on a lognormal trait of raw-scale skewness 1.75:

| scale | previously reported (20-rep MC) | exact | true dispersion effect |
|---|---|---|---|
| raw (lognormal) | 28.1 SD | **33.2** | **0.2341 — real, not spurious** |
| log | 2.7 SD | **3.12** | 0 (exactly) |
| RINT | 3.2 SD | **3.09** | -1.7e-4 |

The old figures came from 20 Monte Carlo replicates, which carry +-4.6 on 28.1 and
+-0.49 on 2.7 and 3.2. All three are consistent with the exact values; only the raw
one was materially imprecise.

**Correction 1: RINT does remove the mean-variance mechanism, completely.** Repeat
the calculation at MAF 0.50, where the genotype skewness `(1-2p)/sqrt(2p(1-p))` is
exactly zero:

| scale | MAF 0.30 | MAF 0.50 |
|---|---|---|
| raw | 33.2 | 34.4 |
| log | 3.12 | **2.8e-13** |
| RINT | 3.09 | **6.4e-14** |

Zero to machine precision. Everything that survives transformation is proportional to
the *genotype* skewness and disappears with it.

**Correction 2: the residual is not `W`-misspecification.** `W` does not know the
allele frequency, so it is equally misspecified in both columns above. The residual
is the genotype-skewness-times-density-curvature term derived independently on the
`sim/dispersion` branch, and their guard `G = 0.6855 * gskew * Z_mu^2 / sqrt(N)`
predicts it to within 0.4 per cent (3.110 predicted against 3.122 exact).

**Correction 3: the raw-scale number is not a bias.** A lognormal's conditional SD is
proportional to its conditional mean, so the true raw-scale variance effect really is
0.2341. Stage 2 returns 0.1405 — it *understates* the true raw-scale effect by 40 per
cent. The problem is the scale of analysis, not the estimator.

**The number that matters.** At N = 461,320, on a raw untransformed trait of skewness
1.75, a pure mean effect of **|Z_mu| = 9.5, about 0.02 SD per allele**, is enough to
push the dispersion parameter past 5e-8. At the milder skewness 0.61 that gives
Young et al.'s UK Biobank BMI coupling of `r_av` = 0.135, the threshold is
**|Z_mu| = 19.7**. Both are ordinary GWAS hits. After log or RINT the threshold moves
out to |Z_mu| = 43-93 and the genotype-skewness mechanism takes over.

Stage 2's GLS orthogonalises `theta_sigma2` against `theta_mu` *within* a variant.
That is a sampling-level statement and is unrelated to the across-variant population
coupling Young et al. remove with `d_l = alpha_vl - r_av * alpha_l`. Since this
project analyses raw `bmi` and `weight` alongside `rint_bmi`, that residualisation is
implemented in `stage2/R/dispersion-guards.R` as an extra output column, alongside
the genotype-skewness guard, which it does not and cannot fix.

## Context

QUAIL (Miao et al. 2022, PNAS) states this as an open problem: whether a mixed model
accounting for relatedness is compatible with a quantile-integral phenotype. The rank-one
argument above answers it - the constant direction *is* the quantile-integral direction.

Separately, Stage 1's sample-level block jackknife is a sandwich estimator, so it is robust
to exactly the heteroskedasticity a vQTL test looks for - more robust than REGENIE's
model-based SE. Iachine et al. 2010 use precisely a clustered Huber-White sandwich for this.

Working files: `/user/work/fh6520/scratch_rel/` (`exact.R`, `rel.R`, `pc.R`, `prof.R`).
