# Power bought by the LOCO change, and the conservatism that is not there

Two questions the project had been asserting without measuring.

1. The 0.3.0 Stage 1 default moved the REGENIE LOCO prediction out of the
   phenotype and into the residualising basis
   ([`LOCO_AS_COVARIATE.md`](LOCO_AS_COVARIATE.md)). Median SE against a fair
   sex-stratified LOCO'd REGENIE went 1.144 -> 1.044 (`beta_left`) and
   1.113 -> 0.970 (`beta_right`). **Did that convert into discoveries?**
2. A permuted-phenotype null rejected 4.2% at alpha = 0.05 and 0.04% at
   alpha = 1e-3, read as "about 2.5x conservative at the stricter threshold".
   **Is G-MoQ leaving power on the table at 5e-8?**

Short answers. **(1) Yes, but in correlated variants, not in loci.** Every
discovery-status change on chr21 and chr22 is a gain -- 35 variants gained and
0 lost at 5e-8, 135 gained and 5 lost at 1e-5 -- and recall against the
stratified GWAS at 5e-8 goes 0.49 -> 1.00 (`beta_left`) and 0.82 -> 0.97
(`beta_right`). At the locus level the same evidence is one extra locus at
5e-8 and six at 1e-5. **(2) No.** At the production configuration the pooled/OLS
parameter test rejects at 0.997 to 1.026 of nominal at 5e-8 (+/- 0.003),
measured by saddlepoint-tilted importance sampling on the real influence
functions -- and where it departs it is *anti*-conservative, at rare variants.
The permuted-null "2.5x" is Monte Carlo noise: with LD, 17,550 chr22 HapMap3
variants carry only 1,200-6,000 independent tests, and the same permuted run
is simultaneously 3.4x *anti*-conservative on `beta_right` at 1e-4.

Everything below is testosterone, N = 399,710, 45 taus on (0.02, 0.98), skew
weights, pooled covariance, OLS projection -- the production configuration.
Run log: [`experiments/power-and-conservatism-20260916-01`](../experiments/power-and-conservatism-20260916-01/README.md).

---

## Part 1 -- power

### The comparator

Sex-stratified REGENIE with its own within-stratum LOCO: female N = 199,572,
male N = 200,138, sum 399,710 = G-MoQ's N, the same people. `beta_left` is read
against the female scan, `beta_right` against the male, `gamma` against a
REGENIE `sex` GWAS on the union with `gamma` negated (REGENIE codes Male = 1).
Alleles harmonised against REGENIE `ALLELE1`; on these inputs all 34,677
variants were already aligned, so no sign flips were needed and sign
concordance at every REGENIE hit is 1.000 in every cell.

### What changed in the estimates

| chr | param | SE ratio, subtract | SE ratio, covariate | median &#124;dbeta&#124; in SE | max &#124;dbeta&#124; in SE |
|---|---|---|---|---|---|
| 22 | `beta_left` | 1.1439 | 1.0445 | 0.291 | 1.80 |
| 22 | `beta_right` | 1.1132 | 0.9705 | 0.342 | 2.06 |
| 22 | `gamma` | 0.9819 | 0.9654 | 0.125 | 0.81 |
| 21 | `beta_left` | 1.1442 | 1.0446 | 0.294 | 1.74 |
| 21 | `beta_right` | 1.1144 | 0.9711 | 0.351 | 1.95 |
| 21 | `gamma` | 0.9818 | 0.9653 | 0.132 | 0.69 |

chr21 reproduces chr22's SE ratios to three decimals, which is worth more than
it looks: the SE effect is a property of the estimator, not of one chromosome's
signal. The point estimates moved at the same time -- a third of an SE at the
median, over two SEs at the maximum -- so nothing below can be read as a pure
SE rescaling.

### Discovery table

`gmoq` and `regenie` are variant counts below the threshold; `loci` is a greedy
500 kb clump, the crudest usable independence proxy.

**chr22, 17,550 HapMap3 variants**

| param | thr | mode | G-MoQ | REGENIE | overlap | precision | recall | G-MoQ loci | REGENIE loci |
|---|---|---|---|---|---|---|---|---|---|
| `beta_left` | 5e-8 | subtract | 20 | 41 | 20 | 1.000 | 0.488 | 1 | 1 |
| `beta_left` | 5e-8 | **covariate** | **42** | 41 | 41 | 0.976 | **1.000** | 1 | 1 |
| `beta_left` | 1e-5 | subtract | 46 | 85 | 45 | 0.978 | 0.529 | 2 | 5 |
| `beta_left` | 1e-5 | **covariate** | **107** | 85 | 76 | 0.710 | **0.894** | **4** | 5 |
| `beta_right` | 5e-8 | subtract | 27 | 33 | 27 | 1.000 | 0.818 | 3 | 3 |
| `beta_right` | 5e-8 | **covariate** | **39** | 33 | 32 | 0.821 | **0.970** | 3 | 3 |
| `beta_right` | 1e-5 | subtract | 59 | 82 | 42 | 0.712 | 0.512 | 7 | 8 |
| `beta_right` | 1e-5 | **covariate** | **116** | 82 | 66 | 0.569 | **0.805** | **9** | 8 |
| `gamma` | 5e-8 | both | 0 | 0 | 0 | -- | -- | 0 | 0 |
| `gamma` | 1e-5 | subtract | 3 | 1 | 0 | 0.000 | 0.000 | 1 | 1 |
| `gamma` | 1e-5 | covariate | 0 | 1 | 0 | -- | 0.000 | 0 | 1 |

**chr21, 17,127 variants**

| param | thr | mode | G-MoQ | REGENIE | overlap | precision | recall | G-MoQ loci | REGENIE loci |
|---|---|---|---|---|---|---|---|---|---|
| `beta_left` | 5e-8 | subtract | 0 | 5 | 0 | -- | 0.000 | 0 | 1 |
| `beta_left` | 5e-8 | **covariate** | **1** | 5 | 1 | 1.000 | **0.200** | **1** | 1 |
| `beta_left` | 1e-5 | subtract | 4 | 30 | 4 | 1.000 | 0.133 | 1 | 4 |
| `beta_left` | 1e-5 | **covariate** | **19** | 30 | 18 | 0.947 | **0.600** | **3** | 4 |
| `beta_right` | 5e-8 / 1e-5 | both | 0 | 0 | 0 | -- | -- | 0 | 0 |
| `gamma` | 1e-5 | both | 1 | 0 | 0 | 0.000 | -- | 1 | 0 |

### How many variants change discovery status

| chr | param | thr | subtract | covariate | gained | lost | McNemar p | SE-only | beta-only | gains that are REGENIE hits |
|---|---|---|---|---|---|---|---|---|---|---|
| 22 | `beta_left` | 5e-8 | 20 | 42 | 22 | 0 | 4.8e-7 | 36 | 35 | 21 / 22 |
| 22 | `beta_left` | 1e-5 | 46 | 107 | 61 | 0 | 8.7e-19 | 86 | 53 | 31 / 61 |
| 22 | `beta_right` | 5e-8 | 27 | 39 | 12 | 0 | 4.9e-4 | 39 | 29 | 5 / 12 |
| 22 | `beta_right` | 1e-5 | 59 | 116 | 59 | 2 | 1.6e-15 | 137 | 52 | 24 / 59 |
| 22 | `gamma` | 1e-5 | 3 | 0 | 0 | 3 | 0.25 | 3 | 0 | 0 / 0 |
| 21 | `beta_left` | 5e-8 | 0 | 1 | 1 | 0 | 1.00 | 0 | 0 | 1 / 1 |
| 21 | `beta_left` | 1e-5 | 4 | 19 | 15 | 0 | 6.1e-5 | 29 | 10 | 14 / 15 |

Totals across both chromosomes and all three parameters: **35 gained and 0 lost
at 5e-8; 135 gained and 5 lost at 1e-5.** The five losses are two `beta_right`
variants and the three `gamma` variants, and the three `gamma` losses are a
gain in disguise: none of them overlaps the single chr22 REGENIE `sex` hit at
1e-5, so covariate mode is dropping three non-replicating calls. The two
"counterfactual" columns hold one half of the change fixed:
`SE-only` recomputes p from the subtract-mode point estimate with the
covariate-mode standard error, `beta-only` does the reverse.

**Was it the SEs or the betas?** Both, and they are not additive.

- `beta_left` at 5e-8 on chr22: baseline 20. SE alone would give 36 (+16); the
  betas alone would give 35 (+15); together they give 42 (+22). The two effects
  overlap heavily -- the same variants are helped by each -- so the joint gain
  is well under the sum.
- `beta_right` at 1e-5 on chr22 runs the other way: the SE improvement alone
  would have delivered **137** hits, and the simultaneous movement of the point
  estimates pulled the realised count back to 116. Roughly a sixth of the SE
  gain was absorbed by the betas moving.

So the answer to "did the SE improvement convert into discoveries, or was it
absorbed by the point estimates moving" is: mostly converted, but not cleanly,
and the two channels contribute comparable amounts.

### Is the precision drop inflation?

No. Every G-MoQ-only hit is a variant REGENIE also likes, just below its own
threshold.

| chr | param | thr | G-MoQ-only hits | their median REGENIE p | fraction with REGENIE p < 1e-3 |
|---|---|---|---|---|---|
| 22 | `beta_left` | 5e-8 | 1 | 3.9e-7 | 1.00 |
| 22 | `beta_left` | 1e-5 | 31 | 3.0e-5 | 1.00 |
| 22 | `beta_right` | 5e-8 | 7 | 5.1e-7 | 1.00 |
| 22 | `beta_right` | 1e-5 | 50 | 4.2e-5 | 1.00 |
| 21 | `beta_left` | 1e-5 | 1 | 1.5e-5 | 1.00 |

Under the null 0.05 of these would sit below REGENIE p = 1e-3; 1.00 do. The
precision column in the discovery table falls because the reference set is a
thresholded GWAS, not because G-MoQ started producing junk.

Genomic control on the real scans points the same way: covariate mode tracks
the stratified GWAS, subtract mode sits below it.

| chr | param | lambda subtract | lambda covariate | lambda REGENIE |
|---|---|---|---|---|
| 22 | `beta_left` | 1.353 | 1.420 | 1.444 |
| 22 | `beta_right` | 1.360 | 1.486 | 1.515 |
| 22 | `gamma` | 1.267 | 1.224 | 1.233 |
| 21 | `beta_left` | 1.174 | 1.177 | 1.195 |
| 21 | `beta_right` | 1.015 | 1.173 | 1.197 |
| 21 | `gamma` | 1.069 | 1.058 | 1.059 |

(These are real scans with real signal, so they are power indices, not null
calibration.)

### How weak is this instrument

Weak, and it must be said plainly.

- A block bootstrap over contiguous 250-SNP blocks puts a 95% interval on the
  variant gain of **[0, 88]** for `beta_left` at 5e-8 on chr22 (point 22),
  **[0, 33]** for `beta_right` (point 12), **[12, 128]** and **[10, 120]** at
  1e-5, **[0, 3]** and **[0, 42]** on chr21. The *paired* evidence (96 gains,
  5 losses, McNemar p down to 9e-19) is far stronger than the count intervals,
  because the pairing removes the between-locus variance that dominates them.
- At the locus level the gain at 5e-8 is **one locus** (chr21 `beta_left`,
  0 -> 1; chr22 gains no locus in either component). At 1e-5 it is six
  (`beta_left` 2 -> 4 on chr22 and 1 -> 3 on chr21, `beta_right` 7 -> 9 on
  chr22). The 2.1x variant count at 5e-8 is one already-detected region
  pushing more of its correlated variants over the line.
- Two chromosomes of one trait. `gamma` has no signal at all here (the REGENIE
  `sex` scan finds nothing on either chromosome at 5e-8), so the change is
  untested for that component.

### Recall against what the SE ratio alone predicts

Take each REGENIE hit's z, divide by the median SE ratio, and ask whether it
still clears the threshold. That is the recall a pure SE rescaling would buy.

| chr | param | thr | n REGENIE | pred. subtract | obs. subtract | pred. covariate | obs. covariate |
|---|---|---|---|---|---|---|---|
| 22 | `beta_left` | 5e-8 | 41 | 0.634 | 0.488 | 0.854 | **1.000** |
| 22 | `beta_left` | 1e-5 | 85 | 0.529 | 0.529 | 0.765 | **0.894** |
| 22 | `beta_right` | 5e-8 | 33 | 0.788 | 0.818 | 1.000 | 0.970 |
| 22 | `beta_right` | 1e-5 | 82 | 0.561 | 0.512 | 1.000 | **0.805** |
| 21 | `beta_left` | 5e-8 | 5 | 0.000 | 0.000 | 0.200 | 0.200 |
| 21 | `beta_left` | 1e-5 | 30 | 0.267 | 0.133 | 0.700 | 0.600 |

The SE ratio *under*-predicts `beta_left` and *over*-predicts `beta_right` at
1e-5. Confirmation that the estimand moved as well as the precision, in
different directions for the two components.

---

## Part 2 -- conservatism

### Why real data cannot answer this

The existing permuted-phenotype run (testosterone shuffled across 421,405
people, chr22, N = 399,710, no LOCO -- `/user/work/fh6520/gmoq_perm`) gives:

| param | lambda (+/- block) | a = 0.05 | a = 1e-3 | a = 1e-4 | implied independent tests |
|---|---|---|---|---|---|
| `gamma` | 0.999 +/- 0.045 | 0.909 +/- 0.104 | 0.342 +/- 0.174 | 0.00 | 1,700-33,000 |
| `beta_left` | 0.904 +/- 0.045 | 0.836 +/- 0.080 | 0.399 +/- 0.171 | 0.570 +/- 0.576 | 3,000-34,000 |
| `beta_right` | 0.993 +/- 0.048 | 0.912 +/- 0.104 | 1.140 +/- 0.615 | **3.419 +/- 3.479** | 830-2,600 |

Cells are observed / expected. Errors are a 2,000-replicate block bootstrap over
100 contiguous genomic blocks, which is what turns "17,550 variants" into the
1,200-6,000 effectively independent tests LD actually leaves. Three things
follow:

1. Every deficit is inside 2.1 block-bootstrap SEs of 1.0. Nothing here is
   significant.
2. The three parameters are not independent evidence -- they are three
   near-collinear contrasts of the same 45 Stage 1 slopes on the same
   genotypes.
3. The run is conservative on two parameters and **anti**-conservative on the
   third (`beta_right`, 20 observed against 17.6 at 1e-3, 6 against 1.8 at
   1e-4). A single number "2.5x conservative" averages over opposite signs.

Extrapolating the median-based scale factor to 5e-8 gives a rejection-rate
factor of 0.98 (`gamma`), 0.20 (`beta_left`), 0.89 (`beta_right`) -- a 5x
spread. That is the measure of how little the permuted null constrains the tail.

### What the test actually is

`validate2.R` reproduces the shipped Stage 2 pooled/OLS output from
`(A, Sigma_0, b, d_resid)` to **max relative error 5e-15** on `theta`, `se` and
`p` for all three parameters over all 17,550 chr22 variants. With
`b = Psi_r' g / D` and `A = (W'W)^-1 W'`, that reduction is exact:

```
t_k  =  (u_k' g) / ( sqrt(D) * sd_u_k ),
u_k  =  Psi_r a_k,      sd_u_k^2 = a_k' Sigma_0 a_k = sum_i u_ki^2 / N
```

A single weighted genotype sum. `Psi_r` was pulled out of the production Stage 1
code path by hooking `compute_chr_rif_residuals` on the real inputs, and the
rebuilt `Sigma_0 = Psi_r' Psi_r / N` matches the shipped `sigma0.tsv.gz` to
4e-16, so `u` is the real per-individual contrast score, not a model of one.

### Simulation

Real `Psi_r`, real `W`, real MAF spectrum, N = 399,710, production tau grid and
projection. Genotypes iid HWE, `D` computed from the realised genotype
residualised on the same basis `Q` Stage 1 uses.

**Plain Monte Carlo, M = 5.0e6 independent null variants.** Ratios are
observed / nominal; errors are Poisson.

| param | Var(t) | mean(t) | kurt(t) | a = 0.05 | a = 1e-3 | a = 1e-5 |
|---|---|---|---|---|---|---|
| `gamma` | 1.00016 | 0.0005 | 2.998 | 1.0009 +/- 0.0020 | 1.0004 +/- 0.0141 | 1.02 +/- 0.14 |
| `beta_left` | 1.00016 | -0.0008 | 3.000 | 1.0010 +/- 0.0020 | 1.0020 +/- 0.0142 | 0.98 +/- 0.14 |
| `beta_right` | 1.00000 | 0.0002 | 3.001 | 1.0009 +/- 0.0020 | 1.0092 +/- 0.0142 | 1.04 +/- 0.14 |

5.0e6 replicates expect 0.25 events at 5e-8, so plain MC cannot reach it. (It
saw 2 for `beta_left`; P(>=2 | 0.25) = 2.6%, and it is the *anti*-conservative
direction.)

**Exponentially tilted importance sampling**, M = 5e5 tilted draws per cell,
tilt at the saddlepoint of the exact CGF. This is a direct, unbiased Monte
Carlo estimate at the real threshold -- no extrapolation.

| param | MAF | MAC | a | rate / nominal (IS) | +/- 1 SE | rate / nominal (saddlepoint) |
|---|---|---|---|---|---|---|
| `gamma` | 0.0098 | 7,805 | 5e-8 | 1.0059 | 0.0025 | 1.0045 |
| `gamma` | 0.2365 | 189,072 | 5e-8 | 0.9973 | 0.0025 | 0.9999 |
| `beta_left` | 0.0098 | 7,805 | 5e-8 | 1.0167 | 0.0025 | 1.0192 |
| `beta_left` | 0.2365 | 189,072 | 5e-8 | 1.0008 | 0.0025 | 0.9998 |
| `beta_right` | 0.0098 | 7,805 | 5e-8 | 1.0257 | 0.0026 | 1.0216 |
| `beta_right` | 0.2365 | 189,072 | 5e-8 | 1.0049 | 0.0025 | 0.9998 |
| `gamma` | 0.0098 | 7,805 | 1e-5 | 0.9987 | 0.0022 | 1.0020 |
| `gamma` | 0.2365 | 189,072 | 1e-5 | 1.0007 | 0.0022 | 1.0001 |
| `beta_left` | 0.0098 | 7,805 | 1e-5 | 1.0088 | 0.0023 | 1.0078 |
| `beta_left` | 0.2365 | 189,072 | 1e-5 | 1.0003 | 0.0022 | 1.0000 |
| `beta_right` | 0.0098 | 7,805 | 1e-5 | 1.0063 | 0.0023 | 1.0086 |
| `beta_right` | 0.2365 | 189,072 | 1e-5 | 0.9996 | 0.0022 | 1.0000 |

Importance sampling and saddlepoint agree to within 0.5% in all twelve cells
(ratio 0.997 to 1.005), which validates the saddlepoint as the instrument for
the full MAF sweep below.

**Exact saddlepoint over the whole real chr22 MAF spectrum** (Lugannani-Rice on
the exact CGF, both tails separately, weighted by the real per-bin variant
counts). This carries no Monte Carlo error at all; its only error is the
saddlepoint's own, which the IS column above bounds at about 2%.

| param | a = 0.05 | a = 1e-3 | a = 1e-5 | a = 5e-8 | worst MAF bin at 5e-8 |
|---|---|---|---|---|---|
| `beta_left` | 1.00007 | 1.00023 | 1.00055 | **1.0011** | 1.0193 at MAF 0.0098 |
| `beta_right` | 1.00007 | 1.00024 | 1.00061 | **1.0013** | 1.0217 at MAF 0.0098 |
| `gamma` | 1.00006 | 1.00014 | 1.00019 | **1.0002** | 1.0046 at MAF 0.0098 |

**The rejection-rate table asked for**, best estimate with its uncertainty:

| alpha | method | rejection rate / nominal | MC error |
|---|---|---|---|
| 0.05 | plain MC, M = 5.0e6 | 1.0009 - 1.0010 | +/- 0.0020 |
| 1e-3 | plain MC, M = 5.0e6 | 1.000 - 1.009 | +/- 0.014 |
| 1e-5 | plain MC, M = 5.0e6 | 0.98 - 1.04 | +/- 0.14 |
| 1e-5 | tilted IS, M = 5e5/cell | 0.999 - 1.009 | +/- 0.002 |
| **5e-8** | **tilted IS, M = 5e5/cell** | **0.997 - 1.026** | **+/- 0.003** |
| 5e-8 | exact saddlepoint, MAF-weighted | 1.0002 - 1.0013 | none (0.5% method error) |

The test is calibrated. Where it departs from nominal at 5e-8 it is
**anti**-conservative, by up to 2.6% at MAC 7,805 and by 0.02-0.13% averaged
over the real chr22 MAF spectrum. There is no 2.5x, and no 1.2x, and nothing
in the conservative direction larger than 0.3%.

### Attribution

Where could the conservatism have come from, and how much does each source
actually contribute at 5e-8?

| candidate | contribution at 5e-8 | why |
|---|---|---|
| **Reference distribution** | **1.0006x conservative** | The pooled path records `ref_df_source = pooled_sigma0_n_samples`, i.e. `nu = 399,710` -- *not* the per-SNP Satterthwaite law. `t_399710` against a normal is a factor 0.99941. Negligible. |
| **Sigma_0 estimation** | **exactly 1.0000x** (1.002x if the phenotype is also treated as random) | Conditional on the realised phenotype -- the same frame a permutation null uses -- `sd_u^2 = a' Sigma_0 a = sum_i u_i^2 / N` and `Var(u'g) = Var(g) * sum_i u_i^2` identically. `Sigma_0` is not an estimate of anything here; it is the exact second moment. Unconditionally it *is* an estimate, with relative SE `sqrt((kurt_u + 2)/4N)` = 0.21% at `beta_right`'s `kurt_u` = 4.77, worth +0.2% at 5e-8 -- and in the **anti**-conservative direction. |
| **Per-SNP scale `D`** | **1.0000 +/- 0.0004** | The package divides by the realised `d_resid`, not `E[D]`. Plain MC measures `Var(t)` = 1.00015 / 1.00026 / 0.99981. |
| **Weight matrix and OLS projection** | **sets the size of the remaining term, which is 0.1%** | `W` enters only through the shape of `u`. `beta_left` / `beta_right` scores have excess kurtosis +1.38 / +1.77 and skew -0.17 / -0.23; `gamma`'s is -1.89 (platykurtic, near two-point). That is exactly why `gamma` is 5x closer to nominal than the two betas -- 0.02% against 0.13%. |
| **Non-normality of the weighted genotype sum** | **1.0002 - 1.0013x anti-conservative** | The only non-zero term. Driven by the induced cumulants of `S`: at MAF 0.0098 the excess kurtosis of `S` is 5.3e-4 (`beta_left`) and its skew -2.0e-3; both fall by two orders of magnitude by MAF 0.5. |

Two further notes on attribution.

- **The conservatism the project was looking for exists, on the path it no
  longer uses.** On the *jackknife* covariance path the per-SNP Satterthwaite
  law gives `nu` around 800-900, and `t_802` against a normal at 5e-8 is a
  factor **0.74** -- 26% conservative, rising to 23% at `nu = 911` and still 8%
  at `nu = 3,000`. That is a real and large effect, and it is precisely the
  reason `cov-source pooled` is the default. Any measurement of conservatism
  taken on the jackknife path does not transfer to the production one.
- **Imputed dosages do not change this.** Under the standard INFO model a
  dosage is an affine map of the hard genotype (`2f + r(g - 2f)`, `r^2 = INFO`),
  and `t` is invariant to an affine map of `g`, so the null law is identical.
  Real dosages add within-genotype spread, which moves the summand further
  toward Gaussian, i.e. shrinks the only non-zero term above. Consistently, the
  real chr22 `d_resid / 2f(1-f)N` has median 0.995 -- a property of the
  genotype data, not a miscalibration, because the package's `D` and the
  permutation null's variance are the same realised quantity.

### What conservatism would have cost, if it existed

On the real chr21+chr22 covariate-mode scans, forcing the critical value up to
the point where the test rejects at `5e-8 / k`:

| k | critical &#124;z&#124; | `beta_left` hits | lost | recall vs REGENIE | `beta_right` hits | lost | recall |
|---|---|---|---|---|---|---|---|
| 1.0 (measured) | 5.451 | 43 | 0 | 0.913 | 39 | 0 | 0.970 |
| 1.5 | 5.523 | 42 | 1 | 0.891 | 36 | 3 | 0.970 |
| 2.5 (the claim) | 5.612 | 39 | 4 | 0.848 | 35 | 4 | 0.970 |
| 5.0 | 5.731 | 37 | 6 | 0.804 | 33 | 6 | 0.970 |

Even the claimed 2.5x would have cost 4 variants per component and 6.5 points
of `beta_left` recall on two chromosomes. The measured factor is 1.00, so the
cost is zero.

**So what does explain the recall gap?** Nothing, any more -- the gap the brief
quotes (0.41-0.47 female, 0.73-0.78 male) is the *pre-fix* number. After the
LOCO change, recall at 5e-8 on chr22 is 1.000 (`beta_left`) and 0.970
(`beta_right`); pooled over chr21+chr22 it is 0.913 and 0.970. The entire
residual gap is 5 variants: 4 of them in one 5 kb window on chr21
(rs2039246 / rs2824140 / rs2824141 / rs2824143, 18.297-18.302 Mb), where G-MoQ
reaches |z| 5.24-5.30 against REGENIE's 5.62-5.68 at an SE ratio of 1.0466 --
which is the SE ratio and nothing else -- plus rs2294917 on chr22
(|z| 5.43 against 5.83). **100% of that residual gap is the remaining SE ratio
(1.045 on `beta_left`) and the estimand difference, and 0% of it is
conservatism.**

---

## What was not measured

- **One trait, two chromosomes.** Everything in Part 1 is testosterone on
  chr21 and chr22 HapMap3 variants. The SE result replicates across the two
  chromosomes; the discovery result rests on 5 independent loci at 5e-8.
- **`gamma` has no signal here.** Neither chromosome has a REGENIE `sex` hit,
  so the LOCO change's effect on the membership parameter's *power* is
  untested. Its SE and its null calibration are measured.
- **Part 2 is conditional on the realised phenotype.** Genotypes are random,
  `Psi_r` is the real one. That is the same frame as a permutation test and the
  right one for a GWAS null, but it is not a statement about repeated sampling
  of the phenotype.
- **Genotypes in Part 2 are iid HWE, with no LD and no HWE departure.** LD does
  not change any single variant's null law, only the Monte Carlo error of an
  empirical estimate across variants -- which is the whole point of Part 2's
  first section. HWE departure is not modelled; a real-dosage permutation null
  at 5e-8 was not run.
- **The jackknife covariance path was not simulated.** The 0.74 factor quoted
  for `nu = 802` is arithmetic on the t distribution, not a measurement of that
  path end to end.
- **`Q` behaviour at extreme MAF.** The rarest chr22 HapMap3 variant has
  MAC 771; below that the saddlepoint numbers here do not apply.
