# Paired density/RIF diagnosis

These runs isolate the Stage 1 RIF input while retaining the shipped
residualisation, C++ block-score, and covariance kernels. The three paired
arms use the same genotype, phenotype, block IDs, and random seed for each
dataset:

- `A_empirical_kde`: empirical type-8 quantile and the production Silverman
  KDE density;
- `B_empirical_exact_f`: the same empirical quantile and exact three-component
  Normal-mixture density evaluated at that quantile;
- `C_population_exact`: exact population quantile and exact density.

The artifacts were then passed through the installed Stage 2 wrapper with
oracle Normal mean/variance W, both pooled+OLS and jackknife+OLS. Pooled and
jackknife point estimates are identical, as expected; their covariance source
only changes standard errors and reference calibration.

## Immutable runs

Both runs used `bp1-compute053`, 45 taus, 1,000 blocks, and two OpenMP
threads. Each diagnostic dataset contains one common-MAF SNP so the RIF-arm
contrast is not diluted by additional variants. The independent dataset is
the Monte Carlo unit.

`density-bias-20260915-01` used source harness commit
`47ec83d777934c4e675fedaaf8ac29ba5817e534`, N=10,000, 100 datasets, and
script SHA-256
`713588c092abb5f8ac0b629acabe4795368e04249351de08848a9302bfcb91a8`.

`density-bias-20260915-02` used harness commit
`eadf3fa7c8146ff9d7ed99d30377728aac6791e4`, N=100,000, 50 datasets, v=.24
only, and script SHA-256
`13e7f1bb86a1f66cb2b2ddb9c5def205206771b1d71775d13486acaf702d2732`.

The exact Stage 1 command for the first run was:

```text
/Users/fergushamilton/bin/bprun --cpus 2 --mem 8G --time 01:00:00 \
  env OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 \
  /user/home/fh6520/miniforge3/envs/fungwas-stage1/bin/python \
  /user/work/fh6520/gmoq_finish_20260915/validity/repo/tests/validity_pilot/density_probe.py \
  --repo /user/work/fh6520/gmoq_finish_20260915/validity/repo \
  --out /user/work/fh6520/gmoq_finish_20260915/validity/density-bias-20260915-01 \
  --run-id density-bias-20260915-01 --reps 100 --N 10000 --blocks 1000 --threads 2 \
  --source-commit 47ec83d777934c4e675fedaaf8ac29ba5817e534
```

The second run used the same command shape with output ID
`density-bias-20260915-02`, `--cell v24`, `--reps 50`, `--N 100000`, and
source commit `eadf3fa7c8146ff9d7ed99d30377728aac6791e4`.

## Results

For the strong variance cell (`v=.24`), the paired difference in the fitted
`beta_sigma2` is:

| N | A − B | B − C |
|---:|---:|---:|
| 10,000 | −0.003989 ± 0.000324 | +0.000368 ± 0.000288 |
| 100,000 | −0.001834 ± 0.000146 | +0.0000037 ± 0.000100 |

The values are identical for pooled and jackknife point estimates. The
corresponding bias against the exact projection target was:

| N | A | B | C |
|---:|---:|---:|---:|
| 10,000 | −0.003707 ± 0.002549 | +0.000282 ± 0.002617 | −0.000086 ± 0.002569 |
| 100,000 | −0.001909 ± 0.001016 | −0.000074 ± 0.001066 | −0.000078 ± 0.001037 |

At N=10,000 the mean `fhat / f_exact(q_emp)` was 0.999669 ± 0.000143;
at N=100,000 it was 0.999636 ± 0.000067. The average ratio hides tau-specific
smoothing, but the paired A−B effect is direct evidence of its impact on the
variance projection. The empirical-to-population quantile RMS fell from
0.01559 ± 0.00048 to 0.00447 ± 0.00017.

The weaker variance cells show the same direction: A−B is
−0.000932 ± 0.000105 at v=.06 and approximately zero at v=0. Mean-effect
contrasts are near zero. B−C is small and vanishes at larger N, while C is
essentially unbiased against the exact projection target.

The A−B contrast falls by roughly the N^-0.4 rate expected from the squared
Silverman bandwidth, whereas B−C follows the finite-quantile error. The
strong-cell discrepancy is therefore primarily production KDE smoothing in
the RIF density denominator. Finite-effect projection still explains the
remaining difference between the scientific target v and the exact
projection target, but it does not explain the production A-arm bias.

The absolute A-arm bias is not expected to equal the earlier six-SNP pilot
bit-for-bit because this diagnosis deliberately uses one SNP and new frozen
draws. The paired A−B and B−C contrasts are the within-run causal
comparisons.

The compact CSVs next to this report contain the complete arm, paired, and
density-ratio summaries for both immutable runs.
