# Fourth-order density repair prototype

This is a test-only, paired diagnostic. It does not modify production
`compute_rif` and does not support a release or repair claim. Each frozen
dataset was sent through the shipped Stage 1 residualisation, C++ block-score,
and covariance kernels twice:

- `A_production_kde`: the production type-8 quantile and Gaussian Silverman KDE;
- `D_fourth_order_kde`: the same empirical quantile and bandwidth, using
  `K(u) = (1.5 - 0.5*u^2) * phi(u)` and a `1e-8` denominator floor.

The exact Normal-mixture projection and scientific targets are held fixed
between arms. The installed Stage 2 wrapper then fits both pooled+OLS and
jackknife+OLS, with oracle and estimated Normal mean/variance W. Coverage uses
the shipped `param_ref_df` reference degrees of freedom; the summary counts
every attempted dataset in the coverage denominator.

## Provenance and commands

Run ID: `density-repair-20260915-01`

Remote output:

```text
/user/work/fh6520/gmoq_finish_20260915/validity/density-repair-20260915-01
```

Stage 1 ran on `bp1-compute019` under Slurm allocation `18966618`, with
N=10,000, 100 datasets per cell/arm, three cells (`v=0`, `.06`, `.24`), 45
taus, 1,000 blocks, two threads, seed `20260915`, and the requested 8 GB
wrapper budget within the held 32 GB allocation. The executed harness source was
commit `28c5316d55679ab07da943e2700e246869953760`; the metadata records script
SHA-256
`ec2153785c2993945d4fe885a279ef68fe71a72b3283c1ee8495fb5c724b3f7f`.

The exact Stage 1 command was:

```text
/Users/fergushamilton/bin/bprun --cpus 2 --mem 8G --time 01:00:00 \
  env OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 \
  /user/home/fh6520/miniforge3/envs/fungwas-stage1/bin/python \
  /user/work/fh6520/gmoq_finish_20260915/validity/repo/tests/validity_pilot/density_repair_probe.py \
  --repo /user/work/fh6520/gmoq_finish_20260915/validity/repo \
  --out /user/work/fh6520/gmoq_finish_20260915/validity/density-repair-20260915-01 \
  --run-id density-repair-20260915-01 --reps 100 --N 10000 --blocks 1000 \
  --threads 2 --seed 20260915 \
  --source-commit 28c5316d55679ab07da943e2700e246869953760 \
  --source-status density-repair-20260915-01-local-clean-28c5316
```

The installed Stage 2 command used R 4.5.1 from
`/software/local/languages/miniforge3/envs/r-4.5.1/bin/Rscript`, the existing
`pilot_stage2.R` wrapper, and the configured `fungwasStage2` test/dependency
libraries (`GMOQ_TEST_LIB=/user/work/fh6520/gmoq_finish_20260915/validity/lib`,
`GMOQ_DEP_LIB=/user/work/fh6520/r_packages`). It wrote 2,400 Stage 2
dataset-arm fits and 4,800 parameter rows from the 600 Stage 1 dataset/arm
combinations. The failure-aware postprocessor is committed in `bc09b45` and
`34439b9`.

The exact Stage 2 command was:

```text
/Users/fergushamilton/bin/bprun --cpus 2 --mem 8G --time 00:30:00 \
  env GMOQ_TEST_LIB=/user/work/fh6520/gmoq_finish_20260915/validity/lib \
  GMOQ_DEP_LIB=/user/work/fh6520/r_packages \
  /software/local/languages/miniforge3/envs/r-4.5.1/bin/Rscript --vanilla \
  /user/work/fh6520/gmoq_finish_20260915/validity/repo/tests/validity_pilot/pilot_stage2.R \
  /user/work/fh6520/gmoq_finish_20260915/validity/density-repair-20260915-01 \
  /user/work/fh6520/gmoq_finish_20260915/validity/density-repair-20260915-01/pilot_rows.csv
```

## Results

All 600 attempted Stage 1 dataset/arm combinations produced artifacts, and
all 2,400 Stage 2 fits succeeded: 100/100 per cell/arm/source/W, with zero
nonfinite estimates, standard errors, or
reference z values. The fourth-order raw density had no nonpositive or floored
tau values in any of the 300 repaired datasets. Its minimum density averaged
0.0458 for `v=.24`; the production minimum averaged 0.0480. The complete
per-dataset and aggregate failure tables are stored beside this report.

For the strong variance cell (`v=.24`, oracle W), the projection-target bias
and coverage were:

| covariance | arm | projection bias ± MCSE | scientific bias ± MCSE | projection coverage ± MCSE | scientific coverage ± MCSE | mean SE / empirical SD |
|---|---|---:|---:|---:|---:|---:|
| jackknife | production | −0.003739 ± 0.002654 | −0.006502 ± 0.002654 | 0.92 ± 0.027 | 0.91 ± 0.029 | 0.894 |
| jackknife | fourth-order | +0.000022 ± 0.002705 | −0.002741 ± 0.002705 | 0.92 ± 0.027 | 0.93 ± 0.026 | 0.891 |
| pooled | production | −0.003739 ± 0.002654 | −0.006502 ± 0.002654 | 0.92 ± 0.027 | 0.90 ± 0.030 | 0.848 |
| pooled | fourth-order | +0.000022 ± 0.002705 | −0.002741 ± 0.002705 | 0.91 ± 0.029 | 0.92 ± 0.027 | 0.845 |

Estimated W gave the same conclusion. The paired fourth-order minus
production shift in `beta_sigma2` was `+0.003761 ± 0.000101` for both
covariance sources; the location shift was negligible (`−0.000010 ±
0.000025`). At `v=.06`, the shift was `+0.000924 ± 0.000048`; at `v=0`, it
was effectively zero (`−0.000030 ± 0.000044`).

The weaker cells had 0.97–0.98 pooled/jackknife coverage for the variance
parameter and approximately 0.96 coverage at v=0, with binomial MCSEs from
0.014 to 0.020. The v=0 arm's finite-pilot point estimate was about
−0.0040 ± 0.0022 in both kernels; this is not evidence of a nonzero null bias
at this precision. The fourth-order change therefore does not establish broad
validity.

The prototype removes the observed strong-cell point bias against the exact
projection target, but its 100-dataset coverage remains too imprecise and the
pooled SE/empirical-SD ratio remains below one. The current jackknife uses a
fixed RIF and omits density-estimation nuisance uncertainty. These results
therefore identify a promising bias diagnostic while leaving covariance
calibration unresolved; they do not justify a production repair.

The per-arm elapsed times were about 1.9–2.0 seconds for 100 datasets. Those
times include the extra raw-density and exact-reference diagnostics in the
harness and are not a fair production-versus-repair speed benchmark.
