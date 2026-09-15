# Bounded biobank covariance audit

This audit replays the preserved Stage 1 binary payload through the current
Stage 2 package. It is deliberately limited to the shipped OLS projection,
with pooled covariance as the production path and one paired jackknife run as
a diagnostic. It does not regenerate Stage 1, simulate new phenotypes, or
alter the package.

The payload is at
`/user/work/fh6520/covsource_run/raw` on BluePebble. Each selected cell has
20 complete filesets, `cellNNN_chunk0000` through `cellNNN_chunk0019`, with
750 replicates per fileset (15,000 replicates per cell), 45 taus and 1,000
blocks. The selected cells are:

| index | cell | MAF | purpose |
|---:|---|---:|---|
| 0 | `A_null_maf0.3_normal` | 0.30 | common-frequency null |
| 1 | `A_loc_maf0.3_normal_Z16` | 0.30 | common-frequency location signal |
| 2 | `A_disp_maf0.3_normal_Z16` | 0.30 | common-frequency dispersion signal |
| 4 | `A_null_maf0.05_normal` | 0.05 | low-frequency null |
| 5 | `A_loc_maf0.05_normal_Z16` | 0.05 | low-frequency location signal |
| 6 | `A_disp_maf0.05_normal_Z16` | 0.05 | low-frequency dispersion signal |

The driver uses the checked-in `snapshot/cells.json` and `snapshot/truth_cs.R`
by default; an explicit snapshot directory can be supplied as the eighth
argument. The source DGM is `Y = a*g + (1 + c*g)*E`, with standard Normal `E` and
counted allele dosage `g`. The dispersion coefficient `c` is a standard
deviation multiplier. The audit compares each OLS result with the exact
finite-mixture RIF projection target in `truth_cs.R`; `Ztarget=16` refers to
that projection target. Thus the reported `beta_sigma2` coverage is coverage
of the OLS projection estimand, on squared-outcome units per counted allele,
and is not silently treated as a variance coefficient from the generating
formula.

The preserved Stage 1 implementation files match source commit `24eef02`
(`stage1` in the validity repository); only scratch packaging metadata differs.
The current Stage 2 package was installed from this branch at commit
`b1be4d7` into `/user/work/fh6520/gmoq_finish_20260915/inference/lib`.
The separate ACTUAL100 validity pilot has provenance `f1b4a2c`; it is not the
source of these N=461,320 replays.

Run on a BP compute node with R 4.5.1 and the isolated library:

```bash
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export PATH=/software/local/languages/miniforge3/envs/r-4.5.1/bin:$PATH
export R_LIBS_USER=/user/work/fh6520/gmoq_finish_20260915/inference/lib:/user/work/fh6520/r_packages

R --vanilla --slave --file=stage2_biobank_audit.R --args \
  0 /user/work/fh6520/covsource_run/raw OUT 2500 \
  /user/work/fh6520/covsource_run \
  /user/work/fh6520/gmoq_finish_20260915/inference/lib
R --vanilla --slave --file=summarise_biobank_audit.R --args OUT \
  /user/work/fh6520/covsource_run \
  /user/work/fh6520/gmoq_finish_20260915/inference/lib
```

Repeat the first command for cell indices `1`, `2`, `4`, `5`, and `6` in the
same `OUT` directory. The default driver is pooled-only. For the paired
jackknife diagnostic, use the optional seventh argument
`pooled_jackknife` for cell `6` and place that RDS file in a separate output
directory before running the summariser.

The checked-in aggregate results contain no raw Stage 1 data. In the pooled
OLS replay, null type-I error and projection coverage are near nominal for
both MAFs and for common-frequency location/dispersion. At MAF 0.05 under the
dispersion signal, pooled `beta_sigma2` coverage is `0.9285` (MCSE `0.0021`)
with SE/empirical-SD `0.9206`; the paired jackknife is `0.9447` (MCSE
`0.0019`) with ratio `0.9860`. This is a bounded signal-level limitation of
the shared `Sigma_0 / d_resid` working-independence covariance. It does not
justify an unvalidated covariance rewrite: retain pooled as primary, retain
explicit jackknife confirmation, and validate a future genotype-weighted
influence covariance against exact matched projection targets.
