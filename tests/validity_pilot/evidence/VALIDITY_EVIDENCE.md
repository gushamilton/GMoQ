# Compact validity evidence

The authoritative main run is `pilot_linear` on `bp1-compute001`:

```text
/user/work/fh6520/gmoq_finish_20260915/validity/pilot_linear
```

Its metadata records source `48171e3`, N=10,000, 100 independent datasets,
45 taus, 1,000 blocks, and two compute threads. The independent dataset is
the Monte Carlo unit. The aggregate CSVs in this directory are copied from
that run; no large result payloads are committed.

The exact compute commands were:

```bash
/Users/fergushamilton/bin/bprun --cpus 2 --mem 32G --time 00:30:00 \
  env OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 \
  /user/home/fh6520/miniforge3/envs/fungwas-stage1/bin/python \
  /user/work/fh6520/gmoq_finish_20260915/validity/repo/tests/validity_pilot/pilot_stage1.py \
  --repo /user/work/fh6520/gmoq_finish_20260915/validity/repo \
  --source-commit 48171e3 --source-status work-pooled-validity-48171e3 \
  --out /user/work/fh6520/gmoq_finish_20260915/validity/pilot_linear \
  --reps 100 --N 10000 --blocks 1000 --threads 2 --nbg 100000 --Nref 200000
```

The Stage 2 command used the installed `fungwasStage2` library and ran
`pilot_stage2.R`, `summarise_pilot.R`, and `report_pilot.R` on the same output.
The compute-node wrapper was:

```bash
module load languages/R/4.5.1
export GMOQ_TEST_LIB=/user/work/fh6520/gmoq_finish_20260915/validity/lib
export GMOQ_DEP_LIB=/user/work/fh6520/r_packages
Rscript --vanilla /user/work/fh6520/gmoq_finish_20260915/validity/repo/tests/validity_pilot/pilot_stage2.R \
  /user/work/fh6520/gmoq_finish_20260915/validity/pilot_linear \
  /user/work/fh6520/gmoq_finish_20260915/validity/pilot_linear/pilot_rows.csv
Rscript --vanilla /user/work/fh6520/gmoq_finish_20260915/validity/repo/tests/validity_pilot/summarise_pilot.R \
  /user/work/fh6520/gmoq_finish_20260915/validity/pilot_linear
Rscript --vanilla /user/work/fh6520/gmoq_finish_20260915/validity/repo/tests/validity_pilot/report_pilot.R \
  /user/work/fh6520/gmoq_finish_20260915/validity/pilot_linear
```

The richer-basis Stage 1 wrapper used the same BP environment and called
`covariate_basis_probe.py --source-commit be70f0a --reps 100 --N 10000
--blocks 1000 --threads 2 --Nref 200000`; its R wrapper then called
`covariate_basis_probe.R` and wrote `probe_rows_summary.csv`.
The per-cell provenance table records the source commit in every Stage 1
metadata file. The `pilot_linear` directory was removed before this Stage 1
command. Existing common/rare payloads happen to hash identically to the
earlier deterministic `pilot_final` payloads because the DGP and seeds for
those cells did not change; their metadata still records the recomputed
`48171e3` run. The four Normal variance-linear cells have no prior payload.

The richer-basis diagnostic is in
`/user/work/fh6520/gmoq_finish_20260915/validity/cov_probe` and compares the
same shipped production pipeline with C versus C,C²,C³, using the same seeds,
N, reps, taus, blocks, covariance sources, and W arms. It was run from source
`be70f0a`; its summary includes reported-SE/empirical-SD ratios as well as
point bias, rejection, and Q.

Interpretation limits: the strong Normal variance cell has beta_sigma2 bias
`-0.00912 +/- 0.00218` under JK/oracle (about 3.8% of v=0.24 but 4.2 MCSE),
so it is statistically distinguishable from zero at this pilot precision. Its
scientific coverage is 0.96 (JK) and 0.94 (pooled). The richer basis removes
the large linear-basis point bias, but pooled rare-SNP rejection remains about
0.30 for beta_sigma2 versus 0.09 for JK, so covariance calibration remains an
open supported-domain limitation.
