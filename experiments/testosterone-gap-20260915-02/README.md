# Corrected testosterone density-interface diagnostic

This is the corrected rerun of `testosterone-gap-20260915-01`. The first
diagnostic centered dosage at the population MAF, while Stage 1's intercept
projection centers the realized dosage; its failed identity output is retained
under the `-01` BP evidence directory. This run uses the same design after
fixing that diagnostic-only error.

Stable ID: `testosterone-gap-20260915-02`

The run uses 60 paired datasets, N=10,000, separation 8, upper-membership
gamma .08, MAF .30, 45 taus (0.02--0.98), seed 20260915, one bounded Python
worker on a BluePebble compute node. It compares historical `-dF/f_model`
weights with proposed same-density `-dF/f_stage1` weights and the equivalent
CDF-scale projection.

Run command:

```bash
bprun --cpus 2 --mem 4G --time 00:30:00 \
  python repo/tests/validity_pilot/testosterone_gap_diagnostic.py \
  --reps 60 --n 10000 --seed 20260915 --out \
  repo/experiments/testosterone-gap-20260915-02/diagnostic.json
```

Source commit: `b45a93d` (runner fix plus BP syntax-check support). The
diagnostic ran as Slurm job `18970538` on `bp1-compute015` through `bprun`,
which invokes `srun --mpi=pmi2`; the output is also retained at the BP path
`/user/work/fh6520/gmoq_finish_20260915/validity/testosterone-gap-20260915-02/diagnostic.json`.
The corrected runner's R syntax check passed as job `18970555` on the same
compute node.

Results: the historical mismatched arm estimated mean gamma `0.035918` (truth
`0.08`), while the same-density arm estimated `0.079556`; the mean model-PDF /
Stage-1-KDE ratio was `1.13293` (mean minimum `0.43787`, mean maximum
`1.22954`). The RIF/CDF identity error was at most `1.79e-13`. These results
support the density correction as the production cause of the attenuation in
this isolated model. They do not establish nuisance-fit uncertainty,
genome-wide calibration, or the actual UK Biobank testosterone result.
