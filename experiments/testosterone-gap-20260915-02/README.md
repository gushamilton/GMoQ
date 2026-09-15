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
