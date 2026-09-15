# Testosterone density-interface diagnostic

This bounded diagnostic asks whether the historical testosterone Stage 2
construction uses a density inconsistent with the production Stage 1 RIF.
Stage 1 uses a Silverman Gaussian KDE at empirical type-8 quantiles. The
historical two-component mixture weights use the fitted-mixture PDF at those
quantiles. The script uses the production Stage 1 KDE calculation and an oracle
two-normal mixture (`sep=8`, upper membership `gamma=.08`, MAF `.30`) so this
isolates the density interface without UK Biobank inputs or mixture fitting.

Stable ID: `testosterone-gap-20260915-01`

The diagnostic source is committed before execution. The run uses 60 paired
datasets, N=10,000, 45 taus (0.02--0.98), seed 20260915, one bounded Python
worker on a BluePebble compute node. It compares the historical mismatched
weights `-dF/f_model` with the proposed same-density weights `-dF/f_stage1` and
the equivalent CDF-scale projection.

Run command:

```bash
bprun --cpus 2 --mem 4G --time 00:30:00 \
  python repo/tests/validity_pilot/testosterone_gap_diagnostic.py \
  --reps 60 --n 10000 --seed 20260915 --out \
  repo/experiments/testosterone-gap-20260915-01/diagnostic.json
```

The output JSON is aggregate evidence only; no participant-level inputs are
used. The expected decision rule is narrow: if the mismatched arm shows the
archive's attenuation of the membership effect and the matched arm agrees with
the CDF arm to floating-point precision, implement the same-density correction
in the main-repo testosterone weight builder. This does not validate nuisance
fit uncertainty, genome-wide calibration, or the full UK Biobank result.
