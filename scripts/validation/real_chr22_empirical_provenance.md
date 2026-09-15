# Bounded empirical log_bmi provenance

The matched exploratory fit was run on `bp1-compute053` with the 24 SNP list
from the known-truth smoke. Participant-level inputs and Stage 1 outputs remain
on BluePebble under
`/user/work/fh6520/gmoq_finish_20260915/inputs/real_chr22_smoke/empirical/`.

The baseline builder calls the same Stage 1 input functions as the production
CLI: `load_analysis_inputs` with `keep_log_bmi` and
`age,sex,PC1..PC10`; `load_loco_predictions` and
`resolve_loco_scale("regression")`; then `get_variants_bgen` and
`match_samples`. It therefore uses natural-log BMI, the chr22 LOCO prediction,
and the exact matched BGEN sample mask before calculating the baseline mean,
SD, and type-8 (`median_unbiased`) quantiles.

The resulting aggregate baseline is:

| quantity | value |
| --- | ---: |
| phenotype rows after LOCO matching | 484,460 |
| BGEN matched rows used for baseline | 461,320 |
| LOCO regression scale | 0.1606963867041275 |
| matched baseline mean | 3.296314455837662 |
| matched baseline SD | 0.1556538337922891 |
| tau count | 45 |

The Stage 1 covariance sidecar independently reports 461,320 samples, the same
12 covariates, the same LOCO path and scale, chromosome 22, and 45 taus. The
serialized Stage 1 tau labels are rounded to three decimals; both the baseline
builder and Stage 2 use those labels, while the RIF implementation internally
generates the corresponding 45-point grid before serializing it.

The final pooled OLS mean-plus-variance output is
`log_bmi_matched_both_pooled.tsv`. For `rs9617528`, it reports
`beta_mu = -0.0005071266731` (SE `0.0006064544180`, p `0.4030340`) and
`beta_sigma2 = 0.0001195505962` (SE `0.0001490502422`, p `0.4225056`;
Q `24.4983972`, Q p `0.1817302`). These are descriptive results for a
bounded empirical SNP set and are not treated as known-truth validation.
