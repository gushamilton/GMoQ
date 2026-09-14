# Independent verification of the tail-calibration fix (v0.3.0)

Run by the coordinating session on 2026-09-14, before merging
`fix/tail-calibration-t-reference`, and independently of the agent that wrote
the code. Recorded here because the acceptance test
(`fungwas_paper/docs/BP1_acceptance_null_power.md` §11) made it a merge
condition, and a verification that exists only in a transcript is not a record.

Method: the branch was cloned and built fresh (`954bc6d`), the Stage 2 package
installed to a scratch library with the loaded path asserted (`~/.Rprofile`
silently overrides `R_LIBS_USER`), and `param_gwas()` called on synthetic Stage
1 input with a positive-definite covariance built as
`diag(se) R_BB diag(se)`. 20 taus, 300 SNPs, `tau_cov_estimator = "ols"`.

## Results

| check | result |
|---|---|
| `max abs(emitted p - 2*pt(-abs(z), nu))`, location parameter | **0.000e+00** |
| `nu` resolved per SNP and per parameter at B=1000, n_eff=4e5 | 996.8 (location), 992.5 (scale) |
| contrast classifier assigns `c` | 0.88 flat / 2.62 scale — correct |
| `ref_dist` / `ref_df_source` recorded in output | `t` / `per_snp_satterthwaite` |
| explicit `nu = 999` honoured | yes, `user_supplied` |
| convergence to normal at `nu = 1e8` | `max abs(p_t - p_norm)` = 1.71e-09 |
| covariance with no metadata | **refuses**, emits no p-value, warns explicitly |
| regression against v0.2.0 on identical inputs | none — see below |

## A non-finding worth recording

On the synthetic covariance used here, the scale parameter's standard error is
`NA` for 143 of 300 SNPs, so its p-value is `NA` too. This is **not** a
regression: the same inputs through the v0.2.0 package produce an identical
count (157/300 finite) and identical values. It is `compute_calibrated_se`'s
variance guard (`isfinite(var) && var > 0`) correctly refusing to emit a
standard error where the synthetic sandwich is degenerate. The guard predates
this work and is behaving as intended.

Recorded because "the new version returns NA" is exactly the kind of
observation that becomes a false regression report if the comparison against
the old version is not made and written down.

## Scope

Synthetic input only; the end-to-end Stage 1 -> Stage 2 check on real BGEN data
was performed separately by the implementing agent and is reported in its own
account. Both were run; only this one was independent of the author.
