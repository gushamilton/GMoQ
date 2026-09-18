# The OLS full-covariance model score

Status: implementation and unit tests, 2026-09-18. Additive: the default
`param_gwas_multi()` path is unchanged.

## The gap this fills

Stage 2 fits `b ~ W theta`, where `b` is the vector of `T` per-quantile effects
for one SNP, `Sigma` is their `T x T` covariance, and `W` (`T x k`) encodes a
candidate parametric shape. Two loss/estimator pairings existed:

| path | fit | residual loss | penalty |
|---|---|---|---|
| `param_gwas(se_mode = "tau_cov", tau_cov_estimator = "ols")` | OLS | diagonal, `r' D^-1 r` | none returned |
| `fit_multi_models` | DWLS | diagonal | `Q + 2k` |
| `fit_multi_models_tau_cov` | GLS | full, `r' Sigma^-1 r` | `Q + 2k` |

The missing combination was an **OLS point estimate** scored with the
**per-SNP full covariance**, which is the natural comparison criterion when
OLS is the primary estimator (`docs/OLS_VS_GLS.md`). That is what
`residual_loss = "full_cov"` provides.

## The mathematics

For an OLS fit,

```
theta_ols = (W'W)^-1 W'b,     H = W (W'W)^-1 W',     r = b - W theta_ols.
```

The full-covariance residual statistic is

```
Q_full = r' Sigma^-1 r.
```

Fitting `k` parameters is optimistic by the expected spurious reduction in the
loss. For a `Sigma^-1`-weighted quadratic loss the expectation is exact:

```
E[spurious reduction] = 2 tr(Sigma^-1 H Sigma) = 2 tr(H) = 2 rank(W).
```

So the coherent OLS model score is

```
AIC = Q_full + 2 rank(W)
BIC = Q_full + rank(W) log(N).
```

`rank(W)` is used, not `ncol(W)`: a duplicated or collinear column must not buy
a full penalty. The rank is reported per model so a rank-deficient design is
visible rather than silently scored.

### The diagonal contrast

With a diagonal loss `r' D^-1 r` the correct optimism is

```
2 tr(D^-1 H Sigma),
```

which is **not** `2k`, because the per-tau errors are strongly correlated
(adjacent taus correlate above 0.99). The existing diagonal-`Q` + `2k` pairing
is therefore not a valid comparison criterion; the full-covariance score above
is. When `Sigma = s^2 I` the two coincide, `Q_full = Q_diagonal` and the
optimism collapses to `2k`.

## API

`residual_loss` is a new argument on `param_gwas_multi()` and
`param_gwas_multi_from_file()`, defaulting to `"auto"` (the previous
behaviour, byte-for-byte). Set it to `"full_cov"` to use the OLS fit and the
full-covariance score. `return_params = TRUE` also returns the fitted OLS
coefficient matrices. `se_mode` must be `"tau_cov"` and `cov_data` must be
supplied (the file wrapper supplies it from the pooled `Sigma_0` sidecar or the
jackknife payload, exactly as for the other tau-covariance paths).

With `residual_loss = "full_cov"` the call returns a list:

| element | shape | meaning |
|---|---|---|
| `Q_full` | `P x n_models` | `r' Sigma^-1 r` per SNP and model |
| `rank` | `n_models` | `rank(W)`, the penalty actually applied |
| `AIC` | `P x n_models` | `Q_full + 2 rank(W)` |
| `BIC` | `P x n_models` | `Q_full + rank(W) log(N)` |
| `reason` | `P x n_models` | `"ok"`, or why the entry is `NA` |
| `rank_deficient` | `n_models` | `TRUE` when `rank(W) < ncol(W)` |
| `params` | list (`return_params = TRUE`) | OLS point estimates, `k x P` |

### Example: score several candidate `W` for the same SNP, and pick the winner

```r
library(fungwasStage2)

## One trait's Stage 1 file; pooled Sigma_0 is the default covariance source.
res <- param_gwas_multi_from_file(
  stage1_file = "testosterone_chr22.stage1.tsv.gz",
  W_list = list(
    additive       = W_add,    # T x 1
    location_scale = W_ls,     # T x 2
    skew           = W_skew    # T x 3
  ),
  N = 461320,
  se_mode = "tau_cov",
  cov_source = "pooled",
  residual_loss = "full_cov"
)

## Lower AIC wins. res$AIC is P x n_models.
best <- max.col(-res$AIC, ties.method = "first")
data.frame(snp = rownames(res$AIC), model = colnames(res$AIC)[best])

## Inspect the applied penalties and any refusals.
res$rank
table(res$reason)
```

The same call against an in-memory Stage 1 object (for example one SNP read
from a scan) uses `param_gwas_multi(stage1, W_list, N, se_mode = "tau_cov",
cov_data = cov_data, residual_loss = "full_cov")`.

## Edge cases and refusal

`reason` records why an entry is `NA` rather than returning a wrong number:

| `reason` | meaning |
|---|---|
| `ok` | score computed |
| `beta_nonfinite` | Stage 1 effect is not finite |
| `offset_bad` | covariance byte offset falls outside the payload |
| `cov_scale_invalid` | per-SNP `cov_scales` value is non-finite or `<= 0` |
| `cov_nonfinite` | covariance entry is not finite |
| `sigma_not_pd` | `Sigma` is not positive definite (no ridge is applied) |
| `sigma_ill_conditioned` | reciprocal condition number below `1e-12` |
| `sigma_solve_fail` | `Sigma^-1 b` or `Sigma^-1 W` failed |
| `theta_solve_fail` | the OLS solve failed |
| `q_invalid` | `Q_full` is not finite, or materially negative |

A rank-deficient `W` is *not* refused: the residual (hence `Q_full`) is
well-defined because it is a projection residual, so the score is returned,
but `rank_deficient` is `TRUE` and the penalty uses `rank(W)`. Where `W` is
singular the returned `params` are a minimum-norm least-squares solution; the
score does not depend on that choice.

Covariance is read through the existing plumbing, so both sources work
unchanged: the pooled route (`Sigma_0` from the `.sigma0.tsv.gz` sidecar with
its per-SNP `1/d_resid` scale) and the jackknife route (`cov_vec` / `offsets` /
`cov_scales`). The `plugin_cor` shared-correlation approximation is a
different covariance definition and is not used.
