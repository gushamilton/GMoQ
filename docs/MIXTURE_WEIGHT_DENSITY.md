# Which density divides a mixture weight column

`make_weights_normal_mixture()`, `make_weights_normal_mixture_vqtl()` and
`make_weights_generic()` take `density = "model"`, `"hybrid"` or `"stage1"`.
**The default is `"hybrid"` whenever the weight matrix has a membership
(`gamma`) column, and `"model"` otherwise.** A membership W therefore needs
`stage1_density` (see below); without it the builder stops and says so.

The 1000 Genomes example (`docs/examples/1kg_chr22/`) shows why: in a
two-component mixture separated by 4 SDs, `"model"` recovers about 70% of the
true membership effect and `"hybrid"` about 90%, within its standard errors.

## Why the choice exists

A weight column is the quantile derivative

    W_j(tau) = dq_tau / dtheta_j = -(dF/dtheta_j)(q_tau) / f(q_tau).

Stage 1 does not return a density-free quantity. Its RIF-OLS slope is

    b(tau) = -(dF/dtheta)(q_tau) / fhat(q_tau),

where `fhat` is the Stage 1 kernel density at the empirical quantile: a
Gaussian kernel with Silverman bandwidth `1.06 sd(Y) N^-1/5`. `b` carries
`fhat`'s scaling. The weights divide by some density `f`, so the question is
which one.

- **Location columns.** For a mixture these are the component shares
  `P(k | q_tau) = pi_k f_k(q) / f(q)`. They only make sense with the model's
  own `f`. Then they sum to exactly 1, so a shift of every component moves
  every quantile by the same amount. The additive column of a shifted family is
  likewise exactly 1.
- **The membership column.** `-dF/dgamma = -p(1-p)(F_1 - F_2)` is a CDF-level
  quantity with no density in it. Dividing it by the model density gives the
  right scale only if the fitted mixture density equals `fhat` at `q_tau`. In
  a well-separated mixture the fitted density is often badly wrong in the gap
  between the components, even when the fitted CDF is right: the fitted
  density can be far below the kernel density there while the fitted CDF
  matches the empirical one closely. Most of the gamma information comes from
  that gap.

## The three modes

| mode | location / scale columns | membership (`gamma`) column | identities | use |
|---|---|---|---|---|
| `"model"` (default without a membership column) | model pdf | model pdf | hold: `W_1 + W_2 = 1`, additive column = 1 | unimodal and overlapping mixtures; families without a membership column |
| `"hybrid"` (default with a membership column) | model pdf (unchanged) | **Stage 1 density** `fhat` | hold (checked; a violation stops) | two-component mixtures with a membership parameter, especially separated ones |
| `"stage1"` | Stage 1 density | Stage 1 density | **broken**, and the builder warns | sensitivity analysis only |

With no membership column, `"hybrid"` gives the same matrix as `"model"`. For
`make_weights_generic()`, state that explicitly with `stage1_columns =
character(0)`. The default `stage1_columns = "gamma"` refuses rather than
silently doing nothing. For `make_weights_normal_mixture_vqtl()`, only `gamma`
changes. The component sd columns are component share times a scale score,
so they keep the model pdf.

## Supplying the Stage 1 density

`"hybrid"` and `"stage1"` need `stage1_density`, the value of `fhat` at each
tau. You can pass any of:

- a numeric vector, one value per tau;
- a list with `fhat`, and optionally `taus` and `q_tau`. When present, these
  are checked against the builder's grid, and a mismatch stops the build;
- the path of a JSON file with arrays `fhat`, `taus` and `q_tau`, read by
  `read_stage1_density()`.

`stage1_kde_density(y, q_tau)` recomputes `fhat` from the phenotype. It uses
the same kernel, the same bandwidth rule (population sd) and the same floor as
Stage 1 (`gmoq.core.compute_rif()` in Python). Use the exact phenotype
vector Stage 1 analysed (after any LOCO subtraction, if you used
`--loco-mode subtract`) and Stage 1's `q_tau` (median-unbiased, R `type = 8`).
Build the weights on that same `q_tau`.

## Provenance

A non-default matrix carries these attributes: `density_mode`,
`density_weight_mode`, `stage1_density_columns`, `stage1_density`,
`model_density` and, for a JSON input, `stage1_density_source`. The provenance
labels are `parametric_model_density`, `hybrid_stage1_membership_density` and
`stage1_rif_density`.

A `"model"` matrix carries no attributes. `weight_density_mode(W)` therefore reports `"model"` for any
matrix without them. It also reads a weight-object list's
`provenance$density_weight_mode`.

`param_gwas_from_file()` and `param_gwas_multi_from_file()` take
`require_density_mode`. When it is set, the run stops unless the weights were
built in that mode.

## Evidence

In simulation (10 mixture scenarios, 400 replicates each, N from 100,000 to
400,000, pooled covariance), for a variant that changes membership only:

- **Separated mixtures where the fitted density differs from the Stage 1
  kernel density in the gap** (including a correctly specified mixture):
  `"model"` recovered 35–74% of the true gamma, `"hybrid"` 98–100%. Coverage
  of 95% intervals rose from at most 13% to 94–97%. Under `"model"`, a third to
  two thirds of gamma also leaked into the lower component's location;
  under `"hybrid"` the leak was below 1%.
- **Heavily overlapping mixtures:** 90–92% under `"model"`, 97–98% under
  `"hybrid"`.
- **Null calibration** was unchanged.

Caveats:
- Gamma's absolute scale is only as good as `fhat` in the gap, which is
  typically one or two quantiles. Its sign and z statistic are robust.
- `"hybrid"` leaves a small location bias (4–9% on the narrower component in
  separated mixtures), because `P(k | q)` ignores the kernel density's own
  smoothing bias. `"stage1"` removes that bias but breaks both identities.

The test suite (`stage2/tests/testthat/test-weights.R`) reproduces the
attenuation and its removal on a small Laplace-component mixture.
