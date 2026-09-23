# fungwasStage2 0.5.0

* The mixture and generic weight builders (`make_weights_normal_mixture()`,
  `make_weights_normal_mixture_vqtl()`, `make_weights_generic()`) take
  `density = c("model", "hybrid", "stage1")`. With `"hybrid"` the
  component-share columns stay on the model density and still sum to 1. The
  membership (`gamma`) column is divided by the Stage 1 density, because the
  Stage 1 coefficient carries that density's scaling. The default `"model"` is
  unchanged: its output is identical to 0.4.0. `"stage1"` reproduces the old
  all-KDE weights and warns. See `docs/MIXTURE_WEIGHT_DENSITY.md`.
* New: `stage1_kde_density()`, which is the Stage 1 Silverman KDE at `q_tau`;
  `read_stage1_density()`, which reads an `fhat`/`taus`/`q_tau` JSON sidecar;
  and `weight_density_mode()`, which reads back the mode a weight matrix was
  built with.
* `param_gwas_from_file()` and `param_gwas_multi_from_file()` gain
  `require_density_mode`, which refuses weights built in another mode.
* `make_weights_generic()` gains `unit_sum` for checking identities such as
  `W_left + W_right = 1` or an additive column of 1.
* The package now has testthat tests (`tests/testthat/`). Tier 0 of the suite
  runs them.
