# gmoq 0.8.1

* Smaller Stage 1 output. `beta_tau_*` is written to 0.001 of its standard error and `se_tau_*`
  to four significant digits (about a third smaller). With `--cov-mode pooled-only`, where the
  standard errors follow exactly from Sigma_0 and `d_resid`, Stage 1 writes `z_tau_*` (to three
  decimals) instead of `beta_tau_*` and `se_tau_*`, about three times smaller; Stage 2 rebuilds
  the effects and standard errors from the `.sigma0` sidecar. Estimates change by at most
  0.0005 standard errors.

# gmoq 0.8.0

* Stage 1 is several times faster end to end. With a `.bgi` index and `--snps`, BGEN records
  are read and decoded in parallel; the per-quantile cross-products no longer grow with the
  number of quantiles; Stage 2 calibrates Q once per pooled covariance rather than
  once per SNP.
* Dosages are now the exact expected allele count of the stored probabilities (previously
  rounded to single precision), so Stage 1 estimates differ from 0.7.0 in about the seventh
  significant digit.
* New `--jobs FILE`: several phenotypes, each with its own samples, covariates and LOCO file,
  analysed in one pass over the genotypes; each output is identical to a separate run.
* Building Stage 1 now needs zlib headers; libdeflate is used when available.

# gmoq 0.7.0

* Renamed from fungwas_stage1 / fungwasStage2 to gmoq. Releases up to 0.6.2 used the old names; results are unchanged.
  The Stage 1 Python package and distribution are now `gmoq` and its command
  is `gmoq-stage1` (`python -m gmoq.cli` also works); the R package is `gmoq`;
  the conda environment is `gmoq`; the Stage 1 build variables are
  `GMOQ_STAGE1_ARCH` and `GMOQ_STAGE1_OPENMP`.
* The `.cov.meta.json` sidecar records a plain-language provenance string for
  its degrees-of-freedom model. No numerical output changes.
* Documentation rewritten for users: `docs/METHOD.md` (how the method works and
  its known limitations) and `docs/INPUT_PREPARATION.md` replace the earlier
  development notes.
* One compact test suite per package, covering the behaviour users rely on.
  The on-demand CI workflow `example.yml` runs the 1000 Genomes example end to
  end.

# gmoq 0.6.2

Stage 1 only; no change to any estimate.

* Installs without conda: `setup.py` finds Armadillo via `ARMADILLO_DIR`,
  conda, Homebrew or the system prefixes, and stops with install instructions
  if it is missing.
* Stage 1 reports how many phenotype rows it dropped and why, warns when a
  phenotype or covariate column contains `-9`, and warns when the extension was
  built without OpenMP. Which rows are analysed is unchanged.

# gmoq 0.6.1

* MIT licence. Development records removed from the repository.

# gmoq 0.6.0

Both packages now carry the same version.

* **The mixture weight builders default to `density = "hybrid"` whenever the
  weight matrix has a membership (`gamma`) column**, and to `"model"`
  otherwise. Hybrid needs the Stage 1 density, so these calls require
  `stage1_density` (for example `stage1_kde_density(y, q_tau)`). Pass
  `density = "model"` to reproduce 0.5.0.
* `param_gwas(transform = "two_normal")` passes `density` and `stage1_density`
  through to the weight builder.
* Stage 1 installs on macOS, including Apple Silicon, probing for working
  OpenMP flags; it no longer needs a `.bgi` index for a whole-file scan.
* `param_gwas_from_file(rtau_file = "*.gz")` reads gzip without `R.utils`.
* New worked example on 1000 Genomes chr22 genotypes (`docs/examples/1kg_chr22/`).

# gmoq 0.5.0

* The mixture and generic weight builders take
  `density = c("model", "hybrid", "stage1")`; see `docs/MIXTURE_WEIGHT_DENSITY.md`.
* New: `stage1_kde_density()`, `read_stage1_density()` and
  `weight_density_mode()`.
* `param_gwas_from_file()` and `param_gwas_multi_from_file()` gain
  `require_density_mode`.
* `make_weights_generic()` gains `unit_sum` for checking identities such as
  `W_left + W_right = 1`.
