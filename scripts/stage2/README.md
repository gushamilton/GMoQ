# Fresh paper Stage 2 preparation

These scripts implement the current-pin paper Stage 2 contract under
`/user/work/fh6520/Gmoq_paper`. They are preparation and runner scripts only;
they do not submit Slurm jobs or allocate BP resources. Run them inside the
single coordinated BP compute allocation.

`prepare_matched_weights.R` is run once for each trait/chromosome. It reads
the Stage 1 header to recover the serialized 45-tau grid, applies the same
complete-case and LOCO alias matching as Stage 1, computes the regression LOCO
scale over all chromosomes before the BGEN sample mask, applies that mask, and
builds W from the resulting phenotype. Stage 1 constructs the RIF on this
LOCO-adjusted phenotype and then residualizes each RIF column on covariates, so
the W baseline is the matched LOCO-adjusted phenotype itself; residualizing the
phenotype before building W would change the estimand.

The output RDS files contain only W, tau values, aggregate baseline parameters,
mixture components where applicable, and non-identifying provenance. No
phenotype, covariate, LOCO, or participant-ID arrays are written. Existing
output files are refused; use a new run directory for every campaign.

For the log-BMI chr22 smoke, supply the frozen aggregate gate values from the
current empirical check:

`--expected-n 461320 --expected-mean 3.296314455837662
--expected-sd 0.1556538337922891 --expected-scale 0.1606963867041275`.
The script stops if any of these values disagree.

Supported families are:

* BMI: `additive`, `shifted`, `both`, and `logscale_vqtl` (shifted
  lognormal weights).
* log-BMI: `additive`, `shifted`, and `both`, with independently generated
  matched log-BMI W. Raw-BMI W is rejected by the provenance/tau checks.
* weight: `additive`, `shifted`, and `both`.
* testosterone and the three shrink phenotypes: `skew`, using the historical
  two-component skew-normal model with components ordered by location.
* `jcv_log_all` and `mcv_log_all`: `mixture`, using the historical
  two-component skew-normal mixture weight construction.

`run_stage2_chr.R` invokes
`fungwasStage2::param_gwas_from_file()` for one chromosome. The primary
source is explicit OLS plus pooled Sigma_0; the sensitivity source is explicit
OLS plus the per-SNP jackknife covariance and IDs sidecar. Pooled and jackknife
outputs carry distinct `cov_source` columns. Two-parameter and mixture fits
retain the package's parameter covariance stream.

`run_stage2_all.sh` prepares W when requested, runs all 22 chromosomes for
one trait/model and one or both covariance sources, and calls
`combine_stage2_results.sh`. The concatenator checks all 22 gzipped chunks
and identical headers, then concatenates the Stage 2 TSV and parameter
covariance gzip members. It never attempts to combine Stage 1 covariance.

A complete main-trait invocation (from the fresh BP checkout) looks like:

```bash
case "$(hostname)" in bp1-compute*) ;; *) exit 2 ;; esac
cd /user/work/fh6520/Gmoq_paper
source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh
conda activate freshRlib
bash review/paper-prep/stage2/run_stage2_all.sh \
  --trait log_bmi --model both \
  --stage1-dir .local/runs/<stage-id>/log_bmi \
  --loco results/paper/regenie_step1/full/log_bmi/log_bmi_step1_pred.list \
  --pheno data/derived/paper/phenotypes_paper.txt \
  --covar data/derived/paper/covars_with_sex_paper.txt \
  --keep data/derived/paper/keep_files/keep_log_bmi.txt \
  --r-lib .local/freshRlib \
  --weights-root .local/runs/<run-id>/weights \
  --out-root .local/runs/<run-id>/stage2 \
  --cov-source pooled,jackknife --prepare-weights \
  --threads "${SLURM_CPUS_PER_TASK:-4}"
```

The Stage 1 directory must contain
`<trait-with-dots-replaced-by-underscores>_chr01.stage1.tsv.gz` through
`chr22`, with sibling `.sigma0.tsv.gz`, `.sigma0.meta.json`,
`.cov.gz`, `.cov.ids.tsv.gz`, `.cov.scale.gz` where applicable, and
`.cov.meta.json`. For antibody traits use the antibody phenotype/covariate/
keep files and the matching antibody LOCO list; pass its Stage 1 directory
explicitly. The antibody all-eligible and seropositive branches are outside
this core runner.

The four requested main families can be run independently by changing only
`--trait` and `--model`:

```
bmi:       additive shifted both logscale_vqtl
log_bmi:   additive shifted both
weight:    additive shifted both
testosterone, testosterone_shrink_0.2, testosterone_shrink_0.5, testosterone_shrink_0.7: skew
jcv_log_all, mcv_log_all: mixture
```

The matched-baseline gate was run on 2026-09-15 through held BP allocation
18968780 (`bp1-compute053`) with the code-only copy at
`.local/stage2-validation-source/prepare_matched_weights.R`:

```bash
module load languages/R/4.4.3
cd /user/work/fh6520/Gmoq_paper
Rscript --vanilla .local/stage2-validation-source/prepare_matched_weights.R \
  --trait log_bmi --chr 22 \
  --stage1 .local/paper-empirical-smoke-20260915-01/log_bmi_chr22.stage1.tsv.gz \
  --pheno /user/work/fh6520/fungwas_paper/data/derived/paper/phenotypes_paper.txt \
  --covar /user/work/fh6520/fungwas_paper/data/derived/paper/covars_with_sex_paper.txt \
  --keep /user/work/fh6520/fungwas_paper/data/derived/paper/keep_files/keep_log_bmi.txt \
  --loco /user/work/fh6520/fungwas_paper/results/paper/regenie_step1/full/log_bmi/log_bmi_step1_pred.list \
  --bgen-sample /bp1/mrcieu1/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr1-22.sample \
  --family shifted_lognormal \
  --out-dir .local/stage2-validation-20260915-01/weights/log_bmi/chr22 \
  --loco-scale regression --r-lib .local/Rlib \
  --expected-n 461320 --expected-mean 3.296314455837662 \
  --expected-sd 0.1556538337922891 --expected-scale 0.1606963867041275
```

It passed and wrote `weights_additive.rds`, `weights_shifted.rds`, and
`weights_both.rds` plus `weights_manifest.tsv`. The manifest records
`N=461320`, factor `0.160696386704128`, mean `3.29631445583766`, and SD
`0.155653833792289`. The first run exposed that the BP R 4.5.1 module lacks
`R.utils`; the recorded command uses the existing R 4.4.3 module, which has
the required dependency. No package or R library was modified.

Stage 1 array outputs are nested by trait and chromosome. Prefer the explicit
root form, which resolves exactly one file per chromosome and errors if a
second candidate also exists:

```bash
--stage1-root .local/runs/<stage1-id>
```

The resolver looks for
`<root>/<trait-tag>/chrNN/<trait-tag>_chrNN.stage1.tsv.gz`, then the equivalent
literal trait directory and trait-specific `chrNN` root. The legacy
`--stage1-dir` form remains available for a flat trait directory.

The testosterone `skew` builder uses the exact Stage 1 Silverman KDE density at
the empirical quantiles when constructing W, matching the density that scaled
the Stage 1 RIF. When the canonical 45-level grid is identified from the
three-decimal TSV labels, the builder also recovers the full-precision grid
used by the production runner; custom grids remain at their serialized
precision. The correction and tau-grid provenance are recorded in the weight
RDS. Keep full logs and all individual-level inputs on BP; only aggregate W,
baseline summaries and Stage 2 result files belong in any exportable evidence
bundle.
