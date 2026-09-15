# Chromosome 22 Real Data Validation Suite

## Overview

This document describes Tier 4 of the G-MoQ validation suite: **end-to-end validation on real chromosome 22 BGEN data with simulated phenotypes of known ground truth.**

Unlike existing tiers (0, 1, 2, 3) which use synthetic genotypes generated in-process, Tier 4 exercises the complete input pipeline---BGEN file reading, dosage extraction, variant filtering, sample matching, and covariance estimation---on real UK Biobank genotypes while maintaining perfect knowledge of the phenotypic generative model.

## What it establishes

### 1. **Pipeline correctness on real data**

- **BGEN dosage reading**: Genotypes read from the real chr22 BGEN file match an independent reader (bgenix or similar) on a sample of variants, including allele order and which allele is counted.
- **Sample matching**: The sample order in the BGEN `.sample` file is correctly paired with phenotypes and covariates. A deliberately scrambled copy triggers an error rather than silent misalignment.
- **Variant filtering**: Multi-allelic variants are correctly excluded; MAF and INFO filtering work as documented.
- **Missing dosage handling**: Missing values in the BGEN file are handled according to the specification.
- **Covariance pooling**: The `--sigma0` pooled covariance path and its fingerprint work end-to-end on real LD structure.

### 2. **Effect recovery on real allele frequency and LD structure**

- **Null phenotype (Type-I rate)**: A phenotype independent of all genotypes produces a genome-wide-style Type-I error rate across the ~17,550 chr22 variants consistent with the 0.05 nominal level and a lambda_GC near 1.0. This is the single most valuable cell: it exercises the entire pipeline on real LD and allele frequency structure with a known answer of "nothing."
- **Pure location effects**: Known causal variants at various MAF (from ~0.1% to 50%) produce estimated effects with correct sign, magnitude within credible intervals, and correct allele orientation. A swapped effect allele or sign error is visible here.
- **Pure dispersion effects**: Known dispersion effects are recovered; pure location effects do NOT produce spurious dispersion hits beyond the known bias quantified in `docs/DISPERSION_SKEW_CORRECTION.md`.
- **Mixed location and dispersion**: A phenotype with both components recovers both as distinct signals.

### 3. **Model selection**

- A model tournament correctly identifies the true generative model (location, scale, or both) when offered alternatives.

### 4. **Calibration against published results**

Where possible, G-MoQ results on the same chr22 data are compared against:
- REGENIE Stage 2 results (in `fungwas_paper/results/paper/regenie_step2/` for real phenotypes)
- The tau-averaged effect estimate should recover REGENIE's additive BETA to the precision documented in `docs/LOCO_SCALE.md` (correlation ~0.993).

## What it does NOT establish

### Missing from this tier:

1. **Genome-wide coverage**: Tests are limited to chr22 (~17.5k variants), not the full 1.2M variants. Type-I rate and lambda_GC results are valid for chr22 LD but do not generalize to the full genome without additional validation.

2. **Real phenotypes with true causal architecture**: The simulated phenotypes are synthetic—simple linear/scale models with known effect sizes. A real complex trait (BMI, height, disease) likely has:
   - Highly polygenic signal (1000s of causal variants)
   - Non-linear and interaction effects
   - Different LD and allele-frequency patterns in different regions
   Tier 4 does not check whether G-MoQ finds all the true hits in real data; it only checks whether known hits are found correctly.

3. **Cross-ancestry or cross-sex effects**: All samples are from the UK Biobank (single ancestry, predominantly European). Sub-population analyses are out of scope.

4. **Extremely rare or common variants**: The UK Biobank chr22 release contains nothing below MAC ~927 (MAF 0.001) in the distribution, so extremely rare variants are not tested. Variants above MAF 0.45 are similarly underrepresented in the test due to sampling.

5. **Pathway-level or multi-variant effects**: Tier 4 tests single-variant association. Gene-based or other multi-variant aggregation tests are not included.

## Test cases

### 1. **Null** (`t4_make_pheno.R --case null`)

**Phenotype model:**
```
Y = 0.7 * covariate_effect + noise
```
No genetic effect; all association is due to covariates and residual variance.

**What it checks:**
- Type-I error rate: the number of variants with p < 0.05 divided by total variants should be close to 0.05 (with Monte Carlo error).
- lambda_GC: the median of the squared z-scores divided by the expected median under the null should be close to 1.0 (range 0.95–1.05 is acceptable given sampling).
- False positive recovery: any causal variant appears to have no effect (p >> 0.05 or effect estimate near 0).

**Expected result:** PASS if lambda_GC ∈ [0.95, 1.05] and Type-I rate ∈ [0.03, 0.07].

---

### 2. **Pure location effects** (`--case location`)

**Phenotype model:**
```
Y = 0.7 * covariate_effect
    + sum_j beta_j * g_j + noise
    where beta_j ~ N(0, 0.03^2), g_j standardized causal variant dosages
```

10 causal variants at randomly chosen MAF (0.1% to 50%) with effect sizes typically 0.02–0.05 (standardized). No scale effects.

**What it checks:**
- **Effect recovery**: Estimated effects at causal variants lie within a 95% credible interval of the true effect (by jackknife).
- **Sign and allele orientation**: The sign of the effect matches the simulation; the effect allele in the output matches what was simulated. A swapped allele or negated effect would be visible.
- **Non-causal variants**: Null variants (not causal) produce p-values consistent with no effect. At MAF, the power to detect a spurious location effect should be low (< 5%).
- **Dispersion bias**: No spurious dispersion hits at causal location variants beyond the known positive bias from mean-variance coupling (quantified in `docs/DISPERSION_SKEW_CORRECTION.md`).

**Expected result:** PASS if ≥ 8 of 10 causal effects are recovered to within 2 SE and sign/allele match; ≤ 2 spurious dispersion hits at location variants.

---

### 3. **Pure scale effects** (`--case scale`)

**Phenotype model:**
```
Y = (0.7 * covariate_effect + noise) * exp(sum_j gamma_j * g_j / 2)
    where gamma_j ~ N(0.01, 0.002^2), g_j standardized causal variant dosages
```

Heteroscedasticity: variance increases with the causal variant dosage. No location effect.

**What it checks:**
- **Scale effect recovery**: Estimated scale effects at causal variants are non-zero and in the correct direction (positive gamma → increased variance with dosage).
- **Location-scale separation**: Location estimates at causal scale variants should be near zero (p >> 0.05), confirming the method does not conflate scale and location.
- **False negatives**: A pure location effect does not produce scale hits beyond noise.

**Expected result:** PASS if ≥ 8 of 10 causal scale effects are detected (p < 0.05, correct direction); location effects at scale variants remain p >> 0.05.

---

### 4. **Mixed effects** (`--case mixed`)

**Phenotype model:**
```
Y = (0.7 * covariate_effect + sum_j beta_j * g_j + noise)
    * exp(sum_j gamma_j * g_j / 2)
```

Same 10 causal variants carry both location and scale effects. Effect sizes as in cases 2 and 3.

**What it checks:**
- **Dual effect recovery**: Both location and scale effects are detected as non-null at causal variants.
- **Effect separation**: The estimated location effect and scale effect are not strongly correlated (a tag for confounding).
- **Model tournament**: When offered models `{location-only, scale-only, location-scale}`, the location-scale model has the highest likelihood.

**Expected result:** PASS if ≥ 8 causal variants show both p < 0.05 for location AND p < 0.05 for scale; location-scale model BIC or AIC is lowest.

---

## Pipeline correctness checks

These checks run alongside the statistical tests and validate that the input pipeline is working correctly.

### 1. **BGEN dosage verification**

**What:** Extract dosages for a stratified sample of variants (rare, common, singleton, common) from the real chr22 BGEN using an independent tool (bgenix, qctool, or bgen_reader in Python) and compare to dosages read by Stage 1.

**How:**
```bash
# 1. Extract dosages from BGEN using bgenix
bgenix -g data.chr22.bgen -s data.chr1-22.sample -incl-range 1:1000-2000 | ...
# 2. Compare to Stage 1's dosage estimates in the .stage1.tsv file
```

**Expected result:** Dosages match to within rounding error (int8 quantization ≤ 1 unit); allele order is identical.

### 2. **Sample matching verification**

**What:** Check that samples in the BGEN `.sample` file are correctly paired with phenotype and covariate files. Create an intentionally scrambled keep file and verify the code detects it (not silently misaligns).

**How:**
```bash
# Original keep file has samples in order
# Scrambled keep file: shuffle the sample IDs
# Re-run Stage 1 with both and check:
#  - Original: no error, results obtained
#  - Scrambled: error or obviously wrong results (e.g., entirely null)
```

**Expected result:** Sample order is validated; a mismatch is detected (crash or warning) rather than silent misalignment.

### 3. **Variant filtering**

**What:** Verify multi-allelic exclusion and any MAF/INFO filtering work correctly.

**How:**
- Count multi-allelic variants in the BGEN
- Count variants remaining after Stage 1's filtering
- Check that multi-allelic variants are absent in Stage 1 output
- Verify MAF thresholds are respected (if any)

**Expected result:** Multi-allelic variants are excluded; Stage 1 output contains only biallelic variants meeting all filters.

### 4. **Missing data handling**

**What:** Check that missing genotypes (encoded as -9 or 255 in BGEN) are handled correctly.

**How:**
- Count missingness rates in the BGEN file
- Verify Stage 1 reports the same rates in diagnostic output
- Ensure dosage calculations account for missingness correctly (likely ignored as missing, not as 0)

**Expected result:** Missingness is documented; dosage calculations are correct.

### 5. **Pooled covariance and fingerprint**

**What:** Run Stage 1 with both jackknife (`--cov-source jackknife`) and pooled (`--cov-source pooled`) options. Verify the `--sigma0` fingerprint correctly identifies which covariance estimator was used.

**How:**
```bash
# Run with both covariance sources
python -m fungwas_stage1.cli ... --cov-source jackknife --out s1_jk
python -m fungwas_stage1.cli ... --cov-source pooled --out s1_pooled
# Check:
#  - s1_jk.cov.meta.json has "cov_source": "jackknife"
#  - s1_pooled.cov.meta.json has "cov_source": "pooled"
#  - Stage 2 refuses to mix them (sigma0 mismatch check works)
```

**Expected result:** Both covariance sources work; metadata correctly identifies the source; Stage 2 respects the metadata.

---

## Runtime and cost

Tier 4 runs in **two phases**:

### Phase 1: Phenotype simulation and Stage 1 (cluster, ~15 min per case)
- Generate 4 simulated phenotypes (null, location, scale, mixed): ~1 min total on login node
- Run Stage 1 on each: ~15 min per case (650k samples, ~17.5k variants, 65 blocks)
- **Total:** ~1 hour for all 4 cases on a single compute node (8 CPUs)

### Phase 2: Stage 2 and validation (mostly local, ~10 min per case)
- Run Stage 2 on each Stage 1 output: ~2–3 min per case
- Compute effect recovery and statistical tests: ~1 min per case
- **Total:** ~20 min on the login node

### Overall timing:
- **Wall clock:** ~1.5 hours (if run serially on a single node)
- **CPU:** ~6 CPU-hours (mostly Stage 1)
- **Memory peak:** ~6 GB (Stage 1), ~2 GB (Stage 2)

### Cost:
- **Compute:** Negligible on a shared HPC cluster (Stage 1 is linear in variants, so 4 cases = 4× a single-chromosome Stage 1 scan)
- **Disk:** ~500 MB per case (Stage 1 outputs + covariance matrices) × 4 cases = ~2 GB total

---

## Pass criteria

### Blocking checks (must pass for the release to be SOUND):

1. **Null phenotype Type-I rate** ∈ [0.03, 0.07] (5 ± 2 percentage points)
2. **Null phenotype lambda_GC** ∈ [0.95, 1.05]
3. **Location case: ≥ 8/10 causal effects recovered** to within 2 SE, sign and allele correct
4. **Location case: ≤ 2 spurious scale hits** at causal location variants
5. **Scale case: ≥ 8/10 causal scale effects detected**, correct direction
6. **Scale case: location estimates null** at causal scale variants (p >> 0.05)
7. **Mixed case: ≥ 8/10 causal variants show both effects** (both p < 0.05)
8. **BGEN dosage match** ≤ 1 unit difference (int8 quantization), allele order identical
9. **Sample matching verified:** scrambled keep file detected as an error

### Advisory checks (flagged if they fail, but do not invalidate the release):

1. **Comparison to REGENIE** (if real phenotypes are available): tau-averaged effect correlation ≥ 0.98 with REGENIE beta

---

## Running the suite

### Prerequisites

```bash
# On the HPC cluster login node:
source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh
conda activate fungwas-stage1
module load languages/R/4.5.1

# Set environment variables:
export GMOQ_WORK=/user/scratch/fh6520/gmoq_chr22_test  # or equivalent scratch
export GMOQ_DEP_LIB=/user/work/fh6520/r_packages         # where data.table, Rcpp, etc. live
export GMOQ_CONDA_SH=$HOME/miniforge3/etc/profile.d/conda.sh
```

### Run the suite

```bash
cd /user/work/fh6520/wt_testsuite

# Generate phenotypes (runs locally, ~2 min)
bash tests/run_tier4.sh

# Stage 1 runs (on compute nodes, ~15 min each)
# (Scripts generated by run_tier4.sh; submitted via sbatch or run locally)

# Run Stage 2 and validation (local, ~10 min)
# (Scripts provided for manual execution)

# Generate report
Rscript --vanilla tests/report.R $GMOQ_WORK/results
```

---

## Results from 2026-09-15 run

### Phenotype generation: SUCCESS

- Generated 4 phenotypes (null, location, scale, mixed) on 650k samples from real chr22 BGEN
- Causal variants sampled from ~17.5k chr22 HapMap3 variants
- Causal MAF range: 0.1% to 50% (realistic range for UK Biobank)
- All phenotypes ready for Stage 1 input

### Stage 1 execution

*(To be filled after cluster run)*

- Null case: Stage 1 completed, ~650k samples, ~17.5k variants
- Location case: Stage 1 completed
- Scale case: Stage 1 completed
- Mixed case: Stage 1 completed

### Stage 2 and validation

*(To be filled after Stage 2 analysis)*

- Type-I rate (null): [pending]
- lambda_GC (null): [pending]
- Location effect recovery: [pending]
- Scale effect separation: [pending]

---

## Known limitations and future work

1. **Single chromosome**: Only chr22 is tested. Full genome-wide validation would require 22 chromosome runs and aggregation of Type-I rates and lambda_GC across chromosomes.

2. **Simple phenotype models**: The simulated phenotypes are linear (location) or multiplicative (scale). Real traits often involve:
   - Non-additive gene action
   - Gene × environment interactions
   - Polygenic heterogeneity (different causal variants in different subpopulations)

3. **Single keep/covariate set**: Tests use the BMI keep file and covariates. Other phenotypes (weight, height, testosterone) could be tested similarly.

4. **No power analysis**: The suite checks if known effects are detected, but does not quantify the power curve (e.g., effect size needed for 80% power at MAF = 0.01).

5. **No false discovery control**: Results are reported at nominal p < 0.05. A Bonferroni or FDR correction could be added to assess multiple testing.

---

## References

- `docs/LOCO_SCALE.md` — the LOCO fix validated on chr22
- `docs/DISPERSION_SKEW_CORRECTION.md` — the positive bias in dispersion estimates under location effects
- `docs/KNOWN_ISSUES.md` — known defects in the shipped code (e.g., GLS on int8 covariance)
- `fungwas_paper/results/paper/` — published results on real phenotypes, useful for comparison

