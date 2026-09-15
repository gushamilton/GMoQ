# Tier 4: Real Chromosome 22 Validation

Historical validation scaffold, retained for reference. These scripts are not a
completed scientific validation: LD-linked noncausal variants need not have null
marginal effects, and generated scale effects must be converted to the fitted
parameter's units before judging recovery. The criteria and timings below are
historical proposals, not verified release evidence.

On BluePebble, run every analysis, interpreter invocation and environment check
inside a compute-node allocation. Login nodes are for inspection, file work and
job control only. Reuse one held allocation and cap threads to its resources.

## Overview

Tiers 0-3 use synthetic genotypes generated in-process, so they never touch the BGEN reader, dosage extraction, variant filtering, or sample matching code paths. Tier 4 exercises the complete input pipeline on real UK Biobank data while maintaining perfect knowledge of the phenotypic model.

## What's included

### Scripts

- **`t4_make_pheno.R`** — Simulate phenotypes on real chr22 data
  - Four test cases: null, location, scale, mixed
  - Output: phenotype files in Stage 1 input format
  - Metadata: ground truth effects and variant lists saved as JSON

- **`t4_stage1.sh`** — Run Stage 1 on each simulated phenotype
  - Inputs: phenotype file, real chr22 BGEN, real samples/covariates
  - Outputs: Stage 1 results (dosages, covariance, etc.)
  - Requires: Stage 1 conda environment

- **`t4_stage2.R`** — Run Stage 2 on Stage 1 outputs
  - Quantile GWAS for each tau
  - Tabular results for validation

- **`t4_analyze.R`** — Validate results against ground truth
  - Effect recovery (sign, magnitude, SE coverage)
  - Type-I rate and lambda_GC
  - Model selection and comparison

## Quick start

### Prerequisites

```bash
# Inside a BluePebble compute-node allocation:
source /user/home/fh6520/miniforge3/etc/profile.d/conda.sh
conda activate fungwas-stage1
module load languages/R/4.5.1

export GMOQ_WORK=/user/scratch/fh6520/gmoq_chr22  # or your scratch
export GMOQ_DEP_LIB=/user/work/fh6520/r_packages
```

### Run phenotype generation

```bash
cd /user/work/fh6520/wt_testsuite

# Generate all 4 phenotypes (null, location, scale, mixed)
bash tests/run_tier4.sh
```

This creates phenotype files in `$GMOQ_WORK/`:
- `pheno_null.txt`
- `pheno_location.txt`
- `pheno_scale.txt`
- `pheno_mixed.txt`

Each also has a `.meta.json` file with ground truth details.

### Run Stage 1 (cluster compute)

For each phenotype case:

```bash
bash tests/tier4/t4_stage1.sh <case> $GMOQ_WORK $GMOQ_WORK/stage1_<case>
```

Example:
```bash
for case in null location scale mixed; do
  sbatch --account=sscm013902 --partition=short --cpus-per-task=8 \
    -o $GMOQ_WORK/stage1_${case}.log \
    -e $GMOQ_WORK/stage1_${case}.err \
    bash tests/tier4/t4_stage1.sh $case $GMOQ_WORK $GMOQ_WORK/stage1_$case
done
```

Stage 1 takes ~15 minutes per case on 8 CPUs.

### Run Stage 2 and validation

After all Stage 1 jobs complete:

```bash
for case in null location scale mixed; do
  Rscript --vanilla tests/tier4/t4_stage2.R \
    --stage1-prefix $GMOQ_WORK/stage1_${case} \
    --case $case \
    --out $GMOQ_WORK/results_${case}.tsv
done
```

Then analyze:

```bash
Rscript --vanilla tests/tier4/t4_analyze.R \
  --work-dir $GMOQ_WORK \
  --out $GMOQ_WORK/validation_results.json
```

## Test cases explained

### 1. **Null** — phenotype ~ covariates only

**Ground truth:** No genetic effects; all association should be random at p = 0.05 level.

**Expected results:**
- Type-I error rate: 0.03–0.07 (5% ± 2 pp)
- lambda_GC: 0.95–1.05
- All z-scores distributed as standard normal

### 2. **Location** — phenotype ~ covariates + genetic location effects

**Ground truth:** 10 causal variants with effect sizes β ~ N(0, 0.03²). No scale effects.

**Expected results:**
- ≥ 8 of 10 causal variants detected with p < 0.05
- Estimated effects within 2 SE of truth
- Sign and allele orientation correct
- No spurious scale hits beyond documented bias

### 3. **Scale** — heteroscedastic phenotype with genetic scale effects

**Ground truth:** Same 10 causal variants with scale effects γ ~ N(0.01, 0.002²). No location effects.

**Expected results:**
- ≥ 8 causal scale effects detected with p < 0.05
- Location estimates near 0 at causal variants (p >> 0.05)
- Scale effect magnitudes and directions match truth

### 4. **Mixed** — both location and scale effects

**Ground truth:** Same 10 causal variants with both effects simultaneously.

**Expected results:**
- Both location and scale detected at causal variants
- Effect sizes match individual cases (location and scale)
- Model selection prefers location-scale over either alone

## Input data (real)

All from UK Biobank, accessible on the HPC cluster:

- **BGEN file:** `/bp1/mrcieu1/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr22.bgen`
  - ~650k samples, ~17.5M variants (biallelic)

- **Sample file:** `.../data.chr1-22.sample`
  - Sample IDs and metadata

- **SNP list:** `/user/work/fh6520/fungwas_paper/inputs/data/w_hm3_chr22_only.txt`
  - ~17.5k HapMap3 variants on chr22 (filtering done by paper analysis)

- **Covariates:** `/user/work/fh6520/fungwas_paper/data/derived/paper/covars_with_sex_paper.txt`
  - Real covariates (age, sex, principal components, etc.)

- **Keep file:** `/user/work/fh6520/fungwas_paper/data/derived/paper/keep_files/keep_bmi.txt`
  - Sample subset for this analysis

## Output structure

After a full run, `$GMOQ_WORK/` contains:

```
$GMOQ_WORK/
├── pheno_null.txt                  # Simulated phenotypes
├── pheno_null.meta.json            # Ground truth metadata
├── pheno_location.txt
├── pheno_location.meta.json
├── pheno_scale.txt
├── pheno_scale.meta.json
├── pheno_mixed.txt
├── pheno_mixed.meta.json
├── stage1_null.*                   # Stage 1 outputs (large)
├── stage1_location.*
├── stage1_scale.*
├── stage1_mixed.*
├── results_null.tsv                # Stage 2 outputs
├── results_location.tsv
├── results_scale.tsv
├── results_mixed.tsv
└── validation_results.json         # Final analysis results
```

## Performance

- **Phenotype simulation:** ~2 min (all 4 cases, compute node)
- **Stage 1 per case:** ~15 min (8 CPUs, linear in variant count)
- **Stage 1 all cases:** ~1 hour (if serial) or ~15 min (if parallel on 4 nodes)
- **Stage 2 per case:** ~3 min (compute node)
- **Analysis:** ~1 min (compute node)

**Total wall time:** ~1.5–2 hours; **CPU hours:** ~6

## Pass criteria (blocking)

1. **Null Type-I rate:** 0.03–0.07
2. **Null lambda_GC:** 0.95–1.05
3. **Location: ≥8/10 effects recovered**, correct sign/allele
4. **Location: ≤2 spurious scale hits**
5. **Scale: ≥8/10 scale effects detected**
6. **Scale: location null at causal variants**
7. **Mixed: both effects at ≥8 causal variants**

## Known limitations

1. **Single chromosome:** Only chr22, not genome-wide. Type-I and lambda_GC are valid for chr22 LD; genome-wide validation requires all 22 chromosomes.

2. **Simple phenotype models:** Linear (location) and multiplicative (scale). Real traits often have:
   - Non-additive gene action
   - Gene × environment interactions
   - Polygenicity (1000s of weak causal variants)

3. **Simulation vs. real genetics:** Simulated effects are applied to simulated dosages (from a beta distribution of MAF), not to real BGEN dosages. This is a practical compromise: reading the full real BGEN in R is complex, and the simulation preserves the MAF distribution and LD structure of real data.

4. **No baseline comparison:** Results are not compared to REGENIE or other tools on the same data. Such a comparison is valuable but would require real phenotypes (not simulated).

## Future enhancements

1. **Read real BGEN dosages:** Use bgen_reader or bgenix to actually read real genotypes, then apply simulated effects on top. This would make the validation completely end-to-end.

2. **Genome-wide validation:** Run the same pipeline for all 22 chromosomes and aggregate Type-I rates and lambda_GC.

3. **Compare to REGENIE:** Use real phenotypes from the paper and compare G-MoQ tau-averaged effects to REGENIE BETA, as done in `docs/LOCO_SCALE.md`.

4. **Pipeline correctness checks:** Explicit tests for:
   - BGEN dosage agreement with independent reader
   - Sample matching with intentional scrambling
   - Variant filtering (multi-allelic, MAF, INFO)
   - Missing data handling
   - Covariance fingerprint validation

5. **Power analysis:** Quantify the effect-size-to-power curve for various MAF and model settings.

## References

- **`docs/CHR22_VALIDATION.md`** — full technical specification
- **`docs/LOCO_SCALE.md`** — LOCO fix validated on chr22; comparison to REGENIE
- **`docs/DISPERSION_SKEW_CORRECTION.md`** — positive bias in dispersion estimates
- **`docs/KNOWN_ISSUES.md`** — known defects in shipped code
