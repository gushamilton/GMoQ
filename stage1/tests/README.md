# Stage 1 Validation Chain

This folder contains sequential proofs that the jackknife standard errors
match analytic OLS results under progressively more complex models.

`test_bgen_fixture.py` exercises the real Stage 1 CLI on a small known-truth
layout-2, unphased, diploid, 8-bit BGEN fixture.  It checks externally supplied
probabilities with an independent decoder before checking CLI dosages, sample
permutation, missingness, allele labels, and output alignment.  Its scope is the
production 8-bit input boundary.  The file also includes a direct 32-bit regression
for the high-bit dosage fallback, with independent probability truth.

## Order of checks

1) `test_jk_ols_equivalence.py`
   - Standard OLS: `Y ~ G` using block jackknife vs analytic OLS.
   - Uses the same preprocessing as the C++ kernel (imputation, MAC/min-var,
     and covariate residualization via QR).

2) `test_jk_equivalence.py`
   - Stage 1 RIF/quantile GWAS: jackknife SEs vs analytic RIF OLS.
   - Uses the same preprocessing as the C++ kernel.

3) Stage 2 comparison (in the R package)
   - See `stage2/tests/proof_stage2_se.R` for Stage 2 SE mode comparisons.

## Running

From the repo root:

```bash
PYTHONPATH=stage1 python stage1/tests/test_jk_ols_equivalence.py
PYTHONPATH=stage1 python stage1/tests/test_jk_equivalence.py
Rscript stage2/tests/proof_stage2_se.R
```
