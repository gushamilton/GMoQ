# Effect alleles and comparison with REGENIE

Verified on the saved chr22 log-BMI results on 2026-09-15. **Stage 1 and
REGENIE report the same effect allele in all 17,550 matched variants in this
comparison. No sign flip is needed.** The previous version of this document
incorrectly claimed the conventions necessarily oppose one another.

## Stage 1 convention

`stage1/fungwas_stage1/cli.py::extract_dosage()` computes the expected count of
the first BGEN allele, `alleles[0]`. For diploid biallelic records at up to
8-bit probability precision it uses `2 - variant.alt_dosage`. Higher-bit
records use `2*P(AA) + P(Aa)` from decoded probabilities to avoid an observed
dosage error in the installed `bgen` library. Bit depth is checked per record,
including files with mixed depths. That allele is written to `effect_allele`;
the other is written to `other_allele`.

## Harmonisation rule

Compare allele **labels**, not their column positions or names across formats.
For a matched variant, compare Stage 1 `(effect_allele, other_allele)` with
REGENIE `(ALLELE1, ALLELE0)`:

- Same pair: retain the effect sign.
- Reversed pair: reverse the effect sign when aligning to REGENIE.
- Neither: investigate or exclude the incompatible pair; do not assume a flip.

REGENIE column ordering can depend on its input/settings, including `--ref-first`.
Do not extrapolate the orientation of one dataset to all REGENIE files. Strand
ambiguity and variant identity still require deliberate handling in comparisons
across datasets; the check below concerns the same chr22 data.

## Reproduced result

| Check | Result |
|---|---:|
| Matched variants with finite effects | 17,550 |
| Same effect and other allele | 17,550 |
| Reversed pair | 0 |
| Incompatible pair | 0 |
| Correlation of mean quantile beta with REGENIE BETA before alignment | 0.9928536466 |
| Correlation after label-based alignment | 0.9928536466 |

The earlier `x[ea != ra, mb := -mb]` harmonisation line is a no-op for these
results. Omitting it does **not** produce a correlation of -0.993.

The comparison used a chr22 log-BMI Stage 1 scan and the matching REGENIE step 2
output. It can be reproduced with `scripts/validation/check_allele_alignment.py
--stage1 ... --regenie ... --out allele_alignment.json` from the tag
`archive/research-records-20261016`; the script reports aggregates and rejects
duplicate variant IDs.

This establishes alignment of the saved results, not independent correctness
of all dosage decoding or unbiasedness of the estimator. The public BGEN fixture
tests should independently establish decoding from known genotype probabilities.
Generating phenotypes with the same reader used in analysis cannot by itself
detect a consistent allele reversal in that reader.
