# Effect allele convention

## Stage 1

Stage 1 counts the **first allele of each BGEN record** (`alleles[0]`). The genotype value for a
person is the expected number of copies of that allele, `2 P(AA) + P(Aa)`, computed from the
genotype probabilities in the BGEN file. Every reported effect is per copy of that allele, which
the output names in `effect_allele`; the other allele is `other_allele`.

Bit depth is checked record by record, so files mixing 8-bit and higher-precision records are
decoded correctly. Missing calls are replaced by the variant's mean dosage.

Which allele comes first is decided when the BGEN file is written. With PLINK 2,
`--export bgen-1.2 ref-first` puts the reference allele first, so the effect allele is the
reference allele (the 1000 Genomes example does this). Check the `effect_allele` column rather
than assuming.

The test suite checks the convention on BGEN files written from known genotype probabilities.

## Comparing with other software

Compare allele **labels**, not column positions, because column order differs between programs
and can depend on their settings (REGENIE's `--ref-first`, for example). For a variant present in
both results, compare the G-MoQ pair `(effect_allele, other_allele)` with the other program's
`(tested allele, other allele)` (`ALLELE1`, `ALLELE0` in REGENIE):

- same pair: keep the sign;
- reversed pair: reverse the G-MoQ sign;
- anything else: investigate or exclude the variant; do not assume a flip.

Strand-ambiguous variants (A/T, C/G) and variant identity across datasets still need deliberate
handling when the two analyses did not use the same genotype files.
