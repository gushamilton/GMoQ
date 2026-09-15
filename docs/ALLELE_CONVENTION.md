# Which allele G-MoQ counts, and how it lines up with REGENIE

Status: measured 2026-09-15 on all 17,550 matched chr22 variants.
**Stage 1 and REGENIE report the same effect allele. No flip is needed when merging.**

An earlier version of this document said the opposite. It was wrong, and wrong in the
harmful direction - following it would have introduced a global sign error into a merge
that is already correct. The correction is recorded in full at the end.

## What Stage 1 does

`stage1/fungwas_stage1/cli.py`, `extract_dosage()`:

```python
# alt_dosage is the expected count of allele[1]
# For diploid, dosage of allele[0] = 2 - alt_dosage
```

Stage 1 counts **`alleles[0]`** - the first allele as written in the BGEN - and reports
that same allele in its `effect_allele` column. Counted and reported agree by
construction, and the docstring says so.

## What REGENIE does

REGENIE reports **`ALLELE1`** as its effect allele, with `A1FREQ` its frequency.

**Without `--ref-first`, REGENIE's `ALLELE1` is the first allele in the BGEN** - the same
one Stage 1 counts. The step-2 command used for this project does not pass `--ref-first`.

## Measured

chr22, `log_bmi`, every variant present in both outputs:

```
matched variants:                          17550
Stage 1 effect_allele == REGENIE ALLELE1 : 17550  (1.0000)
Stage 1 effect_allele == REGENIE ALLELE0 :     0  (0.0000)
```

Example rows - note the pair is the same, printed in the opposite order:

| snp_id | Stage 1 effect | Stage 1 other | REGENIE ALLELE0 | REGENIE ALLELE1 |
|---|---|---|---|---|
| rs1000427 | G | A | A | G |
| rs1000815 | A | G | G | A |
| rs1001022 | T | C | C | T |

**The difference is column order, not allele identity.** Stage 1 writes
`(effect_allele, other_allele)`; REGENIE writes `(ALLELE0, ALLELE1)` with the effect
allele second. A naive column-by-column comparison shows the pair reversed on every
variant, which is what an earlier reading mistook for opposite conventions.

## Merging

Align on the effect allele, as you would with any two summary-statistic files:

```r
merge(gmoq, regenie, by.x = c("snp_id", "effect_allele"),
                     by.y = c("ID", "ALLELE1"))
```

**Do not apply a sign flip.** On this data `effect_allele == ALLELE1` for every variant,
so a flip conditioned on their disagreement is a no-op, and a flip applied unconditionally
reverses every effect.

The condition under which this would change is `--ref-first` in REGENIE step 2, which
swaps which allele REGENIE calls `ALLELE1`. If a future scan uses it, re-run the count
above before merging rather than assuming either convention.

## Evidence the reader itself is right

Realised allele frequency from Stage 1's extracted dosages against REGENIE's `A1FREQ`,
for the allele Stage 1 reports as the effect allele, 10 chr22 variants, n = 461,320:
**maximum difference 0.000000**. Exact frequency agreement is not compatible with counting
the wrong allele. Example: `rs9617528`, Stage 1 AF = 0.745686, REGENIE `A1FREQ` = 0.745686.

## The correction, in full

The earlier version of this document claimed Stage 1 and REGENIE "count the same physical
allele" but label it oppositely, and that merging without a flip "gives systematically
sign-flipped effects". Both halves were wrong: they report the same effect allele, and no
flip is needed.

It also cited the `LOCO_SCALE.md` chr22 comparison, which contains

```r
x[ea != ra, mb := -mb]     # align to REGENIE's effect allele
```

and presented that line as the reason the comparison recovered a correlation of 0.993
rather than -0.993. Since `ea == ra` on every variant, **the line flips nothing**. The
0.993 is correct and unaffected; the explanation attached to it was not.

## What has not been verified

Per-variant agreement of the *sign* of Stage 1's tau-averaged effect against REGENIE's
`BETA`. `LOCO_SCALE.md` reports a correlation of 0.993 over these 17,550 variants, which
implies aggregate agreement. Given exact frequency agreement and a reader that counts and
reports the same allele, a systematic sign error is close to excluded, but the per-variant
check has not been run.
