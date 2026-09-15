# Which allele G-MoQ counts, and how it differs from REGENIE

Status: verified from source and against real chr22 data, 2026-09-15.
**Stage 1's reader is correct. Its allele convention is the opposite of REGENIE's.**

## What Stage 1 does

`stage1/fungwas_stage1/cli.py`, `extract_dosage()`:

```python
# alt_dosage is the expected count of allele[1]
# For diploid, dosage of allele[0] = 2 - alt_dosage
```

Stage 1 counts **`alleles[0]`** - the first allele as written in the BGEN - and reports
that same allele in its `effect_allele` output column. Counted and reported agree by
construction, and the docstring says so.

## What REGENIE does

REGENIE reports **`ALLELE1`**, the second allele, as its effect allele, with `A1FREQ`
its frequency.

## The consequence

Both conventions are internally consistent. They are opposite. On 10 chr22 variants
checked, Stage 1's reported allele pair is flipped relative to REGENIE's in 10 of 10
cases - `T/C` against `C/T` - while both count the same physical allele.

**Merging G-MoQ and REGENIE summary statistics on allele columns without harmonising
gives systematically sign-flipped effects.** The comparison in `LOCO_SCALE.md` that
recovered a correlation of 0.993 harmonised explicitly:

```r
x[ea != ra, mb := -mb]     # align to REGENIE's effect allele
```

Without that line the same comparison returns approximately -0.993, which would look
like a catastrophic disagreement rather than a labelling difference.

## Evidence the reader is right

Realised allele frequency computed from Stage 1's extracted dosages, against REGENIE's
`A1FREQ` for the allele Stage 1 reports as the effect allele, 10 chr22 variants,
n = 461,320: **maximum difference 0.000000**. Exact frequency agreement is not
compatible with counting the wrong allele, so the reader counts what it reports.

Example: `rs9617528`, Stage 1 AF[alleles[0]] = 0.745686, REGENIE A1FREQ = 0.745686.

## What has NOT been verified

Per-variant agreement of the *sign* of Stage 1's tau-averaged effect against REGENIE's
`BETA` after harmonising alleles. `LOCO_SCALE.md` reports a correlation of 0.993 over
17,550 chr22 variants, which implies aggregate sign agreement, but the per-variant check
was not run. Given exact frequency agreement and a reader that counts and reports the
same allele, a systematic sign error is close to excluded; a per-variant confirmation
would close it entirely.

## Note for the validation suite

The tier-4 chr22 suite builds its phenotype from real dosages. If it extracts those
dosages with **Stage 1's own reader**, an allele flip in that reader would cancel: the
phenotype would be built flipped and recovered flipped, and the test would pass. So the
suite validates Stage 1's *analysis* pipeline, not its reader.

That is the right trade. Use Stage 1's reader - it is the tested one, it seeks with the
`.bgi` index, and a parallel implementation is a second thing to get wrong - and keep the
reader honest with the external frequency check above, which is a one-off comparison
rather than a duplicate code path.

`tests/tier4/extract_dosages.py` currently reimplements extraction with the Python `bgen`
package and no index seeking: 69 seconds for 10 variants at n = 461,320, consistent with a
linear scan of the whole file. `fungwas_stage1.cli.extract_dosage` is directly importable
and should replace it.
