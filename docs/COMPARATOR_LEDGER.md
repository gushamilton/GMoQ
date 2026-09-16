# Ledger of G-MoQ vs REGENIE comparisons

**Status: comparator rebuilt 2026-09-16. Run ID
`stratified-loco-comparator-20260916-01`.**

This document is the register of *every* G-MoQ-versus-REGENIE estimate produced
in this project. Its purpose is that a reader can see at a glance which
comparisons are like-for-like and which are not.

Several tables elsewhere in this repository were computed against a REGENIE that
was run with **no LOCO adjustment at all** (`--ignore-pred`). Those tables are
not deleted; they are marked here and at their source, and the numbers against
the rebuilt comparator are given beside them.

---

## 1. The asymmetry, which governs every design choice

G-MoQ is given **only the phenotype**. It must infer sex membership from the
shape of the distribution. Sex-stratified REGENIE is **handed the ground-truth
labels** and told which stratum every participant belongs to.

The fair comparison gives each method its best configuration given what it is
allowed to know. The comparator must therefore not be handicapped: it gets a
proper **within-stratum LOCO**, which is what anyone holding the labels would
do. REGENIE has strictly more information and should win; the question this
ledger answers is *by how much*.

Fitting the ridge **within** a sex stratum also sidesteps the pathology in
[`LOCO_MIXTURE_HAZARD.md`](LOCO_MIXTURE_HAZARD.md): a within-stratum step 1
cannot learn sex, because sex is constant inside the stratum. That is precisely
why this is the comparator's best case rather than a repeat of the broken arm.

---

## 2. Comparator (REGENIE) configurations

All step 2 runs are chr22, HapMap3 variants (`w_hm3_rsids_only.txt`), REGENIE
v4.1.2, `--bsize 1000`, covariates age + 10 PCs (`covars_no_sex.txt`, **no sex
term**), phenotype file `phenotypes_final.txt`.

| ID | Description | Step 1 | Step 2 flag | N | Like-for-like with G-MoQ? |
|---|---|---|---|---|---|
| **R0** | Archived sex-stratified runs used by every earlier table | none | `--ignore-pred` | 199,572 (F) / 200,138 (M) | N-matched, but **no polygenic adjustment** — handicapped |
| **R0s** | Archived sex GWAS used for `gamma` | none | `--ignore-pred`, `--bt` | **462,869** | **No** — 15.8% more people than G-MoQ's 399,710 |
| **R1** | Rerun of R0 from this project's own invocation | none | `--ignore-pred` | 199,572 / 200,138 | yes; byte-identical to R0 (control) |
| **R1s** | Sex GWAS restricted to the G-MoQ analysis set | none | `--ignore-pred`, `--bt` | **399,710** | yes — N-matched replacement for R0s |
| **R2** | **The fair comparator.** Within-stratum whole-genome ridge | female-only / male-only step 1 | `--pred` | 199,572 / 200,138 | yes |
| **R2s** | Sex GWAS with its own LOCO, N-matched | whole-cohort step 1 on `sex`, restricted to the 399,710 | `--pred`, `--bt` | 399,710 | yes |

### Step 1 provenance

The two stratified step-1 runs copy the whole-cohort invocation logged at
`results/paper/regenie_step1/full/testosterone/testosterone_step1.log` verbatim,
changing only the sample restriction and the output path:

- same `--bed .../merged_chr1-22/chr1-22_merged`
- same `--extract .../testosterone_qc_pass.snplist` (517,143 variants)
- same `--bsize 1000 --lowmem`

**Where sex came from.** Sex is not in `covars_no_sex_paper.txt`. It is column 4
of `data/derived/paper/covars_with_sex_paper.txt` and column 5 of
`data/processed/phenotypes_final.txt`, coded **0 = female (264,687), 1 = male
(223,547)**. No new keep files needed to be derived from it: the pipeline
already contains participant-level sex splits at
`data/processed/samples_{female,male}.txt`, which agree exactly with that coding
(264,687 / 223,547). The analysis-set keep files used here are those files
intersected with the bgen sample list, non-missing testosterone, and complete
covariates. They are participant-level and live only in scratch
(`/user/work/fh6520/gmoq_fair_comparator/`); they are **not** in Git.

### The sample sets are exactly the right ones

Reconstructing the step-2 analysis sets independently reproduces REGENIE's own
reported N to the person, their union is exactly G-MoQ's N, and **every one of
them is present in the array `.fam`**, so adding `--pred` drops zero samples.
The R1-vs-R2 contrast therefore isolates the LOCO adjustment and nothing else.

| | female | male | union |
|---|---|---|---|
| reconstructed analysis set | 199,572 | 200,138 | **399,710** |
| REGENIE's own reported N | 199,572 | 200,138 | — |
| also present in array `.fam` (step 1 eligible) | 199,572 | 200,138 | 399,710 |
| **lost by adding `--pred`** | **0** | **0** | **0** |
| G-MoQ N | — | — | **399,710** |

---

## 3. G-MoQ configurations

| ID | Description | Path (BluePebble) |
|---|---|---|
| **G1** | no-LOCO, pooled + OLS | `gmoq_noloco/grid/<trait>/pooled_ols/chr22.tsv.gz` |
| **G2** | two-pass blind LOCO, pooled + OLS (main phenotype only) | `gmoq_twopass/stage2/pooled_ols_resid_mean/chr22.tsv.gz` |
| **G3** | whole-cohort ("broken") LOCO, pooled + OLS — kept for context | `Gmoq_paper/.local/runs/paper-stage2-testosterone-grid-all-20260916-01/<trait>/chr22/pooled_ols/chr22.tsv.gz` |

Mapping and harmonisation: `beta_left` <-> female, `beta_right` <-> male,
`gamma` <-> sex with **`gamma` negated**, because REGENIE codes Male = 1.
G-MoQ `effect_allele`/`other_allele` are matched against REGENIE
`ALLELE1`/`ALLELE0`; the sign is flipped when swapped and variants matching
neither are dropped. In practice **0 variants** were dropped for allele
mismatch and 17,550 of 17,637 merged.

---

## 4. Control evidence

Two things had to be true before any number below could be trusted.

### 4.1 The `--ignore-pred` reference is reproducible (R1 == R0)

Re-running step 2 from this project's own invocation reproduces the archived
results **byte-for-byte**:

| file | `cmp` vs archive |
|---|---|
| `testosterone_female_chr22_testosterone{,_shrink_0.2,_shrink_0.5,_shrink_0.7}.regenie` | identical (4/4) |
| `testosterone_male_chr22_testosterone{,_shrink_0.2,_shrink_0.5,_shrink_0.7}.regenie` | identical (4/4) |
| `sex_chr22_sex.regenie` | identical |

Restricting the keep file from the full sex split to the exact analysis set also
changes nothing (identical, 8/8), confirming the reconstruction. **The
comparison is controlled.**

### 4.2 One step-1 per stratum serves all four phenotypes

The claim is that the shrink series shifts male values by a constant within the
male stratum, so a within-stratum regression slope cannot change. Verified
directly on the phenotype file rather than assumed:

| phenotype | female `diff` from `testosterone` | male `diff` | within-stratum SD ratio |
|---|---|---|---|
| `shrink_0.2` | exactly 0 | constant −0.4867700593 (sd 4.9e-15) | 1.0000000000 |
| `shrink_0.5` | exactly 0 | constant −1.216925148 (sd 4.6e-15) | 1.0000000000 |
| `shrink_0.7` | exactly 0 | constant −1.703695207 (sd 3.1e-15) | 1.0000000000 |

Female values are bit-identical; male values differ by a single constant to
floating-point rounding. Within-stratum correlation with `testosterone` is 1 to
twelve decimals. A constant shift is absorbed by the intercept, so the ridge
fit, the LOCO prediction, the slope and the SE are all invariant. Two step-1
runs (plus one for `sex`) therefore legitimately serve all four phenotypes —
and the step-2 outputs confirm it, being byte-identical across the four
phenotypes within each stratum in both the `--ignore-pred` and `--pred` arms.

---

## 5. The ledger

<!--LEDGER-TABLE-->

---

## 6. Earlier tables in this repository

A full sweep of the tracked tree for numeric G-MoQ-vs-REGENIE comparisons. These
are **marked, not deleted**. "Handicapped" means the REGENIE side had no LOCO.

### 6.1 Against the handicapped comparator (R0 / R0s)

| Location | Quantity | REGENIE config | Status |
|---|---|---|---|
| `LOCO_MIXTURE_HAZARD.md` "Measured cost" | SE ratios 1.1209/1.0677 (LOCO) vs **1.0236/0.9311** (no-LOCO); cor 0.8981/0.8165 vs 0.9699/0.9188 | R0, states `--ignore-pred` | **Handicapped.** Still valid for its own question (what G-MoQ's broken LOCO cost) since both arms share the reference. Banner added. |
| `EFFICIENCY_BOUND.md` "Decomposing the observed 12%" | SE-factor decomposition ending at observed 1.121 | R0, silent | **Handicapped.** Banner added. Bound and decomposition unaffected; the observed ratio compared against moves. |
| `EFFICIENCY_BOUND.md` "### gamma" | `sqrt(462869/399710)=1.0761` correction, ratio 0.9994 | **R0s, N=462,869** | **Handicapped and not N-matched.** The `1.0761` fixes N by rescaling, but the extra 63,159 people are a *different sample*, not a superset restriction. Superseded by R1s/R2s, which restrict to the same 399,710. |
| `EFFICIENCY_BOUND.md` "Predicted versus observed" (12 cells) | obs/predicted 0.972–1.021 | R0/R0s, silent | **Handicapped** in the "obs" numerator. |
| `EFFICIENCY_BOUND.md` EM cross-check | EM 1.127/1.074 vs G-MoQ 1.121/1.068 | R0, silent | **Handicapped**, but it is an EM-vs-G-MoQ contrast against a shared reference, so the comparison it makes survives. |
| `KNOWN_ISSUES.md` "The LOCO workflow can destroy..." | "8–10% in SE, 0.07–0.10 in correlation" | R0, silent | **Handicapped.** Pointer added. |
| `LOCO_MIXTURE_HAZARD.md` recall caveat | recall 0.41–0.47 on the female component | R0, silent | **Handicapped.** |

### 6.2 A different comparator — not part of this rebuild

These compare **whole-cohort log-BMI / BMI / weight** against REGENIE runs drawn
from a `regenie_step1` + `regenie_step2` pair, i.e. a different question
(does G-MoQ recover REGENIE's additive mean effect?) with a different phenotype
family and no sex stratification. They are listed for completeness and are not
restated against R2.

| Location | Quantity | Notes |
|---|---|---|
| `LOCO_SCALE.md` | cor(mean beta, REGENIE BETA) 0.500 / 0.993 / 0.993 across `--loco-scale` modes; SD recovery to 0.02%; N=461,320 | The section's premise is matching REGENIE's `scale_Y`, which only arises when LOCO is applied, so LOCO appears to be in use on both sides. **Not re-verified against the step-2 log in this pass.** `BACKLOG.md` C6 flags the allele harmonisation as UNVERIFIED. |
| `ALLELE_CONVENTION.md` | 17,550 matched, 0 reversed, cor 0.99285 | log-BMI. Allele convention only. |
| `VALIDATION_2026-09-15.md` line 53 | restatement of cor 0.99285 | as above |
| `CHR22_VALIDATION.md` | acceptance criteria cor ~0.993 / >=0.98 | log-BMI acceptance thresholds |
| `PAPER_FIGURES_AND_TABLES.md` Figure 1 | planned use of r=.99285 | already carries its own "match samples, phenotype, covariates, alleles and LOCO scale" caveat |
| `RELATEDNESS.md` | `A_mu = rho`; lambda_mu−1 = 0.0071 | analytic, not a run; flagged UNVERIFIED at `BACKLOG.md` C4 |

### 6.3 Gaps this rebuild closes

1. **Only one file in the repo ever stated the comparator's `--pred` status**
   (`LOCO_MIXTURE_HAZARD.md`). Every other REGENIE comparison was silent, so the
   handicap did not propagate to readers of the downstream tables.
2. **Two different REGENIE comparators ran under one narrative** — a
   sex-stratified `--ignore-pred` pair at ~200k per arm for
   `beta_left`/`beta_right`, and a whole-cohort sex GWAS at N=462,869 for
   `gamma`. Section 5 puts all three components on N-matched references.
3. **Neither headline table had an experiment-register row.** There was no
   recorded REGENIE command line to check the `--ignore-pred` claim against.
   This run is registered as `stratified-loco-comparator-20260916-01`, and
   section 4.1 supplies the byte-identical reproduction that was missing.
4. **Three different N values for testosterone** appear across the docs —
   399,710, 419,890 (`LOCO_SCALE.md`, the whole-cohort step-1 fit) and 421,405
   (the permutation keep file). Section 2 pins the comparison to 399,710 and
   shows how it is reached.

---

## 7. Reproducing this

Scratch tree (participant-level inputs, never committed):
`/user/work/fh6520/gmoq_fair_comparator/`

```
scripts/step1_{female,male,sex}.sh     # REGENIE step 1, within stratum
scripts/step2_{repro,noloco,loco}_*.sh # REGENIE step 2, chr22
scripts/shift.R                        # constant-shift verification (S4.2)
scripts/compare.R                      # harmonisation + all ledger statistics
scripts/report.R                       # renders the ledger table
analysis/comparison_table.tsv          # full machine-readable results
```
