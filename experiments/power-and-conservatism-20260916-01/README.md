# `power-and-conservatism-20260916-01`

Did the 0.3.0 LOCO change buy discoveries, and is the Stage 2 pooled/OLS
parameter test conservative at genome-wide significance? Full statement of the
questions, the results and the reasoning:
[`docs/POWER_AND_CONSERVATISM.md`](../../docs/POWER_AND_CONSERVATISM.md).

No default was changed. Nothing in `stage1/` or `stage2/` was touched.

## What ran

BluePebble compute nodes only. Exclusive scratch root
`/user/work/fh6520/gmoq_power_conserv`; nothing was written into `G-MoQ`,
`Gmoq_paper` or any existing `gmoq_*` tree.

| File | Slurm | What |
|---|---|---|
| `scripts/gmoq_chr21.sh` | 18983160 | Stage 1 chr21 in both LOCO modes plus matched weights and Stage 2 (skew, pooled, OLS), the chr22 configuration moved one chromosome over. Extends the discovery instrument from one chromosome to two. |
| `scripts/rg21_{female,male,sex}.sh` | 18983137-9 | The fair sex-stratified LOCO'd REGENIE comparator for chr21, same keep files and same step-1 predictions as `gmoq_fair_comparator`'s chr22 arm. |
| `scripts/discovery.R` | interactive | Allele harmonisation, discovery counts, precision/recall, 500 kb locus clumping, status-change counts and the two SE/beta counterfactuals. Run once per chromosome. |
| `scripts/supp.R` | interactive | Genomic control on the real scans, what the G-MoQ-only hits look like in REGENIE, recall predicted from the SE ratio alone, and a block bootstrap on the discovery gain. |
| `scripts/perm_null.R` | interactive | What the existing permuted-phenotype run can and cannot establish: rejection rates with a 100-block, 2,000-replicate block bootstrap, and the pure-scale extrapolation to 5e-8. |
| `scripts/validate2.R` | interactive | Reproduces the shipped Stage 2 pooled/OLS `theta`, `se` and `p` from `(A, Sigma_0, b, d_resid)`. Max relative error 5e-15. This is what licenses the one-line reduction of the test used below. |
| `scripts/dump_psi.py` + `run_dump.sh` | 18983170, 18983236 | Pulls the residualised RIF matrix `Psi_r` and the covariate basis `Q` out of the production Stage 1 code path by hooking `compute_chr_rif_residuals`, so the null simulation uses the real influence functions rather than a model of them. |
| `scripts/nullcal.py` | -- | Exact CGF, saddlepoint solver and Lugannani-Rice tail for `S = sum_i u_i g_i`. |
| `scripts/run_nullcal.py`, `run_nullcal2.py` | interactive | Exact rejection rates over the real chr22 MAF spectrum. The second is the self-consistent version and the one the doc quotes; the first sets the threshold from the real `d_resid` and is kept because the difference between them is itself the genotype-variance diagnostic. |
| `scripts/mc_null.py` + `mc.sbatch` | 18983281 (plain, array 0-49, M = 1e5 each), 18983282 then 18983517 (IS, array 0-11, M = 1e5 then 5e5 each) | Plain Monte Carlo at M = 5.0e6 null variants, and exponentially tilted importance sampling at the 5e-8 and 1e-5 thresholds. |
| `scripts/agg_mc.py`, `attrib.py`, `cost.R`, `check.R` | interactive | Aggregation, the attribution arithmetic, what a hypothetical conservatism would cost on the real scans, and the harmonisation / missed-hit checks. |

Stage 1 chr21 ran at 8 threads, N = 399,710, 17,127 HapMap3 chr21 variants,
45 taus, 1,000 jackknife blocks, seed 42 -- identical settings to the chr22
arm in `loco-as-covariate-20260916-01`.

## Checks that had to pass first

- **Harmonisation.** 34,677 variants merged across chr21 and chr22; all 34,677
  aligned to REGENIE `ALLELE1` with 0 flips, so no sign flip was applied and
  sign concordance at every REGENIE hit is 1.000.
- **Stage 2 reduction.** `validate2.R` reproduces the shipped chr22
  covariate-mode output to 5e-15 relative on `theta`, `se` and `p`.
- **`Psi_r` fidelity.** `Sigma_0 = Psi_r' Psi_r / N` rebuilt from the dumped
  matrix matches the shipped `sigma0.tsv.gz` to 4e-16.
- **Simulation calibration.** Plain MC gives `Var(t)` = 1.00016 / 1.00016 /
  1.00000 and `kurt(t)` = 2.998 / 3.000 / 3.001 at M = 5.0e6.
- **Two independent tail estimators agree.** Tilted importance sampling and
  the saddlepoint agree to within 0.5% in all twelve (param, MAF, alpha) cells
  at M = 5e5 tilted draws each.

## Results

Part 1, power:

- `discovery_chr22.tsv`, `discovery_chr21.tsv`, `discovery_pooled.tsv`
  -- hits, overlap, precision, recall, 500 kb locus counts, sign concordance.
- `statuschange_chr22.tsv`, `statuschange_chr21.tsv` -- gained / lost /
  McNemar, plus the `SE-only` and `beta-only` counterfactuals.
- `seratio_chr22.tsv`, `seratio_chr21.tsv` -- SE ratios and how far the
  point estimates moved, in SE units.
- `gain_bootstrap.tsv` -- 95% block-bootstrap interval on the gain.
- `gmoq_only_hits.tsv` -- where the G-MoQ-only hits sit in REGENIE.
- `inflation_real.tsv`, `recall_pred.tsv`.

Part 2, conservatism:

- `perm_null_rates.tsv` -- the permuted null with block-bootstrap errors
  and the implied number of independent tests.
- `spa_scan_selfconsistent.tsv`, `spa_bins_selfconsistent.tsv` -- exact
  saddlepoint rejection rates, MAF-weighted and per MAF bin.
- `mc_plain.tsv`, `mc_is.tsv` -- Monte Carlo and importance sampling.
- `attribution.json`, `cost_of_conservatism.tsv`.
- `A_contrasts.tsv` -- the three OLS contrast vectors over taus, the
  only object from the weights that Part 2 needs.

All committed output is aggregate over SNPs or over simulation replicates. No
participant-level data left BluePebble.

## Load-bearing caveats

- The discovery instrument is two chromosomes of one trait, and at 5e-8 it is
  **five independent loci**. The paired status-change evidence (35 gains and
  0 losses at 5e-8, 135 and 5 at 1e-5) is much stronger than the count
  intervals, which mostly include zero.
- `gamma` has no signal on either chromosome, so the LOCO change's effect on
  the membership parameter's power is untested.
- Part 2's null holds the realised phenotype fixed and randomises the genotype
  -- the frame a permutation test uses. Genotypes are iid HWE: LD is absent by
  construction, which is deliberate, because LD is exactly what makes the real
  permuted null unable to answer the question.
- The residual departure from nominal at 5e-8 is in the **anti**-conservative
  direction and is largest at the rarest MAF bin (MAC 7,805): +1.7% to +2.6%
  there, +0.02% to +0.13% averaged over the real MAF spectrum. Nothing in the
  conservative direction exceeds 0.3%.
