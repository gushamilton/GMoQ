# Simulation inventory — what each campaign asked, what it concluded, and whether that is still true

Written 2026-09-17 against `main` at `9c47e7b`, version 0.4.0.

This is an audit of every simulation campaign accumulated in this repository and
in `gushamilton/Gmoq_paper`. Its purpose is to say, campaign by campaign, whether
the conclusion survives the 0.4.0 defaults, and which campaigns a live document
still leans on. It is also the record of what was archived and where it went.

**Frozen archive ref:** tag `archive/sim-campaigns-20260917` and branch
`archive/sim/campaigns-20260917`, both at `9c47e7b`. They hold every campaign
below in full, including the ones removed from `main` by this change. Nothing was
deleted. To read an archived campaign:

    git show archive/sim-campaigns-20260917:sim/density/tables/02_coverage_typeI.txt
    git checkout archive/sim-campaigns-20260917 -- sim/out          # into a scratch worktree

---

## The defaults these campaigns are being judged against

| | 0.4.0 default | was |
|---|---|---|
| Stage 2 `cov_source` | `pooled` (`stage2/R/gwas.R:851`) | unchanged |
| Stage 2 `tau_cov_estimator` | `ols` (`stage2/R/gwas.R:253`) | unchanged |
| Stage 1 `--loco-mode` | `covariate` when `--loco-preds` given, else `subtract` | `subtract` always |
| Stage 1 kernel | restructured "fast" kernel; `--exact-kernel` is the old path | exact only |
| Stage 1 `--blocks` | 1000, 45 taus in [0.02, 0.98] | unchanged |
| Stage 1 `--cov-mode` | `{full, pooled, pooled-only}` | did not exist |

Two of those movements turn out not to reach any campaign below, and it is worth
saying why up front:

* **The LOCO flip is inert here.** Every campaign in this repository simulates
  genotypes and phenotypes directly and calls `core.run_stage1()` with no LOCO
  predictor at all (`docs/BACKLOG.md` item A1 makes the same observation from the
  other direction, as a *testing* gap). `sim/projeff` is the only one that reads
  anything from a LOCO-era run, and what it reads is a fitted weight matrix
  (`gmoq_noloco/.../weights_skew.rds`) whose column geometry is the object of
  study. No conclusion below was measured under `--loco-mode subtract` in a way
  that the flip to `covariate` disturbs.
* **The kernel restructure is inert here too.** `docs/STAGE1_FAST_KERNEL.md`
  measures the fast kernel against the exact kernel at reassociation level:
  `max abs(delta beta_tau) / SE = 5.35e-14`, `d_resid` max relative `1.01e-14`,
  identical accept/reject on all 176 batches, zero discovery-status changes on
  real chr22. Every coverage and Type-I number below is quoted to four decimals
  at best. The kernel cannot move any of them.

So the live threats to these conclusions are narrower than the version bump
suggests: the *W* a campaign was run against, and the `cov_source` /
`tau_cov_estimator` arm it assumed.

---

## Inventory

### 1. `sim/covsource` — the 780,000-scan covariance-source head-to-head

**Question.** Which covariance source should Stage 2 default to, measured on the
comparison that had never been run — pooled+OLS against jackknife+OLS?

**Conclusion.** `docs/FINDINGS_covsource_default.md`. 52 cells x 15,000 paired
replicates at N = 461,320, `--blocks 1000`. Under OLS the point estimate never
touches `Sigma`, so `sd(estimate)` is 1.00000 between arms in all 44 main cells
and the two arms differ only in the SE they report. Pooled decays monotonically
with rarity under a **true dispersion** effect: `beta_sigma2` coverage
0.9467 / 0.9285 / 0.9027 / 0.8475 / 0.8742 / 0.7919 across MAF 0.30–0.005, against
jackknife 0.944–0.951 throughout. Extrapolated to 5e-8, `beta_mu` under a true
dispersion effect is inflated 1.26x / 295x / 2086x (pooled) against 0.94x / 1.08x /
0.99x (jackknife). At the null the two are indistinguishable, 0.047–0.054 — which
is why null simulations never caught it.

**Still valid under 0.4.0?** **VALID-BUT-NARROWED.** Every arm it measured is
still a shipped arm; the pooled degrees-of-freedom convention it used
(`nu = n_samples`, `nu_source = "pooled_sigma0_n_samples"`) is still what
`stage2/R/gwas.R:1011` sets. What narrows it is the *W*: this campaign ran the
**2-column** `make_weights_vqtl` W (`beta_mu`, `beta_sigma2`). The paper and the
package ship a **3-column** mixture W with a membership column. Campaign 2 exists
precisely because that difference was not assumable in either direction, and it
found the defect **concentrates** rather than disappears. Read this campaign as
the mechanism and campaign 2 as the shipped-W magnitude.

Note also that its own recommendation — jackknife+OLS — was **overruled**.
`BACKLOG.md` B1 is CLOSED as "pooled, with a documented boundary". The campaign is
not wrong; the default is a judgement made on top of it.

**Load-bearing?** Yes, twice over. `docs/FINDINGS_covsource_default.md:6` names
`sim/covsource/tables/` as its aggregates, and `Gmoq_paper` consumes the *same
bytes* — `outputs/simulations/biobank_covariance/full_audit/{cells.json,
exact_table.tsv, paired_cs.tsv, power_cs.tsv, summary_cs.tsv, tables.md}` are
byte-identical to `sim/covsource/tables/`, and drive paper Figure 3
(`scripts/figures/render_paper_figures.py:129`) and the supplementary tables
(`scripts/tables/build_supplementary_tables.py:50`).

**Promotable?** **Already promoted.** The aggregates are in the paper repo and a
referee is already looking at them. Nothing further to do.

**Disposition: ARCHIVED** to `archive/sim-campaigns-20260917`, and
`docs/FINDINGS_covsource_default.md` updated to name that ref and the paper-repo
copy. The aggregates a referee needs continue to live in `Gmoq_paper`, so removing
the duplicate from `main` loses no evidence.

---

### 2. `sim/covsource_mixtureW` — does the defect transfer to the W we actually ship?

**Question.** `sim/covsource` decided `cov_source` on a 2-column means-only W. Does
the pooled-covariance defect look the same on the 3-column mixture/skew W
(`gamma`, `beta_left`, `beta_right`) that production builds?

**Conclusion.** No — it **concentrates**. Exact asymptotic surface over 158 cells
plus 19 Monte Carlo cells at 5,000 reps. At MAF 0.005 / Z = 30 with a true
dispersion effect, `gamma` is *protected* (pooled SE ratio **1.0136**, 5e-8
inflation 0.658x) because it takes essentially all its variance on the
tau-covariance's leading eigendirection, while `beta_left`/`beta_right` sit at 90
degrees to it with zero share and are inflated **3.51e+03x**; at MAF 0.01 / Z = 30
the same pair is inflated **236x**. The onset of a 5% SE shortfall on the component
means is |Z| = **6.9** at MAF 0.005, **9.8** at 0.01, **23.2** at 0.05 and **87.9**
at 0.30. Under the null the two sources agree exactly, provably — the discrepancy
term `Sigma_true - Sigma_pooled = diag(1/f) E diag(1/f)` vanishes. At the real
chr22-testosterone operating point (MAF 0.05, |Z| ≤ 11) the pooled SE ratio on
`beta_1` is 0.978 and on `gamma` 1.0009.

Seven controls pass before any number is drawn, including the one that licenses
the whole exercise: `make_weights_normal_mixture(include_membership = TRUE)`
equals production's `build_skew_mixture()` W at `alpha = 0` to max rel diff
3.9e-10.

**Still valid under 0.4.0?** **VALID.** This is the campaign run against the W and
the arms the package ships today. It is the most current piece of simulation
evidence in the tree.

**Load-bearing?** **Yes — and more exposed than anything else here.**
`README.md:58` states the 3510x / 236x figures, the `gamma` protection and the
whole onset ladder as the package's headline caveat on its own default, and cites
**`sim/covsource_mixtureW/` the directory**. `BACKLOG.md` B1 — the CLOSED decision
to keep `pooled` — quotes the same numbers.

**There is no `docs/FINDINGS_covsource_mixtureW.md`.** The campaign's own README
promises one and it was never written. So the *only* readable statement of this
result outside a paragraph of `README.md` is
`sim/covsource_mixtureW/tables/tables_exact.txt`. Archiving this directory would
leave the package's central caveat about its own default with no evidence on
`main` at all.

**Promotable?** **Yes — the strongest candidate on this list, and the one a
referee will actually ask for.** It is the direct answer to "you ship pooled; show
me what pooled costs on the W you ship", it is exact rather than Monte Carlo over
most of its surface, and it carries its own verification suite. It needs the
findings document it was promised.

**Disposition: KEPT ON `main`, unarchived.** See "Open gap" below.

---

### 3. `sim/projeff` — why is the OLS-vs-GLS penalty 57% in simulation and 0.8% on real data?

**Question.** BACKLOG H5. The projection efficiency loss from using OLS instead of
GLS is enormous in simulation and negligible on real testosterone. Which is wrong?

**Conclusion.** `docs/EFFICIENCY_BOUND.md` (register
`projection-efficiency-20260916-01`). Neither — it is the **membership column**,
not the distribution family. The simulation used a 2-column means-only W;
production's `build_skew_mixture()` always emits `gamma`. At matched separation and
K, skew and normal agree (S = 4: 1.593 vs 1.568/1.732); deleting `gamma` from the
real testosterone W moves the penalty 1.013 -> 2.13/2.51. `W_gamma` sits
0.011–0.16 degrees from `Sigma`'s leading eigenvector, so `col(W)` is
`Sigma`-invariant and OLS is already efficient. Diagnostic `rho_max` has Spearman
0.993 against the realised penalty across 38 configurations and the bound is never
violated. Every shipped W sits at `rho_max` 0.03–0.29 (penalty < 4%); the
means-only W reaches 0.67–0.99. GLS under a 16% density error in `Sigma` reports
0.40–0.82 of its true SE while OLS stays ≥ 1.

**Still valid under 0.4.0?** **VALID.** It explains the shipped default (`ols`) and
was measured on the shipped W family. Nothing moved underneath it.

**Load-bearing?** Yes. `docs/EFFICIENCY_BOUND.md:286` says outright: "Full
evidence, the mixture-family prediction and the `Sigma`-misspecification fragility
result are in `sim/projeff/src`, with aggregate tables in `sim/projeff/tables`."
`BACKLOG.md` H5 is marked RESOLVED against it.

**Promotable?** The `rho_max` diagnostic and the `Sigma`-misspecification result
are worth a supplementary paragraph — they are the honest answer to "why not GLS?",
which a referee will ask. The 113-row exact grid itself is development detail.

**Disposition: KEPT ON `main`.** 12 files, 88K; cited by path from a live document;
archiving it would buy nothing and cost a working link.

---

### 4. `sim/nscale` — does calibration measured at N = 30,000 survive to N = 461,320?

**Question.** `OLS_VS_GLS.md` established coverage at N = 30,000 only. Production N
is 461,320. Does the "sag to ~0.94 under a true scale effect" hold, shrink, or
grow?

**Conclusion.** `docs/N_SCALING.md`. 36 cells x 20,000 replicates at N = 30,000 /
100,000 / 461,320 (720,000 replicates, ~470 CPU-hours). Under the null, calibration
is exact at every tau, MAF and MAC tested. Under a true scale effect it is **not**
a constant of the method: it is misspecification bias measured in standard errors
that shrink. The absolute misfit `max|beta_tau - W theta|` is 0.12162 at all three
rungs to five decimals; lack of fit *in SE units* grows like `sqrt(N_eff)`,
measured 10.92 / 19.93 / 42.81 against predicted ratios 1.83 and 3.92 — exact.
Coverage of `beta_sigma2` under a true scale effect at MAF 0.30 follows it down:
OLS 0.9440 / 0.9434 / 0.9413; **GLS 0.9316 / 0.9064 / 0.8437**.

**Still valid under 0.4.0?** **VALID for the default; the alarming half is on the
opt-in path.** The column that collapses with N is GLS, which is not the default.
The shipped OLS path is flat at 0.941–0.944 across a 15x range of N. This is a
general property of projection onto a misspecified `W` and is not sensitive to any
0.4.0 change.

**Load-bearing?** Yes, mildly. `docs/N_SCALING.md:364` names `sim/nscale/src/` as
the reproduction recipe, script by script. It is a live document (docs index,
"How does this behave as N grows?").

**Promotable?** The `sqrt(N_eff)` scaling law with its exact predicted-vs-measured
ratios is a clean result and belongs in a supplement if the paper claims
large-N validity. The per-cell grid is development history.

**Disposition: KEPT ON `main`.** 16 files, 88K, scripts and `cells.json` only — no
bulk. It is the reproduction section of a live document.

---

### 5. `sim/density` — does G-MoQ hold up on phenotypes that are not smooth?

**Question.** `BP1_density_estimate` showed the Silverman bandwidth is 3.2–3.4x too
wide for bimodal traits and distorts `f_hat` by up to 32%. Nobody had translated
that into coverage or Type-I error. Does it reach inference?

**Conclusion.** `docs/BP2_lumpy_phenotypes.md`. 16 cells across five families of
non-smooth phenotype (ties/rounding, floor and ceiling censoring, bimodality,
count-like, heavy skew), simulated genotypes at MAF 0.05/0.20/0.50, no real data.
**The failure is not where the density doc pointed.** The KDE error is real and
large — up to 2.0x at a bimodal median trough, up to 38x below a censoring limit —
but `b` and `se` each scale as `1/f_hat`, so it **cancels** in every single-quantile
statistic: `z` is invariant to 5e-14 under a 0.5x/2x bandwidth perturbation. Stage 1
Type-I at 0.05 over 225 cell-by-tau combinations has median 0.0501 and max 0.0575,
with no tau anti-conservative anywhere. What breaks is Stage 2 **GLS**, for an
unrelated reason: a tie block spanning two requested taus gives them the same
empirical quantile, their RIF columns coincide, `Sigma` goes singular (condition
number 1e2 -> 1e13), ridge repair fires on 100% of SNPs and inflates SEs by
1e2–1e3x. OLS held 0.9491–0.9508 across every shape tested, including a four-valued
count trait. The boundary is not a tie percentage: a tau is safe iff its tie
fraction is smaller than the gap to its nearest neighbouring tau.

**Still valid under 0.4.0?** **VALID-BUT-NARROWED, and it was already narrowed at
birth.** The commit that brought it in (`56f1059`) records the qualification: the
study's own summary calls GLS the default and so overstates the severity, because
it was written against a branch predating the roxygen correction in `24eef02`. The
broken path is opt-in and has been throughout 0.4.0. The OLS result — which is the
shipped path — stands as measured.

Caveat that is *not* resolved: `BACKLOG.md` C3 lists two of this document's claims
as **UNVERIFIED** — "KDE error cancels exactly in z" (5e-14) and "ceiling censoring
270x at 1e-6", whose extrapolation to 5e-8 is unmeasured and carries ~14% Monte
Carlo error. Anyone closing C3 will need this campaign back out of the archive.

**Load-bearing?** Indirectly. `docs/BP2_lumpy_phenotypes.md` is a live document
(docs index, "My phenotype is lumpy / has a point mass") but it cites its code
only relatively — `src/shapes.py`, `src/keepset.py` — and never names
`sim/density`. It reproduces its own tables inline. `tests/lib/gmoq_sim.py:57`
mentions `sim/density/src` in a comment about how Stage 1 is driven.

**Promotable?** No. This is a robustness investigation whose useful output is the
*rule* ("a tau is safe iff its tie fraction is smaller than the gap to its nearest
neighbour") and the tie-drop guard, both already in the document. The 52 per-cell
JSON summaries are development history.

**Disposition: ARCHIVED** to `archive/sim-campaigns-20260917`. 2.4M and 82 files,
the largest single item on `main`. `docs/BP2_lumpy_phenotypes.md` gains a
provenance line naming the ref — which it previously lacked entirely.

---

### 6. `sim/out` + `sim/src` — the pooled-covariance dispersion study

**Question.** Is the pooled covariance safe as G-MoQ's default, and what actually
produces false dispersion findings?

**Conclusion.** `docs/FINDINGS_pooled_dispersion.md`, plus the exact engine behind
`docs/DISPERSION_SKEW_CORRECTION.md`. "Yes for the scan, no for the loci the scan
finds — but the reason is not the one the estimator has been charged with." The
pooled working-independence error is real, always anti-conservative, bounded and
small at common MAF — and is **not** what produces false vQTLs. Those come from a
second, separate mechanism: a pure location effect at a genotype with non-zero
skewness leaves a second-order residue
`-(delta^2/2) * gskew * f0'(q) / (f0(q) sqrt(Var g))` which, for a normal baseline,
is exactly proportional to the dispersion weight column, and is absorbed by
`beta_sigma2` in full. It scales as `Z_mu^2`, **grows** as `sqrt(N)`, and vanishes
identically at MAF 0.5 where genotype skewness is zero — the MAF 0.50 / Z_mu 60
cell showing no inflation is the proof of mechanism. Exact false-dispersion rates
at 5e-8, N = 461,320, zero true dispersion: 10.6x nominal at MAF 0.30 / Z_mu 30;
6.5e4x at MAF 0.05 / Z_mu 30; **0.912 at MAF 0.01 / Z_mu 30** — 91% of such variants
falsely call a dispersion effect. No covariance estimator touches it, which is why
every SE-based check missed it.

**Still valid under 0.4.0?** **SUPERSEDED IN PART — and the document says so
itself, in two separate ways.**

* Section 8 item 1, "keep pooled as the scan default", is **withdrawn** by
  `FINDINGS_covsource_default.md` §9.
* Its claim that 5e-8 is "exactly nominal for the jackknife" is estimator-specific
  and does not say so: true for jackknife+**OLS**, false for jackknife+**GLS**
  (3.45–4.95x at MAF 0.30, 457x at MAF 0.005/Z30, and inflated at the null before
  any effect exists).
* Section 5 proposes *flagging* dispersion findings at `G >= 0.2`. That artefact
  has since been **corrected** rather than flagged, by `correct_dispersion_skew()`.
  After correction the 5e-8 rate is within 1.03x nominal everywhere in this
  document's own sweep. The flag's meaning has inverted into a power warning.
* Its lead cell, MAF 0.005 / Z 60, requires the minor homozygote to carry about
  **12x the phenotypic SD** and is not a real variant.

`docs/README.md` already files this document under "Audit trail — point-in-time
records, **not** current guidance". What survives untouched is the diagnosis of the
false-vQTL artefact as an estimand gap, and that survives in the two documents that
took it over.

**Load-bearing?** No live claim rests on the directory. `sim/README.md` describes
it; `FINDINGS_pooled_dispersion.md` is explicitly retired guidance. One open item
touches it: `BACKLOG.md` C1 asks of this document "which estimator produced
5.00e-08?" and is UNVERIFIED — answerable from the archive ref.

**Promotable?** Only the exact false-dispersion engine (`sim/src/fdr_exact.py`),
and only because `DISPERSION_SKEW_CORRECTION.md` reports that it reproduces Table 1
"to every published digit". The Monte Carlo was salvaged from a `../poolcal` run
whose driving agent died mid-campaign and was never re-run — provenance a referee
would not enjoy. Development history.

**Disposition: ARCHIVED** to `archive/sim-campaigns-20260917`. This is the clearest
case on the list: the study behind an explicitly retired document, whose central
recommendation has been withdrawn and whose proposed remedy has been replaced by a
shipped correction.

---

### 7. `sim_meanvar` — how a pure mean effect manufactures a vQTL

**Question.** On a raw, untransformed, right-skewed trait, how large a pure mean
effect is needed before the dispersion parameter goes past 5e-8 with no true
dispersion effect at all?

**Conclusion.** `docs/MEAN_VARIANCE_ARTEFACT.md`. Exact population calculation,
cross-validated against an independent implementation. **|Z_mu| = 9.5 — about 0.02
SD per allele — is enough.** That is an ordinary GWAS hit, not an outlier. Log or
rank-inverse-normal transformation removes this mechanism *exactly* and pushes the
threshold out to |Z_mu| = 43–93, at which point a second, unrelated mechanism takes
over that no transformation touches. It also corrects `RELATEDNESS.md`'s numbers
for the same configuration (28.1 SD -> 33.2 raw; 2.7 -> 3.12 log; 3.2 -> 3.09 RINT)
and shows the raw-scale case carries a **real** dispersion effect of 0.2341, not a
spurious one.

**Still valid under 0.4.0?** **VALID.** It is an exact population calculation about
the estimand, not about an estimator arm. No default touches it. It is also the
origin of `stage2/R/dispersion-guards.R`, which ships.

**Load-bearing?** Yes, by path, three times: `docs/MEAN_VARIANCE_ARTEFACT.md`
lines 4, 361 and 455 cite `sim_meanvar/src/`. It is a live document (docs index,
"I found a variance effect — is it real?").

**Promotable?** The headline — |Z_mu| = 9.5 on a raw skewed trait, and that a
transformation removes it exactly — is a result the paper should state, because it
is a warning to anyone running a vQTL scan and not a detail of this implementation.
The working code is development history.

**One thing found while auditing that is worth fixing separately.**
`sim_meanvar/src/dispersion-guards.R` is a **340-line snapshot of what is now the
1,317-line shipped `stage2/R/dispersion-guards.R`**, sitting alongside its own
`test-dispersion-guards.R`. It has drifted by roughly a thousand lines and nothing
marks it as a copy. That is a live trap for anyone reading the tree for the guard's
behaviour, and archiving it removes the trap.

**Disposition: ARCHIVED** to `archive/sim-campaigns-20260917`, with all three path
citations in `docs/MEAN_VARIANCE_ARTEFACT.md` updated to name the ref.

---

### 8. `sim_skewcorr` — validating the genotype-skewness correction

**Question.** The false vQTL from campaign 6 is bias, not a variance error, so it
can be subtracted rather than only flagged. Does subtracting it work, and what does
it cost?

**Conclusion.** `docs/DISPERSION_SKEW_CORRECTION.md`. On a Normal baseline the
correction is exact to O(delta^4) and removes the artefact completely: the worst
false-positive rate in the sweep falls from **0.91 to 1.5e-9** at genome-wide
significance, at no cost on null variants and none on genuine dispersion effects —
except at a rare variant carrying a colossal location effect, where it costs about
26% of honest power. Three conditions are enforced in code, all of which make the
correction actively harmful rather than merely useless if broken: a *signed* effect
allele frequency, `make_weights_vqtl`'s *numeric* design rather than merely its
column names, and a baseline within the skewness envelope for the chosen estimator
(about 0.65 under GLS, about 1.0 under OLS).

**Still valid under 0.4.0?** **VALID.** `correct_dispersion_skew()` ships, its three
guards ship, and the campaign measured the estimand gap rather than an estimator
arm. Its envelope is stated per-estimator, so the OLS default is the wider of the
two.

**Load-bearing?** Weakly. `docs/DISPERSION_SKEW_CORRECTION.md` — a live document —
does **not** cite the directory by path at all; it is self-contained. The single
path citation anywhere is `docs/EMPIRICAL_W.md:50`, which points at
`sim_skewcorr/src/exact.py` for one algebraic identity.

**Promotable?** The 0.91 -> 1.5e-9 result and the 26% power cost belong in the
paper wherever `correct_dispersion_skew()` is described — that is the correction's
entire justification. The `a1`/`a2`/`a3` sweep tables are development history.

**Disposition: ARCHIVED** to `archive/sim-campaigns-20260917`, with
`docs/EMPIRICAL_W.md:50` updated to name the ref.

---

### 9. `tests/` (tier0–tier4) — **not a campaign, not audited, not touched**

The tiered suite with `criteria.json` is the package's test suite, not a paper
campaign. It is governed by `docs/TEST_SUITE.md`. Left exactly as it was.
`tests/lib/gmoq_sim.py` carries two comments referring to archived campaigns as
examples of how Stage 1 is driven; those comments are updated to name the archive
ref, and nothing executable changes.

Worth restating from `BACKLOG.md` A1, because this audit is otherwise a list of
things that *were* measured: `extract_dosage`, `match_samples`,
`load_analysis_inputs`, `read_keep_ids` and `write_stats_batch` have **zero test
references**, and all six simulation campaigns call `core.run_stage1()` directly,
bypassing the input pipeline entirely. No campaign in this inventory is evidence
about reading a BGEN. Tier 4 is built and unrun.

---

### 10. `experiments/` — the run register. Not a campaign, not touched.

`experiments/README.md` is a dated, immutable ledger of individual runs, and is the
thing that makes the rest of this auditable. Three of its entries are the
provenance records for campaigns above (`projection-efficiency-20260916-01`,
`biobank-covariance-20260915-01`, `loco-as-covariate-20260916-01`). Left alone.

---

## Paper repo (`gushamilton/Gmoq_paper`) — all live, nothing archived

Checked and deliberately untouched. Every file below is consumed by a script that
builds a paper figure or table, so none of it is dead:

| path | consumed by |
|---|---|
| `outputs/simulations/biobank_covariance/summary_pooled_ols.tsv`, `summary_jk_compare.tsv` | Figure 2 (`render_paper_figures.py:86-87`), supplementary tables (`build_supplementary_tables.py:48-49`) |
| `outputs/simulations/biobank_covariance/full_audit/` | Figure 3 power grid (`render_paper_figures.py:129`), supplementary design (`build_supplementary_tables.py:50-51`), `HANDOFF.md:46` |
| `outputs/simulations/validity_pilot/` | Figure S3 (`render_paper_figures.py:193-195`) |
| `outputs/stage2-benchmark-20260915-05/timing.tsv` | Figure 4 (`render_paper_figures.py:167`) |

**Data check.** No participant-level content. The largest tracked file in
`outputs/simulations/` is 871 lines (`full_audit/tables.md`); every table is keyed
by simulation cell, MAF, shape and parameter, and `validity_pilot/cell_provenance.csv`
is one row per cell. `outputs/simulations/README.md` states the exclusion explicitly
and the files match the statement. Nothing was moved or copied.

Also noted and deliberately undisturbed: the remote branches
`archive/work/pooled-validity` and `archive/evidence/pilot-f1b4a2c`, which carry a
`tests/validity_pilot/` tree.

---

## Summary table

| campaign | question | conclusion | still valid at 0.4.0? | cited by | promotable? | disposition |
|---|---|---|---|---|---|---|
| `sim/covsource` | pooled vs jackknife, which default? | jackknife+OLS nominal; pooled 2086x inflated at 5e-8 under true dispersion at MAF 0.005; estimates identical | VALID-BUT-NARROWED — measured on the 2-column W, not the shipped 3-column one | `FINDINGS_covsource_default.md:6`; **paper Fig 3 + supp tables** | already promoted (byte-identical in `Gmoq_paper`) | **ARCHIVED**, doc updated |
| `sim/covsource_mixtureW` | does it transfer to the shipped W? | it concentrates: `gamma` protected (SE ratio 1.0136), `beta_left/right` inflated 3510x / 236x; onset \|Z\| 6.9/9.8/23.2/87.9 | **VALID** — the shipped W, the shipped arms | **`README.md:58`**, `BACKLOG.md` B1 | **yes — strongest candidate** | **KEPT** |
| `sim/projeff` | why is the OLS/GLS penalty 57% vs 0.8%? | the membership column; `col(W)` is `Sigma`-invariant; `rho_max` Spearman 0.993 | VALID | `EFFICIENCY_BOUND.md:286`, `BACKLOG.md` H5 | `rho_max` + GLS fragility, as a supplement | **KEPT** |
| `sim/nscale` | does N=30k calibration survive to N=461k? | null exact at all N; under a true scale effect misfit is constant and grows as `sqrt(N_eff)`; GLS coverage 0.9316 -> 0.8437, OLS flat 0.944 -> 0.941 | VALID for OLS; the collapsing arm is opt-in GLS | `N_SCALING.md:364` | the `sqrt(N_eff)` law | **KEPT** |
| `sim/density` | does the KDE bandwidth error reach inference? | no — it cancels in `z` to 5e-14. GLS breaks instead, on tie-induced singular `Sigma` | VALID-BUT-NARROWED — the broken path is opt-in; C3 claims unverified | `BP2_lumpy_phenotypes.md` (relative paths only) | no | **ARCHIVED**, doc gains provenance |
| `sim/out` + `sim/src` | is pooled safe as the default? | false vQTLs are an estimand gap, not a covariance error; 91% false-call rate at MAF 0.01 / Z_mu 30 | **SUPERSEDED IN PART** — §8.1 withdrawn, §5's flag replaced by `correct_dispersion_skew()`, lead cell unphysical | `sim/README.md`; parent doc is retired guidance | no | **ARCHIVED** |
| `sim_meanvar` | how small a mean effect fakes a vQTL? | \|Z_mu\| = 9.5 on a raw skewed trait; log/RINT removes it exactly | VALID — exact, about the estimand | `MEAN_VARIANCE_ARTEFACT.md` x3 | the \|Z_mu\| = 9.5 headline | **ARCHIVED**, doc updated |
| `sim_skewcorr` | does subtracting the skew bias work? | worst false-positive rate 0.91 -> 1.5e-9; ~26% power cost in one corner | VALID | `EMPIRICAL_W.md:50` only | the 0.91 -> 1.5e-9 result | **ARCHIVED**, doc updated |
| `tests/` tier0–4 | (test suite) | — | — | `TEST_SUITE.md` | — | **untouched** |
| `experiments/` | (run register) | — | — | everywhere | — | **untouched** |
| `Gmoq_paper/outputs/` | (paper evidence) | — | — | figure and table scripts | — | **untouched** |

---

## Open gap this audit found

**`docs/FINDINGS_covsource_mixtureW.md` does not exist.** `sim/covsource_mixtureW/README.md`
names it as the campaign's findings document and it was never written.

This matters more than a missing file usually would. `README.md:58` — the
package's own statement of the boundary on its own shipped default — quotes that
campaign's numbers directly (3510x, 236x, the `gamma` SE ratio 1.0136, the onset
ladder), and `BACKLOG.md` B1 closes the `cov_source` decision on them. Every one
of those figures traces to `sim/covsource_mixtureW/tables/tables_exact.txt`, a raw
output dump with no prose around it.

That is why this campaign was kept on `main` rather than archived: archiving it
would have left the package asserting a numerical caveat about its own default with
no evidence and no document behind it. Writing that findings document is the
highest-value follow-up from this audit.

---

## Recovering an archived campaign

    git fetch origin 'refs/tags/archive/sim-campaigns-20260917:refs/tags/archive/sim-campaigns-20260917'
    git show archive/sim-campaigns-20260917:sim/out/tables.txt
    git worktree add /tmp/simarchive archive/sim-campaigns-20260917   # whole tree, read-only use
