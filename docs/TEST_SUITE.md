# The G-MoQ validation suite

Historical suite design retained from the parallel validation branch. The
current reporter labels its outcomes as regression checks, not a certificate
of universal statistical validity. See [the current validation report](VALIDATION_2026-09-15.md)
and [experiment register](../experiments/README.md) for verified results.
All BluePebble execution belongs on compute nodes, including interpreter and
environment checks.

What it establishes, what it does not, what it costs, and how to run one tier
without the others.

The suite lives in `tests/`. Every tier has a single entry point and leaves a
machine-readable `results.jsonl`; `tests/report.R` turns those into the one
artefact a human reads to decide whether a release is sound.

```
tests/
  run_tier0.sh     unit          seconds      no cluster, no Stage 1
  run_tier1.sh     integration   ~1 minute    no cluster; needs the Stage 1 conda env
  run_tier2.sh     validity      cluster      preregistered grid, counted rates
  run_tier3.sh     scale/timing  cluster      production N on a synthetic BGEN
  report.R         the release table
  criteria.json    the preregistered pass criteria for tier 2
  lib/             harness, model catalogue, simulator, environment
  tier0/ tier1/ tier2/ tier3/
```

---

## Running it

```bash
export GMOQ_WORK=/some/scratch            # never inside the repo
export GMOQ_DEP_LIB=/path/to/cran/library # where data.table, Rcpp etc. live
export GMOQ_CONDA_SH=$HOME/miniforge3/etc/profile.d/conda.sh   # tiers 1-3

bash tests/run_tier0.sh                   # installs stage2/ and runs the unit tier
bash tests/run_tier1.sh
bash tests/run_tier2.sh submit            # Slurm array
bash tests/run_tier2.sh status
bash tests/run_tier2.sh collate
bash tests/run_tier2.sh quick             # 3 cells, no Slurm, a few minutes
bash tests/run_tier3.sh 461320 5000
Rscript --vanilla tests/report.R $GMOQ_WORK/results
```

Each tier is independent: tier 2 does not require tier 1 to have run, and the
report says plainly which tiers did not run rather than treating them as
passes.

**Two environment traps the suite handles for you.** `~/.Rprofile` on the
analysis cluster contains a bare `.libPaths(...)`, which *replaces* the search
path and makes `R_LIBS_USER` silently ineffective. Every script therefore runs
under `Rscript --vanilla`, sets the path from `GMOQ_TEST_LIB`, and then asserts
that `fungwasStage2` really loaded from there. It also asserts that the package's
shared object is loaded, because a `NAMESPACE` that has lost `useDynLib`
installs cleanly and fails only at the first C++ call.

---

## Tier 0 -- unit

**Entry point:** `bash tests/run_tier0.sh`. **Runtime: 67 s** including the
package build. No cluster, no Stage 1, no network, no simulation tolerance that
was tuned to make a test pass.

**What it establishes.** That the arithmetic is right, against closed forms:

| group | what is checked against what |
|---|---|
| `t0_rif.R` | the RIF against its analytic influence function; the two-valued step; exact location-scale equivariance (`RIF(aY+b) = a RIF(Y) + b`), which is *why* `b/se` is free of the density estimate; the exact sample identities for `Var[RIF]` and the cross-tau covariance; convergence on the Brownian-bridge form |
| `t0_weights.R` | every weight builder against a **central difference of the actual quantile function** -- not against a re-typed copy of the algebra, so a sign error in the published derivation fails here. Covers `make_weights_vqtl`, both mixture builders (against a root-found mixture quantile), `make_weights_generic`, `make_fd_grad`, and the suite's own sinh-arcsinh truth columns at both a symmetric and a skewed baseline |
| `t0_covariance.R` | `compute_calibrated_se`, `compute_param_cov`, `compute_gls_param_gwas` against the linear algebra written out in R on a 5-tau grid; the **wire format** (row-major packed upper triangle, byte offsets); the pooled identity `Sigma_v = Sigma_0 / d_resid` |
| `t0_qpval.R` | the `Q` reference distribution. The moments against an exact eigen-decomposition; the closed-form `chi2_{T-k}` special case; the tail against `CompQuadForm::davies` (exact numerical inversion). Also pins that the naive `chi2` reference really is wrong, so a silent reversion to it is visible |
| `t0_tailref.R` | the degrees-of-freedom law `1/nu = 1/(B-1) + c/N_eff`, the contrast classifier, every rung of the fallback ladder, and each refusal. Including: there is **no** path that emits a normal-reference p-value without recording that it did |
| `t0_guards.R` | both dispersion guards against their published closed forms, including reproducing the documented `G` values, and `estimate_rav` recovering a known coupling with each of its two moment corrections isolated |
| `t0_package.R` | the package installs; `useDynLib` survives in source *and* in the installed `NAMESPACE`; the shared object loads; a compiled kernel actually executes; every man page parses (`tools::checkRd`); the expected exports exist |

It also re-runs eight pre-existing `stage2/tests/*.R` scripts unchanged.

**What tier 0 does NOT establish.** Nothing statistical. Every check here would
still pass if the estimator were systematically biased, if its intervals had
80% coverage, or if the method could not distinguish any two models. It tests
that the formulas are implemented as written, not that the formulas are right
for the job.

---

## Tier 1 -- integration

**Entry point:** `bash tests/run_tier1.sh`. **Runtime: 13 s of checks on top of
~40 s of fixture building** (which runs the shipped Stage 1 nine times). No
Slurm. Needs the Stage 1 conda environment because Stage 1 is Python and C++.

Fixtures are built by `tests/tier1/t1_make_fixtures.py` and
`t1_families.py`, which run **the shipped Stage 1** (`fungwas_stage1.core.run_stage1`,
`pooled_cov.compute_sigma0`) on simulated data and write the exact artefact set
the Stage 1 CLI writes -- including both JSON sidecars, built by Stage 1's own
`cov_meta` and `pooled_cov` modules rather than hand-written, so a schema change
breaks the suite loudly instead of silently.

**What it establishes.**

1. **Recovery of a known theta, end to end.** `param_gwas_from_file()` over a
   real artefact set recovers the exact population projection `theta_pop = A b(tau)`
   to within 5 standard errors on every causal variant, under jackknife and
   pooled covariance and under OLS and GLS. The int8 covariance reader agrees
   with the float32 payload to within 2% on the standard errors.
2. **The suite's own truth machinery is correct.** `T1.SIM.*` checks that the
   exactly-computed population estimand `b(tau)` agrees, per tau, with what the
   shipped Stage 1 estimates. Without this, every coverage number in tiers 1-3
   would be measuring an unverified target.
3. **Every refusal path.** Each is matched on its *message*, not merely on the
   fact of an error, so a refusal that starts firing for a different reason is a
   failure:

   | refusal | why it must stop rather than warn |
   |---|---|
   | covariance IDs mismatched / missing | rows silently misaligned |
   | `Sigma_0` tau grid differs | `Sigma_0` is indexed by tau |
   | `Sigma_0` covariate set differs | measured cost +1.7% on bmi, **+25.6% on weight**, same direction for every variant |
   | `Sigma_0` sample count differs | wrong scale, invisible downstream |
   | **`Sigma_0` LOCO scale mode differs** | the `docs/LOCO_SCALE.md` defect: shape effects attenuated 2.2x and SEs inflated 2.2x, with nothing downstream able to see it |
   | no `analysis` block in the metadata | nothing to check the sidecar against |
   | missing `d_resid` | no per-variant scale exists |
   | future schema version | the package refuses to guess |
   | `se_mode="tau_cov"` with no covariance | no silent downgrade to `diagonal` |
   | `nrow(W)` wrong / rank-deficient `W` | no recycling, no nonsense |

4. **The defaults are the recommended ones.** `tau_cov_estimator` defaults to
   `"ols"` and `cov_source` to `"jackknife"`, which is what `docs/OLS_VS_GLS.md`
   and `docs/KNOWN_ISSUES.md` section 10 concluded they must be.
5. **All eight model families**, each with known truth: null, pure location,
   pure scale, location+scale, skew, tail weight (kurtosis), a two-component
   mixture where the SNP moves the mixing weight, and one where it moves a
   component mean. For each: the true model recovers `theta_pop`, the true
   model is not rejected by `Q`, and the models that cannot represent the truth
   are.

**What tier 1 does NOT establish.** Rates. Every number in tier 1 comes from a
single replicate, so "recovers within 5 SE" is a smoke test of correctness, not
a measurement of coverage. Nothing here would catch a coverage of 0.90 instead
of 0.95.

---

## Tier 2 -- statistical validity

**Entry points:** `bash tests/run_tier2.sh {submit,status,collate,quick}`.
Cluster. This is the tier whose numbers decide whether a release is sound.

**The grid is preregistered** in `tests/tier2/cells.py` and the **pass criteria
in `tests/criteria.json`**, both committed before the run they judge. Every band
in `criteria.json` cites the measured result in `docs/` that it came from, and
is widened only by the Monte Carlo error of the tier 2 design. Nothing is
computed from the run it judges.

**Design.** 33 cells.

- 28 main cells at N = 50,000 on the production 45-tau grid, 1,000 jackknife
  blocks, 300 replicates each. Six sinh-arcsinh families
  (`null, location, scale, locscale, skew, kurtosis`) x {MAF 0.30, 0.01} x
  {normal, skewed baseline}, plus two mixture families x {MAF 0.30, 0.01}.
- Each replicate carries 190 null variants and 10 causal ones at |Z| = 8, so a cell yields
  **57,000 null observations** (Monte Carlo error on a 0.05 rate: 0.0009) and
  3,000 causal ones (error on a 0.95 coverage: 0.004).
- Every cell is analysed under **four estimator arms at no extra Stage 1 cost**:
  {jackknife, pooled} x {OLS, GLS}, each fitting both the location-scale model a
  real scan would use and the family's own true model, plus the models that
  cannot represent the truth.
- 5 pin cells reproduce established failure modes (below).

**Coverage is of `theta_pop`,** the projection of Stage 1's own exact population
estimand, computed from the conditional CDF rather than by simulation. That is
what the estimator is entitled to be judged on. Coverage of the *generative*
parameter is a different number under a misspecified `W`, and is reported
separately rather than conflated with it.

**The blocking criteria**, in full, from `criteria.json`:

| | OLS | GLS | from |
|---|---|---|---|
| null coverage of 95% intervals | 0.940-0.958 | 0.934-0.960 | BP2 s2 measures 0.9491-0.9508 (OLS) and 0.9433-0.9438 (jackknife GLS) |
| null Type-I at 0.05 | 0.042-0.058 | 0.040-0.062 | MEAN_VARIANCE s6: 0.0488-0.0503 over six cells |
| null Type-I at 1e-3, over nominal | 0.4-2.2 | same | KNOWN_ISSUES s7: 1.042 / 1.074 after the `t` fix |
| `|sd(z) - 1|` under the null | <= 0.03 | same | MEAN_VARIANCE s6: 0.9955-1.0016 |
| `Q` rejection at 0.05, correctly specified | 0.02-0.12 | same | calibrated reference |
| coverage under a true effect, normal baseline | >= 0.925 | >= 0.900 | OLS_VS_GLS: 0.9409 / 0.9373 under a scale effect |
| coverage under a true effect, skewed baseline | >= 0.880 | >= 0.850 | FINDINGS s3: down to 0.9114 at MAF 0.05, Z=30 |
| `|bias| / empirical sd` | <= 0.15 | <= 0.30 | OLS_VS_GLS binned by lack of fit |
| reported SE / empirical sd | 0.92-1.08 | same | FINDINGS s4: 0.898-1.003 across arms |
| `Q` rejection of a non-spanning model | >= 0.50 | same | no prior measurement; see the caveat below |

**The pins.** Four established failure modes. For each, the *blocking* half is
that the default, recommended path stays sound; the presence of the defect on a
non-default path is recorded, and raised as advisory where the evidence base
says it must be fixed.

| pin | blocking | reported |
|---|---|---|
| **false vQTL** from a pure location effect at skewed MAF | the `G` guard does not MISS: where `G < 0.05` the false-dispersion rate is 0.035-0.070, and at MAF 0.50 it is nominal for a `Z_mu = 60` location effect | the rate in the `G >= 0.2` bucket, i.e. that the artefact is still reproducible |
| **tie / GLS singularity** on a rounded phenotype | OLS keeps its nominal size | the GLS collapse, and that ridge repair fires |
| **ceiling-censoring anti-conservatism** | OLS is not anti-conservative at `alpha = 1e-4` | the GLS far-tail ratio; **advisory** if it is above 2x, since nothing in the shipped code warns a user |
| **LOCO scale mismatch** | tested in tier 1 as a refusal | -- |

**What tier 2 does NOT establish.**

- **Rates at 5e-8.** The deepest level any of this reaches by counting is
  `1e-4`, and only in the pin cell built for it. Genome-wide significance is
  unreachable by simulation (about `1e10` replicates per cell). The existing
  evidence base extrapolates to 5e-8 analytically; this suite does not, and does
  not check the extrapolation.
- **Real genotypes.** Every genotype is an independent 3-point HWE draw. No LD,
  no genotyping error, no Hardy-Weinberg departure, no relatedness, no
  population structure.
- **LOCO residualisation.** The simulated cells have no LOCO predictor, so the
  interaction between LOCO and the rest of the pipeline is untested beyond the
  refusal in tier 1.
- **Real phenotypes.** Nothing here has been checked against an actual trait.
- **Covariates with real effects.** The cells carry two covariates that are
  residualised out but do not affect Y.
- **`Q` on shape families.** Measured in tier 1 and reported, not asserted:
  the lack-of-fit statistic rejects a model missing a **location** or **scale**
  dimension every time, and does **not** reject one missing only a **shape**
  (skew or tail-weight) dimension, at any effect size where the linearisation
  still holds. The reason is in tier 0: the shape columns retain roughly 1-5%
  of their norm after projecting out location and scale, so they are weakly
  identified. A criterion asserting shape separation would be asserting
  something that is not true, so `wrong_shape` models are reported rather than
  required to be rejected, and the BIC selection accuracy for those families is
  reported rather than floored.
- **Model selection accuracy in general.** `docs/KNOWN_ISSUES.md` section 4
  measures raw-`Q` BIC recovering an additive truth 59% of the time. The suite
  reports selection accuracy per family; the only blocking criterion is that a
  **spanning** model is selected, and only for the families where separation is
  possible at all.

---

## Tier 3 -- realistic scale and timing

**Entry point:** `bash tests/run_tier3.sh [n_samples] [n_variants]`. Cluster.

**What it does.** Generates a **synthetic** BGEN with `plink2 --dummy` at
production sample size -- no UK Biobank, no ALSPAC, nothing real -- and runs the
Stage 1 CLI and then Stage 2 over it, measuring:

- wall clock, CPU time and peak resident set size, separately for each stage,
  from `/usr/bin/time -v`;
- a phase breakdown: BGEN read, RIF build, `Sigma_0`, kernel, write for Stage 1;
  TSV read, covariance read, kernel (OLS and GLS), write for Stage 2. Stage 1's
  own log already instruments compute and write, so the BGEN read phase is what
  remains of the batch loop after those two -- obtained without editing the
  shipped code;
- how Stage 2's cost grows with the variant count, by repeating the covariance
  payload to 10x and 40x. That part is synthetic and is labelled as such.

The hardware is recorded alongside every number (`hardware.txt`: CPU model,
socket and core counts, threads used, total memory, Slurm node), because a
timing with no machine attached to it is not reproducible and should not be
quoted.

### Measured, 2026-09-14

N = 461,320, 5,000 variants, the production 45-tau grid, 1,000 jackknife blocks,
int8 covariance, 8 threads on an **Intel Xeon Gold 6126 @ 2.60 GHz** (2 sockets
x 12 cores, 187 GB), BluePebble node `bp1-compute010`. The BGEN is 513 MB of
`plink2 --dummy` output.

| | wall | CPU | peak RSS |
|---|---|---|---|
| **Stage 1** | 68.3 s | 248.3 s | 3.09 GB |
| **Stage 2** | 21.7 s | 105.7 s | 2.23 GB |

Stage 1 by phase: **BGEN read 13.1 s**, kernel 27.5 s, write 15.4 s, RIF build
0.49 s, `Sigma_0` 0.07 s. Throughput 118.8 variants/s. The write phase is over
half the cost of the kernel, which is worth knowing before optimising the
kernel.

Stage 2 by phase (of 21.7 s wall, most of which is R start-up and package load):
TSV read 0.21 s, covariance read 0.48 s, OLS kernel + Q calibration 0.43 s,
write 0.05 s.

Stage 2 scaling, covariance payload tiled (synthetic):

| variants | kernel seconds | `cov_vec` in memory |
|---|---|---|
| 5,000 | 0.38 | 0.04 GB |
| 50,000 | 3.72 | 0.39 GB |
| 200,000 | 14.46 | 1.54 GB |

Linear in the variant count. Extrapolated to 1.2M variants: **Stage 1 about 16.6
CPU hours**, Stage 2 about 1.4 wall hours and roughly 9 GB of covariance in
memory -- consistent with the ~10 GB that `docs/KNOWN_ISSUES.md` section 10
quotes.

### One finding tier 3 produced

**GLS solved 0 of 5,000 variants.** The shipped default covariance is a block
jackknife stored as int8, and quantisation destroys positive semi-definiteness,
so there is nothing for the GLS solve to invert. This reproduces
`docs/KNOWN_ISSUES.md` section 10 exactly ("GLS solved, 1,986 chr22 variants:
jackknife **0**, pooled 1,986). `tau_cov_estimator = "gls"` is therefore
unusable as shipped unless `cov_source = "pooled"` or `--cov-dtype float32` is
also chosen. Raised as **advisory**, since OLS is the default and is unaffected:
`T3.S2.gls_unusable_on_int8_jackknife`.

**What tier 3 does NOT establish.** Anything about correctness -- the phenotype
is a null standard normal. It also does not measure a genome-wide scan: the
per-variant costs are multiplied out to 1.2M variants and that number is
labelled an extrapolation, defensible for time (Stage 1 is linear in the variant
count by construction) and *not* used for memory, which the scaling rows measure
directly instead. And it measures one machine: a different node, a different
filesystem or a different thread count will give different numbers, which is why
the hardware is recorded with them.

---

## The release report

```bash
Rscript --vanilla tests/report.R $GMOQ_WORK/results
```

writes `RELEASE_REPORT.json` and `RELEASE_REPORT.md`. Every check carries a
stable id, a status, a severity, its measured statistic, its Monte Carlo error
where it has one, the criterion it was judged against, and one sentence saying
what a pass would establish.

**The verdict rule, stated before any numbers are seen:**

- **SOUND** -- every blocking check passed, and tiers 0, 1 and 2 all ran.
- **SOUND (PARTIAL)** -- every blocking check that ran passed, but a tier did
  not run. The report names which, and what it would have established.
- **NOT SOUND** -- at least one blocking check failed.

Three severities:

- **blocking** -- a failure means the release is not sound.
- **advisory** -- a real defect that does not by itself invalidate results (a
  documentation or hygiene problem, or a defect confined to a non-default code
  path). Listed; does not change the verdict.
- **report** -- a measurement with no pass/fail at all, recorded so the number
  is on the record. Never counted as a pass.

A tier that did not run is never counted as a pass.

---

## Known gaps in the shipped code that the suite detects but does not block on

These are recorded here rather than buried in a log:

1. **`param_gwas_from_file()` is exported but undocumented.** Its roxygen block
   is attached to the internal `.validate_cov_ids()` instead, so `man/` has a
   `dot-validate_cov_ids.Rd` carrying the public function's documentation.
   `R CMD check` would flag it. `T0.PKG.exports_documented`, advisory.
2. **No tie guard.** `docs/BP2_lumpy_phenotypes.md` section 5.5 proposes
   refusing GLS when two taus share a sample quantile. Nothing in the shipped
   code does. `T1.TIE.no_shipped_tie_guard`, reported.
3. **No warning for GLS on a censored trait.** The one anti-conservative
   failure in the evidence base. `T2.PIN.ceiling_gls_still_unguarded`, advisory.
4. **The `G >= 0.2` flag is hard-coded** in `dispersion_report()`.
   `docs/MEAN_VARIANCE_ARTEFACT.md` states it is calibrated for OLS on a normal
   phenotype and that GLS on a skewed trait needs `G >= 0.06`; there is no way
   to set it. `T0.G.threshold_is_not_gls_calibrated`, reported.
5. **The R `quantile_gwas()` bandwidth is not the Stage 1 bandwidth.** R uses
   `stats::bw.nrd0`, Stage 1 uses `1.06 sd N^-0.2`. By the equivariance
   identity this cancels from every per-tau Wald test, but not from `Q`.
   `T0.RIF.bandwidth_rules_differ`, reported.

---

## Two notes on the simulator, because they change what the results mean

**The scale effect is multiplicative.** `sd = s0 exp(b.g / s0)`, not `s0 + b.g`.
An additive scale coefficient large enough to be interesting drives the
conditional standard deviation negative for some genotypes, at which point the
"conditional distribution" is a reflected normal and the simulation no longer
has the truth it claims. This was found by the suite's own
`T1.SIM.*` check during development: the exact estimand disagreed with Stage 1
by 24 standard errors.

**The shape families use a sinh-arcsinh baseline, not a Cornish-Fisher
expansion.** A Cornish-Fisher cubic stops being monotone once the shape
coefficients grow, which caps the largest simulable kurtosis effect at
`|Z| ~ 1.5` at N = 50,000. That cap is a property of the *simulation*, and in
the results it would have been indistinguishable from a property of the
*method*. The sinh-arcsinh family is monotone for every skewness and every
positive tail weight, and its inverse is closed form, so no root-finding enters
the truth at all.
