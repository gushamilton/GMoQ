# Stage 1 fast panel kernel

Measured on BluePebble, 2026-09-16, against `main` at `11007fa` (v0.3.0).
Everything below is chr22 testosterone in the accepted gate's configuration:
N = 399,710 samples, 17,550 variants, T = 45 taus, K = 13 covariate basis
columns, `--blocks 1000`, `--batch-size 100`, `--threads 4`.

## What was wrong

The exact panel kernel (`panel_block_stats`, `src/stage1_pybind.cpp`) is a pure
re-nesting of the original per-SNP loops, so it is bit-identical to them and
inherits their access pattern. It was 139.9 s of a 197 s run, at about
0.9 GFLOP/s per core - roughly 7% of the SSE2 baseline peak. Three causes:

| cause | effect |
| --- | --- |
| `G_ptr[i + cols[p]*N]` strides down the long axis, and is read twice per (sample, SNP) | every genotype access is a separate cache line, fetched twice |
| the K = 13 covariate projection `proj = sum_k qrow[k]*qt[k]` is recomputed per (sample, SNP) | a serial reduction chain in the innermost loop; nothing pipelines |
| every sample accumulates into the block tile its block id selects | with i.i.d. uniform block ids that is an unrelated 5.8 kB tile per sample; the write-hot set is `n_blocks * P * (1+T) * 8` bytes *per thread* and leaves L3 as soon as a few threads run |

The last one is why the exact kernel stops scaling: on a 2 x 14-core Xeon Gold
5120 (19.25 MB L3) it is the same speed at 16 threads as at 8.

## What the fast kernel does

Write `r_ip = g_ip - Q_i . q_p`, `q_p = Q^T g_p`. The jackknife needs

    Out[b][p][0]   = sum_{i in b} r_ip^2
    Out[b][p][1+t] = sum_{i in b} r_ip * RIF[i,t]

Expanding `r` and collecting the SNP-independent parts gives

    Out[b][p][1+t] = W_b[p][t] - sum_k q_p[k] S_b[k][t]
    Out[b][p][0]   = s2_b[p] - 2 sum_k q_p[k] u_b[k][p] + q_p^T A_b q_p

with `S_b = Q_b^T RIF_b` and `A_b = Q_b^T Q_b` depending only on Q, RIF and the
block assignment. Those are built **once per Stage 1 run** (`BlockPlan`), not
once per (sample, SNP): the covariate projection that used to run `N * P` times
per panel now runs `n_blocks * P` times, 400x less often at production N.

The samples are counting-sorted by block id once, so each block owns a
contiguous row range. The three genotype-dependent primitives - `W_b`, `u_b`,
`s2_b` - are then one fused streaming sweep over that range, sequential in every
operand, one genotype read per (sample, SNP). The block tiles are *written*
once each instead of being scattered into, so the write-hot set disappears.

Implementation notes that turned out to matter:

* **Panel-major repacking.** The panel's genotype columns are packed once into an
  `(N, LD)` row-major tile in block-sorted order, `LD = round_up(P, 8)` so tile
  rows are 64-byte aligned. The gather reads down G's contiguous axis; the
  scatter moves whole rows. Written with non-temporal stores, which removes the
  read-for-ownership traffic (the tile is far bigger than L3 and is read exactly
  once afterwards): 48 ms -> 24 ms per batch at 4 threads.
* **Register blocking.** The sweep keeps a `CBLK x PBLK = 6 x 8` accumulator tile
  in 12 ymm registers across the whole sample loop: two vector loads, six
  broadcasts, twelve FMAs per sample and no other memory traffic. GCC does not
  find this shape on its own - it sinks the accumulator to the stack - and that
  one change was worth 2.7x on the sweep (0.196 s -> 0.073 s per batch at 4
  threads, against the same loop auto-vectorised with AVX2 available).
* **ISA.** The hot leaves are compiled twice, once at the translation unit's
  baseline and once under `#pragma GCC target("avx,avx2,fma")`, and dispatched on
  `__builtin_cpu_supports`. That gives the fast kernel AVX2+FMA without raising
  the ISA of the exact kernel, whose arithmetic must not change. (GCC's
  `target_clones("arch=haswell", ...)` is **not** a substitute: it dispatches on
  `__builtin_cpu_is("haswell")`, i.e. the CPU model, so a Skylake-SP or Cascade
  Lake node silently takes the baseline clone.)
* **Scratch reuse.** The per-batch buffers live on the `BlockPlan` and are
  allocated uninitialised and once. Zero-filling a 333 MB tile that is about to
  be overwritten, and paying its first-touch page faults, was 20% of the kernel.
* **Centring.** The tile holds `g - gbar` rather than `g`. See Numerics.

## Speed

Production-shaped batch (M = 100 SNPs, one call), best of 3, exclusive node,
2 x 14-core Xeon Gold 5120 @ 2.2 GHz, `OMP_PLACES=cores`. Thread count equals
**physical** cores throughout.

| physical cores | exact (s) | fast (s) | speedup | exact scaling | fast scaling |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 2.041 | 0.522 | 3.91x | 1.00 | 1.00 |
| 2 | 1.318 | 0.272 | 4.84x | 1.55 | 1.92 |
| 4 | 0.724 | 0.146 | 4.96x | 2.82 | 3.57 |
| 8 | 0.426 | 0.078 | 5.45x | 4.79 | 6.68 |
| 12 | 0.427 | 0.062 | 6.85x | 4.78 | 8.37 |
| 16 | 0.425 | 0.048 | 8.78x | 4.80 | 10.77 |

The exact kernel is flat from 8 cores on; the fast kernel is still gaining at 16.
GFLOP/s here counts only the dominant `2*N*M*T` term: the fast kernel reaches
74 GFLOP/s at 16 cores against the exact kernel's 8.5.

Full chr22, same node class, 4 threads, gate configuration:

| | exact (`--exact-kernel`) | fast (default) |
| --- | ---: | ---: |
| Stage 1 wall clock | 190 s | 104 s |
| of which kernel ("batch compute") | 138.7 s | 33.0 s |
| kernel speedup | | **4.20x** |

What each change bought, per batch at 4 threads (cumulative, in the order applied):

| step | kernel (s) | note |
| --- | ---: | --- |
| exact panel kernel | 0.659 | baseline |
| decomposed + sorted + panel-major repack, baseline ISA | 0.412 | the restructuring alone, no SIMD work |
| + AVX2/FMA clone of the hot leaves, + stop zero-filling the reused scratch | 0.284 | the zero-fill alone was 0.082 of the 0.412 |
| + register-blocked sweep, AVX2 intrinsics (CBLK = 4) | 0.160 | sweep 0.196 -> 0.073; GCC will not find this shape |
| + non-temporal stores in the tile packer, scratch kept across batches | 0.130 | packer 0.048 -> 0.024 |
| + CBLK 4 -> 6, + mean-centring and its exact correction | 0.131 | sweep -10%, centring gives it back; kept for the numerics |

Per-pass breakdown of the final kernel (4 physical cores, per 100-SNP batch):

| pass | s | % |
| --- | ---: | ---: |
| MAC filter (one sequential read of G) | 0.013 | 10 |
| pack the block-sorted genotype tile | 0.027 | 20 |
| streaming sweep (W, u, s2) | 0.067 | 50 |
| q = sum_b u_b | 0.001 | 0.5 |
| covariate correction + jackknife | 0.027 | 20 |

## Numerics

The fast kernel is **not** bit-identical to the exact one: the same quantity is
formed from a different, algebraically equal expression. Summation order is
otherwise preserved - within a block the samples are visited in ascending
original row index, which is exactly the order the exact kernel adds them in,
because the counting sort is stable.

The error is governed by the cancellation ratio of the decomposition,

    kappa = ||g||^2 / D,   D = sum_i r_i^2

i.e. by how much of the genotype's second moment the intercept and the covariate
basis absorb. Measured on a synthetic SNP swept across kappa (uncentred form):

| sd(dosage) | kappa | max abs error / max abs value | eps * kappa |
| ---: | ---: | ---: | ---: |
| 0.5 | 5.4e0 | 2.6e-14 | 1.2e-15 |
| 0.1 | 1.0e2 | 2.2e-13 | 2.3e-14 |
| 0.01 | 1.0e4 | 2.5e-12 | 2.2e-12 |
| 0.001 | 1.0e6 | 2.6e-11 | 2.3e-10 |

At production this bites on variants whose effect-allele frequency is close to 0
or 1, where kappa reaches ~1e4 and the block statistics lose four digits. The
kernel therefore works on the **mean-centred** genotype, which removes the
intercept's contribution to kappa. Centring is exact only when the intercept lies
in col(Q), so the plan also carries `e = 1 - Q (Q^T 1)` and the exact correction
terms `e_b^T RIF_b`, `Q_b^T e_b`, `e_b^T e_b`; `e` is zero to rounding in every
configuration the CLI builds, and `Q^T 1` is accumulated with compensated
summation because a naive sum of the intercept column carries an O(N * eps)
error that would otherwise inject more noise than the centring removes.

### (a) fast vs exact, full chr22, float64

Both kernels evaluated on every one of the 176 real genotype batches in the same
process (`bench/dualcheck.py`), worst disagreement anywhere on the chromosome:

| quantity | max abs difference | scale | ratio |
| --- | ---: | ---: | ---: |
| block statistics | 2.31e-9 | 1.28e3 | 1.80e-12 of the array max |
| `beta_tau_*` | 8.44e-15 | 1.11 | **5.35e-14** = max abs(delta beta)/SE |
| `se_tau_*` | - | - | **8.65e-13** max relative |
| packed covariance | 5.25e-13 | 0.387 | 1.35e-12 max relative |
| `d_resid` | 1.83e-9 | - | **1.01e-14** max relative |
| NaN pattern (accept/reject) | identical on all 176 batches | | |

The SE is the one number above pure reassociation level, and that is the
jackknife's own conditioning, not the kernel's: the SE is the spread of 1000
delete-one estimates that all sit within `SE/sqrt(B-1)` of each other, so a
relative perturbation eps in the block statistics appears in the SE amplified by
roughly `sqrt(B-1) * abs(beta)/SE`. The measured amplification, 8.65e-13 against
beta's 5.35e-14, is 16x - which is what that formula gives.

Before centring these numbers were 8.70e-11 (beta/SE), 3.36e-11 (SE) and
3.33e-11 (d_resid); centring bought between 39x and 3300x.

### (a) at output resolution

Stage 1 writes `%.6g`. Across all 17,550 x 96 fields of `stage1.tsv.gz`, the fast
kernel's output is **byte-identical** to the exact kernel's - and therefore to the
accepted gate, which `--exact-kernel` reproduces byte for byte. `cov.gz` is
float32 and differs in 36 of 18,164,250 elements, each by exactly one float32 ULP
(max relative 1.18e-7 = `float32 eps`); the largest implied change in any SE read
off that file is 4.9e-8 relative.

### (b) Stage 2

Stage 2 skew family, `pooled` covariance, `ols` tau-covariance estimator, reusing
the accepted gate's own `weights_skew.rds` and its pooled Sigma_0 sidecar so that
Stage 1's table is the only thing that varies:

| Stage 1 input | Stage 2 `chr22.tsv.gz` vs the accepted Stage 2 gate |
| --- | --- |
| the gate's own Stage 1 (harness check) | IDENTICAL |
| `--exact-kernel` | IDENTICAL |
| fast kernel (default) | **IDENTICAL** |

Zero difference in any column, and therefore zero variants changing discovery
status at 5e-8 or 1e-5 (reference hit counts: `beta_left_p` 20 / 46,
`beta_right_p` 27 / 59, `gamma_p` 0 / 3, `Q_pval` 0 / 0 at 5e-8 / 1e-5).

For completeness, the *pre-centring* fast kernel did move Stage 2 slightly - all
point estimates still bit-identical, 3 of 17,550 SNPs differing in SE/p columns
with worst relative difference 7.1e-9, and still zero discovery-status changes.
Those three rows were exactly the three where `d_resid` flipped a `%.6g` digit.

### (c) the bit-exact path

`--exact-kernel` selects the unmodified panel kernel. Run in the gate's own
configuration it reproduces the accepted gate byte for byte:

| file | `--exact-kernel` vs gate |
| --- | --- |
| `stage1.tsv.gz` | IDENTICAL |
| `cov.gz` | IDENTICAL |
| `cov.ids.tsv.gz` | IDENTICAL |
| `sigma0.tsv.gz` | DIFFERS - see below |

`sigma0.tsv.gz` differs for both kernels identically: the sidecar's schema went
from v1 to v2 between the gate's commit (`ccc2d97`) and `main` at `11007fa`. That
is pre-existing drift, not this branch. Similarly, `main`'s `--loco-mode` now
defaults to `covariate` where the gate ran `subtract`; with the default the whole
table moves by up to 4.4 SE regardless of which kernel is used, which is why every
comparison above pins `--loco-mode subtract`.

## Invariants the rewrite relies on, and where they are pinned

`stage1/tests/test_fast_kernel.py` (24 tests):

* agreement with the exact kernel to reassociation level, on block statistics,
  beta, SE and covariance;
* identical accept/reject decisions (the MAC filter is computed in the original
  order and is bit-identical);
* **bit-exact** independence of thread count (1, 2, 3, 4);
* **bit-exact** independence of panel width (0, 1, 2, 7, 8, 16, 64);
* **bit-exact** independence of how SNPs are split into batches (1, 5, 11, 23);
* the fused path agrees with finalising the raw statistics;
* out-of-range block ids, `T = 1`, `K = 1`, and panels that are not a multiple of
  the register-block width;
* inputs are not mutated.

Full Stage 1 suite: 146 tests before, 170 after, all passing.

## Where the time is now, and what would be needed to go further

Full chr22 at 4 threads, fast kernel: 104 s wall, of which

| | s |
| --- | ---: |
| setup (phenotype, covariates, RIF, LOCO, Sigma_0) | 18 |
| kernel ("batch compute") | 33 |
| write path (NumPy jackknife + TSV/cov writing) | 30 |
| BGEN decode not hidden behind the kernel, pipeline sync | 23 |

**The kernel is no longer the critical path.** It is now roughly tied with the
NumPy write path, and the decode that the double-buffered pipeline used to hide
behind a 140 s kernel no longer fits behind a 33 s one.

To go further, in order of expected payoff:

1. **The write path (30 s).** The per-SNP jackknife finalisation still runs in
   NumPy. The C++ fused finalisation already exists (`--fused-finalize`) and does
   that work in 0.027 s per batch inside the OpenMP loop - but on that path
   `d_resid` comes from `residual_sum_squares`, which builds two 320 MB NumPy
   temporaries per batch and cost 142 s of write time in a measured run (total
   206 s, i.e. *slower* than the default path). Returning `d_resid` from the
   kernel, which already has it, would make `--fused-finalize` a clear win. That
   is a Python-side change and is out of scope here; it is recorded as a negative
   result for the fused path as it stands.
2. **The decode (23 s).** Single-threaded Python per variant; now exposed.
3. **The kernel's remaining 33 s.** The sweep is half of it and runs at about
   21 GFLOP/s per core, roughly 6 flops/cycle against an AVX2 issue limit near 16.
   It is not bandwidth-bound (1 -> 4 cores is 97% efficient) and not L1-capacity
   bound (varying the sample-slice size from 4 kB to 1 MB moves it by under 10%),
   so the remaining gap is issue/latency in the sweep loop and would need
   profiling at the counter level rather than more restructuring. The MAC filter
   (10%) could be folded into the packing pass, but only by accumulating the
   genotype sums out of order, which would change the imputation mean and lose the
   bit-identical filter.
4. **AVX-512.** The nodes have it; the sweep is written for 4-wide ymm. Doubling
   the register block would halve the sweep's issue count. Not attempted.

## Negative results

* `--fused-finalize` with the fast kernel: 206 s wall against 104 s for the
  default path, because of `residual_sum_squares`'s NumPy temporaries (above).
* GCC's auto-vectoriser on the register-blocked sweep: it sinks the accumulator
  tile to the stack and produces a scalar-add reduction chain; 1.5 flops/cycle
  against 6 for the intrinsics version.
* `target_clones("arch=haswell", "default")`: compiles and produces an ifunc, but
  dispatches on CPU *model*, so it never selects the AVX2 clone on Skylake-SP or
  Cascade Lake. Silent 35% loss.
* Sample-slice (L1) blocking of the sweep: expected to matter, measured not to
  (under 10% across a 256x range of slice sizes). Kept, because it costs nothing
  and does help slightly.
* Previously recorded, not revisited: whole-batch NumPy BLAS (0.667 s) and
  sample-chunked BLAS (0.713 s) both lost to the old hand kernel (0.434 s) on
  their 320 MB temporaries. The formulation here avoids materialising `Gi` and
  `R` at all, which is why it wins.
