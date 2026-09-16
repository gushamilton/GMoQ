# Can Stage 1 fuse several phenotypes into one pass over the genotypes?

**Yes, exactly, and it is a 1.8-2.4x loss. Do not do it.**

The implementation is real, bit-identical and tested; it is kept on
`archive/perf/multi-pheno` rather than on `main`, because shipping a slower
path means every future kernel change has to keep it working. This document is
the result; the branch is the evidence.

## The premise was wrong

The idea was that decoding the BGEN once and reusing it across P phenotypes
would amortise a large fixed cost. Under the default `--pipeline` there is no
such cost to amortise: the single-phenotype run spends 145.5 s of a 163.8 s
batch loop inside the worker, so the *unhidden* decode is about 18 s and setup
about 21 s. Sharing both across four phenotypes caps the prize at **~6%**.

The 74 s decode quoted in the old docstring was a `--no-pipeline` attribution.

## What the existing multi path already did

It was already fused as well as it can be: one `_decode_variants` loop, one
`G_batch`, one `_process_batch` per batch; one shared `Q`, one projection, one
`g_resid`. What repeats per phenotype is only the `N x T` accumulation
`stats[base+1+t] += g_i * rif_p[i + t*N]`, plus the per-phenotype writes.

That is the part that *cannot* be shared. Per variant the kernel is one `N x K`
projection (K=13, common to all phenotypes) and one `N x T` RIF accumulation
(T=45, per phenotype by definition).

## Measured (full chr22, 17,550 variants, N=399,710, T=45, blocks=1000, `--cov-mode full`, 4 physical cores)

Each configuration ran its runs sequentially on one node; nodes differ by up to
35%, so rows from different blocks are not comparable.

No LOCO, existing `--pheno-cols` path (bp1-compute198):

| run | P | wall | vs separate | peak RSS |
|---|---|---|---|---|
| single x4 | 1 | 167/160/159/164 s | - | 2.8-2.9 GB |
| fused | 2 | 655 s | 327 s separate = **2.00x** | 3.3 GB |
| fused | 4 | 1,156 s | 650 s separate = **1.78x** | 4.1 GB |

`--loco-mode covariate`, new `--fuse-loco-multi` (bp1-compute154):

| run | P | wall | vs separate | peak RSS |
|---|---|---|---|---|
| single x4 | 1 | 291/290/287/283 s | - | ~2.87 GB |
| fused | 2 | 1,526 s | 581 s separate = **2.63x** | 3.4 GB |
| fused | 4 | 2,756 s | 1,151 s separate = **2.39x** | 4.1 GB |

Genome-wide projection for the four-phenotype testosterone series: 19.2 h run
separately against ~46 h fused, i.e. **27 hours lost**. Even if fusion were
free it would have saved 0.84 h (4.4%).

## It is now worse than these numbers say

The fused path refuses `--panel-kernel` and silently turns it off for P > 1.
Since the restructured kernel landed, `--panel-kernel` **is** the default and
is 4.2x faster on chr22. So fusion today costs the measured 1.8-2.4x *and* the
fast kernel on top of it. The measurements above were taken with both arms on
the old kernel and are the generous case.

## What was proved on the way

- **Missingness is not a blocker** for this series: all four phenotypes have
  exactly the same 421,405 observed samples under the shared keep+covar files,
  symmetric difference 0. Covariates are not a blocker either.
- **LOCO is only bookkeeping.** The bases `[X, L_p(chr)]` differ in one column,
  and Householder QR builds reflector *j* from columns 1..*j* only, so the
  leading K columns of `qr([X, L_p])` are the leading K columns of `qr(X)`,
  bit for bit. Verified rather than assumed, with a Gram-Schmidt fallback.
  That is what let the fused kernel share `proj_base` per variant and add only
  each phenotype's own `u'g` term.
- **Memory is not the constraint.** Peak RSS 2.9/3.4/4.1 GB at P=1/2/4, about
  0.4 GB per extra phenotype. Roughly 150 phenotypes would fit in a 64 GB
  allocation. Speed is the constraint.

## Exactness

Fused output was **byte-identical** to the separate runs -- `stage1.tsv.gz`,
`sigma0.tsv.gz`, `cov.gz` and `cov.ids.tsv.gz`, on decompressed content, full
chr22, every column class, at P=2 and P=4, with and without LOCO. Max absolute
difference 0, so no Stage 2 comparison was needed. A single-phenotype run on
the modified source was byte-identical to the unmodified source.

## Not attributed

The real runs lose more than a kernel-only benchmark predicts (that benchmark
puts fused at 0.87-1.22x of separate, a few percent either way). The leading
candidate is that P resident RIF matrices of 144 MB each, streamed per variant,
saturate the memory bandwidth the decoder thread also needs. This was not
pinned down: the microbenchmark's absolute throughput was 3x slower than the
same kernel inside the real run on the same node, unexplained, so it is trusted
only for ratios at fixed conditions.

## Provenance note, unrelated to fusion

A fresh single run differs from the published
`paper-core-genome-20260915-02/testosterone/chr22` Stage 1, max
|delta beta_tau_0.500| = 1.2. The unmodified-source control was byte-identical
to the modified one, so this is the `--loco-mode` default flip in 0.4.0, not a
regression. Someone should confirm the published artefacts are meant to be
`subtract`-era.

## Recommendation

Run phenotypes separately. Revisit only if the kernel stops being the critical
path -- `archive/perf/multi-pheno` holds the mechanism, and it is correct.
