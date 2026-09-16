"""Full-precision fast-vs-exact kernel comparison on the real chr22 batches.

The Stage 1 TSV is written at %.6g, which hides differences of the size this
restructuring makes.  This runs the real CLI pipeline but, for every genotype
batch, evaluates BOTH kernels in float64 and accumulates the worst disagreement
seen anywhere on the chromosome.

Usage: dualcheck.py -- <all the normal fungwas_stage1.cli arguments>
"""
import sys
import numpy as np

from fungwas_stage1 import core, cli   # noqa: E402

STATE = {
    "batches": 0, "snps": 0,
    "raw_abs": 0.0, "raw_scale": 0.0,
    "d_abs": 0.0, "d_rel": 0.0,
    "beta_abs": 0.0, "beta_scale": 0.0, "beta_over_se": 0.0,
    "se_abs": 0.0, "se_rel": 0.0,
    "cov_abs": 0.0, "cov_scale": 0.0, "cov_rel_diag": 0.0,
    "nan_mismatch": 0,
}


def _acc(key_abs, key_scale, a, b):
    m = np.isfinite(a) & np.isfinite(b)
    if not m.any():
        return
    STATE[key_abs] = max(STATE[key_abs], float(np.max(np.abs(a[m] - b[m]))))
    STATE[key_scale] = max(STATE[key_scale], float(np.max(np.abs(a[m]))))


_orig = core.compute_block_scores_panel_fast


def patched(G, Y, Q, block_ids, n_blocks, min_mac=5.0, min_var=1e-8,
            n_threads=1, panel_width=0, plan=None):
    fast = _orig(G, Y, Q, block_ids, n_blocks, min_mac, min_var, n_threads,
                 panel_width, plan)
    exact = core.compute_block_scores_panel(G, Y, Q, block_ids, n_blocks,
                                            min_mac, min_var, n_threads, panel_size=16)
    if not np.array_equal(np.isnan(fast[:, 0, 0]), np.isnan(exact[:, 0, 0])):
        STATE["nan_mismatch"] += 1
    _acc("raw_abs", "raw_scale", exact, fast)

    # d_resid as the writer forms it: the block D column summed over blocks.
    df = np.nansum(fast[:, :, 0], axis=1)
    de = np.nansum(exact[:, :, 0], axis=1)
    ok = np.isfinite(de) & (de > 0)
    if ok.any():
        STATE["d_abs"] = max(STATE["d_abs"], float(np.max(np.abs(df[ok] - de[ok]))))
        STATE["d_rel"] = max(STATE["d_rel"],
                             float(np.max(np.abs(df[ok] - de[ok]) / de[ok])))

    # beta / SE / covariance, both kernels run through the same finalisation.
    bf, sf, cf = core.compute_block_finalize_panel_fast(
        G, Y, Q, block_ids, n_blocks, min_mac, min_var, n_threads, panel_width, plan)
    be, se_, ce = core.compute_block_finalize_panel(
        G, Y, Q, block_ids, n_blocks, min_mac, min_var, n_threads, panel_size=16)
    _acc("beta_abs", "beta_scale", be, bf)
    m = np.isfinite(se_) & np.isfinite(sf) & (se_ > 0)
    if m.any():
        STATE["se_abs"] = max(STATE["se_abs"], float(np.max(np.abs(se_[m] - sf[m]))))
        STATE["se_rel"] = max(STATE["se_rel"],
                              float(np.max(np.abs(se_[m] - sf[m]) / se_[m])))
        STATE["beta_over_se"] = max(
            STATE["beta_over_se"], float(np.max(np.abs(be[m] - bf[m]) / se_[m])))
    _acc("cov_abs", "cov_scale", ce, cf)

    STATE["batches"] += 1
    STATE["snps"] += G.shape[1]
    return fast


core.compute_block_scores_panel_fast = patched
cli.compute_block_scores_panel_fast = patched

if __name__ == "__main__":
    argv = sys.argv[1:]
    if argv and argv[0] == "--":
        argv = argv[1:]
    sys.argv = ["fungwas-stage1"] + argv
    try:
        cli.main()
    finally:
        print("\n=== full-precision fast vs exact kernel, real chr22 batches ===")
        print(f"batches={STATE['batches']} snps={STATE['snps']} "
              f"nan_pattern_mismatched_batches={STATE['nan_mismatch']}")
        rows = [
            ("block stats", STATE["raw_abs"], STATE["raw_scale"], None),
            ("beta_tau_*", STATE["beta_abs"], STATE["beta_scale"], STATE["beta_over_se"]),
            ("se_tau_*", STATE["se_abs"], None, STATE["se_rel"]),
            ("cov (packed)", STATE["cov_abs"], STATE["cov_scale"], None),
            ("d_resid", STATE["d_abs"], None, STATE["d_rel"]),
        ]
        print(f"{'quantity':<14} {'max |diff|':>12} {'max |value|':>12} "
              f"{'ratio':>12}")
        for name, a, sc, r in rows:
            scs = f"{sc:.4e}" if sc is not None else "-"
            rs = f"{r:.4e}" if r is not None else (f"{a/sc:.4e}" if sc else "-")
            print(f"{name:<14} {a:>12.4e} {scs:>12} {rs:>12}")
        print("ratio = max relative difference, except beta_tau_* where it is "
              "max |delta beta| / SE.")
