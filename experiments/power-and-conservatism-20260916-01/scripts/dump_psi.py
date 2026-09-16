"""Dump the residualised RIF matrix Psi_r that Stage 1 builds, via the real code path.

Runs the production CLI on a one-SNP extract list and intercepts
compute_chr_rif_residuals, so Psi_r is exactly what the shipped run used.
"""
import os, sys
import numpy as np
import fungwas_stage1.cli as cli

OUTDIR = os.environ["PSI_OUT"]
os.makedirs(OUTDIR, exist_ok=True)

_orig = cli.compute_chr_rif_residuals

def hooked(pheno_cols, y_values_full, X_full, taus, loco_predictions, chrom,
           valid_mask, loco_scales=None, **kw):
    out = _orig(pheno_cols, y_values_full, X_full, taus, loco_predictions, chrom,
                valid_mask, loco_scales=loco_scales, **kw)
    rif_resid_list, q_chr, q_tau_list = out
    P = np.asarray(rif_resid_list[0], dtype=np.float64)
    tag = str(chrom)
    np.save(os.path.join(OUTDIR, f"psi_chr{tag}.npy"), P)
    # The residualising basis the kernels also apply to the genotype.
    if loco_predictions is not None and kw.get("loco_mode") == "covariate":
        loco_col = cli._loco_column_for_chr(loco_predictions, pheno_cols,
                                            cli.normalize_chr_label(chrom), valid_mask)
        X_chr = np.column_stack([X_full[valid_mask], loco_col])
    else:
        X_chr = X_full[valid_mask]
    Q, _ = np.linalg.qr(X_chr)
    np.save(os.path.join(OUTDIR, f"Q_chr{tag}.npy"), np.asarray(Q, dtype=np.float64))
    sys.stderr.write(f"[hook] Q {Q.shape} max|Q'Q-I| "
                     f"{np.abs(Q.T @ Q - np.eye(Q.shape[1])).max():.3e}\n")
    np.save(os.path.join(OUTDIR, f"qtau_chr{tag}.npy"), np.asarray(q_tau_list[0]))
    np.save(os.path.join(OUTDIR, f"taus.npy"), np.asarray(taus))
    np.save(os.path.join(OUTDIR, f"y_chr{tag}.npy"),
            np.asarray(y_values_full[pheno_cols[0]], dtype=np.float64)[valid_mask])
    with open(os.path.join(OUTDIR, f"meta_chr{tag}.txt"), "w") as fh:
        fh.write(f"n={P.shape[0]}\nT={P.shape[1]}\nchrom={tag}\nvalid={int(valid_mask.sum())}\n")
        fh.write(f"X_full_cols={X_full.shape[1]}\nloco_mode={kw.get('loco_mode')}\n")
    sys.stderr.write(f"[hook] dumped Psi_r {P.shape} for chr {tag}\n")
    return out

cli.compute_chr_rif_residuals = hooked
sys.argv = ["fungwas_stage1.cli"] + sys.argv[1:]
cli.main()
