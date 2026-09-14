#!/usr/bin/env python
"""Collapse per-chunk npz into flat binaries Stage 2 (R) can readBin."""
import argparse, glob, json, os, sys
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    files = sorted(glob.glob(os.path.join(a.dir, "raw_*.npz")))
    assert files, a.dir
    metas = sorted(glob.glob(os.path.join(a.dir, "meta_*.json")))
    meta = json.load(open(metas[0]))
    T = len(meta["taus"]); nel = T * (T + 1) // 2
    os.makedirs(a.out, exist_ok=True)

    fb = open(os.path.join(a.out, "beta.bin"), "wb")
    fs = open(os.path.join(a.out, "se.bin"), "wb")
    fc = open(os.path.join(a.out, "cov.bin"), "wb")
    fp = open(os.path.join(a.out, "pooled.bin"), "wb")
    fd = open(os.path.join(a.out, "meta_snp.bin"), "wb")   # rep, snp_idx, maf, d_resid
    diag = dict(fhat=[], tie=[], fgap=[], spacing=[], hbw=[], q_tau=[],
                ndistinct=[], maxatom=[], ndist_h=[])
    nsnp = 0
    rep_off = 0
    for f in files:
        z = np.load(f)
        B, S, C, DR, SG = z["beta"], z["se"], z["cov"], z["d_resid"], z["sigma0"]
        R, M, _ = B.shape
        fb.write(np.ascontiguousarray(B.reshape(R * M, T), dtype=np.float64).tobytes())
        fs.write(np.ascontiguousarray(S.reshape(R * M, T), dtype=np.float64).tobytes())
        fc.write(np.ascontiguousarray(C.reshape(R * M, nel), dtype=np.float64).tobytes())
        pooled = SG[:, None, :] / DR[:, :, None]
        fp.write(np.ascontiguousarray(pooled.reshape(R * M, nel), dtype=np.float64).tobytes())
        mm = np.empty((R * M, 4))
        mm[:, 0] = np.repeat(np.arange(R) + rep_off, M)
        mm[:, 1] = np.tile(np.arange(M), R)
        mm[:, 2] = np.tile(z["mafs"], R)
        mm[:, 3] = DR.reshape(-1)
        fd.write(np.ascontiguousarray(mm, dtype=np.float64).tobytes())
        for k in diag:
            diag[k].append(z[k])
        nsnp += R * M
        rep_off += R
    for f_ in (fb, fs, fc, fp, fd):
        f_.close()
    out = {k: np.concatenate(v).tolist() for k, v in diag.items()}
    meta["n_snp"] = nsnp
    meta["n_rep"] = rep_off
    meta["M"] = M
    with open(os.path.join(a.out, "meta.json"), "w") as f:
        json.dump(meta, f)
    np.savez_compressed(os.path.join(a.out, "diag.npz"), **{k: np.concatenate(v) for k, v in diag.items()})
    print(f"{a.dir}: {rep_off} reps x {M} snps = {nsnp}")


if __name__ == "__main__":
    main()
