#!/usr/bin/env python
"""Truth-free tau guard #2: drop a tau whose sample quantile collides with the
previously kept tau's sample quantile.

This is strictly weaker than the tie-fraction guard (it removes only the exact
rank-deficiency, keeping every tau that carries independent information) and is
computable from the Stage 1 output alone -- q_tau is already returned.
"""
import argparse, json, os
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--dir", required=True)
ap.add_argument("--collide-max", type=float, default=0.5,
                help="drop tau if it collides with the last kept tau in more "
                     "than this fraction of replicates")
a = ap.parse_args()

z = np.load(os.path.join(a.dir, "diag.npz"))
Q = np.asarray(z["q_tau"])          # reps x T
R, T = Q.shape

keep = [0]
coll = [0.0]
for j in range(1, T):
    p = float(np.mean(Q[:, j] == Q[:, keep[-1]]))
    coll.append(p)
    if p <= a.collide_max:
        keep.append(j)

out = {"keep": [int(i + 1) for i in keep], "collide": coll,
       "collide_max": a.collide_max, "n_keep": len(keep), "T": T}
p = os.path.join(a.dir, "keep_q.json")
if len(keep) >= 3 and len(keep) < T:
    json.dump(out, open(p, "w"))
    print("qdedup keep", len(keep), "of", T, "->", out["keep"])
else:
    if os.path.exists(p):
        os.remove(p)
    print("qdedup not applicable:", len(keep), "of", T)
