#!/usr/bin/env python
"""Truth-free tau guard: drop taus whose sample tie fraction exceeds a threshold."""
import argparse, json, os
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--dir", required=True)
ap.add_argument("--tie-max", type=float, default=0.01)
a = ap.parse_args()
z = np.load(os.path.join(a.dir, "diag.npz"))
tie = np.asarray(z["tie"]).mean(axis=0)
keep = [int(i + 1) for i in range(len(tie)) if tie[i] <= a.tie_max]
out = {"keep": keep, "tie": tie.tolist(), "tie_max": a.tie_max}
if len(keep) >= 3 and len(keep) < len(tie):
    json.dump(out, open(os.path.join(a.dir, "keep.json"), "w"))
    print("keep", len(keep), "of", len(tie))
else:
    p = os.path.join(a.dir, "keep.json")
    if os.path.exists(p):
        os.remove(p)
    print("no tie-drop applicable:", len(keep), "of", len(tie))
