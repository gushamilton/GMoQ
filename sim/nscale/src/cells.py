"""Cell design for the N-scaling coverage study.

Three rungs of the N ladder: simulation N (30,000), an intermediate (100,000)
and production N (461,320 -- the bmi sample actually matched into the BGEN).

Three arms, so that N, MAF and MAC can be separated:

  A  MAF 0.30            common; MAC grows with N (18k -> 60k -> 277k)
  B  MAF 0.01            constant MAF; MAC grows with N (600 -> 2000 -> 9226)
  C  MAF = 922/(2N)      constant MAC ~ 922, matching the rarest chr22 variant
                         at production (MAC 927).  MAF falls with N.

Arm B and arm C coincide at no rung, so a difference between them at fixed N
is a MAF effect, and a change along arm C is a pure N effect at fixed MAC.
"""
import json
import sys

NS = [30000, 100000, 461320]
TARGET_MAC = 922.0
TAG = {30000: "N30k", 100000: "N100k", 461320: "N461k"}


def build_cells():
    cells = []
    for N in NS:
        tag = TAG[N]
        # Uniform 20,000 replicates at every rung, so that the Monte Carlo
        # error on a 0.95 coverage rate is 0.0015 everywhere and rungs are
        # directly comparable.  Measured cost at N = 461,320 is 1.46 s/rep.
        big = 20000
        small = 20000
        arms = []
        # arm A: common variant, all three DGPs, both shapes
        for shape in ("normal", "gamma2"):
            for dgp, a, c in (("null", 0.0, 0.0), ("scale", 0.0, 0.05),
                              ("mixed", 0.05, 0.05)):
                arms.append(("A", 0.30, "maf0.3", shape, dgp, a, c, big))
        # arm B: constant MAF 0.01
        for shape in ("normal", "gamma2"):
            for dgp, a, c in (("null", 0.0, 0.0), ("scale", 0.0, 0.05)):
                arms.append(("B", 0.01, "maf0.01", shape, dgp, a, c, small))
        # arm C: constant MAC ~ 922
        maf_c = TARGET_MAC / (2.0 * N)
        for dgp, a, c in (("null", 0.0, 0.0), ("scale", 0.0, 0.05)):
            arms.append(("C", maf_c, "mac922", "normal", dgp, a, c, small))
        for arm, maf, mtag, shape, dgp, a, c, reps in arms:
            nm = dgp + "_" + mtag + "_" + shape + "_" + tag
            cells.append(dict(name=nm, arm=arm, dgp=dgp, maf=round(maf, 8),
                              shape=shape, a=a, c=c, N=N, reps=reps,
                              n_blocks=1000))
    return cells


if __name__ == "__main__":
    cells = build_cells()
    with open(sys.argv[1], "w") as fh:
        json.dump(cells, fh, indent=1)
    tot = 0
    for i, c in enumerate(cells):
        mac = 2 * c["N"] * c["maf"]
        tot += c["reps"]
        print("%3d %-34s arm%s N=%7d maf=%.6f MAC=%9.0f reps=%6d"
              % (i, c["name"], c["arm"], c["N"], c["maf"], mac, c["reps"]))
    print("cells:", len(cells), " total reps:", tot)
