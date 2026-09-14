import json, sys

# N = 30000 throughout: the same scale gls_coverage used, and deliberately a
# HARDER test of the second-order expansion than a real scan, because reaching a
# given Z_mu at N = 30k needs a much larger delta/sigma than at N = 461k.
N, REPS = 30000, 20000

C = []
def add(name, maf, shape, delta, gamma):
    C.append(dict(name=name, maf=maf, shape=shape, delta=delta, gamma=gamma,
                  N=N, reps=REPS))

# --- null: the correction must not break a clean variant -------------------
add("null_maf30_normal",  0.30, "normal", 0.00, 0.00)
add("null_maf01_normal",  0.01, "normal", 0.00, 0.00)
# --- pure location, zero true dispersion: the headline ---------------------
add("loc10_maf30_normal", 0.30, "normal", 0.10, 0.00)
add("loc25_maf30_normal", 0.30, "normal", 0.25, 0.00)
add("loc15_maf01_normal", 0.01, "normal", 0.15, 0.00)
add("loc35_maf01_normal", 0.01, "normal", 0.35, 0.00)
add("loc35_maf50_normal", 0.50, "normal", 0.35, 0.00)   # control: exact no-op
# --- genuine dispersion, no location: power must not degrade ---------------
add("disp_maf30_normal",  0.30, "normal", 0.00, 0.05)
add("disp_maf01_normal",  0.01, "normal", 0.00, 0.12)
# --- mixed -----------------------------------------------------------------
add("mix_maf30_normal",   0.30, "normal", 0.20, 0.05)
add("mix_maf01_normal",   0.01, "normal", 0.25, 0.10)
# --- non-normal baselines --------------------------------------------------
add("loc25_maf30_skew",   0.30, "skew",   0.25, 0.00)
add("loc35_maf01_skew",   0.01, "skew",   0.35, 0.00)
add("loc25_maf30_heavy",  0.30, "heavy",  0.25, 0.00)
add("loc35_maf01_heavy",  0.01, "heavy",  0.35, 0.00)
add("mix_maf01_skew",     0.01, "skew",   0.25, 0.10)

if __name__ == "__main__":
    json.dump(C, open(sys.argv[1], "w"), indent=1)
    for i, c in enumerate(C):
        print(i, c["name"], c["N"], c["reps"])
