#!/usr/bin/env python
"""Build the headline tables from the exact sweep."""
import glob, csv, sys
import numpy as np

rows = []
for f in glob.glob("/user/work/fh6520/wt_dispersion/sim/out/exact/*.tsv"):
    for r in csv.DictReader(open(f), delimiter="\t"):
        rows.append(r)
def fl(r, k):
    try: return float(r[k])
    except Exception: return float("nan")

N0 = 461320
sel = [r for r in rows if int(float(r["N"])) == N0 and r["shape"] == "normal"]

print(f"### FALSE DISPERSION RATE - pure location effect, zero true dispersion")
print(f"### normal outcome, N = {N0}.  Exact: no Monte Carlo error.\n")
hdr = ("maf", "Zloc", "gskew", "z_sig", "r_pool", "pool@0.05", "jack@0.05",
       "pool@5e-8", "jack@5e-8", "SEonly@5e-8", "infl@5e-8")
print("%-7s %5s %7s %8s %8s %10s %10s %11s %11s %11s %10s" % hdr)
for maf in ["0.005", "0.01", "0.05", "0.2", "0.3", "0.5"]:
    rs = sorted([r for r in sel if r["maf"] == maf and r["kind"] == "loc"],
                key=lambda r: float(r["target_z"]))
    for r in rs:
        b = "both"
        print("%-7s %5.0f %7.3f %8.4f %8.5f %10.4g %10.4g %11.4g %11.4g %11.4g %10.4g" % (
            maf, fl(r, "target_z"), fl(r, "geno_skew"),
            fl(r, f"{b}.pooled.ols.beta_sigma2.z_true"),
            fl(r, f"{b}.pooled.ols.beta_sigma2.se_ratio"),
            fl(r, f"{b}.pooled.ols.beta_sigma2.rate_5e-2"),
            fl(r, f"{b}.jack.ols.beta_sigma2.rate_5e-2"),
            fl(r, f"{b}.pooled.ols.beta_sigma2.rate_5e-8"),
            fl(r, f"{b}.jack.ols.beta_sigma2.rate_5e-8"),
            fl(r, f"{b}.pooled.ols.beta_sigma2.rate_se_only_5e-8"),
            fl(r, f"{b}.pooled.ols.beta_sigma2.rate_5e-8") / 5e-8))
    print()

print("\n### SE RATIO (reported/true) - direction of error")
print("### r < 1 => SE too small => ANTI-CONSERVATIVE\n")
print("%-7s %-6s %5s %10s %10s %10s %10s %10s" % (
    "maf", "kind", "Z", "mu_pool", "sig_pool", "tau_min", "tau_max", "worst_lo"))
for kind in ("loc", "disp"):
    for maf in ["0.005", "0.05", "0.3"]:
        rs = sorted([r for r in sel if r["maf"] == maf and r["kind"] == kind],
                    key=lambda r: float(r["target_z"]))
        for r in rs:
            b = "both"
            print("%-7s %-6s %5.0f %10.5f %10.5f %10.5f %10.5f %10.5f" % (
                maf, kind, fl(r, "target_z"),
                fl(r, f"{b}.pooled.ols.beta_mu.se_ratio"),
                fl(r, f"{b}.pooled.ols.beta_sigma2.se_ratio"),
                fl(r, "tau_se_ratio_min"), fl(r, "tau_se_ratio_max"),
                fl(r, "worst_se_ratio_lo")))
        print()

print("\n### TRUE DISPERSION EFFECT: Type-I/coverage of the LOCATION parameter")
print("### and SE ratio of the dispersion parameter itself\n")
print("%-7s %5s %10s %10s %11s %11s" % ("maf", "Zdisp", "r_sig_pool", "r_mu_pool",
                                        "muFP@5e-8", "muFP_jack"))
for maf in ["0.005", "0.05", "0.3", "0.5"]:
    rs = sorted([r for r in sel if r["maf"] == maf and r["kind"] == "disp"],
                key=lambda r: float(r["target_z"]))
    for r in rs:
        b = "both"
        print("%-7s %5.0f %10.5f %10.5f %11.4g %11.4g" % (
            maf, fl(r, "target_z"),
            fl(r, f"{b}.pooled.ols.beta_sigma2.se_ratio"),
            fl(r, f"{b}.pooled.ols.beta_mu.se_ratio"),
            fl(r, f"{b}.pooled.ols.beta_mu.rate_5e-8"),
            fl(r, f"{b}.jack.ols.beta_mu.rate_5e-8")))
    print()
