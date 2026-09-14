#!/usr/bin/env python
"""Aggregate the per-cell summaries into the report tables."""
import glob, json, os, sys
import numpy as np
from scipy import stats as sps

SUMM = "/user/work/fh6520/kde_lumpy/summ"
GW = 5e-8

ORDER = ["normal_smooth", "lognorm_smooth", "bimodal_smooth", "bimodal_uneven",
         "round_0p10", "round_0p25", "round_0p50", "round_1p00",
         "floor_05", "floor_20", "floor_40", "ceil_20",
         "count_8", "count_4", "bimodal_round"]


def load():
    out = {}
    for f in sorted(glob.glob(os.path.join(SUMM, "*.json"))):
        nm = os.path.basename(f)[:-5]
        if nm == "TEST":
            continue
        try:
            out[nm] = json.load(open(f))
        except Exception as e:
            print("SKIP", nm, e, file=sys.stderr)
    return out


def mc(p, n):
    return float(np.sqrt(max(p * (1 - p), 0.0) / max(n, 1)))


def gw_factor(sd_z):
    """Type-I at 5e-8 implied by a z that is N(0, sd_z), as a multiple of nominal."""
    if not np.isfinite(sd_z) or sd_z <= 0:
        return np.nan
    zc = sps.norm.ppf(1 - GW / 2)
    return float(2 * sps.norm.sf(zc / sd_z) / GW)


def key(nm):
    c = nm.split("__")[-1]
    return (ORDER.index(c) if c in ORDER else 99, nm)


def sec(t):
    print()
    print("=" * 118)
    print(t)
    print("=" * 118)


def main():
    S = load()
    runA = {k: v for k, v in S.items() if k.startswith("runA__")}
    print("cells summarised: %d  (runA %d)" % (len(S), len(runA)))

    taus = np.array(S[list(S)[0]]["taus"])
    T = len(taus)

    # ---------------------------------------------------------------- Stage 2
    for tag in ("full", "tiedrop"):
        sec("STAGE 2 -- NULL SNPs, %s tau set.  Coverage of theta=(shift,scale), "
            "Type-I at 0.05, and genome-wide inflation from sd(z)=1/se_ratio" % tag.upper())
        print("%-16s %-18s %8s | %-24s | %-24s | %9s" %
              ("cell", "estimator", "n_null", "cov95 [shift, scale]",
               "typeI@.05 [shift,scale]", "GW x nom"))
        for nm in sorted(runA, key=key):
            d = S[nm]
            cell = d["cell"]
            for est in ("jackknife_ols", "jackknife_gls", "pooled_ols", "pooled_gls"):
                k = "%s|%s" % (tag, est)
                if k not in d["s2"]:
                    continue
                r = d["s2"][k]["null"]
                n = r["n"]
                c95 = r["cov95"]; t1 = r["rej05"]
                sdz = [1.0 / x if x else np.nan for x in r["se_ratio"]]
                gwf = [gw_factor(z) for z in sdz]
                print("%-16s %-18s %8d | %6.4f(%.4f) %6.4f(%.4f) | "
                      "%6.4f(%.4f) %6.4f(%.4f) | %4.2f %4.2f" %
                      (cell, est, n, c95[0], mc(c95[0], n), c95[1], mc(c95[1], n),
                       t1[0], mc(t1[0], n), t1[1], mc(t1[1], n), gwf[0], gwf[1]))

    sec("STAGE 2 -- CAUSAL SNP, FULL tau set. Coverage of the TRUE theta=(0, s0), "
        "bias in empirical sd, SE ratio")
    print("%-16s %-18s %8s | %-17s | %-17s | %-15s" %
          ("cell", "estimator", "n", "cov95[sh,sc]", "bias_sd[sh,sc]", "se_ratio[sh,sc]"))
    for nm in sorted(runA, key=key):
        d = S[nm]
        for est in ("jackknife_gls", "pooled_gls"):
            k = "full|%s" % est
            if k not in d["s2"]:
                continue
            r = d["s2"][k]["causal"]
            n = r["n"]
            print("%-16s %-18s %8d | %7.4f %7.4f | %7.3f %7.3f | %7.3f %7.3f" %
                  (d["cell"], est, n, r["cov95"][0], r["cov95"][1],
                   r["bias_sd"][0], r["bias_sd"][1],
                   r["se_ratio"][0], r["se_ratio"][1]))

    sec("STAGE 2 -- cov99 (null) and Q lack-of-fit rejection, FULL vs TIEDROP, GLS jackknife")
    print("%-16s | %-17s %-17s | %-13s %-13s" %
          ("cell", "cov99 full", "cov99 tiedrop", "Qrej05 full", "Qrej05 tied"))
    for nm in sorted(runA, key=key):
        d = S[nm]
        a = d["s2"].get("full|jackknife_gls", {}).get("null")
        b = d["s2"].get("tiedrop|jackknife_gls", {}).get("null")
        if a is None:
            continue
        bb = (["%7.4f" % x for x in b["cov99"]], "%.4f" % b["Q_rej05"]) if b else (["   --  "] * 2, "  --  ")
        print("%-16s | %7.4f %7.4f  %s %s | %11.4f %11s" %
              (d["cell"], a["cov99"][0], a["cov99"][1], bb[0][0], bb[0][1],
               a["Q_rej05"], bb[1]))

    # ---------------------------------------------------------------- Stage 1
    sec("STAGE 1 -- PER-TAU, NULL SNPs: SE ratio (mean SE / empirical sd) and Type-I at 0.05")
    print("tau".ljust(18) + "".join("%7.2f" % t for t in taus))
    for nm in sorted(runA, key=key):
        d = S[nm]
        r = d["s1_null_all"]
        print(("%-12s se/sd" % d["cell"][:12]).ljust(18) +
              "".join("%7.3f" % x for x in r["se_ratio"]))
    print()
    print("tau".ljust(18) + "".join("%7.2f" % t for t in taus))
    for nm in sorted(runA, key=key):
        d = S[nm]
        r = d["s1_null_all"]
        print(("%-12s t1.05" % d["cell"][:12]).ljust(18) +
              "".join("%7.3f" % x for x in r["t1_05"]))

    sec("STAGE 1 -- PER-TAU, CAUSAL SNP: coverage of b_ship (the shipped estimand) "
        "and of b_latent (what the user wants)")
    print("tau".ljust(20) + "".join("%7.2f" % t for t in taus))
    for nm in sorted(runA, key=key):
        d = S[nm]
        r = d["s1_causal"]
        print(("%-13s cov_ship" % d["cell"][:13]).ljust(20) +
              "".join("%7.3f" % x for x in r["cov95_ship"]))
    print()
    print("tau".ljust(20) + "".join("%7.2f" % t for t in taus))
    for nm in sorted(runA, key=key):
        d = S[nm]
        r = d["s1_causal"]
        print(("%-13s cov_lat" % d["cell"][:13]).ljust(20) +
              "".join("%7.3f" % x for x in r["cov95_latent"]))

    # ------------------------------------------------------ observable guards
    sec("OBSERVABLE DIAGNOSTICS (no truth needed): sample tie fraction at q_tau")
    print("tau".ljust(20) + "".join("%7.2f" % t for t in taus))
    for nm in sorted(runA, key=key):
        d = S[nm]
        print(("%-13s tie" % d["cell"][:13]).ljust(20) +
              "".join("%7.3f" % x for x in d["diag"]["tie"]))
    print()
    print("bandwidth-to-local-spacing ratio  h / spacing  (>1 means the kernel spans distinct values)")
    print("tau".ljust(20) + "".join("%7.2f" % t for t in taus))
    for nm in sorted(runA, key=key):
        d = S[nm]
        v = d["diag"]["h_over_spacing"]
        print(("%-13s h/sp" % d["cell"][:13]).ljust(20) +
              "".join("%7.2f" % x for x in v))

    sec("GUARD CALIBRATION: does the observable tie fraction predict the damage?")
    print("%-16s %6s %6s %8s %9s %9s %9s" %
          ("cell", "maxtie", "ndist", "dup_qtau", "cov95_sh", "cov95_sc", "T1_scale"))
    for nm in sorted(runA, key=key):
        d = S[nm]
        g = d["s2"].get("full|jackknife_gls", {}).get("null")
        if g is None:
            continue
        print("%-16s %6.3f %6.0f %8.2f %9.4f %9.4f %9.4f" %
              (d["cell"], max(d["diag"]["tie"]), d["diag"]["ndistinct"],
               d["diag"]["dup_qtau"], g["cov95"][0], g["cov95"][1], g["rej05"][1]))

    # ------------------------------------------------- per-tau vs tie scatter
    sec("LOCALITY: per-tau null SE ratio against the per-tau tie fraction, pooled over runA cells")
    x, y = [], []
    for nm in runA:
        d = S[nm]
        x += list(d["diag"]["tie"])
        y += list(d["s1_null_all"]["se_ratio"])
    x = np.array(x); y = np.array(y)
    for lo, hi in [(0, .001), (.001, .01), (.01, .05), (.05, .10),
                   (.10, .20), (.20, .50)]:
        s = (x >= lo) & (x < hi)
        if s.sum():
            print("tie in [%.3f,%.3f): n=%4d  se_ratio mean %.4f  min %.4f  max %.4f" %
                  (lo, hi, s.sum(), y[s].mean(), y[s].min(), y[s].max()))

    # ----------------------------------------------------- bandwidth variants
    sec("BANDWIDTH RULES (runB, 500 reps): Stage 2 null coverage and Type-I, jackknife GLS, FULL")
    print("%-18s %-11s %8s | %-17s | %-17s" %
          ("cell", "bw", "n", "cov95[sh,sc]", "rej05[sh,sc]"))
    for nm in sorted([k for k in S if k.startswith("runB__")], key=key):
        d = S[nm]
        r = d["s2"].get("full|jackknife_gls", {}).get("null")
        if r is None:
            continue
        print("%-18s %-11s %8d | %7.4f %7.4f | %7.4f %7.4f" %
              (d["cell"], d["bw_rule"], r["n"], r["cov95"][0], r["cov95"][1],
               r["rej05"][0], r["rej05"][1]))

    sec("SAMPLE SIZE (runC, N=80000) and EFFECT SIZE (runD): jackknife GLS, FULL, null")
    print("%-28s %8s | %-17s | %-17s" % ("run/cell", "n", "cov95[sh,sc]", "rej05[sh,sc]"))
    for nm in sorted([k for k in S if k.startswith(("runC__", "runD__"))], key=key):
        d = S[nm]
        r = d["s2"].get("full|jackknife_gls", {}).get("null")
        if r is None:
            continue
        print("%-28s %8d | %7.4f %7.4f | %7.4f %7.4f" %
              (nm.replace("run", "").replace("__", "/"), r["n"],
               r["cov95"][0], r["cov95"][1], r["rej05"][0], r["rej05"][1]))

    sec("COVARIANCE CONDITIONING (Stage 2 tau-covariance), full vs tiedrop")
    print("%-18s %-8s %12s %10s | %-8s %12s %10s" %
          ("cell", "tag", "cond_med", "dup_pairs", "tag", "cond_med", "dup_pairs"))
    for nm in sorted(runA, key=key):
        d = S[nm]
        a = d["cond"].get("full"); b = d["cond"].get("tiedrop")
        if not a:
            continue
        bs = ("%12.3g %10.1f" % (b["cond_median"], b["dup_pairs_median"])) if b else "%12s %10s" % ("--", "--")
        print("%-18s %-8s %12.3g %10.1f | %-8s %s" %
              (d["cell"], "full", a["cond_median"], a["dup_pairs_median"],
               "tiedrop", bs))


if __name__ == "__main__":
    main()
