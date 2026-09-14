#!/usr/bin/env python
"""Coverage / Type-I / bias summaries for one collated cell."""
from __future__ import annotations
import argparse, json, os
import numpy as np
from scipy import stats as sps

KEYS = ["jackknife_ols", "jackknife_gls", "pooled_ols", "pooled_gls"]
TAGS = ["full", "tiedrop", "qdedup"]


def mc(p, n):
    return float(np.sqrt(max(p * (1 - p), 0) / max(n, 1)))


def load(d):
    meta = json.load(open(os.path.join(d, "meta.json")))
    T = len(meta["taus"]); nel = T * (T + 1) // 2; P = meta["n_snp"]
    beta = np.fromfile(os.path.join(d, "beta.bin")).reshape(P, T)
    se = np.fromfile(os.path.join(d, "se.bin")).reshape(P, T)
    snpm = np.fromfile(os.path.join(d, "meta_snp.bin")).reshape(P, 4)
    diag = dict(np.load(os.path.join(d, "diag.npz")))
    s2 = {}
    for tg in TAGS:
        for k in KEYS:
            f = os.path.join(d, f"s2_{tg}_{k}.bin")
            if os.path.exists(f):
                s2[f"{tg}|{k}"] = np.fromfile(f).reshape(P, 10)
    cond = {}
    for tg in TAGS:
        f = os.path.join(d, f"cond_{tg}.json")
        if os.path.exists(f):
            cond[tg] = json.load(open(f))
    return meta, beta, se, snpm, diag, s2, cond


def summarise(d):
    meta, beta, se, snpm, diag, s2, cond = load(d)
    taus = np.array(meta["taus"]); T = len(taus)
    b_lat = np.array(meta["b_latent"]); b_ship = np.array(meta["b_ship"])
    b_obs = np.array([np.nan if v is None else v for v in meta["b_obs"]], dtype=float)
    f_ship = np.array(meta["f_ship"])
    f_true = np.array([np.nan if v is None else v for v in meta["f_obs_true"]], dtype=float)
    atom = np.array(meta["atom"])
    causal = snpm[:, 1] == 0
    nullm = ~causal
    nu = meta["blocks"] - 1  # per-tau reference (jackknife)
    out = {"cell": meta["cell"], "s0": meta["s0"], "N": meta["N"],
           "n_rep": meta["n_rep"], "M": meta["M"],
           "taus": taus.tolist(), "atom": atom.tolist(),
           "f_ship": f_ship.tolist(), "f_true": f_true.tolist(),
           "b_latent": b_lat.tolist(), "b_ship": b_ship.tolist(),
           "b_obs": b_obs.tolist(),
           "bw_rule": meta.get("bw_rule", "silverman"),
           "h_realised": meta.get("h_realised"),
           "h_pop": meta.get("h_pop"), "sd_obs": meta.get("sd_obs")}

    # ---------------- Stage 1, causal SNP, per tau
    bc, sc = beta[causal], se[causal]
    n = bc.shape[0]
    emp_sd = bc.std(axis=0, ddof=1)
    out["s1_causal"] = dict(
        n=n,
        mean_b=bc.mean(axis=0).tolist(),
        emp_sd=emp_sd.tolist(),
        mean_se=sc.mean(axis=0).tolist(),
        se_ratio=(sc.mean(axis=0) / emp_sd).tolist(),
        bias_vs_latent=((bc.mean(axis=0) - b_lat) / emp_sd).tolist(),
        bias_vs_ship=((bc.mean(axis=0) - b_ship) / emp_sd).tolist(),
        cov95_latent=[float(np.mean(np.abs(bc[:, t] - b_lat[t]) <= sps.t.ppf(.975, nu) * sc[:, t])) for t in range(T)],
        cov99_latent=[float(np.mean(np.abs(bc[:, t] - b_lat[t]) <= sps.t.ppf(.995, nu) * sc[:, t])) for t in range(T)],
        cov95_ship=[float(np.mean(np.abs(bc[:, t] - b_ship[t]) <= sps.t.ppf(.975, nu) * sc[:, t])) for t in range(T)],
        z_power=[float(np.mean(np.abs(bc[:, t] / sc[:, t]) > sps.t.ppf(.975, nu))) for t in range(T)],
    )
    # ---------------- Stage 1, null SNPs, per tau (type I + SE calibration)
    for tag, sel in (("maf05", nullm & (snpm[:, 2] == 0.05)),
                     ("maf20", nullm & (snpm[:, 2] == 0.20)),
                     ("maf50", nullm & (snpm[:, 2] == 0.50)),
                     ("all", nullm)):
        bn, sn = beta[sel], se[sel]
        if bn.shape[0] == 0:
            continue
        z = bn / sn
        nn = bn.shape[0]
        d_ = dict(n=nn,
                  emp_sd=bn.std(axis=0, ddof=1).tolist(),
                  mean_se=sn.mean(axis=0).tolist(),
                  se_ratio=(sn.mean(axis=0) / bn.std(axis=0, ddof=1)).tolist(),
                  mean_b=bn.mean(axis=0).tolist())
        for al, lab in ((0.05, "t1_05"), (0.01, "t1_01"), (1e-3, "t1_1e3"), (1e-4, "t1_1e4")):
            crit = sps.t.ppf(1 - al / 2, nu)
            r = [float(np.mean(np.abs(z[:, t]) > crit)) for t in range(T)]
            d_[lab] = r
            d_[lab + "_mc"] = [mc(x, nn) for x in r]
        # tail behaviour of |z| against the t_nu reference
        az = np.abs(z)
        d_["z_q"] = {}
        for lev in (0.99, 0.999, 0.9999, 0.99999):
            emp = np.quantile(az, lev, axis=0)
            ref = sps.t.ppf(0.5 + lev / 2, nu)
            d_["z_q"][str(lev)] = dict(emp=emp.tolist(), ref=float(ref),
                                       ratio=(emp / ref).tolist())
        d_["z_max"] = az.max(axis=0).tolist()
        out[f"s1_null_{tag}"] = d_

    # ---------------- Stage 2
    th = np.array(meta["theta_true"])
    s2out = {}
    for k, A in s2.items():
        p1, p2 = A[:, 0], A[:, 1]
        e1, e2 = A[:, 2], A[:, 3]
        pv1, pv2 = A[:, 4], A[:, 5]
        df1, df2 = A[:, 6], A[:, 7]
        Q, Qp = A[:, 8], A[:, 9]
        r = {}
        for tag, sel in (("causal", causal), ("null", nullm)):
            nn = int(sel.sum())
            tt = th if tag == "causal" else np.zeros(2)
            c95 = [float(np.mean(np.abs(p1[sel] - tt[0]) <= sps.t.ppf(.975, df1[sel]) * e1[sel])),
                   float(np.mean(np.abs(p2[sel] - tt[1]) <= sps.t.ppf(.975, df2[sel]) * e2[sel]))]
            c99 = [float(np.mean(np.abs(p1[sel] - tt[0]) <= sps.t.ppf(.995, df1[sel]) * e1[sel])),
                   float(np.mean(np.abs(p2[sel] - tt[1]) <= sps.t.ppf(.995, df2[sel]) * e2[sel]))]
            r[tag] = dict(
                n=nn,
                mean=[float(p1[sel].mean()), float(p2[sel].mean())],
                emp_sd=[float(p1[sel].std(ddof=1)), float(p2[sel].std(ddof=1))],
                mean_se=[float(e1[sel].mean()), float(e2[sel].mean())],
                se_ratio=[float(e1[sel].mean() / p1[sel].std(ddof=1)),
                          float(e2[sel].mean() / p2[sel].std(ddof=1))],
                bias_sd=[float((p1[sel].mean() - tt[0]) / p1[sel].std(ddof=1)),
                         float((p2[sel].mean() - tt[1]) / p2[sel].std(ddof=1))],
                cov95=c95, cov95_mc=[mc(c95[0], nn), mc(c95[1], nn)],
                cov99=c99, cov99_mc=[mc(c99[0], nn), mc(c99[1], nn)],
                rej05=[float(np.mean(pv1[sel] < 0.05)), float(np.mean(pv2[sel] < 0.05))],
                rej_1e3=[float(np.mean(pv1[sel] < 1e-3)), float(np.mean(pv2[sel] < 1e-3))],
                rej_1e4=[float(np.mean(pv1[sel] < 1e-4)), float(np.mean(pv2[sel] < 1e-4))],
                rej_1e6=[float(np.mean(pv1[sel] < 1e-6)), float(np.mean(pv2[sel] < 1e-6))],
                Q_rej05=float(np.mean(Qp[sel] < 0.05)),
                Q_rej_1e3=float(np.mean(Qp[sel] < 1e-3)),
                Q_median=float(np.nanmedian(Q[sel])),
                p_q={str(l): [float(np.nanquantile(pv1[sel], 1 - l)),
                              float(np.nanquantile(pv2[sel], 1 - l))]
                     for l in (0.99, 0.999, 0.9999)},
                p_min=[float(np.nanmin(pv1[sel])), float(np.nanmin(pv2[sel]))],
            )
        s2out[k] = r
    out["s2"] = s2out
    out["cond"] = cond

    # ---------------- observable diagnostics (sample, no truth)
    out["diag"] = dict(
        tie=np.asarray(diag["tie"]).mean(axis=0).tolist(),
        fgap=np.abs(np.asarray(diag["fgap"])).mean(axis=0).tolist(),
        h_over_spacing=(np.asarray(diag["hbw"])[:, None] /
                        np.asarray(diag["spacing"])).mean(axis=0).tolist(),
        ndist_h=(np.asarray(diag["ndist_h"]).mean(axis=0).tolist() if "ndist_h" in diag else None),
        fhat_mean=np.asarray(diag["fhat"]).mean(axis=0).tolist(),
        fhat_sd=np.asarray(diag["fhat"]).std(axis=0, ddof=1).tolist(),
        q_tau_sd=np.asarray(diag["q_tau"]).std(axis=0, ddof=1).tolist(),
        q_tau_mean=np.asarray(diag["q_tau"]).mean(axis=0).tolist(),
        dup_qtau=float(np.mean([len(np.asarray(r)) - len(np.unique(np.asarray(r)))
                                for r in np.asarray(diag["q_tau"])])),
        ndistinct=float(np.mean(diag["ndistinct"])),
        maxatom=float(np.mean(diag["maxatom"])),
        hbw=float(np.mean(diag["hbw"])),
    )
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    json.dump(summarise(a.dir), open(a.out, "w"))
    print("wrote", a.out)
