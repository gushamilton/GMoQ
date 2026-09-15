#!/usr/bin/env python3
"""Assemble tier 3 into one JSON: phases, CPU time, peak RSS, hardware.

Stage 1's own log already instruments the phases; /usr/bin/time -v supplies
wall clock, user and system CPU time and maximum resident set size. The BGEN
read phase is the part of the batch loop that is neither compute nor write,
which is how it is obtained without modifying the shipped code.
"""
import argparse, json, os, re


def parse_time_v(path):
    out = {}
    if not os.path.exists(path):
        return out
    for line in open(path):
        line = line.strip()
        if line.startswith("Elapsed (wall clock) time"):
            v = line.split(": ", 1)[1]
            parts = [float(x) for x in v.split(":")]
            sec = 0.0
            for p in parts:
                sec = sec * 60 + p
            out["wall_sec"] = sec
        elif line.startswith("User time (seconds)"):
            out["user_sec"] = float(line.split(":")[1])
        elif line.startswith("System time (seconds)"):
            out["system_sec"] = float(line.split(":")[1])
        elif line.startswith("Maximum resident set size"):
            out["max_rss_gb"] = float(line.split(":")[1]) / 2 ** 20
        elif line.startswith("Percent of CPU this job got"):
            out["cpu_percent"] = float(line.split(":")[1].strip().rstrip("%"))
    out["cpu_sec"] = out.get("user_sec", 0.0) + out.get("system_sec", 0.0)
    return out


def parse_stage1_log(path):
    txt = open(path).read() if os.path.exists(path) else ""
    def g(pat, default=None):
        m = re.search(pat, txt)
        return float(m.group(1)) if m else default
    d = {
        "pheno_load_sec": g(r"Phenotype load time: ([0-9.]+)s"),
        "rif_build_sec": g(r"RIF build for [^:]+: ([0-9.]+)s"),
        "sigma0_sec": g(r"Sigma_0 compute time: ([0-9.]+)s"),
        "batch_compute_sec": g(r"Batch compute time: ([0-9.]+)s"),
        "batch_write_sec": g(r"write time: ([0-9.]+)s"),
        "batch_total_sec": g(r"batch total: ([0-9.]+)s"),
        "total_sec": g(r"Done! \d+ SNPs in ([0-9.]+)s"),
        "n_snps": g(r"Done! (\d+) SNPs"),
        "snps_per_sec": g(r"\(([0-9.]+) SNPs/s\)"),
    }
    if d["batch_total_sec"] is not None:
        d["bgen_read_sec"] = max(
            d["batch_total_sec"] - (d["batch_compute_sec"] or 0.0)
            - (d["batch_write_sec"] or 0.0), 0.0)
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    hw = {}
    hp = os.path.join(a.dir, "hardware.txt")
    if os.path.exists(hp):
        for line in open(hp):
            if "=" in line:
                k, v = line.strip().split("=", 1)
                hw[k] = v

    s1 = parse_stage1_log(os.path.join(a.dir, "stage1.log"))
    s1.update({f"proc_{k}": v for k, v in parse_time_v(os.path.join(a.dir, "stage1.time")).items()})
    s2t = parse_time_v(os.path.join(a.dir, "stage2.time"))
    s2 = {}
    sp = os.path.join(a.dir, "stage2.json")
    if os.path.exists(sp):
        s2 = json.load(open(sp))
    s2.update({f"proc_{k}": v for k, v in s2t.items()})

    n = s1.get("n_snps") or float(hw.get("n_variants", 0) or 0)
    per_variant = {}
    if n:
        per_variant["stage1_cpu_sec_per_variant"] = s1.get("proc_cpu_sec", 0.0) / n
        per_variant["stage1_wall_sec_per_variant"] = s1.get("proc_wall_sec", 0.0) / n
        per_variant["stage2_wall_sec_per_variant"] = s2.get("proc_wall_sec", 0.0) / n
        for gw in (1_200_000,):
            per_variant[f"stage1_cpu_hours_at_{gw}"] = \
                s1.get("proc_cpu_sec", 0.0) / n * gw / 3600
            per_variant[f"stage2_wall_hours_at_{gw}"] = \
                s2.get("proc_wall_sec", 0.0) / n * gw / 3600

    payload = dict(hardware=hw, stage1=s1, stage2=s2, extrapolation=per_variant,
                   extrapolation_note=(
                       "Per-variant costs multiplied out to 1.2M variants. Stage 1 is "
                       "linear in the variant count by construction (a fixed RIF, then "
                       "a per-variant kernel), so this is a defensible extrapolation "
                       "for TIME. It is NOT a measurement, and it says nothing about "
                       "Stage 2 memory, which grows with the stored covariance and is "
                       "measured directly by the scaling rows instead."))
    json.dump(payload, open(a.out, "w"), indent=2)
    print(json.dumps(payload, indent=2)[:4000])


if __name__ == "__main__":
    main()
