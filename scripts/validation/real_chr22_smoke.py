#!/usr/bin/env python3
"""Run a small real-chr22 Stage 1 smoke with a known dosage-linked effect.

The script selects a bounded set of biallelic records, makes one phenotype from
the first selected record's real BGEN dosage, and invokes the shipped Stage 1
CLI on those records.  It writes only the selected IDs, phenotype, and compact
truth metadata; the full scan outputs stay in the requested work directory.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np


def read_snp_ids(path: Path) -> list[str]:
    ids = []
    seen = set()
    with path.open() as fh:
        for line in fh:
            fields = line.split()
            if not fields:
                continue
            rsid = fields[0]
            if rsid not in seen:
                ids.append(rsid)
                seen.add(rsid)
    return ids


def select_real_dosages(bgen_path: Path, sample_path: Path, candidates: list[str], n_variants: int):
    from bgen import BgenReader

    reader = BgenReader(str(bgen_path), sample_path=str(sample_path), delay_parsing=True)
    sample_ids = list(reader.samples)
    wanted = set(candidates)
    selected = []
    found = set()
    for variant in reader:
        if variant.rsid not in wanted or variant.rsid in found:
            continue
        if len(variant.alleles) != 2:
            continue
        dosage = 2.0 - np.asarray(variant.alt_dosage, dtype=np.float64)
        called = np.isfinite(dosage)
        if called.sum() < max(100, int(0.99 * len(dosage))):
            continue
        mean = float(np.nanmean(dosage))
        maf = min(mean / 2.0, 1.0 - mean / 2.0)
        if not 0.05 <= maf <= 0.50:
            continue
        selected.append({
            "rsid": variant.rsid,
            "alleles": list(variant.alleles),
            "chrom": str(variant.chrom),
            "bp": int(variant.pos),
            "dosage": dosage,
            "maf": maf,
        })
        found.add(variant.rsid)
        if len(selected) >= n_variants:
            break
    reader.close()
    if len(selected) < n_variants:
        raise RuntimeError(
            f"found only {len(selected)} usable biallelic variants from the requested list; "
            f"need {n_variants}"
        )
    return sample_ids, selected


def write_phenotype(path: Path, sample_ids: list[str], dosage: np.ndarray, effect: float, seed: int):
    rng = np.random.default_rng(seed)
    imputed = np.where(np.isfinite(dosage), dosage, np.nanmean(dosage))
    centered = imputed - np.mean(imputed)
    y = effect * centered + rng.normal(0.0, 1.0, size=len(imputed))
    with path.open("w") as fh:
        fh.write("FID IID trait\n")
        for sample_id, value in zip(sample_ids, y):
            fh.write(f"{sample_id} {sample_id} {value:.12g}\n")
    return y


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bgen", type=Path, required=True)
    ap.add_argument("--sample", type=Path, required=True)
    ap.add_argument("--snps", type=Path, required=True, help="candidate SNP IDs, one per line")
    ap.add_argument("--out", type=Path, required=True, help="smoke output directory")
    ap.add_argument("--n-variants", type=int, default=32)
    ap.add_argument("--effect", type=float, default=0.25)
    ap.add_argument("--seed", type=int, default=20260915)
    ap.add_argument("--blocks", type=int, default=65)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--batch-size", type=int, default=16)
    args = ap.parse_args()
    if args.n_variants < 2:
        ap.error("--n-variants must be at least 2")
    if args.blocks <= 45:
        ap.error("--blocks must exceed the shipped 45-tau grid")

    args.out.mkdir(parents=True, exist_ok=True)
    snp_ids = read_snp_ids(args.snps)
    sample_ids, selected = select_real_dosages(args.bgen, args.sample, snp_ids, args.n_variants)
    causal = selected[0]
    pheno_path = args.out / "phenotype.txt"
    write_phenotype(pheno_path, sample_ids, causal["dosage"], args.effect, args.seed)
    selected_path = args.out / "selected_snps.txt"
    selected_path.write_text("".join(f"{item['rsid']}\n" for item in selected))

    truth = {
        "causal_snp": causal["rsid"],
        "effect_per_allele": args.effect,
        "causal_effect_alleles": causal["alleles"],
        "causal_maf": causal["maf"],
        "n_samples": len(sample_ids),
        "n_variants": len(selected),
        "seed": args.seed,
        "selected_snps": [item["rsid"] for item in selected],
        "expected_beta_mu": args.effect,
        "note": "Non-causal rows are descriptive; LD means they are not asserted null.",
    }
    truth_path = args.out / "truth.json"
    truth_path.write_text(json.dumps(truth, indent=2) + "\n")

    stage1_prefix = args.out / "stage1" / "chr22_smoke"
    stage1_prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, "-m", "fungwas_stage1.cli",
        "--bgen", str(args.bgen), "--sample", str(args.sample),
        "--pheno", str(pheno_path), "--pheno-col", "trait",
        "--snps", str(selected_path), "--blocks", str(args.blocks),
        "--batch-size", str(args.batch_size), "--threads", str(args.threads),
        "--cov-dtype", "float32", "--out", str(stage1_prefix),
    ]
    print("running:", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)
    print(json.dumps({"stage1_prefix": str(stage1_prefix), **truth}, indent=2), flush=True)


if __name__ == "__main__":
    main()
