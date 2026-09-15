#!/usr/bin/env python3
"""Record the matched LOCO-adjusted log_bmi baseline for Stage 2 weights."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bgen", type=Path, required=True)
    ap.add_argument("--sample", type=Path, required=True)
    ap.add_argument("--snps", type=Path, required=True)
    ap.add_argument("--pheno", type=Path, required=True)
    ap.add_argument("--covar", type=Path, required=True)
    ap.add_argument("--keep", type=Path, required=True)
    ap.add_argument("--loco", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    from fungwas_stage1 import cli

    covars = ["age", "sex", *[f"PC{i}" for i in range(1, 11)]]
    sample_ids, sample_aliases, y_raw, _, _ = cli.load_analysis_inputs(
        str(args.pheno), ["log_bmi"], str(args.covar),
        keep_file=str(args.keep), covar_cols=covars,
    )
    loco_files = cli.resolve_loco_files(str(args.loco), ["log_bmi"])
    preds, loco_mask = cli.load_loco_predictions(
        str(args.loco), ["log_bmi"], sample_aliases,
    )
    y = np.asarray(y_raw["log_bmi"], dtype=np.float64)[loco_mask]
    sample_ids = [sid for sid, keep in zip(sample_ids, loco_mask) if keep]
    scale = cli.resolve_loco_scale("regression", y, preds["log_bmi"], "log_bmi")
    y_loco = y - scale * preds["log_bmi"]["22"]
    snps = [line.split()[0] for line in args.snps.read_text().splitlines() if line.split()]
    _, bgen_samples, handle = cli.get_variants_bgen(str(args.bgen), str(args.sample), snps)
    valid_mask, _ = cli.match_samples(bgen_samples, sample_ids)
    handle.close()
    y_match = y_loco[valid_mask]
    taus = np.linspace(0.02, 0.98, 45)
    result = {
        "phenotype": "log_bmi",
        "definition": "log(bmi), LOCO chr22 subtraction, then BGEN sample match",
        "loco_file": loco_files["log_bmi"],
        "loco_scale": float(scale),
        "n_loco_rows": int(y.size),
        "n_matched_rows": int(y_match.size),
        "mean": float(np.mean(y_match)),
        "sd": float(np.std(y_match, ddof=1)),
        "taus": taus.tolist(),
        "q_tau": np.quantile(y_match, taus, method="median_unbiased").tolist(),
    }
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
