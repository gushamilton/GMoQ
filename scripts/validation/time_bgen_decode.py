#!/usr/bin/env python3
"""Time native dosage extraction against the checked per-record BGEN path."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bgen", type=Path, required=True)
    ap.add_argument("--sample", type=Path, required=True)
    ap.add_argument("--snps", type=Path, required=True)
    args = ap.parse_args()

    import numpy as np
    from fungwas_stage1 import cli

    snps = [line.split()[0] for line in args.snps.read_text().splitlines() if line.split()]

    def load():
        gen, _, handle = cli.get_variants_bgen(str(args.bgen), str(args.sample), snps)
        variants = list(gen)
        return variants, handle

    variants, handle = load()
    t0 = time.perf_counter()
    native = [2.0 - np.asarray(v.alt_dosage, dtype=np.float64) for v in variants]
    native_s = time.perf_counter() - t0
    handle.close()

    variants, handle = load()
    t0 = time.perf_counter()
    checked = [cli.extract_dosage(v) for v in variants]
    checked_s = time.perf_counter() - t0
    handle.close()

    max_diff = max(float(np.max(np.abs(a - b))) for a, b in zip(native, checked))
    print(json.dumps({
        "n_variants": len(variants),
        "native_alt_to_allele0_seconds": native_s,
        "checked_extract_dosage_seconds": checked_s,
        "checked_minus_native_seconds": checked_s - native_s,
        "overhead_percent": 100.0 * (checked_s / native_s - 1.0),
        "max_abs_dosage_difference": max_diff,
    }, indent=2))


if __name__ == "__main__":
    main()
