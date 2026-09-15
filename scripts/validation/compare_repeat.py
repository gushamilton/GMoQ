#!/usr/bin/env python3
"""Compare decompressed Stage 1 outputs from two repeated runs."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with gzip.open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--first", type=Path, required=True)
    ap.add_argument("--second", type=Path, required=True)
    args = ap.parse_args()
    first, second = digest(args.first), digest(args.second)
    print(json.dumps({"first": first, "second": second, "equal": first == second}, indent=2))
    if first != second:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
