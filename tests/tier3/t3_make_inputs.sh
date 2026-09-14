#!/bin/bash
# ---------------------------------------------------------------------------
# Tier 3 inputs: a SYNTHETIC BGEN at production sample size.
#
#   bash tests/tier3/t3_make_inputs.sh <dir> <n_samples> <n_variants>
#
# plink2 --dummy generates the genotypes, so nothing about this touches UK
# Biobank, ALSPAC or any other real data, and nothing produced here may or does
# enter the repository. A real BGEN is needed because the point of tier 3 is to
# measure the BGEN read phase, which no in-memory harness exercises.
#
# The phenotype is a standard normal independent of the genotypes: tier 3
# measures time and memory, not inference, and a null phenotype keeps the
# measurement free of any effect-size choice.
# ---------------------------------------------------------------------------
set -euo pipefail
dir="${1:?output directory}"
N="${2:-461320}"
M="${3:-5000}"
mkdir -p "$dir"

if [ -f "$dir/geno.bgen" ] && [ -f "$dir/pheno.txt" ]; then
  echo "tier 3 inputs already present in $dir"
  exit 0
fi

echo "== generating a synthetic BGEN: $N samples x $M variants"
command -v plink2 >/dev/null 2>&1 || {
  echo "plink2 not on PATH. On BluePebble: module load apps/plink2/2.00a6LM" >&2
  exit 1
}

t0=$(date +%s)
plink2 --dummy "$N" "$M" 0 dosage-freq=0.0 \
       --out "$dir/geno" --export bgen-1.2 bits=8 --threads "${GMOQ_THREADS:-4}" \
  >"$dir/plink2.log" 2>&1
echo "   plink2 took $(( $(date +%s) - t0 ))s"
ls -lh "$dir/geno.bgen" "$dir/geno.sample"

# Stage 1 keys the phenotype and covariate tables on FID, and takes sample IDs
# from the BGEN's own .sample file.
python - "$dir" <<'PY'
import os, sys
import numpy as np
d = sys.argv[1]
with open(os.path.join(d, "geno.sample")) as fh:
    lines = [l.split() for l in fh.read().splitlines() if l.strip()]
# .sample has two header rows (names, then types)
ids = [r[0] for r in lines[2:]]
rng = np.random.default_rng(31415)
n = len(ids)
y = rng.normal(size=n)
C = rng.normal(size=(n, 10))
with open(os.path.join(d, "pheno.txt"), "w") as fh:
    fh.write("FID\tIID\tsimpheno\n")
    for i, s in enumerate(ids):
        fh.write(f"{s}\t{s}\t{y[i]:.6f}\n")
with open(os.path.join(d, "covar.txt"), "w") as fh:
    fh.write("FID\tIID\t" + "\t".join(f"PC{j+1}" for j in range(C.shape[1])) + "\n")
    for i, s in enumerate(ids):
        fh.write(f"{s}\t{s}\t" + "\t".join(f"{v:.6f}" for v in C[i]) + "\n")
print(f"wrote pheno.txt and covar.txt for {n} samples")
PY
echo "tier 3 inputs ready in $dir"
