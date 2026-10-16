#!/usr/bin/env bash
# Download 1000 Genomes phase 3 chr22 (public, ~200 MB) into $GMOQ_1KG_DIR/data.
set -euo pipefail
cd "$(dirname "$0")"; source config.sh
mkdir -p "$(dirname "$VCF")"
B=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502
F=ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
[ -s "$VCF" ] || curl -fL -o "$VCF" "$B/$F"
[ -s "$VCF.tbi" ] || curl -fL -o "$VCF.tbi" "$B/$F.tbi"
echo "Downloaded to $(dirname "$VCF")"
