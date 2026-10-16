#!/usr/bin/env bash
# chr22 biallelic SNPs, MAF > 0.05 -> 5 PCs from LD-pruned SNPs -> 2,000 random SNPs as BGEN.
set -euo pipefail
cd "$(dirname "$0")"; source config.sh; mkdir -p "$WORK"
P="$PLINK2 --threads $THREADS --memory 4000 --seed 2026"
$P --vcf "$VCF" --snps-only just-acgt --max-alleles 2 --maf 0.05 \
   --set-all-var-ids 'chr@:#:$r:$a' --rm-dup exclude-all \
   --make-pgen --out "$WORK/common"
$P --pfile "$WORK/common" --indep-pairwise 500kb 0.2 --out "$WORK/prune"
$P --pfile "$WORK/common" --extract "$WORK/prune.prune.in" --pca 5 --out "$WORK/pca"
$P --pfile "$WORK/common" --thin-count 2000 --write-snplist --out "$WORK/snps"
# BGEN v1.2, 8-bit; ref-first => allele 1 = REF = the effect allele Stage 1 reports.
# Stage 1 needs no .bgi index for a whole-file scan.
$P --pfile "$WORK/common" --extract "$WORK/snps.snplist" \
   --export bgen-1.2 bits=8 ref-first --out "$WORK/geno"
# Dosages for the simulation (02_simulate.R aligns them to allele 1).
$P --pfile "$WORK/common" --extract "$WORK/snps.snplist" --export A --out "$WORK/geno"
rm -f "$WORK"/common.{pgen,pvar,psam} "$WORK/prune.prune.out"
