#!/usr/bin/env python3
"""
extract_dosages.py - Extract real dosages from BGEN file for specific variants and samples.
Iterates through BGEN to find requested variants.
"""
import sys
import argparse
import numpy as np
from pathlib import Path
import json
import logging

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(description='Extract real dosages from BGEN')
    parser.add_argument('--bgen', required=True, help='BGEN file path')
    parser.add_argument('--sample', required=True, help='Sample file path')
    parser.add_argument('--variant-ids', required=True, help='Comma-separated variant IDs or file path')
    parser.add_argument('--sample-ids', required=True, help='Comma-separated sample IDs or file path')
    parser.add_argument('--out', required=True, help='Output TSV file')
    parser.add_argument('--meta-out', help='Optional: output metadata JSON')
    parser.add_argument('--min-samples', type=int, default=100, help='Minimum samples required')
    parser.add_argument('--min-variants', type=int, default=5, help='Minimum variants required')
    return parser.parse_args()


def read_id_list(spec):
    """Parse ID spec: either file path or comma-separated."""
    spec = spec.strip()
    
    # Try as file first
    try:
        with open(spec, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]
            ids = []
            for line in lines:
                # Handle FID IID format - take first column
                parts = line.split()
                if parts:
                    ids.append(parts[0])
            return ids
    except:
        pass
    
    # Fall back to comma/tab separated
    if ',' in spec or '\t' in spec:
        return [s.strip() for s in spec.replace('\t', ',').split(',') if s.strip()]
    
    # Single ID
    return [spec.strip()]


def main():
    args = parse_args()

    logger.info(f"Loading BGEN from {args.bgen}")
    logger.info(f"Sample file: {args.sample}")

    try:
        from bgen import BgenReader
    except ImportError:
        logger.error("bgen package not available")
        sys.exit(1)

    try:
        bfile = BgenReader(args.bgen, sample_path=args.sample, delay_parsing=True)
    except Exception as e:
        logger.error(f"Failed to load BGEN: {e}")
        sys.exit(1)

    # Get sample IDs
    bgen_samples = bfile.samples
    if callable(bgen_samples):
        bgen_samples = list(bgen_samples())
    else:
        bgen_samples = list(bgen_samples)
    logger.info(f"BGEN has {len(bgen_samples)} samples")

    # Parse IDs
    variant_ids = read_id_list(args.variant_ids)
    sample_ids = read_id_list(args.sample_ids)

    logger.info(f"Requesting {len(variant_ids)} variants and {len(sample_ids)} samples")

    # Match samples
    bgen_to_idx = {s: i for i, s in enumerate(bgen_samples)}
    sample_indices = []
    sample_ids_found = []
    
    for sid in sample_ids:
        if sid in bgen_to_idx:
            sample_indices.append(bgen_to_idx[sid])
            sample_ids_found.append(sid)

    if len(sample_indices) < args.min_samples:
        logger.error(f"Only {len(sample_indices)} samples, need {args.min_samples}")
        sys.exit(1)

    logger.info(f"Matched {len(sample_indices)} samples")

    # Iterate through BGEN to extract dosages
    variant_set = set(variant_ids)
    dosages = []
    effect_alleles = []
    allele_pairs = []
    mafs_observed = []
    variants_found = []

    logger.info("Extracting dosages by iterating through BGEN...")
    n_checked = 0
    sample_indices_list = list(sample_indices)

    for var in bfile:
        n_checked += 1
        if n_checked % 10000 == 0:
            logger.info(f"  Checked {n_checked} variants, found {len(dosages)} requested")

        if var.rsid not in variant_set:
            continue

        # Extract dosage for requested samples
        # alt_dosage is allele[1] count
        # dosage for allele[0] = 2 - alt_dosage
        try:
            alt_dosages = var.alt_dosage
            dosage = 2.0 - np.array([alt_dosages[i] for i in sample_indices_list], dtype=np.float32)
        except Exception as e:
            logger.warning(f"Failed to extract {var.rsid}: {e}")
            continue

        dosages.append(dosage)
        effect_alleles.append(var.alleles[0])
        allele_pairs.append((var.alleles[0], var.alleles[1]))
        variants_found.append(var.rsid)

        # Compute realized MAF
        maf = np.mean(dosage) / 2.0
        maf_observed = min(maf, 1.0 - maf)
        mafs_observed.append(float(maf_observed))

        # Stop early if we found all variants
        if len(dosages) == len(variant_set):
            break

    logger.info(f"Checked {n_checked} variants, extracted {len(dosages)}")

    if len(dosages) < args.min_variants:
        logger.error(f"Only {len(dosages)} variants extracted, need {args.min_variants}")
        sys.exit(1)

    # Stack into matrix
    dosage_matrix = np.array(dosages, dtype=np.float32)
    logger.info(f"Dosage matrix shape: {dosage_matrix.shape}")

    # Write output TSV
    logger.info(f"Writing to {args.out}")
    with open(args.out, 'w') as f:
        f.write('variant_id')
        for sid in sample_ids_found:
            f.write(f'\t{sid}')
        f.write('\n')

        for j, rsid in enumerate(variants_found):
            f.write(rsid)
            for i in range(len(sample_ids_found)):
                f.write(f'\t{dosage_matrix[j, i]:.6g}')
            f.write('\n')

    logger.info(f"Wrote {len(dosages)} variants x {len(sample_ids_found)} samples")

    # Write metadata
    if args.meta_out:
        metadata = {
            'n_variants': len(dosages),
            'n_samples': len(sample_ids_found),
            'variants': variants_found,
            'samples': sample_ids_found,
            'effect_alleles': effect_alleles,
            'allele_pairs': [list(ap) for ap in allele_pairs],
            'realized_mafs': mafs_observed,
        }
        with open(args.meta_out, 'w') as f:
            json.dump(metadata, f, indent=2)
        logger.info(f"Wrote metadata to {args.meta_out}")

    logger.info("Done!")


if __name__ == '__main__':
    main()
