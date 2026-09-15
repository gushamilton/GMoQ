"""Compare labelled Stage 1 and REGENIE ADD results; never infer a global flip.

Run on a compute node when inputs are on BluePebble. Produces aggregates only.
"""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np


def rows(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        header = handle.readline().split()
        for line in handle:
            yield dict(zip(header, line.split()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage1', required=True)
    parser.add_argument('--regenie', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    reg = {}
    for row in rows(args.regenie):
        if row.get('TEST') != 'ADD':
            continue
        key = row['ID']
        if key in reg:
            raise ValueError(f'Duplicate REGENIE ADD ID: {key}')
        reg[key] = row
    counts = dict(same_effect_allele=0, reversed_effect_allele=0,
                  incompatible_pair=0, stage1_only=0)
    raw, aligned, reference = [], [], []
    seen = set()
    for row in rows(args.stage1):
        key = row['snp_id']
        if key in seen:
            raise ValueError(f'Duplicate Stage 1 ID: {key}')
        seen.add(key)
        if key not in reg:
            counts['stage1_only'] += 1
            continue
        other = reg[key]
        pair = row['effect_allele'], row['other_allele']
        target = other['ALLELE1'], other['ALLELE0']
        if pair == target:
            sign = 1
            counts['same_effect_allele'] += 1
        elif pair == target[::-1]:
            sign = -1
            counts['reversed_effect_allele'] += 1
        else:
            counts['incompatible_pair'] += 1
            continue
        beta = np.mean([float(v) for k, v in row.items() if k.startswith('beta_tau_')])
        ref = float(other['BETA'])
        if np.isfinite(beta) and np.isfinite(ref):
            raw.append(beta)
            aligned.append(sign * beta)
            reference.append(ref)
    if len(reference) < 2 or np.std(reference) == 0 or np.std(raw) == 0:
        raise ValueError('Need at least two compatible, finite, varying effects for correlation.')
    result = {'stage1': str(Path(args.stage1).resolve()),
              'regenie': str(Path(args.regenie).resolve()), **counts,
              'finite_pairs': len(reference),
              'correlation_before_alignment': float(np.corrcoef(raw, reference)[0, 1]),
              'correlation_after_alignment': float(np.corrcoef(aligned, reference)[0, 1]),
              'interpretation': 'Allele-label and empirical agreement check, not proof of unbiasedness.'}
    Path(args.out).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
