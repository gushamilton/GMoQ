#!/usr/bin/env python3
"""Summarise timing.csv on a BP compute node."""
import csv
import statistics
from collections import defaultdict
import sys

groups = defaultdict(lambda: [[], []])
with open(sys.argv[1], newline="") as fh:
    for row in csv.DictReader(fh):
        key = (row["cell"], row["path"])
        groups[key][0].append(float(row["stage_seconds"]))
        groups[key][1].append(float(row["hc0_seconds"]))
for key in sorted(groups):
    stage, hc0 = groups[key]
    ms = statistics.mean(stage)
    mh = statistics.mean(hc0)
    print(*key, "n", len(stage), "stage_mean", ms,
          "hc0_mean", mh, "hc0_over_stage", mh / ms)
