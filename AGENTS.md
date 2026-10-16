# G-MoQ continuation

`main` holds the two packages, their tests, CI, user docs and the 1000 Genomes
example, and nothing else. Keep it that way: run records, simulation campaigns
and one-off diagnostics do not go on `main`.

Before launching an analysis or simulation, read the experiment register and the
evidence for the same question at tag `archive/research-records-20261016`
(`experiments/README.md`, `sim/`, `docs/`). Reuse completed evidence; do not rerun
an existing experiment merely because it is absent from current conversation
context. Record new runs with distinct IDs, committed source, settings, status
and output paths, outside `main`. Preserve failed results as well as successes.

OLS with pooled covariance is the primary approach. Known calibration limits
are recorded in `docs/KNOWN_ISSUES.md` (and, in full, in
`docs/VALIDATION_2026-09-15.md` at the archive tag); passing implementation tests
is not proof that those statistical limits are repaired.

Use Luna xhigh subagents for the user-authorized implementation work. All
BluePebble execution runs on compute nodes, including small R/Python commands,
tests and builds. Login nodes are for inspection, files, Git and job control.
Coordinate one held allocation and bounded worker resources; never start a
second allocation independently or alter another campaign's jobs.
