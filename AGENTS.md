# G-MoQ continuation

Before launching an analysis or simulation, read `experiments/README.md` and the
linked evidence for the same question. Reuse completed evidence; do not rerun
or overwrite an existing experiment merely because it is absent from current
conversation context. Register new runs with distinct IDs, committed source,
settings, status and output paths. Preserve failed results as well as successes.

OLS with pooled covariance is the primary approach. Known calibration limits
are recorded in `docs/VALIDATION_2026-09-15.md`; passing implementation tests is
not proof that those statistical limits are repaired.

Use Luna xhigh subagents for the user-authorized implementation work. All
BluePebble execution runs on compute nodes, including small R/Python commands,
tests and builds. Login nodes are for inspection, files, Git and job control.
Coordinate one held allocation and bounded worker resources; never start a
second allocation independently or alter another campaign's jobs.
