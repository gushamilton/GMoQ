# `sim/` — simulation campaigns

Each subdirectory is one campaign. The readable result is always the named
document in `docs/`; the directory holds the scripts and the aggregate tables.
None of this is needed to use the package.

`docs/SIMULATION_INVENTORY.md` audits all of them — what each asked, what it
concluded, and whether that conclusion still holds under the 0.4.0 defaults.
Read it before reusing any number from here.

## On `main`

| directory | question | result |
|---|---|---|
| `covsource_mixtureW/` | Does the pooled-covariance defect transfer to the 3-column mixture `W` the paper ships? | `tables/` — **no findings document yet**; `README.md:58` and `docs/BACKLOG.md` B1 quote its numbers |
| `projeff/` | Why is the OLS-vs-GLS projection penalty 57% in simulation and 0.8% on real testosterone? | [`docs/EFFICIENCY_BOUND.md`](../docs/EFFICIENCY_BOUND.md) |
| `nscale/` | Does the calibration measured at N = 30,000 survive to N = 461,320? | [`docs/N_SCALING.md`](../docs/N_SCALING.md) |

## Archived

Removed from `main` and frozen at tag `archive/sim-campaigns-20260917` (branch
`archive/sim/campaigns-20260917`), both at `9c47e7b`. Nothing was deleted.

| was | question | result | why archived |
|---|---|---|---|
| `covsource/` | Which covariance source should Stage 2 default to? | [`docs/FINDINGS_covsource_default.md`](../docs/FINDINGS_covsource_default.md) | its aggregates are byte-identical to `Gmoq_paper/outputs/simulations/biobank_covariance/full_audit/`, which is what the paper figures read |
| `density/` | Does G-MoQ hold up on phenotypes that are not smooth? | [`docs/BP2_lumpy_phenotypes.md`](../docs/BP2_lumpy_phenotypes.md) | the document reproduces its own tables; the failure it found is on the opt-in GLS path |
| `out/`, `src/` | Is the pooled covariance safe as the default? | [`docs/FINDINGS_pooled_dispersion.md`](../docs/FINDINGS_pooled_dispersion.md) | that document is retired guidance: its recommendation is withdrawn and its proposed flag has been replaced by `correct_dispersion_skew()` |

To read one:

    git show archive/sim-campaigns-20260917:sim/out/tables.txt
    git worktree add /tmp/simarchive archive/sim-campaigns-20260917
