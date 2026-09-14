# Pooled-covariance dispersion study (branch `sim/dispersion`)

Findings: `docs/FINDINGS_pooled_dispersion.md`.

## Layout

    src/truth.py           exact population law of Y given g (copied from ../poolcal)
    src/make_truth.py      exact Sigma_0 and exact HC0 sandwich (copied from ../poolcal)
    src/fdr_exact.py       NEW - exact false-dispersion rate engine
    src/validate_exact.py  NEW - exact engine vs the shipped-code Monte Carlo
    src/guard.py           NEW - fits and tests the per-variant guard G
    src/guardmc.py         NEW - independent MC check of the guard
    src/report.py          NEW - builds the headline tables
    src/analyse.R          aggregates poolcal's 790 MC chunks (copied from ../poolcal)
    src/run_analyse.sh     sbatch wrapper for analyse.R
    src/run_exact.sh       sbatch array (36 tasks) for the exact sweep
    out/exact/*.tsv        exact sweep: 2 shapes x 6 MAF x 3 N x Z-ladder 3..60
    out/mc/                salvaged MC: params.tsv, tau.tsv, popexact.tsv
    out/validation.tsv     exact vs MC, per cell
    out/tables.txt         headline tables

## Reproducing

    sbatch src/run_analyse.sh          # salvage the MC in ../poolcal/rds
    sbatch src/run_exact.sh            # exact sweep
    python src/report.py               # tables
    python src/validate_exact.py       # exact vs MC
    python src/guard.py
    python src/guardmc.py

Python needs the fungwas-stage1 conda env; analyse.R needs
"module load languages/R/4.5.1" and "Rscript --vanilla".

## Provenance

The Monte Carlo was NOT re-run here. ../poolcal had completed all 790 chunks
across 61 cells before its driving agent died; what had failed was the R package
rebuild (a stale 00LOCK directory), which left two job arrays in
DependencyNeverSatisfied. Those were cancelled and the finished RDS files
aggregated as-is. The MC drives the shipped Stage 1 (compute_rif, residualize,
compute_sigma0, block jackknife) and the shipped Stage 2 (param_gwas), with the
exact per-individual HC0 sandwich as reference.

The exact engine in fdr_exact.py is new and is the basis for every 5e-8 number,
since that level is unreachable by simulation.
