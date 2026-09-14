#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Tier 3: turn the measured timings into report rows.
#
#   Rscript --vanilla tests/tier3/t3_report.R <tier3.json>
#
# Nothing here is a pass/fail. Tier 3 measures cost, not correctness, so every
# row is severity "report" -- except two structural checks that ARE assertions:
# that the run completed at the requested scale, and that the hardware was
# recorded, because an unattributed timing is not reproducible.
# ---------------------------------------------------------------------------
suppressWarnings(suppressMessages({
  .self <- local({
    a <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", a[grep("^--file=", a)])
    normalizePath(f[1L])
  })
  source(file.path(dirname(.self), "..", "lib", "harness.R"))
}))
gmoq_init_roots(.self)
gmoq_use_lib()
gmoq_begin("tier3", "scale_and_timing")

p <- commandArgs(trailingOnly = TRUE)[1L]
j <- jsonlite::fromJSON(p, simplifyVector = TRUE)
hw <- j$hardware; s1 <- j$stage1; s2 <- j$stage2; ex <- j$extrapolation

num <- function(x) if (is.null(x) || !length(x)) NA_real_ else as.numeric(x)[1]
rep_row <- function(id, what, v, crit) gmoq_check(id, what, TRUE, severity = "report",
                                                  stat = num(v), criterion = crit)

gmoq_check("T3.STRUCT.completed_at_requested_scale",
           "Stage 1 ran to completion at the requested sample and variant count.",
           is.finite(num(s1$n_snps)) && num(s1$n_snps) > 0,
           stat = num(s1$n_snps),
           criterion = sprintf("variants processed > 0 at N = %s", hw$n_samples))
gmoq_check("T3.STRUCT.hardware_recorded",
           "The hardware the timings were measured on is recorded with them.",
           !is.null(hw$cpu_model) && nzchar(hw$cpu_model),
           criterion = "CPU model, core counts, memory and node recorded",
           detail = paste(hw$cpu_model, "|", hw$slurm_node, "|", hw$threads_used, "threads"))

hwtxt <- sprintf("N=%s, M=%s, %s threads on %s", hw$n_samples, hw$n_variants,
                 hw$threads_used, hw$cpu_model)

rep_row("T3.S1.wall_sec", "Stage 1 wall clock.", s1$proc_wall_sec, hwtxt)
rep_row("T3.S1.cpu_sec", "Stage 1 CPU time (user + system).", s1$proc_cpu_sec, hwtxt)
rep_row("T3.S1.peak_rss_gb", "Stage 1 peak resident set size.", s1$proc_max_rss_gb, hwtxt)
rep_row("T3.S1.phase_bgen_read_sec", "Stage 1 phase: BGEN read (batch loop less compute and write).", s1$bgen_read_sec, hwtxt)
rep_row("T3.S1.phase_kernel_sec", "Stage 1 phase: the block-score kernel.", s1$batch_compute_sec, hwtxt)
rep_row("T3.S1.phase_write_sec", "Stage 1 phase: writing the covariance and TSV.", s1$batch_write_sec, hwtxt)
rep_row("T3.S1.phase_rif_sec", "Stage 1 phase: RIF build.", s1$rif_build_sec, hwtxt)
rep_row("T3.S1.phase_sigma0_sec", "Stage 1 phase: the pooled covariance sidecar.", s1$sigma0_sec, hwtxt)
rep_row("T3.S1.snps_per_sec", "Stage 1 throughput.", s1$snps_per_sec, hwtxt)

rep_row("T3.S2.wall_sec", "Stage 2 wall clock (whole process).", s2$proc_wall_sec, hwtxt)
rep_row("T3.S2.cpu_sec", "Stage 2 CPU time.", s2$proc_cpu_sec, hwtxt)
rep_row("T3.S2.peak_rss_gb", "Stage 2 peak resident set size.", s2$proc_max_rss_gb, hwtxt)
rep_row("T3.S2.phase_read_tsv_sec", "Stage 2 phase: reading the Stage 1 table.", s2$stage2_read_tsv_sec, hwtxt)
rep_row("T3.S2.phase_read_cov_sec", "Stage 2 phase: reading and dequantising the covariance.", s2$stage2_read_cov_sec, hwtxt)
rep_row("T3.S2.phase_kernel_ols_sec", "Stage 2 phase: the OLS projection and Q calibration.", s2$stage2_kernel_ols_sec, hwtxt)
rep_row("T3.S2.phase_kernel_gls_sec", "Stage 2 phase: the GLS solve.", s2$stage2_kernel_gls_sec, hwtxt)
rep_row("T3.S2.phase_write_sec", "Stage 2 phase: writing the results table.", s2$stage2_write_sec, hwtxt)
rep_row("T3.S2.gls_solved_fraction", "Fraction of variants the GLS solve completed.", s2$gls_solved_fraction, hwtxt)

rep_row("T3.EXTRAP.stage1_cpu_hours_1.2M",
        "Stage 1 CPU hours extrapolated to 1.2M variants. An extrapolation, not a measurement: Stage 1 is linear in the variant count by construction.",
        ex$stage1_cpu_hours_at_1200000, hwtxt)
rep_row("T3.EXTRAP.stage2_wall_hours_1.2M",
        "Stage 2 wall hours extrapolated to 1.2M variants. Time only: memory is measured directly by the scaling rows.",
        ex$stage2_wall_hours_at_1200000, hwtxt)

if (!is.null(s2$scaling) && length(s2$scaling)) {
  sc <- s2$scaling
  for (i in seq_len(nrow(as.data.frame(sc)))) {
    d <- as.data.frame(sc)[i, ]
    rep_row(sprintf("T3.SCALE.stage2_sec_at_%d", as.integer(d$n_variants)),
            "Stage 2 kernel seconds at this variant count (covariance payload tiled; synthetic).",
            d$seconds, sprintf("%d variants, cov_vec %.2f GB in memory",
                               as.integer(d$n_variants), d$cov_vec_gb))
  }
}

gmoq_finish()
