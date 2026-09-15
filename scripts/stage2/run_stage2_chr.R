#!/usr/bin/env Rscript

## Run one current-pin Stage 2 fit from one current-pin Stage 1 chromosome.
## The caller is expected to invoke this once per chromosome and concatenate
## the resulting TSV/gz streams only after all chromosome checks pass.

args <- commandArgs(trailingOnly = TRUE)
if (!grepl("^bp1-compute", Sys.info()[["nodename"]])) {
  stop("Run Stage 2 on a BP compute node (bp1-compute*).", call. = FALSE)
}
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) args[[i + 1L]] else default
}

r_lib <- get_arg("--r-lib")
if (!is.null(r_lib)) .libPaths(c(r_lib, .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
})
if (!requireNamespace("fungwasStage2", quietly = TRUE)) {
  stop("fungwasStage2 is unavailable; point --r-lib at the current freshRlib.",
       call. = FALSE)
}

trait <- get_arg("--trait")
model <- get_arg("--model")
chr_raw <- get_arg("--chr")
stage1_file <- get_arg("--stage1")
weights_file <- get_arg("--weights")
out_prefix <- get_arg("--out")
cov_source <- get_arg("--cov-source", "pooled")
cov_file <- get_arg("--cov-file")
sigma0_file <- get_arg("--sigma0-file")
cov_meta_file <- get_arg("--cov-meta-file")
n_threads <- suppressWarnings(as.integer(get_arg(
  "--threads", Sys.getenv("SLURM_CPUS_PER_TASK", "1")
)))
if (!is.finite(n_threads) || n_threads < 1L) n_threads <- 1L

if (any(vapply(list(trait, model, chr_raw, stage1_file, weights_file, out_prefix),
              is.null, logical(1)))) {
  stop(paste(
    "Usage: Rscript run_stage2_chr.R --trait <trait> --model <model>",
    "--chr <1..22> --stage1 <file> --weights <file> --out <prefix>",
    "--cov-source <pooled|jackknife> [--sigma0-file <file> | --cov-file <file>]",
    sep = "\n"
  ), call. = FALSE)
}
if (!trait %in% c(
  "bmi", "log_bmi", "weight",
  "testosterone", "testosterone_shrink_0.2", "testosterone_shrink_0.5",
  "testosterone_shrink_0.7", "jcv_log_all", "mcv_log_all"
)) stop("Unsupported paper trait: ", trait, call. = FALSE)
valid_models <- list(
  bmi = c("additive", "shifted", "both", "logscale_vqtl"),
  log_bmi = c("additive", "shifted", "both"),
  weight = c("additive", "shifted", "both"),
  testosterone = "skew",
  testosterone_shrink_0.2 = "skew",
  testosterone_shrink_0.5 = "skew",
  testosterone_shrink_0.7 = "skew",
  jcv_log_all = "mixture",
  mcv_log_all = "mixture"
)
if (!model %in% valid_models[[trait]]) {
  stop("Unsupported model ", model, " for trait ", trait, call. = FALSE)
}
if (!cov_source %in% c("pooled", "jackknife")) {
  stop("--cov-source must be pooled or jackknife.", call. = FALSE)
}
chr_num <- suppressWarnings(as.integer(chr_raw))
if (!is.finite(chr_num) || chr_num < 1L || chr_num > 22L) {
  stop("--chr must be an autosome from 1 through 22.", call. = FALSE)
}
chr <- as.character(chr_num)
for (path in c(stage1_file, weights_file)) {
  if (!file.exists(path)) stop("Missing required file: ", path, call. = FALSE)
}

## Stage 1's header is the authoritative serialised grid.  Do this check here
## too, so a weight file from raw BMI cannot accidentally be supplied to log-BMI.
hdr <- names(fread(stage1_file, nrows = 0L, showProgress = FALSE))
beta_cols <- grep("^beta_tau_", hdr, value = TRUE)
if (length(beta_cols) != 45L) stop("Expected 45 beta_tau_* columns in ", stage1_file,
                                    call. = FALSE)
taus <- as.numeric(sub("^beta_tau_", "", beta_cols))
if (any(!is.finite(taus)) || any(diff(taus) <= 0)) {
  stop("Stage 1 tau columns are not finite and increasing.", call. = FALSE)
}
## Stage 1 serialises the current pin's 45-level grid to three decimals.
## Pin that grid here as well so a legacy or hand-edited Stage 1 file cannot
## silently select a different W indexing scheme.
expected_taus <- round(seq(0.02, 0.98, length.out = 45L), 3L)
if (!isTRUE(all.equal(taus, expected_taus, tolerance = 5e-7))) {
  stop("Stage 1 tau grid is not the pinned 45-level grid; refusing Stage 2.",
       call. = FALSE)
}

wobj <- readRDS(weights_file)
W <- if (is.list(wobj) && !is.null(wobj$W)) wobj$W else wobj
w_taus <- if (is.list(wobj) && !is.null(wobj$taus)) as.numeric(wobj$taus) else NULL
W <- as.matrix(W)
if (!is.null(w_taus)) {
  if (nrow(W) != length(w_taus) && ncol(W) == length(w_taus)) W <- t(W)
  if (length(w_taus) != length(taus) ||
      !isTRUE(all.equal(w_taus, taus, tolerance = 5e-7))) {
    stop("Weight tau grid differs from the serialized Stage 1 tau grid; refusing W.",
         call. = FALSE)
  }
}
if (nrow(W) != length(taus) || any(!is.finite(W))) {
  stop("Weight matrix must be finite with one row per Stage 1 tau.", call. = FALSE)
}
if (is.null(colnames(W))) colnames(W) <- paste0("param", seq_len(ncol(W)))
if (any(!nzchar(colnames(W)))) stop("Weight matrix has unnamed parameters.", call. = FALSE)
if (!is.list(wobj) || is.null(wobj$provenance)) {
  stop("Weight file lacks matched-baseline provenance; regenerate it with ",
       "prepare_matched_weights.R (legacy raw-BMI W is not accepted).", call. = FALSE)
}
p <- wobj$provenance
if (!identical(as.character(p$phenotype), trait)) {
  stop("Weight provenance phenotype is ", p$phenotype, "; expected ", trait,
       call. = FALSE)
}
if (!identical(as.character(p$chromosome), chr)) {
  stop("Weight provenance chromosome is ", p$chromosome, "; expected ", chr,
       call. = FALSE)
}

if (identical(cov_source, "pooled")) {
  if (is.null(sigma0_file)) stop("Pooled Stage 2 requires --sigma0-file explicitly.",
                                 call. = FALSE)
  if (!file.exists(sigma0_file)) stop("Missing Sigma_0 file: ", sigma0_file, call. = FALSE)
  sigma0_meta <- sub("\\.sigma0\\.tsv\\.gz$", ".sigma0.meta.json", sigma0_file)
  if (!file.exists(sigma0_meta)) stop("Missing Sigma_0 metadata: ", sigma0_meta,
                                     call. = FALSE)
  inferred_cov_meta <- sub("\\.sigma0\\.tsv\\.gz$", ".cov.meta.json", sigma0_file)
  if (is.null(cov_meta_file)) cov_meta_file <- inferred_cov_meta
  if (!file.exists(cov_meta_file)) stop("Missing Stage 1 covariance metadata: ",
                                        cov_meta_file, call. = FALSE)
} else {
  if (is.null(cov_file)) stop("Jackknife Stage 2 requires --cov-file explicitly.",
                              call. = FALSE)
  if (!file.exists(cov_file)) stop("Missing jackknife covariance file: ", cov_file,
                                   call. = FALSE)
  cov_ids_file <- sub("\\.cov\\.gz$", ".cov.ids.tsv.gz", cov_file)
  if (!file.exists(cov_ids_file)) stop("Missing jackknife covariance IDs: ", cov_ids_file,
                                       call. = FALSE)
  if (is.null(cov_meta_file)) cov_meta_file <- sub("\\.cov\\.gz$", ".cov.meta.json", cov_file)
  if (!file.exists(cov_meta_file)) stop("Missing jackknife covariance metadata: ",
                                        cov_meta_file, call. = FALSE)
}

return_cov <- if (model %in% c("both", "logscale_vqtl", "skew", "mixture")) {
  "upper"
} else {
  "none"
}
param_cov_file <- if (identical(return_cov, "upper")) paste0(out_prefix, ".param_cov.gz") else NULL
dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
existing_outputs <- c(paste0(out_prefix, ".tsv.gz"),
                      paste0(out_prefix, ".manifest.tsv"), param_cov_file)
existing_outputs <- existing_outputs[file.exists(existing_outputs)]
if (length(existing_outputs)) {
  stop("Refusing to overwrite existing Stage 2 output(s): ",
       paste(existing_outputs, collapse = ", "),
       ". Use a new exclusive run root.", call. = FALSE)
}

message("Stage 2 trait=", trait, " model=", model, " chr=", chr,
        " cov_source=", cov_source, " estimator=ols")
fit <- fungwasStage2::param_gwas_from_file(
  stage1_file = stage1_file,
  W = W,
  se_mode = "tau_cov",
  cov_file = if (identical(cov_source, "jackknife")) cov_file else NULL,
  cov_source = cov_source,
  sigma0_file = if (identical(cov_source, "pooled")) sigma0_file else NULL,
  cov_meta_file = cov_meta_file,
  tau_cov_estimator = "ols",
  return_cov = return_cov,
  param_cov_file = param_cov_file,
  n_threads = n_threads,
  q_method = "saddlepoint"
)
out_dt <- if (is.list(fit) && !is.null(fit$results)) fit$results else fit
if (!is.data.frame(out_dt) || !nrow(out_dt)) stop("Stage 2 returned no result rows.",
                                                    call. = FALSE)
out_path <- paste0(out_prefix, ".tsv.gz")
tmp_path <- paste0(out_path, ".tmp.", Sys.getpid())
fwrite(out_dt, tmp_path, sep = "\t")
if (!file.rename(tmp_path, out_path)) stop("Could not move Stage 2 output into place: ",
                                           out_path, call. = FALSE)

## This sidecar contains no participant-level information.  It lets the
## concatenator verify that all chromosome chunks used identical model/source/W
## choices without loading any covariance payload itself.
manifest <- data.table(
  trait = trait, model = model, chromosome = chr, covariance_source = cov_source,
  tau_cov_estimator = "ols", n_taus = length(taus), n_parameters = ncol(W),
  n_snps = nrow(out_dt), stage1_file = normalizePath(stage1_file, mustWork = TRUE),
  weights_file = normalizePath(weights_file, mustWork = TRUE), output_file = out_path
)
fwrite(manifest, paste0(out_prefix, ".manifest.tsv"), sep = "\t")
message("Wrote ", nrow(out_dt), " SNPs to ", out_path)
