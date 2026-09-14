# ---------------------------------------------------------------------------
# The pooled cross-quantile covariance as an alternative to the block jackknife
#
# Stage 1 can emit one T x T matrix per trait per chromosome,
#
#     Sigma_0 = (1/n) Psi_r' Psi_r
#
# with Psi_r the covariate-residualised recentered influence functions. The
# sampling covariance of the per-tau slopes for a variant is then
#
#     Sigma_v = Sigma_0 / d_resid,     d_resid = sum_i g_resid_i^2
#
# because the cross-quantile correlation is a property of the phenotype and the
# covariates, not of the genotype; only the scale is per-variant, and d_resid is
# a column of the Stage 1 TSV since v0.3.0.
#
# WHY IT IS WORTH HAVING. Measured on the published chr22 Stage 1 inputs
# (N = 461,320) against an exact per-individual HC0 sandwich, 500 variants per
# trait: median diagonal ratio 1.00000 and median cross-tau correlation
# deviation 0.0053 for both bmi and weight, against the 65-block jackknife's
# 5-95% spread of 0.725 to 1.307 on the variance. Stage 2 standard errors agree
# to within 0.03%. It is positive semi-definite by construction, has no block
# count and no seed, and costs 8 KB per chromosome against 229-289 bytes per
# variant.
#
# WHY THE JACKKNIFE PATH STAYS. Every published result used it. cov_source is a
# switch, not a replacement, and "jackknife" remains the default.
#
# WHY THIS REFUSES RATHER THAN WARNS. Sigma_v = Sigma_0 / d_resid holds only
# when Sigma_0 was built from the same phenotype, covariate set, tau grid, LOCO
# predictor and sample as the scan. A mismatched Sigma_0 leaves the point
# estimates untouched and the standard errors plausible: the measured failure
# for a covariance blind to the covariates is a 1.7% variance inflation on bmi
# and 25.6% on weight, in the same direction for every variant, with nothing
# downstream able to see it. So Stage 1 writes a fingerprint into both the
# Sigma_0 sidecar and the covariance metadata, and a disagreement stops the run.
# This is the same precedent as refusing to emit a p-value when the degrees of
# freedom are unknown: a silent fallback is how the normal-reference defect
# survived.
# ---------------------------------------------------------------------------

#' Chromosome key used when Stage 1 ran without a per-chromosome adjustment
#' @keywords internal
SIGMA0_ALL_CHROM <- "ALL"

#' Highest \code{.sigma0.meta.json} schema version this Stage 2 understands
#' @keywords internal
SIGMA0_META_VERSION_MAX <- 1L

#' Normalise a chromosome label the way Stage 1 does
#'
#' Numeric autosomes lose leading zeroes, a \code{chr} prefix is stripped, and
#' non-numeric labels are upper-cased -- so \code{04}, \code{4} and \code{chr04}
#' all key the same matrix.
#'
#' @param x Vector of chromosome labels.
#' @return Character vector of normalised labels.
#' @keywords internal
.normalize_chr <- function(x) {
  s <- trimws(as.character(x))
  s <- sub("^(chr|CHR)", "", s)
  num <- grepl("^[0-9]+$", s)
  s[num] <- as.character(as.integer(s[num]))
  s[!num] <- toupper(s[!num])
  s
}

#' Extract a numeric array from a JSON document with base R
#'
#' Companion to \code{.json_number}, used only on sidecars Stage 1 writes
#' itself, so that reading them never depends on \pkg{jsonlite}.
#'
#' @param txt Whole JSON document as one string.
#' @param key Key name.
#' @return Numeric vector, possibly length zero.
#' @keywords internal
.json_number_array <- function(txt, key) {
  v <- .json_capture(txt, paste0("\"", key, "\"[[:space:]]*:[[:space:]]*\\[([^]]*)\\]"))
  if (is.na(v)) return(numeric(0))
  parts <- trimws(strsplit(v, ",", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts)]
  suppressWarnings(as.numeric(parts))
}

#' Extract a string array from a JSON document with base R
#'
#' @param txt Whole JSON document as one string.
#' @param key Key name.
#' @return Character vector, possibly length zero.
#' @keywords internal
.json_string_array <- function(txt, key) {
  v <- .json_capture(txt, paste0("\"", key, "\"[[:space:]]*:[[:space:]]*\\[([^]]*)\\]"))
  if (is.na(v)) return(character(0))
  m <- gregexpr("\"([^\"]*)\"", v)[[1]]
  if (length(m) == 1L && m[1L] == -1L) return(character(0))
  hits <- regmatches(v, gregexpr("\"[^\"]*\"", v))[[1]]
  gsub("^\"|\"$", "", hits)
}

#' Path of the Sigma_0 metadata sidecar for a matrix file
#'
#' @param sigma0_file Path to a \code{.sigma0.tsv.gz}.
#' @return Path to the matching \code{.sigma0.meta.json}.
#' @keywords internal
.sigma0_meta_file <- function(sigma0_file) {
  sub("\\.sigma0\\.tsv\\.gz$", ".sigma0.meta.json", sigma0_file)
}

#' Path of the Stage 1 covariance metadata for a Sigma_0 matrix file
#'
#' The covariance metadata is what the Sigma_0 sidecar is checked against, so
#' it is required even when the (much larger) \code{.cov.gz} is not read.
#'
#' @param sigma0_file Path to a \code{.sigma0.tsv.gz}.
#' @return Path to the matching \code{.cov.meta.json}.
#' @keywords internal
.sigma0_cov_meta_file <- function(sigma0_file) {
  sub("\\.sigma0\\.tsv\\.gz$", ".cov.meta.json", sigma0_file)
}

#' Read the Sigma_0 metadata sidecar
#'
#' @param path Path to a \code{.sigma0.meta.json}.
#' @return A list of the recorded fields.
#' @keywords internal
.read_sigma0_meta <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stop("Missing pooled-covariance metadata: ", path,
         ". Sigma_0 cannot be used without it: it is what records the covariate ",
         "set, the tau grid and the sample the matrix was built on, and those ",
         "are the only things that make Sigma_0 / d_resid the covariance of ",
         "this scan rather than of some other one.", call. = FALSE)
  }
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    m <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  } else {
    txt <- paste(readLines(path, warn = FALSE), collapse = " ")
    m <- list(
      sigma0_meta_version = .json_number(txt, "sigma0_meta_version"),
      package_version     = .json_string(txt, "package_version"),
      estimator           = .json_string(txt, "estimator"),
      phenotype           = .json_string(txt, "phenotype"),
      n_taus              = .json_number(txt, "n_taus"),
      taus                = .json_number_array(txt, "taus"),
      n_samples           = .json_number(txt, "n_samples"),
      n_covariates        = .json_number(txt, "n_covariates"),
      covariates          = .json_string_array(txt, "covariates"),
      loco_predictor      = .json_string(txt, "loco_predictor"),
      chromosomes         = .json_string_array(txt, "chromosomes"),
      fingerprint         = .json_string(txt, "fingerprint")
    )
  }
  v <- suppressWarnings(as.numeric(m$sigma0_meta_version)[1L])
  if (!is.finite(v) || v > SIGMA0_META_VERSION_MAX) {
    stop("Pooled-covariance sidecar ", path, " declares schema version ",
         m$sigma0_meta_version, ", which this version of fungwasStage2 does not ",
         "understand (it reads up to ", SIGMA0_META_VERSION_MAX, "). Upgrade the ",
         "package rather than reading a schema it may misinterpret.", call. = FALSE)
  }
  m
}

#' Read the Sigma_0 matrices
#'
#' @param path Path to a \code{.sigma0.tsv.gz}.
#' @param n_taus Expected number of tau levels.
#' @return Named list of \code{n_taus x n_taus} matrices, keyed by normalised
#'   chromosome label.
#' @keywords internal
.read_sigma0_matrices <- function(path, n_taus) {
  if (!file.exists(path)) {
    stop("Missing pooled-covariance matrix file: ", path, call. = FALSE)
  }
  # Base R rather than fread: the file is a few hundred rows, the chromosome
  # column must stay character ("X", "MT"), and this keeps the sidecar readable
  # with no dependency at all.
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  dt <- utils::read.table(con, header = TRUE, sep = "\t", quote = "",
                          colClasses = "character", check.names = FALSE)
  n_taus <- as.integer(n_taus)
  if (ncol(dt) != n_taus + 2L) {
    stop(path, " has ", ncol(dt), " columns; expected ", n_taus + 2L,
         " (chr, tau, and one per tau).", call. = FALSE)
  }
  chr <- .normalize_chr(dt[[1L]])
  vals <- matrix(as.numeric(unlist(dt[, 3:ncol(dt)], use.names = FALSE)),
                 nrow = nrow(dt), ncol = n_taus)
  out <- list()
  for (k in unique(chr)) {
    rows <- which(chr == k)
    if (length(rows) != n_taus) {
      stop(path, ": chromosome ", k, " has ", length(rows), " rows, expected ",
           n_taus, ".", call. = FALSE)
    }
    S <- vals[rows, , drop = FALSE]
    dimnames(S) <- NULL
    if (!all(is.finite(S))) {
      stop(path, ": Sigma_0 for chromosome ", k, " is not finite.", call. = FALSE)
    }
    out[[k]] <- 0.5 * (S + t(S))
  }
  out
}

#' Read a Stage 1 pooled covariance sidecar
#'
#' @param sigma0_file Path to a \code{.sigma0.tsv.gz} written by Stage 1.
#' @return List with \code{meta} (the sidecar metadata) and \code{matrices}
#'   (named list of T x T matrices, keyed by normalised chromosome label, or by
#'   \code{"ALL"} when Stage 1 ran without \code{--loco-preds}).
#' @export
read_sigma0 <- function(sigma0_file) {
  meta <- .read_sigma0_meta(.sigma0_meta_file(sigma0_file))
  list(meta = meta,
       matrices = .read_sigma0_matrices(sigma0_file, meta$n_taus))
}

#' Refuse a Sigma_0 that does not belong to the scan being analysed
#'
#' Checks, in order: the tau grid against the Stage 1 column names; the
#' provenance block against the Stage 1 covariance metadata (phenotype,
#' covariate columns, sample count, LOCO predictor, fingerprint); and any
#' expectation the caller stated explicitly.
#'
#' @param s0_meta List from \code{.read_sigma0_meta()}.
#' @param cov_meta List from \code{.read_cov_meta()}, or \code{NULL}.
#' @param taus Tau grid recovered from the Stage 1 file.
#' @param cov_meta_path Path reported in error messages.
#' @param require_covariates Optional character vector the sidecar must match.
#' @param require_n_samples Optional sample count the sidecar must match.
#' @return Invisible \code{TRUE}; stops otherwise.
#' @keywords internal
.check_sigma0_compatible <- function(s0_meta, cov_meta, taus, cov_meta_path,
                                     require_covariates = NULL,
                                     require_n_samples = NULL) {
  refuse <- function(...) stop("Refusing to use the pooled covariance: ", ...,
                               "\nSigma_v = Sigma_0 / d_resid is the sampling covariance of ",
                               "beta_tau only when Sigma_0 was built from the same phenotype, ",
                               "covariate set, tau grid, LOCO predictor and sample as this scan. ",
                               "A mismatched Sigma_0 biases every standard error in the same ",
                               "direction and leaves the point estimates untouched, so nothing ",
                               "downstream can detect it. Re-run Stage 1, or use ",
                               "cov_source = \"jackknife\".", call. = FALSE)

  # 1. the tau grid, against what the Stage 1 file actually carries.
  theirs <- round(as.numeric(s0_meta$taus), 3)
  ours <- round(as.numeric(taus), 3)
  if (length(theirs) != length(ours) || !isTRUE(all.equal(theirs, ours, tolerance = 0))) {
    refuse("the tau grid differs. Sigma_0 was built on ", length(theirs),
           " levels (", paste(theirs[seq_len(min(4L, length(theirs)))], collapse = ", "),
           if (length(theirs) > 4) ", ..." else "", ") and the Stage 1 file carries ",
           length(ours), " (", paste(ours[seq_len(min(4L, length(ours)))], collapse = ", "),
           if (length(ours) > 4) ", ..." else "", "). Sigma_0 is indexed by tau. ")
  }

  # 2. the provenance block written by the same Stage 1 run.
  if (is.null(cov_meta)) {
    refuse("there is no Stage 1 covariance metadata at ", cov_meta_path,
           ", so there is nothing to check this Sigma_0 against. That file is ",
           "written by every Stage 1 run beside the covariance itself. ")
  }
  analysis <- if (is.list(cov_meta)) cov_meta[["analysis"]] else NULL
  if (!is.list(analysis) || is.null(analysis$fingerprint) ||
      all(is.na(analysis$fingerprint))) {
    refuse("the Stage 1 covariance metadata at ", cov_meta_path,
           " carries no 'analysis' block, so there is nothing to check this ",
           "Sigma_0 against. That file was written by a Stage 1 older than the ",
           "pooled covariance. ")
  }
  cmp <- function(field, label, mine, theirs_) {
    a <- mine; b <- theirs_
    if (is.null(a)) a <- NA
    if (is.null(b)) b <- NA
    same <- length(a) == length(b) && all(as.character(a) == as.character(b), na.rm = FALSE)
    if (isTRUE(all(is.na(a))) && isTRUE(all(is.na(b)))) same <- TRUE
    if (!isTRUE(same)) {
      refuse(label, " differs: ", paste(as.character(a), collapse = ", "),
             " in the Sigma_0 sidecar against ", paste(as.character(b), collapse = ", "),
             " in ", cov_meta_path, ". ")
    }
  }
  cmp("phenotype", "the phenotype", s0_meta$phenotype, analysis$phenotype)
  cmp("n_samples", "the analysed sample count",
      as.integer(s0_meta$n_samples), as.integer(analysis$n_samples))
  cmp("covariates", "the covariate set",
      as.character(s0_meta$covariates), as.character(analysis$covariates))
  cmp("loco_predictor", "the LOCO predictor",
      s0_meta$loco_predictor, analysis$loco_predictor)
  cmp("fingerprint", "the analysis fingerprint",
      s0_meta$fingerprint, analysis$fingerprint)

  # 3. anything the caller asserted independently.
  if (!is.null(require_covariates)) {
    if (!identical(as.character(require_covariates), as.character(s0_meta$covariates))) {
      refuse("the covariate set does not match the one the caller required: ",
             paste(as.character(s0_meta$covariates), collapse = ", "),
             " in the sidecar against ",
             paste(as.character(require_covariates), collapse = ", "), " required. ")
    }
  }
  if (!is.null(require_n_samples)) {
    if (as.integer(require_n_samples) != as.integer(s0_meta$n_samples)) {
      refuse("the sample count does not match the one the caller required: ",
             s0_meta$n_samples, " in the sidecar against ", require_n_samples,
             " required. ")
    }
  }
  invisible(TRUE)
}

#' Build tau-covariance data for Stage 2 from Sigma_0 and d_resid
#'
#' Produces exactly the \code{list(cov_vec, offsets)} layout that
#' \code{.read_tau_cov_data()} produces for the jackknife, so every downstream
#' consumer -- the OLS delta-method SEs, the GLS fit, the saddlepoint Q
#' calibration, the parameter covariance output -- runs unchanged. That is
#' deliberate: the covariance source is the only thing that differs between the
#' two paths, and it should be the only thing that differs in the code.
#'
#' @param matrices Named list of T x T matrices from \code{.read_sigma0_matrices()}.
#' @param chr Chromosome label per SNP (Stage 1 \code{chr} column).
#' @param d_resid Numeric vector, Stage 1's \code{d_resid} per SNP.
#' @param n_taus Number of tau levels.
#' @return List with \code{cov_vec} and \code{offsets}, plus \code{n_bad}, the
#'   number of SNPs whose \code{d_resid} was unusable and whose records were
#'   filled with \code{NA} (the downstream kernels skip those).
#' @keywords internal
.pooled_cov_data <- function(matrices, chr, d_resid, n_taus) {
  n_taus <- as.integer(n_taus)
  n_snps <- length(d_resid)
  n_cov_elements <- as.integer(n_taus * (n_taus + 1L) / 2L)
  iu <- which(upper.tri(matrix(0, n_taus, n_taus), diag = TRUE))
  # Row-major upper triangle, the order the Stage 1 writer and the C++ readers
  # both use: (1,1) (1,2) ... (1,T) (2,2) ...
  ord <- order(row(matrix(0, n_taus, n_taus))[iu], col(matrix(0, n_taus, n_taus))[iu])
  iu <- iu[ord]

  keys <- .normalize_chr(chr)
  if (length(matrices) == 1L && names(matrices)[1L] == SIGMA0_ALL_CHROM) {
    keys[] <- SIGMA0_ALL_CHROM
  }
  missing_keys <- setdiff(unique(keys), names(matrices))
  if (length(missing_keys)) {
    stop("Refusing to use the pooled covariance: no Sigma_0 for chromosome(s) ",
         paste(missing_keys, collapse = ", "), "; the sidecar carries ",
         paste(names(matrices), collapse = ", "),
         ". Sigma_0 is per chromosome because the LOCO predictor changes the ",
         "phenotype, so another chromosome's matrix is the wrong object.",
         call. = FALSE)
  }

  upper <- vapply(names(matrices), function(k) as.numeric(matrices[[k]][iu]),
                  numeric(n_cov_elements))
  if (is.null(dim(upper))) upper <- matrix(upper, nrow = n_cov_elements)
  colnames(upper) <- names(matrices)

  d <- as.numeric(d_resid)
  bad <- !is.finite(d) | d <= 0
  d[bad] <- NA_real_

  cov_vec <- numeric(n_cov_elements * n_snps)
  for (k in names(matrices)) {
    idx <- which(keys == k)
    if (!length(idx)) next
    # outer() over (element, SNP): one T(T+1)/2 vector scaled by 1/d per SNP.
    block <- outer(upper[, k], 1 / d[idx])
    pos <- rep((idx - 1L) * n_cov_elements, each = n_cov_elements) +
      rep(seq_len(n_cov_elements), times = length(idx))
    cov_vec[pos] <- as.numeric(block)
  }

  list(cov_vec = cov_vec,
       offsets = seq(0, by = n_cov_elements * 4, length.out = n_snps),
       n_bad = sum(bad))
}

#' Assemble the pooled covariance for a Stage 1 scan
#'
#' Reads the sidecar, refuses on any mismatch with the scan, and returns the
#' \code{cov_data} Stage 2 consumes plus the degrees of freedom that go with it.
#'
#' @param sigma0_file Path to \code{.sigma0.tsv.gz}.
#' @param chr Chromosome label per SNP.
#' @param d_resid Stage 1's \code{d_resid} per SNP.
#' @param taus Tau grid recovered from the Stage 1 file.
#' @param cov_meta_file Path to the Stage 1 \code{.cov.meta.json} to check
#'   against; defaults to the one beside \code{sigma0_file}.
#' @param require_covariates,require_n_samples Optional caller assertions.
#' @return List with \code{cov_data}, \code{meta} and \code{nu}.
#' @keywords internal
.pooled_tau_cov <- function(sigma0_file, chr, d_resid, taus,
                            cov_meta_file = NULL,
                            require_covariates = NULL,
                            require_n_samples = NULL) {
  s0 <- read_sigma0(sigma0_file)
  if (is.null(cov_meta_file)) cov_meta_file <- .sigma0_cov_meta_file(sigma0_file)
  cov_meta <- .read_cov_meta(cov_meta_file)
  .check_sigma0_compatible(s0$meta, cov_meta, taus, cov_meta_file,
                           require_covariates = require_covariates,
                           require_n_samples = require_n_samples)
  if (is.null(d_resid)) {
    stop("Refusing to use the pooled covariance: the Stage 1 file has no ",
         N_EFF_COLUMN, " column. Sigma_v = Sigma_0 / d_resid needs the per-SNP ",
         "scale; without it there is no per-variant covariance to form. That ",
         "column has been emitted since v0.3.0; re-run Stage 1, or use ",
         "cov_source = \"jackknife\".", call. = FALSE)
  }
  cd <- .pooled_cov_data(s0$matrices, chr, d_resid, length(taus))
  if (cd$n_bad > 0) {
    warning(cd$n_bad, " of ", length(d_resid), " SNPs have a non-positive or ",
            "non-finite ", N_EFF_COLUMN, "; their pooled covariance is NA and ",
            "they are skipped downstream.", call. = FALSE)
  }
  list(cov_data = list(cov_vec = cd$cov_vec, offsets = cd$offsets),
       meta = s0$meta,
       nu = as.numeric(s0$meta$n_samples))
}
