#' fungwasStage2: Stage 2 of the Genomic Method of Quantiles
#'
#' Fits `b ~ W theta` across the quantile grid, where `b` holds the per-quantile
#' effects Stage 1 estimated for a variant and `W` encodes a candidate shape for
#' the distributional effect (additive, shifted, both, or a variance QTL).
#'
#' @section Why this file exists:
#' The two directives below were previously hand-written at the foot of
#' `NAMESPACE`. `devtools::document()` regenerates that file from roxygen
#' comments alone, so every documentation pass silently dropped them, and any
#' commit taken while they were missing shipped a package whose compiled kernel
#' never loads -- every tau-covariance path that reaches C++ then fails at run
#' time. Declaring them here means roxygen regenerates them rather than
#' discarding them. Do not move them back into `NAMESPACE` by hand.
#'
#' `stage2/tests/test-namespace-dynlib.R` asserts they survive.
#'
#' @useDynLib fungwasStage2, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats dnorm lm.fit pnorm qnorm quantile
#' @keywords internal
"_PACKAGE"
