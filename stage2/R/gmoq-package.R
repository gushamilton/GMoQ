#' gmoq: Stage 2 of the Genomic Method of Quantiles (G-MoQ)
#'
#' Maps the per-quantile SNP effects estimated by Stage 1 onto the parameters of
#' a distributional model. For each variant it fits `b ~ W theta` across the
#' quantile grid, where `b` holds the quantile effects and `W` encodes the
#' model (for example a mean and a variance, or a two-component mixture), and
#' propagates the Stage 1 covariance to standard errors and tests.
#'
#' The `useDynLib` directive is declared here, rather than written into
#' `NAMESPACE` by hand, so that regenerating the documentation keeps it.
#'
#' @useDynLib gmoq, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats dnorm lm.fit pnorm qnorm quantile
#' @keywords internal
"_PACKAGE"
