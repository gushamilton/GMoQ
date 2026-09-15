#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("pass at least one R source file")
for (path in args) {
  exprs <- parse(file = path)
  message("R_PARSE_OK ", path)
  if (grepl("prepare_matched_weights\\.R$", path)) {
    kde_expr <- exprs[vapply(exprs, function(e) {
      is.call(e) && identical(as.character(e[[1L]]), "<-") &&
        identical(as.character(e[[2L]]), "stage1_kde_at")
    }, logical(1))]
    if (length(kde_expr) != 1L) stop("stage1_kde_at definition not found exactly once")
    eval(kde_expr[[1L]], envir = .GlobalEnv)
    set.seed(7)
    y <- rnorm(300L)
    x <- c(-1.2, 0, 1.2)
    got <- stage1_kde_at(y, x)
    n <- length(y)
    h <- 1.06 * sqrt(sum((y - mean(y))^2) / n) * n^(-0.2)
    want <- vapply(x, function(q) {
      mean(exp(-0.5 * ((y - q) / h)^2)) / (h * sqrt(2 * pi))
    }, numeric(1))
    if (max(abs(got - pmax(want, 1e-8))) > 1e-14) {
      stop("stage1_kde_at does not match the production formula")
    }
    message("STAGE1_KDE_FORMULA_OK")
  }
}
