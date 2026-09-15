#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("pass at least one R source file")
for (path in args) {
  parse(file = path)
  message("R_PARSE_OK ", path)
}
