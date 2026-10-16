# Load the INSTALLED fungwasStage2 (including internal functions) into the global
# environment, so test scripts exercise the package as built rather than
# recompiling src/ with Rcpp::sourceCpp().
suppressPackageStartupMessages(library(fungwasStage2))
local({
  ns <- asNamespace("fungwasStage2")
  for (nm in ls(ns, all.names = TRUE)) assign(nm, get(nm, envir = ns), envir = globalenv())
})
