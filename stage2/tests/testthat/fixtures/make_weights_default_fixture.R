#!/usr/bin/env Rscript
# Record the weight builders' output under a PREVIOUS package install, so the
# current default can be checked to be identical(). Run on a compute node:
#   Rscript --vanilla make_weights_default_fixture.R <old_lib> <out.rds> <commit>
a <- commandArgs(trailingOnly = TRUE)
library(fungwasStage2, lib.loc = a[1])
self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(normalizePath(self)), "..", "helper-weight-cases.R"))
out <- weight_cases_default()
attr(out, "recorded_with") <- list(commit = a[3], lib = a[1],
                                   version = as.character(packageVersion("fungwasStage2", lib.loc = a[1])),
                                   R = R.version.string)
saveRDS(out, a[2])
cat("wrote", a[2], "\n")
