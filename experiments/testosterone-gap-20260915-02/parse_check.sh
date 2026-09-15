#!/usr/bin/env bash
set -euo pipefail
module load languages/R/4.4.3
Rscript --vanilla /user/work/fh6520/gmoq_finish_20260915/validity/testosterone-gap-20260915-02/parse_check.R \
  /user/work/fh6520/gmoq_finish_20260915/validity/testosterone-gap-20260915-02/prepare_matched_weights.R
