# Final package verification

Source: `ccc2d9712ac3211e369fcddec5c14c0f7a526d2b`, preserved by
`experiments/2026-09-15/final-package-source`.

The source archive was built with vignettes, checked, installed into a fresh
library and tested on `bp1-compute019`, allocation `18966618`. The runner
exited successfully. The allocation was then released and its status confirmed
as `no live allocation`.

Results:

- R CMD check: zero errors and notes; one existing GNU Makevars portability
  warning. Examples and vignette rebuilding passed. Tests were run separately.
- `test-multi-from-file-unsafe.R`: passed, including raw covariance-ID
  validation before duplicate filtering, a nonzero retained covariance offset,
  known model scores, and refusal of mismatched IDs.
- `test-pooled-covariance.R`: passed, including the explicit jackknife path,
  pooled/file metadata checks, chromosome matching and caller assertions.
- `test-cov-ids-mandatory.R`: passed.

Exact command, launched through the held allocation:

```sh
bprun --cpus 2 --mem 8G bash \
  /user/work/fh6520/gmoq_finish_20260915/root/final-source/scripts/validation/final_package_check.sh \
  /user/work/fh6520/gmoq_finish_20260915/root/final-source \
  /user/work/fh6520/gmoq_finish_20260915/root/final-package-check-20260915-01 \
  ccc2d9712ac3211e369fcddec5c14c0f7a526d2b \
  test-multi-from-file-unsafe.R test-pooled-covariance.R test-cov-ids-mandatory.R
```

The archived source package and installed library remain in the remote output
directory. The logs and environment are preserved here. This final check
complements the earlier full Stage 1, mathematical, integration, ten-script
Stage 2 and compact/dense comparison results; it does not repeat the simulation
campaign or claim universal statistical calibration.

Earlier build troubleshooting belongs to the pooled-memory run: a check first
referenced an archive in the wrong directory and was skipped, and a check with
vignettes disabled produced two avoidable missing-vignette warnings. Neither
was counted as the final standard build result. The standard build was then
verified, followed by this final archive after the alignment fix.
