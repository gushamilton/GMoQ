# pooled-memory-repair-20260915-01

The registered fixture run for the compact pooled covariance repair. Run
`stage2/tests/test-pooled-memory-repair.R` with the package sources on a BP
compute node after the source commit is frozen. The fixture supplies two
independent SPD Sigma_0 matrices and interleaves chr21/chr22 SNPs, so no Stage 1
reader output is reused as truth. The test asserts equality of the dense
compatibility and compact representations through calibrated OLS SEs, GLS,
Q-calibration, parameter covariance, and model-tournament kernels. It also
checks that one bad `d_resid` is skipped and that the compact payload contains
one upper triangle per chromosome plus one scale and offset per SNP.
