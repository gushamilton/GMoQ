#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Tier 0 -- the reference distribution for the parameter Wald test.
#
# theta/SE is t_nu, not normal, because Sigma is estimated. Reading it against a
# normal inflates Type-I error 17x at p = 5e-8 when B = 65
# (docs/KNOWN_ISSUES.md section 2). Everything that decides nu is closed form
# and is checked here, including the ladder of fallbacks and the refusals --
# a silent fallback is precisely how the original defect survived.
#
# Run:  Rscript --vanilla tests/tier0/t0_tailref.R
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  .self <- local({
    a <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (!length(f)) stop("Run as: Rscript --vanilla tests/tier0/t0_tailref.R")
    normalizePath(f[1L])
  })
  source(file.path(dirname(.self), "..", "lib", "harness.R"))
}))
gmoq_init_roots(.self)
gmoq_load_stage2()
gmoq_begin("tier0", "tailref")

eff_nu    <- gmoq_internal(".effective_nu")
flatness  <- gmoq_internal(".contrast_flatness")
contrastc <- gmoq_internal(".contrast_c")
param_ref <- gmoq_internal(".param_reference")
param_p   <- gmoq_internal(".param_pvalues")
resolve_b <- gmoq_internal(".resolve_block_count")
resolve_m <- gmoq_internal(".resolve_df_model")
resolve_n <- gmoq_internal(".resolve_nu")
read_meta <- gmoq_internal(".read_cov_meta")
jnum      <- gmoq_internal(".json_number")
jstr      <- gmoq_internal(".json_string")
jsarr     <- gmoq_internal(".json_string_array")
jnarr     <- gmoq_internal(".json_number_array")
C_LOCATION_DEFAULT     <- gmoq_internal("C_LOCATION_DEFAULT")
C_SCALE_DEFAULT        <- gmoq_internal("C_SCALE_DEFAULT")
FALLBACK_N_EFF_DEFAULT <- gmoq_internal("FALLBACK_N_EFF_DEFAULT")
N_EFF_COLUMN           <- gmoq_internal("N_EFF_COLUMN")
COV_META_VERSION_MAX   <- gmoq_internal("COV_META_VERSION_MAX")

# =====================================================================
# 1. The df law, written out
#      1/nu = 1/(B - 1) + c / N_eff
# =====================================================================
B <- 1000; ne <- c(9134, 50000, 1e6)
gmoq_close("T0.DF.law_location",
           "The effective df is 1/nu = 1/(B-1) + c/N_eff with c = 0.88 for a location contrast.",
           eff_nu(B, ne, 0.88), 1 / (1 / (B - 1) + 0.88 / ne), 1e-10)
gmoq_close("T0.DF.law_scale",
           "The same law with c = 2.62 for a scale contrast.",
           eff_nu(B, ne, 2.62), 1 / (1 / (B - 1) + 2.62 / ne), 1e-10)
gmoq_check("T0.DF.nu_below_blocks",
           "nu is strictly below B - 1, which is the measured fact that rules out a B-1 default.",
           all(eff_nu(B, ne, 0.88) < B - 1),
           stat = max(eff_nu(B, ne, 0.88)), criterion = "nu < B-1 for every finite N_eff")
gmoq_check("T0.DF.scale_nu_below_location_nu",
           "A scale contrast gets fewer degrees of freedom than a location one, as measured.",
           all(eff_nu(B, ne, 2.62) < eff_nu(B, ne, 0.88)),
           criterion = "nu_scale < nu_location at every N_eff")

# The shipped constants must be the calibrated ones.
gmoq_check("T0.DF.constants",
           "The shipped Satterthwaite coefficients are the calibrated 0.88 / 2.62 and N_eff fallback 9134.",
           C_LOCATION_DEFAULT == 0.88 && C_SCALE_DEFAULT == 2.62 &&
             FALLBACK_N_EFF_DEFAULT == 9134 && N_EFF_COLUMN == "d_resid",
           criterion = "c_location=0.88, c_scale=2.62, fallback_n_eff=9134, column=d_resid")

# =====================================================================
# 2. The 17x claim itself, as arithmetic
# =====================================================================
z58 <- qnorm(5e-8 / 2, lower.tail = FALSE)
ratio65 <- 2 * pt(-z58, df = 64) / 5e-8
gmoq_check("T0.DF.normal_reference_cost",
           "Reading theta/SE against a normal at B = 65 costs ~17x at 5e-8, the documented figure.",
           abs(ratio65 - 17.1) < 1.0, stat = ratio65,
           criterion = "2*pt(-5.4513, 64)/5e-8 within 1.0 of the documented 17.1")

# =====================================================================
# 3. Contrast classification
# =====================================================================
Tt <- 45L; taus <- seq(0.02, 0.98, length.out = Tt)
Wls <- cbind(mu = 1, q = taus - mean(taus))
A <- solve(crossprod(Wls), t(Wls))
gmoq_check("T0.DF.flatness_flat",
           "A constant contrast has flatness 1; a mean-zero one has flatness 0.",
           abs(flatness(rep(1, Tt)) - 1) < 1e-12 &&
             abs(flatness(taus - mean(taus))) < 1e-12,
           criterion = "flatness(1) == 1 and flatness(centred) == 0")
gmoq_check("T0.DF.contrast_c_selection",
           "The location row of A gets c = 0.88 and the scale row gets c = 2.62.",
           contrastc(A[1, ]) == 0.88 && contrastc(A[2, ]) == 2.62,
           criterion = "A[1,] -> c_location, A[2,] -> c_scale")
gmoq_check("T0.DF.ambiguous_contrast_is_conservative",
           "A contrast that is not clearly flat gets the scale coefficient, the conservative one.",
           contrastc(c(1, rep(0, Tt - 1))) == 2.62,
           criterion = "an unclassifiable contrast -> c_scale")

# =====================================================================
# 4. The p-value itself
# =====================================================================
theta <- matrix(c(0.1, -0.2), 2, 1); se <- matrix(c(0.02, 0.05), 2, 1)
ref_t <- list(dist = "t", df = matrix(c(800, 560), 2, 1))
gmoq_close("T0.DF.pvalue_is_two_sided_t",
           "The parameter p-value is 2 pt(-|theta/SE|, nu).",
           param_p(theta, se, ref_t),
           2 * pt(-abs(c(0.1 / 0.02, -0.2 / 0.05)), df = c(800, 560)), 1e-14)
gmoq_check("T0.DF.pvalue_none_returns_null",
           "When no reference could be established the p-value is NULL, not a normal one.",
           is.null(param_p(theta, se, list(dist = "none"))),
           criterion = "dist='none' -> NULL")

# =====================================================================
# 5. The ladder of fallbacks, rung by rung
# =====================================================================
P <- 4L
n_eff <- c(9134, 50000, 200000, NA_real_)

r1 <- param_ref(A, P, nu = 777, nu_source = "pinned")
gmoq_check("T0.DF.ladder1_explicit_nu",
           "Rung 1: an explicit nu wins and its provenance is recorded.",
           r1$dist == "t" && all(r1$df == 777) && r1$source == "pinned",
           criterion = "dist=t, df=777, source recorded")

r2 <- suppressWarnings(param_ref(A, P, n_blocks = 1000, n_eff = n_eff))
gmoq_check("T0.DF.ladder2_per_snp",
           "Rung 2: block count plus per-SNP d_resid gives one nu per SNP per parameter.",
           r2$source == "per_snp_satterthwaite" && dim(r2$df)[1] == 2L && dim(r2$df)[2] == P &&
             abs(r2$df[1, 2] - eff_nu(1000, 50000, 0.88)) < 1e-8 &&
             abs(r2$df[2, 2] - eff_nu(1000, 50000, 2.62)) < 1e-8,
           criterion = "per-SNP nu matches the law, K x P shape")
gmoq_check("T0.DF.ladder2_fills_missing_n_eff",
           "A SNP with no usable d_resid falls back to fallback_n_eff, and that is warned about.",
           abs(r2$df[1, 4] - eff_nu(1000, FALLBACK_N_EFF_DEFAULT, 0.88)) < 1e-8,
           criterion = "missing n_eff -> fallback_n_eff")

r3 <- suppressWarnings(param_ref(A, P, n_blocks = 1000, n_eff = NULL))
gmoq_check("T0.DF.ladder3_flat_fallback",
           "Rung 3: block count alone gives a flat df at fallback_n_eff, recorded as such.",
           r3$source == "flat_fallback_n_eff" &&
             abs(r3$df[1, 1] - eff_nu(1000, FALLBACK_N_EFF_DEFAULT, 0.88)) < 1e-8,
           criterion = "source == flat_fallback_n_eff")

r5 <- suppressWarnings(param_ref(A, P, unknown_nu = "omit"))
gmoq_check("T0.DF.ladder5_omit",
           "With nothing usable the default is to emit NO p-value, not a normal one.",
           r5$dist == "none" && is.null(r5$df),
           criterion = "unknown_nu='omit' -> dist='none'")
r5n <- suppressWarnings(param_ref(A, P, unknown_nu = "normal"))
gmoq_check("T0.DF.ladder5_normal_is_recorded",
           "unknown_nu='normal' is available but records ref_dist='normal' in the output.",
           r5n$dist == "normal",
           criterion = "dist == 'normal'")
gmoq_refuses("T0.DF.ladder5_error",
             "unknown_nu='error' stops rather than guessing a reference.",
             param_ref(A, P, unknown_nu = "error"),
             "Degrees of freedom of the Stage 1 covariance are unknown")

gmoq_check("T0.DF.no_silent_normal",
           "There is no path that emits a normal-reference p-value without saying so.",
           r5$dist != "normal",
           criterion = "the default when nu is unknown is not 'normal'")

# df_coef override
r_ov <- suppressWarnings(param_ref(A, P, n_blocks = 1000, n_eff = n_eff,
                                   df_coef = c(mu = 1.5)))
gmoq_check("T0.DF.df_coef_override",
           "df_coef overrides c per named parameter, for a W the constants were not calibrated on.",
           abs(r_ov$c[["mu"]] - 1.5) < 1e-12 && abs(r_ov$c[["q"]] - 2.62) < 1e-12,
           criterion = "named override applies to that parameter only")

# =====================================================================
# 6. Reading the Stage 1 sidecar
# =====================================================================
tmp <- file.path(tempdir(), paste0("t0meta_", as.integer(runif(1, 1, 1e6))))
dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

meta_json <- function(version = 4, n_blocks = 1000, nonempty = 998) sprintf('{
  "cov_meta_version": %d,
  "cov_dtype": "float32",
  "n_taus": 45,
  "n_cov_elements": 1035,
  "record_bytes": 4140,
  "scale_file": false,
  "jackknife": {
    "method": "delete_one_block",
    "n_blocks": %d,
    "n_blocks_nonempty": %d,
    "df_model": {
      "form": "1/nu = 1/(n_blocks_nonempty - 1) + c / n_eff",
      "n_eff_column": "d_resid",
      "c_location": 0.88,
      "c_scale": 2.62,
      "fallback_n_eff": 9134,
      "source": "test"
    },
    "seed": 42
  }
}', version, n_blocks, nonempty)

p_ok <- file.path(tmp, "x.cov.meta.json")
writeLines(meta_json(), p_ok)
m <- read_meta(p_ok)
gmoq_check("T0.DF.meta_prefers_nonempty_blocks",
           "The block count used is n_blocks_nonempty: an empty block contributes no replicate.",
           resolve_b(m) == 998,
           stat = resolve_b(m), criterion = "n_blocks_nonempty wins over n_blocks")
dm <- resolve_m(m)
gmoq_check("T0.DF.meta_df_model",
           "The df law is read from the sidecar rather than assumed.",
           dm$c_location == 0.88 && dm$c_scale == 2.62 &&
             dm$fallback_n_eff == 9134 && dm$n_eff_column == "d_resid",
           criterion = "coefficients round-trip from JSON")
gmoq_check("T0.DF.meta_no_nu_by_default",
           "A sidecar with no explicit nu yields no constant nu -- there is no B-1 fallback.",
           !is.finite(resolve_n(m)$nu),
           criterion = "resolve_nu returns NA when 'nu' is absent")

p_future <- file.path(tmp, "future.cov.meta.json")
writeLines(meta_json(version = 99), p_future)
gmoq_refuses("T0.DF.meta_future_schema_refused",
             "A sidecar schema newer than this package understands is refused, not guessed at.",
             read_meta(p_future), "does not.*understand|schema version")
gmoq_check("T0.DF.meta_version_max",
           "The declared maximum understood schema version is 4.",
           COV_META_VERSION_MAX == 4L, stat = COV_META_VERSION_MAX,
           criterion = "COV_META_VERSION_MAX == 4")

gmoq_check("T0.DF.meta_missing_is_null",
           "A missing sidecar returns NULL rather than raising, so the ladder can continue.",
           is.null(read_meta(file.path(tmp, "nope.cov.meta.json"))),
           criterion = "absent file -> NULL")

# The base-R JSON helpers, used whenever jsonlite is not installed.
txt <- paste(readLines(p_ok, warn = FALSE), collapse = " ")
gmoq_check("T0.DF.json_helpers",
           "The dependency-free JSON readers recover numbers, strings and arrays correctly.",
           jnum(txt, "n_blocks") == 1000 && jstr(txt, "n_eff_column") == "d_resid" &&
             identical(jsarr('{"covariates": ["age","sex","PC1"]}', "covariates"),
                       c("age", "sex", "PC1")) &&
             identical(jnarr('{"taus": [0.02, 0.5, 0.98]}', "taus"), c(0.02, 0.5, 0.98)),
           criterion = "number, string and array extraction all correct")

gmoq_finish()
