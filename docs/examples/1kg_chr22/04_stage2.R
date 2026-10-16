# Stage 2: map each phenotype's 45 quantile effects onto model parameters.
#   null, shift, variance -> Normal mean/variance W (beta_mu, beta_sigma2)
#   mixture               -> two-component Normal mixture W with membership (gamma, beta_1, beta_2),
#                            fitted twice: density = "hybrid" (the default since 0.6.0, written
#                            as mixture_hybrid) and density = "model" (pre-0.6.0, written as mixture)
# Each fitted with pooled (primary) and jackknife covariance.
suppressPackageStartupMessages({ library(data.table); library(fungwasStage2) })
args <- commandArgs(trailingOnly = TRUE)
work <- if (length(args)) args[1] else "work"

pheno <- fread(file.path(work, "pheno.txt"))
taus <- unlist(jsonlite::fromJSON(file.path(work, "s1_null.sigma0.meta.json"))$taus)
p1 <- as.numeric(sub("p1=", "", readLines(file.path(work, "mixture_baseline.txt"))))

weights_for <- function(fit) {
  ph <- sub("_hybrid", "", fit)
  y <- pheno[[ph]]
  q <- as.numeric(quantile(y, taus, type = 8))
  if (ph == "mixture") {
    # Baseline = the simulated components; component 1 is N(0,1) with share p1.
    # The default ("hybrid") divides the gamma column by the Stage 1 kernel density.
    if (grepl("_hybrid", fit)) {
      make_weights_normal_mixture(taus, q, p1 = p1, mu1 = 0, sd1 = 1, mu2 = 4, sd2 = 1,
                                  include_membership = TRUE,
                                  stage1_density = stage1_kde_density(y, q))
    } else {
      make_weights_normal_mixture(taus, q, p1 = p1, mu1 = 0, sd1 = 1, mu2 = 4, sd2 = 1,
                                  include_membership = TRUE, density = "model")
    }
  } else {
    make_weights_vqtl(taus, q, mu = mean(y), sd = sd(y))
  }
}

out <- list()
for (ph in c("null", "shift", "variance", "mixture", "mixture_hybrid")) {
  W <- weights_for(ph)
  pre <- file.path(work, paste0("s1_", sub("_hybrid", "", ph)))
  for (src in c("pooled", "jackknife")) {
    res <- param_gwas_from_file(
      stage1_file = paste0(pre, ".stage1.tsv.gz"), W = W,
      se_mode = "tau_cov", tau_cov_estimator = "ols", cov_source = src,
      sigma0_file = if (src == "pooled") paste0(pre, ".sigma0.tsv.gz"),
      cov_file = if (src == "jackknife") paste0(pre, ".cov.gz"))
    res <- as.data.table(res)
    fwrite(res, file.path(work, sprintf("s2_%s_%s.tsv.gz", ph, src)), sep = "\t")
    out[[paste(ph, src)]] <- res
    cat(sprintf("stage2 %-8s %-9s %d SNPs\n", ph, src, nrow(res)))
  }
}
