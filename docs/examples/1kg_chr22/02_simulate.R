# Simulate four phenotypes with known genetic effects on the 1000 Genomes chr22 genotypes.
#
#   null     Y = 0.3*PC1 + e                               no SNP effects
#   shift    Y = 0.3*PC1 + b*G + e,           b = 0.25     truth: beta_mu = b,  beta_sigma2 = 0
#   variance Y = 0.3*PC1 + sqrt(1 + d*(G - mean G)) * e, d = 0.30   truth: beta_mu = 0, beta_sigma2 = +/-d
#   mixture  Y = N(0,1) w.p. p1(G), else N(4,1),
#            logit p1 = a + g*(G - mean G), g = 0.50       truth: gamma = g, beta_1 = beta_2 = 0
#
# G = dosage of BGEN allele 1 (= REF, the effect allele Stage 1 reports). e ~ N(0,1).
# Causal SNPs (MAF 0.2-0.5): 5 shift, 3 variance, 5 mixture, all distinct. A SNP counts as null
# for a model when its r^2 with every causal SNP of that model is < 0.01 (written to ld_to_causal.tsv).
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
work <- if (length(args)) args[1] else "work"
seed <- if (length(args) > 1) as.integer(args[2]) else 2026L
set.seed(seed)

raw <- fread(file.path(work, "geno.raw"))
pcs <- fread(file.path(work, "pca.eigenvec"))
setnames(pcs, 1, "IID")
stopifnot(identical(raw$IID, pcs$IID))

# Dosage of allele 1 (REF). Column names are <ID>_<counted allele>; ID is chr:pos:REF:ALT.
gcols <- names(raw)[-(1:6)]
id <- sub("_[ACGT]+$", "", gcols)
counted <- sub("^.*_", "", gcols)
ref <- vapply(strsplit(id, ":"), `[`, "", 3)
G <- as.matrix(raw[, ..gcols])
flip <- counted != ref
G[, flip] <- 2 - G[, flip]
colnames(G) <- id
G[is.na(G)] <- 0  # phase 3 has no missing calls; guard only
maf <- pmin(colMeans(G) / 2, 1 - colMeans(G) / 2)

cand <- which(maf >= 0.2 & maf <= 0.5)
pick <- sample(cand, 13)
causal <- list(shift = pick[1:5], variance = pick[6:8], mixture = pick[9:13])
N <- nrow(G)
pc1 <- as.numeric(scale(pcs$PC1))
gsum <- function(cols, eff) as.numeric(G[, cols, drop = FALSE] %*% rep(eff, length(cols)))

b <- 0.25; d <- 0.30; g <- 0.50
null <- 0.3 * pc1 + rnorm(N)
shift <- 0.3 * pc1 + gsum(causal$shift, b) + rnorm(N)
Gv <- G[, causal$variance, drop = FALSE]
# Variance rises with the MINOR allele (keeps v > 0), so the truth on the allele-1 scale is +/- d.
d_sign <- ifelse(colMeans(Gv) / 2 <= 0.5, 1, -1)
v <- 1 + as.numeric(sweep(Gv, 2, colMeans(Gv)) %*% (d * d_sign))
stopifnot(all(v > 0))
variance <- 0.3 * pc1 + sqrt(v) * rnorm(N)
Gm <- G[, causal$mixture]
eta <- qlogis(0.7) + as.numeric(sweep(Gm, 2, colMeans(Gm)) %*% rep(g, 5))
z1 <- rbinom(N, 1, plogis(eta)) == 1
mixture <- ifelse(z1, rnorm(N, 0, 1), rnorm(N, 4, 1))

# FID must equal the BGEN sample ID (Stage 1 matches column 1 only).
pheno <- data.table(FID = raw$IID, IID = raw$IID, null, shift, variance, mixture)
covar <- data.table(FID = raw$IID, IID = raw$IID, pcs[, paste0("PC", 1:5), with = FALSE])
fwrite(pheno, file.path(work, "pheno.txt"), sep = " ")
fwrite(covar, file.path(work, "covar.txt"), sep = " ")

truth <- rbindlist(list(
  data.table(model = "shift",    snp_id = id[causal$shift],    param = "beta_mu",     truth = b),
  data.table(model = "variance", snp_id = id[causal$variance], param = "beta_sigma2", truth = d * d_sign),
  data.table(model = "mixture",  snp_id = id[causal$mixture],  param = "gamma",       truth = g)
))
truth[, maf := maf[match(snp_id, id)]]
fwrite(truth, file.path(work, "truth.tsv"), sep = "\t")
r2 <- data.table(snp_id = id)
for (m in names(causal)) r2[[m]] <- apply(cor(G, G[, causal[[m]], drop = FALSE])^2, 1, max)
fwrite(r2, file.path(work, "ld_to_causal.tsv"), sep = "\t")
# Baseline mixture parameters used to build W in 04_stage2.R.
writeLines(sprintf("p1=%.4f", mean(z1)), file.path(work, "mixture_baseline.txt"))
cat(sprintf("Simulated N = %d samples, %d SNPs, seed %d; %d causal SNPs written to truth.tsv\n",
            N, ncol(G), seed, nrow(truth)))
