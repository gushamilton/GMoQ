d <- readRDS("/user/work/fh6520/gmoq_projeff/out/diag.rds")
options(width = 200)
k <- d[grepl("normal-mix 2col|normal-mix 3col|skew-mix 2col|skew-mix 3col|REAL", d$model), ]
k <- k[, c("model","param","se_ratio","R2","rho_max","bound_rho")]
k[, 3:6] <- round(k[, 3:6], 4)
print(k, row.names = FALSE)
cat("\n--- rho_max vs penalty, all families pooled (is rho_max a usable rule?) ---\n")
a <- d[d$algebra_check < 1e-6, ]
a$worst <- ave(a$se_ratio, a$model, FUN = max)
u <- unique(a[, c("model","rho_max","worst","bound_rho")])
u <- u[order(u$rho_max), ]
u[, 2:4] <- round(u[, 2:4], 4)
print(u, row.names = FALSE)
cat("\nspearman(rho_max, worst se_ratio) =",
    round(cor(u$rho_max, u$worst, method = "spearman"), 4), "\n")
cat("bound never violated:", all(u$worst <= u$bound_rho + 1e-9), "\n")
