## 补充统计：n=15（剔除 P1r）、R2 的 bootstrap 置信区间、组间对比
main <- read.csv("data/processed/benchmark_R.csv", stringsAsFactors = FALSE)
rev  <- read.csv("data/processed/revision_metrics.csv", stringsAsFactors = FALSE)
m <- merge(main, rev, by = "pair")
write.csv(m, "data/processed/benchmark_merged.csv", row.names = FALSE, fileEncoding = "UTF-8")

r2fun <- function(x, y) summary(lm(y ~ x))$r.squared
boot_r2 <- function(x, y, B = 2000, seed = 1) {
  set.seed(seed); n <- length(x); out <- numeric(B)
  for (b in seq_len(B)) {
    i <- sample.int(n, n, replace = TRUE)
    if (length(unique(x[i])) < 3) next
    out[b] <- r2fun(x[i], y[i])
  }
  quantile(out, c(0.025, 0.975), na.rm = TRUE)
}

vars <- c("mmd_all_norm", "mmd_k05", "mmd_k1", "mmd_k2", "mmd_norm", "cv_auc", "agree_delta")
res <- data.frame()
for (v in vars) {
  x <- m[[v]]; y <- m$drop
  f <- lm(y ~ x); s <- summary(f)
  ci <- boot_r2(x, y)
  res <- rbind(res, data.frame(var = v, n = nrow(m), r2 = round(s$r.squared, 3),
                               p = signif(s$coefficients[2, 4], 3),
                               lo = round(ci[1], 3), hi = round(ci[2], 3)))
}
cat("=== n = 16 ===\n"); print(res, row.names = FALSE)

m15 <- m[m$pair != "P1r", ]
res15 <- data.frame()
for (v in vars) {
  x <- m15[[v]]; y <- m15$drop
  s <- summary(lm(y ~ x)); ci <- boot_r2(x, y)
  res15 <- rbind(res15, data.frame(var = v, n = nrow(m15), r2 = round(s$r.squared, 3),
                                   p = signif(s$coefficients[2, 4], 3),
                                   lo = round(ci[1], 3), hi = round(ci[2], 3)))
}
cat("\n=== n = 15 (P1r removed) ===\n"); print(res15, row.names = FALSE)

cat("\n=== conformal: within-cohort vs cross-cohort ===\n")
cat(sprintf("within-cohort coverage: median %.2f (range %.2f-%.2f); abstention median %.0f%%\n",
            median(m$cov_within), min(m$cov_within), max(m$cov_within), 100 * median(m$abst_within)))
cat(sprintf("cross-cohort coverage: median %.2f (range %.2f-%.2f); abstention median %.0f%%\n",
            median(m$cov90), min(m$cov90), max(m$cov90), 100 * median(m$abstain90)))
cat(sprintf("coverage >= 0.90:  within %d/16, cross %d/16\n",
            sum(m$cov_within >= 0.9), sum(m$cov90 >= 0.9)))
cat(sprintf("coverage >= 0.85:  within %d/16, cross %d/16\n",
            sum(m$cov_within >= 0.85), sum(m$cov90 >= 0.85)))

cat("\n=== group medians (matched vs mismatched) ===\n")
mm <- c("P2","P3","P4","P5","P6","P9","P10","P11","P12","P13","P15")
cat(sprintf("matched    n=%d median drop %.3f\n", sum(m$pair %in% mm), median(m$drop[m$pair %in% mm])))
cat(sprintf("mismatched n=%d median drop %.3f\n", sum(!m$pair %in% mm), median(m$drop[!m$pair %in% mm])))

cat("\n=== external AUC thresholds ===\n")
cat(sprintf("<0.6: %d (%s)\n", sum(m$external_auc < 0.6), paste(m$pair[m$external_auc < 0.6], collapse = ", ")))
cat(sprintf(">0.98: %d (%s)\n", sum(m$external_auc > 0.98), paste(m$pair[m$external_auc > 0.98], collapse = ", ")))

write.csv(res, "data/processed/r2_table.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(res15, "data/processed/r2_table_n15.csv", row.names = FALSE, fileEncoding = "UTF-8")
