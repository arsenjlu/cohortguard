## 审稿意见对应的补充计算：
##   1) 真正的无标签 MMD（全基因集）
##   2) 核宽敏感性（median heuristic ×0.5 / ×1 / ×2）
##   3) 同队列内部 conformal 基线（交换性成立时的覆盖率）
## 结果写入 data/processed/revision_metrics.csv
suppressWarnings(source(file.path("R", "cohortguard.R"), encoding = "UTF-8"))

RAW <- file.path("data", "raw")
PROC <- file.path("data", "processed")

cfg <- read.csv(file.path("data", "datasets.csv"), encoding = "UTF-8",
                stringsAsFactors = FALSE, check.names = FALSE)
cfg$role <- tolower(trimws(cfg$role))

cg_find_gpl <- function(dir, platform = NA) {
  loc <- list.files(dir, pattern = "^GPL", full.names = TRUE)
  if (length(loc) == 1) return(loc)
  if (!is.na(platform)) {
    m <- loc[startsWith(basename(loc), platform)]
    if (length(m)) return(m[1])
    sh <- file.path(dirname(dir), "_platform")
    if (dir.exists(sh)) {
      m2 <- list.files(sh, pattern = paste0("^", platform), full.names = TRUE)
      if (length(m2)) return(m2[1])
    }
  }
  if (length(loc)) return(loc[1])
  stop("no GPL: ", dir)
}

cg_build_dataset <- function(dir, group_field, case_re, ctrl_re,
                             tissue_field = NA, tissue_re = NA, platform = NA) {
  sm_file <- list.files(dir, pattern = "series_matrix.*\\.gz$", full.names = TRUE)
  sm <- cg_read_series_matrix(sm_file[1])
  map <- cg_read_gpl(cg_find_gpl(dir, platform))
  f <- sm$fields
  key <- names(f)[grepl(tolower(group_field), names(f), fixed = TRUE)][1]
  grp <- f[[key]]; n <- length(sm$samples)
  mask <- rep(TRUE, n)
  if (!is.na(tissue_field) && !is.na(tissue_re)) {
    tk <- names(f)[grepl(tolower(tissue_field), names(f), fixed = TRUE)][1]
    if (!is.na(tk)) mask <- mask & grepl(tissue_re, f[[tk]], ignore.case = TRUE)
  }
  lab <- ifelse(grepl(case_re, grp, ignore.case = TRUE), 1L,
                ifelse(grepl(ctrl_re, grp, ignore.case = TRUE), 0L, -1L))
  mask <- mask & lab >= 0
  idx <- which(mask)
  M <- cg_to_gene(sm$mat[, idx, drop = FALSE], sm$ids, map)
  list(genes = rownames(M), M = M, y = lab[idx])
}

## 同队列内部：fit 50% / cal 25% / test 25%，交换性成立
cg_within_conformal <- function(A, y, k = 50, seed = 7, alpha = 0.1) {
  set.seed(seed)
  pick <- function(idx) {
    s <- sample(idx); n <- length(s)
    a <- max(2, floor(n * 0.5)); b <- max(2, floor(n * 0.25))
    list(fit = s[seq_len(a)], cal = s[a + seq_len(b)], te = s[a + b + seq_len(n - a - b)])
  }
  p0 <- pick(which(y == 0)); p1 <- pick(which(y == 1))
  fit <- c(p0$fit, p1$fit); cal <- c(p0$cal, p1$cal); te <- c(p0$te, p1$te)
  sel <- order(-abs(cg_tstat(A[, fit, drop = FALSE], y[fit])))[seq_len(k)]
  b <- cg_fit_logit(t(A[sel, fit, drop = FALSE]), y[fit])
  cf <- cg_conformal(cg_predict(t(A[sel, cal, drop = FALSE]), b), y[cal],
                     cg_predict(t(A[sel, te, drop = FALSE]), b), y[te], alpha)
  c(coverage = cf$coverage, abstention = cf$abstention)
}

pairs <- unique(cfg$pair_id[cfg$role %in% c("train", "test", "validation")])
out <- list()
for (pid in pairs) {
  sub <- cfg[cfg$pair_id == pid, ]
  tr <- sub[sub$role == "train", ][1, ]
  te <- sub[sub$role %in% c("test", "validation"), ][1, ]
  if (is.na(tr$dataset_id) || is.na(te$dataset_id)) next
  cat("\n===", pid, "===\n")
  d1 <- cg_build_dataset(file.path(RAW, tr$dataset_id), tr$group_field, tr$case_re, tr$ctrl_re,
                         tr$tissue_field, tr$tissue_re, tr$platform)
  d2 <- cg_build_dataset(file.path(RAW, te$dataset_id), te$group_field, te$case_re, te$ctrl_re,
                         te$tissue_field, te$tissue_re, te$platform)
  common <- intersect(d1$genes, d2$genes)
  A <- cg_zscore(d1$M[common, , drop = FALSE])
  B <- cg_zscore(d2$M[common, , drop = FALSE])
  ok <- is.finite(rowSums(A)) & is.finite(rowSums(B))
  A <- A[ok, , drop = FALSE]; B <- B[ok, , drop = FALSE]

  sel <- order(-abs(cg_tstat(A, d1$y)))[seq_len(50)]
  Af <- A[sel, , drop = FALSE]; Bf <- B[sel, , drop = FALSE]

  ## label-free: 全部共同基因
  mmd_all <- cg_mmd_topvar(A, B, nv = 2000)
  ## 核宽敏感性（模型特征上）
  mmd_05 <- cg_mmd(Af, Bf, mult = 0.5)
  mmd_1  <- cg_mmd(Af, Bf, mult = 1)
  mmd_2  <- cg_mmd(Af, Bf, mult = 2)
  wc <- cg_within_conformal(A, d1$y, k = 50, seed = 7)

  out[[pid]] <- data.frame(pair = pid,
                           mmd_all = mmd_all, mmd_all_norm = mmd_all / (1 + mmd_all),
                           mmd_k05 = mmd_05, mmd_k1 = mmd_1, mmd_k2 = mmd_2,
                           cov_within = wc["coverage"], abst_within = wc["abstention"],
                           stringsAsFactors = FALSE)
  cat(sprintf("  MMD all-genes %.4f | k0.5 %.4f | k1 %.4f | k2 %.4f | within-cohort coverage %.2f (abstain %.0f%%)\n",
              mmd_all, mmd_05, mmd_1, mmd_2, wc["coverage"], 100 * wc["abstention"]))
}

res <- do.call(rbind, out); rownames(res) <- NULL
write.csv(res, file.path(PROC, "revision_metrics.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## 与主表合并后计算各种 R²
main <- read.csv(file.path(PROC, "benchmark_R.csv"), stringsAsFactors = FALSE)
m <- merge(main, res, by = "pair")
r2 <- function(x, y) { f <- lm(y ~ x); summary(f)$r.squared }
cat("\n--- R2 with drop (n =", nrow(m), ") ---\n")
for (v in c("mmd_all_norm", "mmd_k05", "mmd_k1", "mmd_k2", "mmd_norm", "cv_auc", "agree_delta")) {
  cat(sprintf("  %-12s R2 = %.3f\n", v, r2(m[[v]], m$drop)))
}
m15 <- m[m$pair != "P1r", ]
cat("\n--- excluding P1r (n =", nrow(m15), ") ---\n")
for (v in c("mmd_all_norm", "mmd_norm", "cv_auc", "agree_delta")) {
  f <- lm(drop ~ m15[[v]]); s <- summary(f)
  cat(sprintf("  %-12s R2 = %.3f  p = %.4f\n", v, s$r.squared,
              s$coefficients[2, 4]))
}
write.csv(m, file.path(PROC, "benchmark_merged.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("\nwritten revision_metrics.csv / benchmark_merged.csv\n")
