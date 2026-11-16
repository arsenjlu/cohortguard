## M12 敏感性：只用训练队列的均值/方差去标准化目标队列
## 用法：Rscript scripts_r/m12_sensitivity.R
suppressWarnings(source(file.path("R", "cohortguard.R"), encoding = "UTF-8"))

RAW <- file.path("data", "raw"); PROC <- file.path("data", "processed")

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
  stop("找不到 GPL 注释: ", dir)
}

cg_build_dataset <- function(dir, group_field, case_re, ctrl_re,
                             tissue_field = NA, tissue_re = NA, platform = NA) {
  sm_file <- list.files(dir, pattern = "series_matrix.*\\.gz$", full.names = TRUE)
  if (!length(sm_file)) stop("找不到 series_matrix: ", dir)
  sm <- cg_read_series_matrix(sm_file[1])
  map <- cg_read_gpl(cg_find_gpl(dir, platform))
  f <- sm$fields
  key <- names(f)[grepl(tolower(group_field), names(f), fixed = TRUE)][1]
  if (is.na(key)) stop("找不到分组字段 ", group_field)
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
  list(genes = rownames(M), M = M, y = lab[idx], n_used = length(idx))
}

## 只用训练统计量标准化的 zscore：返回 (A_std, B_std)
train_stats_zscore <- function(Ar, Br) {
  mu <- rowMeans(Ar, na.rm = TRUE)
  sd <- apply(Ar, 1, sd, na.rm = TRUE)
  sd[!is.finite(sd) | sd == 0] <- NA
  A <- (Ar - mu) / sd
  B <- (Br - mu) / sd   # 目标队列用训练队列的 mu/sd
  list(A = A, B = B)
}

pairs <- unique(cfg$pair_id[cfg$role %in% c("train", "test", "validation")])
res <- list()
for (pid in pairs) {
  sub <- cfg[cfg$pair_id == pid, ]
  tr <- sub[sub$role == "train", ][1, ]
  te <- sub[sub$role %in% c("test", "validation"), ][1, ]
  if (is.na(tr$dataset_id) || is.na(te$dataset_id)) next
  d1 <- cg_build_dataset(file.path(RAW, tr$dataset_id), tr$group_field, tr$case_re, tr$ctrl_re,
                         tr$tissue_field, tr$tissue_re, tr$platform)
  d2 <- cg_build_dataset(file.path(RAW, te$dataset_id), te$group_field, te$case_re, te$ctrl_re,
                         te$tissue_field, te$tissue_re, te$platform)
  common <- intersect(d1$genes, d2$genes)
  Ar <- d1$M[common, , drop = FALSE]; Br <- d2$M[common, , drop = FALSE]
  z <- train_stats_zscore(Ar, Br)
  A <- z$A; B <- z$B
  ok <- is.finite(rowSums(A)) & is.finite(rowSums(B))
  A <- A[ok, , drop = FALSE]; B <- B[ok, , drop = FALSE]
  r <- cg_evaluate(A, d1$y, B, d2$y, k = 50, seed = 7)
  mmd_lf <- cg_mmd_topvar(A, B, nv = 2000); mmd_lf <- mmd_lf / (1 + mmd_lf)
  res[[pid]] <- data.frame(pair = pid, n_genes = nrow(A),
                           cv_auc = round(r$cv, 4), external_auc = round(r$ext, 4),
                           drop = round(r$drop, 4),
                           mmd_feature = round(r$mmd_norm, 4), mmd_labelfree = round(mmd_lf, 4),
                           agree_delta = round(r$agree_delta, 3),
                           cov90 = round(r$cov90, 4), abstain90 = round(r$abstain90, 4))
  cat(sprintf("%s  drop=%.3f  ext=%.3f  cv=%.3f  mmd_lf=%.4f\n",
              pid, r$drop, r$ext, r$cv, mmd_lf))
}
out <- do.call(rbind, res); rownames(out) <- NULL
write.csv(out, file.path(PROC, "m12_sensitivity.csv"), row.names = FALSE, fileEncoding = "UTF-8")

## 对比：与标准（队列内各自 zscore）的掉分差异
std <- read.csv(file.path(PROC, "benchmark_R.csv"), stringsAsFactors = FALSE)
m <- merge(std, out, by = "pair", suffixes = c("_std", "_ts"))
m$drop_diff <- m$drop_ts - m$drop_std
cat("\n=== 掉分差异（训练统计量 − 队列内各自）===\n")
print(m[, c("pair", "drop_std", "drop_ts", "drop_diff")], row.names = FALSE)
cat(sprintf("\n掉分差异：均值 %.3f，最大绝对差异 %.3f\n",
            mean(m$drop_diff), max(abs(m$drop_diff))))
cat(sprintf("external~CV (训练统计量): R²=%.3f\n",
            summary(lm(external_auc ~ cv_auc, data = out))$r.squared))
cat(sprintf("label-free MMD ~ drop (训练统计量): R²=%.3f\n",
            summary(lm(drop ~ mmd_labelfree, data = out))$r.squared))
