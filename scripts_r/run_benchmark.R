## 驱动脚本：读 data/datasets.csv，构建所有配对并评估
## 用法：Rscript R/run_benchmark.R
suppressWarnings(source(file.path("R", "cohortguard.R"), encoding = "UTF-8"))

RAW <- file.path("data", "raw")
PROC <- file.path("data", "processed")
dir.create(PROC, showWarnings = FALSE, recursive = TRUE)

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
  if (is.na(key)) stop("找不到分组字段 ", group_field, "；可选: ", paste(names(f), collapse = ", "))
  grp <- f[[key]]; n <- length(sm$samples)
  if (length(grp) != n) stop("分组字段长度不符")

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
  list(genes = rownames(M), M = M, y = lab[idx], n_total = n,
       n_used = length(idx), probe2gene = length(map))
}

pairs <- unique(cfg$pair_id[cfg$role %in% c("train", "test", "validation")])
done <- character(0)
res <- list()
bf <- file.path(PROC, "benchmark_R.csv")
if (file.exists(bf)) {
  old <- read.csv(bf, stringsAsFactors = FALSE)
  done <- old$pair
  for (i in seq_len(nrow(old))) res[[old$pair[i]]] <- old[i, ]
  cat("已有结果:", paste(done, collapse = ", "), "\n")
}
for (pid in pairs) {
  if (pid %in% done) { cat("[已存在]", pid, "\n"); next }
  sub <- cfg[cfg$pair_id == pid, ]
  tr <- sub[sub$role == "train", ][1, ]
  te <- sub[sub$role %in% c("test", "validation"), ][1, ]
  if (is.na(tr$dataset_id) || is.na(te$dataset_id)) { cat("[跳过]", pid, "\n"); next }
  cat("\n===", pid, "===\n")
  d1 <- cg_build_dataset(file.path(RAW, tr$dataset_id), tr$group_field, tr$case_re, tr$ctrl_re,
                         tr$tissue_field, tr$tissue_re, tr$platform)
  cat(sprintf("  train %s: %d 基因 × %d 样本（病例 %d / 对照 %d）\n",
              tr$dataset_id, nrow(d1$M), ncol(d1$M), sum(d1$y == 1), sum(d1$y == 0)))
  d2 <- cg_build_dataset(file.path(RAW, te$dataset_id), te$group_field, te$case_re, te$ctrl_re,
                         te$tissue_field, te$tissue_re, te$platform)
  cat(sprintf("  test  %s: %d 基因 × %d 样本（病例 %d / 对照 %d）\n",
              te$dataset_id, nrow(d2$M), ncol(d2$M), sum(d2$y == 1), sum(d2$y == 0)))

  common <- intersect(d1$genes, d2$genes)
  A <- cg_zscore(d1$M[common, , drop = FALSE])
  B <- cg_zscore(d2$M[common, , drop = FALSE])
  ok <- is.finite(rowSums(A)) & is.finite(rowSums(B))
  A <- A[ok, , drop = FALSE]; B <- B[ok, , drop = FALSE]
  cat("  共同基因", nrow(A), "\n")

  r <- cg_evaluate(A, d1$y, B, d2$y, k = 50, seed = 7)
  cat(sprintf("  CV %.3f → 外部 %.3f（掉分 %.3f）| MMD %.4f | Δ一致性 %.3f | 90%%覆盖 %.2f | 弃权 %.0f%%\n",
              r$cv, r$ext, r$drop, r$mmd, r$agree_delta, r$cov90, r$abstain90 * 100))
  row <- data.frame(pair = pid, train = tr$dataset_id, test = te$dataset_id,
                           n_genes = nrow(A), n_train = ncol(A), n_test = ncol(B),
                           cv_auc = round(r$cv, 4), external_auc = round(r$ext, 4),
                           drop = round(r$drop, 4), mmd = round(r$mmd, 4),
                           mmd_norm = round(r$mmd_norm, 4),
                           agree_delta = round(r$agree_delta, 3),
                           cov90 = round(r$cov90, 4), abstain90 = round(r$abstain90, 4),
                           stringsAsFactors = FALSE)
  res[[pid]] <- row
  out0 <- do.call(rbind, res); rownames(out0) <- NULL
  write.csv(out0, file.path(PROC, "benchmark_R.csv"), row.names = FALSE, fileEncoding = "UTF-8")
}

if (length(res)) {
  out <- do.call(rbind, res); rownames(out) <- NULL
  write.csv(out, file.path(PROC, "benchmark_R.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  cat("\n--- benchmark_R.csv ---\n"); print(out, row.names = FALSE)
  if (nrow(out) >= 3) {
    fit <- lm(drop ~ agree_delta, data = out)
    cat(sprintf("\n拟合：drop = %.3f + %.2f·Δagreement  (R²=%.2f, n=%d)\n",
                coef(fit)[1], coef(fit)[2], summary(fit)$r.squared, nrow(out)))
  }
}
