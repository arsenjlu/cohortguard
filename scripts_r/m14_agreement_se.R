## agreement 漂移 Δ 的 bootstrap 标准误（审稿人第 5 条）
## 每对用 20 个种子重跑 cg_agreement，取 Δ 的 SD 作为 SE；重点报告 n<50 的小队列
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

NSEEDS <- 20
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
  A <- cg_zscore(d1$M[common, , drop = FALSE]); B <- cg_zscore(d2$M[common, , drop = FALSE])
  ok <- is.finite(rowSums(A)) & is.finite(rowSums(B))
  A <- A[ok, , drop = FALSE]; B <- B[ok, , drop = FALSE]

  deltas <- numeric(NSEEDS)
  for (sd in seq_len(NSEEDS)) {
    sr <- cg_split_ref(d1$y, seed = sd)
    ag <- cg_agreement(A[, sr$tr, drop = FALSE], d1$y[sr$tr],
                       A[, sr$ref, drop = FALSE], B, k = 50, seed = sd)
    deltas[sd] <- ag["delta"]
  }
  res[[pid]] <- data.frame(pair = pid, n_train = ncol(A),
                           delta_mean = round(mean(deltas, na.rm = TRUE), 3),
                           delta_se = round(sd(deltas, na.rm = TRUE), 3))
  cat(sprintf("%s  n_train=%d  Δ=%.3f ± %.3f (SE)\n", pid, ncol(A),
              mean(deltas, na.rm = TRUE), sd(deltas, na.rm = TRUE)))
}
out <- do.call(rbind, res); rownames(out) <- NULL
write.csv(out, file.path(PROC, "agreement_se.csv"), row.names = FALSE, fileEncoding = "UTF-8")

small <- out[out$n_train < 50, ]
large <- out[out$n_train >= 50, ]
cat(sprintf("\n=== 汇总 ===\n"))
cat(sprintf("n<50 的队列（%d 个）：Δ 的 SE 均值 %.3f（范围 %.3f–%.3f）\n",
            nrow(small), mean(small$delta_se), min(small$delta_se), max(small$delta_se)))
cat(sprintf("n>=50 的队列（%d 个）：Δ 的 SE 均值 %.3f\n", nrow(large), mean(large$delta_se)))
cat(sprintf("Δ 本身的 SD（跨 16 对）: %.3f\n", sd(out$delta_mean)))
