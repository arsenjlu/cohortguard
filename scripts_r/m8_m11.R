## M8 + M11：一致性漂移连续概率版 + 各对病例比例 + 外部 AUC 的 DeLong 95% CI
## 用法：Rscript scripts_r/m8_m11.R
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

## DeLong AUC 方差
cg_auc_delong <- function(y, s) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 < 2 || n0 < 2) return(list(auc = NA, se = NA, lo = NA, hi = NA))
  s1 <- s[y == 1]; s0 <- s[y == 0]
  Vx <- vapply(s1, function(x) mean(ifelse(s0 < x, 1, ifelse(s0 == x, 0.5, 0))), numeric(1))
  Vy <- vapply(s0, function(yy) mean(ifelse(s1 > yy, 1, ifelse(s1 == yy, 0.5, 0))), numeric(1))
  auc <- mean(Vx)
  v10 <- if (n1 > 1) sum((Vx - auc)^2) / (n1 - 1) else 0
  v01 <- if (n0 > 1) sum((Vy - auc)^2) / (n0 - 1) else 0
  se <- sqrt(v10 / n1 + v01 / n0)
  list(auc = auc, se = se, lo = auc - 1.96 * se, hi = auc + 1.96 * se)
}

## 连续概率版一致性：12 个模型在参考/目标上的两两 Pearson 相关均值
cg_agreement_cont <- function(Xtr, ytr, Xref, Xtgt, k = 50, nboot = 12, seed = 0) {
  set.seed(seed)
  sel <- order(-abs(cg_tstat(Xtr, ytr)))[seq_len(k)]
  n <- ncol(Xtr); Pr <- list(); Pt <- list()
  for (i in seq_len(nboot)) {
    idx <- sample.int(n, n, replace = TRUE)
    if (length(unique(ytr[idx])) < 2) next
    b <- cg_fit_logit(t(Xtr[sel, idx, drop = FALSE]), ytr[idx])
    Pr[[length(Pr) + 1]] <- cg_predict(t(Xref[sel, , drop = FALSE]), b)
    Pt[[length(Pt) + 1]] <- cg_predict(t(Xtgt[sel, , drop = FALSE]), b)
  }
  pair_cor <- function(P) {
    M <- length(P); if (M < 2) return(NA_real_)
    s <- 0; c <- 0
    for (a in 1:(M - 1)) for (bb in (a + 1):M) {
      cc <- suppressWarnings(cor(P[[a]], P[[bb]]))
      if (is.finite(cc)) { s <- s + cc; c <- c + 1 }
    }
    if (c == 0) NA_real_ else s / c
  }
  cr <- pair_cor(Pr); ct <- pair_cor(Pt)
  c(ref = cr, tgt = ct, delta = cr - ct)
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
  A <- cg_zscore(d1$M[common, , drop = FALSE]); B <- cg_zscore(d2$M[common, , drop = FALSE])
  ok <- is.finite(rowSums(A)) & is.finite(rowSums(B))
  A <- A[ok, , drop = FALSE]; B <- B[ok, , drop = FALSE]

  # 外部 AUC + DeLong
  sel <- order(-abs(cg_tstat(A, d1$y)))[seq_len(50)]
  b <- cg_fit_logit(t(A[sel, , drop = FALSE]), d1$y)
  p_ext <- cg_predict(t(B[sel, , drop = FALSE]), b)
  dl <- cg_auc_delong(d2$y, p_ext)

  # 连续概率一致性
  sr <- cg_split_ref(d1$y, seed = 7)
  ag <- cg_agreement_cont(A[, sr$tr, drop = FALSE], d1$y[sr$tr],
                          A[, sr$ref, drop = FALSE], B, k = 50, seed = 7)

  res[[pid]] <- data.frame(pair = pid,
                           prev_train = round(sum(d1$y == 1) / length(d1$y), 3),
                           prev_target = round(sum(d2$y == 1) / length(d2$y), 3),
                           ext_auc = round(dl$auc, 3), ext_se = round(dl$se, 3),
                           ext_lo = round(dl$lo, 3), ext_hi = round(dl$hi, 3),
                           agree_cont_delta = round(ag["delta"], 3),
                           agree_cont_ref = round(ag["ref"], 3), agree_cont_tgt = round(ag["tgt"], 3))
  cat(sprintf("%s  prev_tr=%.2f prev_tg=%.2f  ext=%.3f (%.3f–%.3f)  Δcont=%.3f\n",
              pid, res[[pid]]$prev_train, res[[pid]]$prev_target,
              dl$auc, dl$lo, dl$hi, ag["delta"]))
}
out <- do.call(rbind, res); rownames(out) <- NULL
write.csv(out, file.path(PROC, "m8_m11.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat("\n=== 连续概率一致性 vs 掉分 ===")
bm <- read.csv(file.path(PROC, "benchmark_R.csv"), stringsAsFactors = FALSE)
mm <- merge(bm[, c("pair", "drop")], out, by = "pair")
cat(sprintf("\nΔcont ~ drop: R²=%.3f (二元版为 0.02)\n",
            summary(lm(drop ~ agree_cont_delta, data = mm))$r.squared))
cat(sprintf("\n外部 AUC 95%% CI 宽度：均值 %.3f，最宽 %.3f（%s）\n",
            mean(mm$ext_hi - mm$ext_lo), max(mm$ext_hi - mm$ext_lo),
            mm$pair[which.max(mm$ext_hi - mm$ext_lo)]))
