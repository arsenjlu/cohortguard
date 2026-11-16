## 第 5 项：covariate-shift 加权共形（Tibshirani et al. 2019 思路）作为对照
## 密度比 w(x)=dP_target/dP_train 用 50 个标签选中基因上的岭逻辑回归估计
## 用法：Rscript scripts_r/m13_shift_conformal.R
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

## 密度比：在 50 个标签基因上，岭逻辑回归区分 train(0)/test(1)，w = p/(1-p)
cg_density_ratio <- function(X_train, X_test, lam = 1) {
  X <- cbind(X_train, X_test)
  z <- c(rep(0, ncol(X_train)), rep(1, ncol(X_test)))
  b <- cg_fit_logit(t(X), z, lam = lam)
  p <- cg_predict(t(X), b)
  w <- p / (1 - p)
  w <- w / sum(w)
  list(w_train = w[seq_len(ncol(X_train))], w_test = w[ncol(X_train) + seq_len(ncol(X_test))])
}

## 加权共形：校准分数按 w 加权取 (1-alpha) 分位
cg_conformal_weighted <- function(p_cal, y_cal, p_te, y_te, w_cal, alpha = 0.1) {
  s <- 1 - ifelse(y_cal == 1, p_cal, 1 - p_cal)
  ord <- order(s); s_ord <- s[ord]; w_ord <- w_cal[ord]
  cum <- cumsum(w_ord) / sum(w_ord)
  k <- which(cum >= 1 - alpha)[1]
  if (is.na(k)) k <- length(s_ord)
  q <- s_ord[k]; thr <- 1 - q
  S1 <- p_te >= thr; S0 <- (1 - p_te) >= thr
  size <- as.integer(S1) + as.integer(S0)
  hit <- ifelse(y_te == 1, S1, S0)
  list(coverage = mean(hit), abstention = mean(size == 2))
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

  sel <- order(-abs(cg_tstat(A, d1$y)))[seq_len(50)]
  sr <- cg_split_ref(d1$y, seed = 7)
  Xtr <- A[sel, sr$tr, drop = FALSE]; ytr <- d1$y[sr$tr]
  Xref <- A[sel, sr$ref, drop = FALSE]; yref <- d1$y[sr$ref]
  Xtg <- B[sel, , drop = FALSE]; ytg <- d2$y

  # 密度比：用训练队列全部样本 vs 目标队列全部样本
  ratio <- cg_density_ratio(A[sel, , drop = FALSE], B[sel, , drop = FALSE])
  w_cal <- ratio$w_train[sr$ref]  # 校准样本对应的权重（取自训练队列内部）

  # 与 cg_evaluate 一致：测试预测用全训练模型 b_full，校准预测用 72% 模型 b2
  b_full <- cg_fit_logit(t(A[sel, , drop = FALSE]), d1$y)
  p_te <- cg_predict(t(Xtg), b_full)
  b2 <- cg_fit_logit(t(Xtr), ytr)
  p_ref <- cg_predict(t(Xref), b2)

  cf_std <- cg_conformal(p_ref, yref, p_te, ytg, 0.1)
  cf_w <- cg_conformal_weighted(p_ref, yref, p_te, ytg, w_cal, 0.1)

  res[[pid]] <- data.frame(pair = pid,
                           cov90_std = round(cf_std$coverage, 3),
                           cov90_w = round(cf_w$coverage, 3),
                           abstain_std = round(cf_std$abstention, 3),
                           abstain_w = round(cf_w$abstention, 3))
  cat(sprintf("%s  cov_std=%.2f cov_w=%.2f  abstain_std=%.2f abstain_w=%.2f\n",
              pid, cf_std$coverage, cf_w$coverage, cf_std$abstention, cf_w$abstention))
}
out <- do.call(rbind, res); rownames(out) <- NULL
write.csv(out, file.path(PROC, "m13_shift_conformal.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat(sprintf("\n覆盖率：标准中位 %.2f，加权中位 %.2f\n",
            median(out$cov90_std), median(out$cov90_w)))
cat(sprintf("达 0.90 的比例：标准 %d/16，加权 %d/16\n",
            sum(out$cov90_std >= 0.9), sum(out$cov90_w >= 0.9)))
