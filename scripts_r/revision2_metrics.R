## 补充计算：
##   #6  随机 50 基因 MMD（不用标签选基因）—— 区分"维度"还是"标签"
##   #7  去掉 P8 后的所有回归 R²
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

pairs <- unique(cfg$pair_id[cfg$role %in% c("train", "test", "validation")])
NS <- 5
rand <- matrix(NA_real_, nrow = length(pairs), ncol = NS,
               dimnames = list(pairs, paste0("s", 1:NS)))
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
  ok <- is.finite(rowSums(A)) & is.finite(rowSums(B)); A <- A[ok, ]; B <- B[ok, ]
  for (s in 1:NS) {
    set.seed(1000 + s)
    sel <- sample.int(nrow(A), 50)
    rand[pid, s] <- cg_mmd(A[sel, , drop = FALSE], B[sel, , drop = FALSE])
  }
  cat(pid, "random-50 MMD mean", round(mean(rand[pid, ]), 4), "\n")
}
write.csv(rand, file.path(PROC, "random50_mmd.csv"))

## R² of drop ~ random-50 MMD (per seed)
main <- read.csv(file.path(PROC, "benchmark_R.csv"), stringsAsFactors = FALSE)
d <- merge(main, data.frame(pair = rownames(rand), rand), by = "pair")
r2fun <- function(x, y) summary(lm(y ~ x))$r.squared
cat("\n=== random-50 MMD R2 (label-free, dimensionality control) ===\n")
for (s in 1:NS) cat(sprintf("  seed %d: R2 = %.3f\n", s, r2fun(d[[paste0("s", s)]], d$drop)))
cat(sprintf("  mean across seeds: %.3f\n", mean(sapply(1:NS, function(s) r2fun(d[[paste0("s", s)]], d$drop)))))

## #7 去掉 P8
cat("\n=== exclude P8 (n = 15) ===\n")
d8 <- d[d$pair != "P8", ]
for (v in c("mmd_all_norm", "mmd_k05", "mmd_k1", "mmd_k2", "mmd_norm", "cv_auc", "agree_delta")) {
  cat(sprintf("  %-12s R2 = %.3f\n", v, r2fun(d8[[v]], d8$drop)))
}
