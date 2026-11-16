## cohortGuard —— 纯 R 实现（无外部依赖）
## 与 Python 版数值对齐：AUC 处理并列、逻辑回归带 ridge、共形用 split conformal

## ---------------- 读取 ----------------
cg_open <- function(path) {
  if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
}

## 读 GEO series matrix：返回 mat(基因×样本)、samples、fields(列表)
cg_read_series_matrix <- function(path) {
  con <- cg_open(path); on.exit(close(con))
  fields <- list(); ids <- character(0); vals <- list()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "!series_matrix_table_begin")) {
      hdr <- strsplit(readLines(con, n = 1, warn = FALSE), "\t")[[1]]
      samples <- gsub('"', "", hdr[-1])
      repeat {
        l <- readLines(con, n = 1, warn = FALSE)
        if (length(l) == 0 || startsWith(l, "!series_matrix_table_end")) break
        cc <- strsplit(l, "\t")[[1]]
        if (length(cc) != length(hdr)) next
        ids <- c(ids, gsub('"', "", cc[1]))
        vals[[length(vals) + 1]] <- suppressWarnings(as.numeric(cc[-1]))
      }
      break
    }
    if (startsWith(line, "!Sample_")) {
      cc <- strsplit(line, "\t")[[1]]
      v <- gsub('"', "", cc[-1])
      key <- sub("^!Sample_", "", cc[1])
      if (length(v) > 0 && grepl(":", v[1])) {
        key <- trimws(tolower(sub(":.*$", "", v[1])))
        v <- sub("^[^:]*: ", "", v)
      }
      if (is.null(fields[[key]])) fields[[key]] <- v
    }
  }
  mat <- do.call(rbind, vals)
  list(mat = mat, ids = ids, samples = samples, fields = fields)
}

## 读 GPL 注释：返回 probe -> symbol 的命名向量
cg_read_gpl <- function(path) {
  con <- cg_open(path); on.exit(close(con))
  sym_col <- NA; ga_col <- NA; extra <- NA
  keys <- character(0); syms <- character(0)
  while (TRUE) {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "#")) next
    cc <- strsplit(line, "\t")[[1]]
    low <- tolower(cc)
    if (is.na(sym_col) && is.na(ga_col)) {
      sym_col <- match(TRUE, low %in% c("gene symbol", "genesymbol", "symbol", "gene_symbol", "gene.symbol"))
      ga_col <- match(TRUE, grepl("gene_assignment", low))
      if (is.na(sym_col) && is.na(ga_col)) stop("注释表头无法识别: ", paste(head(cc, 10), collapse = " | "))
      extra <- if (length(low) > 1 && grepl("probeset", low[2])) 2L else NA_integer_
      next
    }
    s <- NA_character_
    if (!is.na(ga_col) && length(cc) >= ga_col) {
      ga <- cc[ga_col]
      if (nzchar(ga) && !startsWith(ga, "---")) {
        first <- strsplit(ga, " /// ", fixed = TRUE)[[1]][1]
        parts <- strsplit(first, " // ", fixed = TRUE)[[1]]
        if (length(parts) >= 2) s <- trimws(parts[2])
      }
    }
    if (is.na(s) && !is.na(sym_col) && length(cc) >= sym_col) {
      cand <- trimws(strsplit(cc[sym_col], " /// ", fixed = TRUE)[[1]][1])
      if (!is.na(cand) && nzchar(cand) && cand != "---") s <- cand
    }
    if (!is.na(s) && nzchar(s) && s != "---") {
      keys <- c(keys, cc[1]); syms <- c(syms, s)
      if (!is.na(extra) && length(cc) >= extra && nzchar(cc[extra]) && cc[extra] != "---") {
        keys <- c(keys, cc[extra]); syms <- c(syms, s)
      }
    }
  }
  setNames(syms, keys)
}

cg_to_gene <- function(mat, ids, map) {
  g <- unname(map[ids])
  keep <- !is.na(g)
  mat <- mat[keep, , drop = FALSE]; g <- g[keep]
  ug <- sort(unique(g))
  out <- matrix(NA_real_, nrow = length(ug), ncol = ncol(mat),
                dimnames = list(ug, colnames(mat)))
  for (i in seq_along(ug)) {
    idx <- which(g == ug[i])
    out[i, ] <- if (length(idx) == 1) mat[idx, ] else colMeans(mat[idx, , drop = FALSE], na.rm = TRUE)
  }
  out
}

cg_zscore <- function(M) {
  mu <- rowMeans(M, na.rm = TRUE)
  sd <- apply(M, 1, sd, na.rm = TRUE); sd[!is.finite(sd) | sd == 0] <- NA
  (M - mu) / sd
}

## ---------------- 统计 ----------------
cg_auc <- function(y, s) {
  r <- rank(s, ties.method = "average")
  n1 <- sum(y == 1); n0 <- length(y) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

cg_tstat <- function(Xg, y) {
  a <- Xg[, y == 1, drop = FALSE]; b <- Xg[, y == 0, drop = FALSE]
  sd <- sqrt(apply(a, 1, var) / ncol(a) + apply(b, 1, var) / ncol(b)) + 1e-12
  (rowMeans(a) - rowMeans(b)) / sd
}

cg_fit_logit <- function(X, y, lam = 1, iters = 60) {
  Xb <- cbind(1, X); p <- ncol(Xb); b <- rep(0, p)
  for (i in seq_len(iters)) {
    eta <- pmin(pmax(drop(Xb %*% b), -30), 30)
    mu <- 1 / (1 + exp(-eta)); w <- pmax(mu * (1 - mu), 1e-6)
    z <- eta + (y - mu) / w
    A <- crossprod(Xb, w * Xb) + lam * diag(p); A[1, 1] <- A[1, 1] - lam / 2
    nb <- tryCatch(solve(A, crossprod(Xb, w * z)), error = function(e) b)
    if (max(abs(nb - b)) < 1e-8) { b <- nb; break }
    b <- nb
  }
  b
}

cg_predict <- function(X, b) {
  eta <- pmin(pmax(drop(cbind(1, X) %*% b), -30), 30)
  1 / (1 + exp(-eta))
}

cg_mmd <- function(A, B, mult = 1) {
  Z <- t(cbind(A, B)); nA <- ncol(A); nB <- ncol(B)
  D <- as.matrix(dist(Z))^2
  med <- median(D[D > 0]); if (!is.finite(med) || med <= 0) med <- 1
  med <- med * mult^2
  K <- exp(-D / med)
  iA <- seq_len(nA); iB <- nA + seq_len(nB)
  mean(K[iA, iA]) + mean(K[iB, iB]) - 2 * mean(K[iA, iB])
}

cg_cv_auc <- function(Xg, y, k = 50, folds = 5, reps = 4, seed = 0) {
  set.seed(seed); X <- t(Xg); out <- numeric(reps)
  for (r in seq_len(reps)) {
    i0 <- which(y == 0); i1 <- which(y == 1)
    i0 <- sample(i0); i1 <- sample(i1)
    oof <- numeric(length(y))
    for (f in seq_len(folds)) {
      te <- c(i0[seq(f, length(i0), by = folds)], i1[seq(f, length(i1), by = folds)])
      if (!length(te)) next
      tr <- setdiff(seq_along(y), te)
      st <- abs(cg_tstat(Xg[, tr, drop = FALSE], y[tr]))
      sel <- order(-st)[seq_len(k)]
      b <- cg_fit_logit(X[tr, sel, drop = FALSE], y[tr])
      oof[te] <- cg_predict(X[te, sel, drop = FALSE], b)
    }
    out[r] <- cg_auc(y, oof)
  }
  mean(out, na.rm = TRUE)
}

cg_agreement <- function(Xtr, ytr, Xref, Xtgt, k = 50, nboot = 12, seed = 0) {
  set.seed(seed)
  sel <- order(-abs(cg_tstat(Xtr, ytr)))[seq_len(k)]
  n <- ncol(Xtr); Pr <- list(); Pt <- list()
  for (i in seq_len(nboot)) {
    idx <- sample.int(n, n, replace = TRUE)
    if (length(unique(ytr[idx])) < 2) next
    b <- cg_fit_logit(t(Xtr[sel, idx, drop = FALSE]), ytr[idx])
    Pr[[length(Pr) + 1]] <- cg_predict(t(Xref[sel, , drop = FALSE]), b) > 0.5
    Pt[[length(Pt) + 1]] <- cg_predict(t(Xtgt[sel, , drop = FALSE]), b) > 0.5
  }
  pair_agree <- function(P) {
    P <- do.call(rbind, P); M <- nrow(P); s <- 0; c <- 0
    for (a in seq_len(M - 1)) for (bb in (a + 1):M) { s <- s + mean(P[a, ] == P[bb, ]); c <- c + 1 }
    if (c == 0) NA_real_ else s / c
  }
  ar <- pair_agree(Pr); at <- pair_agree(Pt)
  c(ref = ar, tgt = at, delta = ar - at)
}

cg_split_ref <- function(y, frac = 0.28, seed = 0) {
  set.seed(seed)
  i0 <- sample(which(y == 0)); i1 <- sample(which(y == 1))
  n0 <- max(2, round(length(i0) * frac)); n1 <- max(2, round(length(i1) * frac))
  ref <- c(i0[seq_len(n0)], i1[seq_len(n1)])
  list(tr = setdiff(seq_along(y), ref), ref = ref)
}

cg_conformal <- function(p_cal, y_cal, p_te, y_te, alpha = 0.1) {
  s <- 1 - ifelse(y_cal == 1, p_cal, 1 - p_cal)
  k <- min(length(s), ceiling((length(s) + 1) * (1 - alpha)))
  q <- sort(s)[k]; thr <- 1 - q
  S1 <- p_te >= thr; S0 <- (1 - p_te) >= thr
  size <- as.integer(S1) + as.integer(S0)
  hit <- ifelse(y_te == 1, S1, S0)
  list(coverage = mean(hit), abstention = mean(size == 2), n_conf = sum(size == 1))
}

## ---------------- 单对评估 ----------------
cg_evaluate <- function(Az, yA, Bz, yB, k = 50, seed = 0) {
  cv <- cg_cv_auc(Az, yA, k = k, seed = seed)
  sel <- order(-abs(cg_tstat(Az, yA)))[seq_len(k)]
  b <- cg_fit_logit(t(Az[sel, , drop = FALSE]), yA)
  p_ext <- cg_predict(t(Bz[sel, , drop = FALSE]), b)
  ext <- cg_auc(yB, p_ext)
  mmd <- cg_mmd(Az[sel, , drop = FALSE], Bz[sel, , drop = FALSE])
  sr <- cg_split_ref(yA, seed = seed)
  ag <- cg_agreement(Az[, sr$tr, drop = FALSE], yA[sr$tr], Az[, sr$ref, drop = FALSE], Bz,
                     k = k, seed = seed)
  b2 <- cg_fit_logit(t(Az[sel, sr$tr, drop = FALSE]), yA[sr$tr])
  cf <- cg_conformal(cg_predict(t(Az[sel, sr$ref, drop = FALSE]), b2), yA[sr$ref], p_ext, yB, 0.1)
  list(cv = cv, ext = ext, drop = cv - ext, mmd = mmd, mmd_norm = mmd / (1 + mmd),
       agree_ref = ag["ref"], agree_tgt = ag["tgt"], agree_delta = ag["delta"],
       cov90 = cf$coverage, abstain90 = cf$abstention)
}


## 无标签 MMD：取合并数据中方差最大的 nv 个基因（不使用标签）
cg_mmd_topvar <- function(A, B, nv = 2000, mult = 1) {
  X <- cbind(A, B)
  v <- apply(X, 1, var)
  idx <- order(-v)[seq_len(min(nv, nrow(X)))]
  cg_mmd(A[idx, , drop = FALSE], B[idx, , drop = FALSE], mult = mult)
}
