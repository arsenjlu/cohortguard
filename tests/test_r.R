## R 单元测试：Rscript tests/test_r.R
source(file.path("R", "cohortguard.R"), encoding = "UTF-8")

ok <- 0; fail <- 0
check <- function(name, cond) {
  if (isTRUE(cond)) { cat("ok   ", name, "\n"); ok <<- ok + 1 }
  else { cat("FAIL ", name, "\n"); fail <<- fail + 1 }
}

check("AUC 完美", abs(cg_auc(c(0, 0, 1, 1), c(.1, .2, .8, .9)) - 1) < 1e-9)
check("AUC 反转", abs(cg_auc(c(0, 0, 1, 1), c(.9, .8, .2, .1)) - 0) < 1e-9)
check("AUC 并列", abs(cg_auc(c(0, 1, 0, 1), c(.5, .5, .5, .5)) - 0.5) < 1e-9)

set.seed(1)
X <- matrix(rnorm(40 * 300), 40, 300)
y <- sample(rep(c(0, 1), each = 150))          # 必须打乱，否则前 150 例只有一类
X[, y == 1] <- X[, y == 1] + 0.6
Z <- cg_zscore(X)
sel <- order(-abs(cg_tstat(Z[, 1:150], y[1:150])))[1:20]
b <- cg_fit_logit(t(Z[sel, 1:150]), y[1:150])
p <- cg_predict(t(Z[sel, 151:300]), b)
check("逻辑回归可分离到较好 AUC", cg_auc(y[151:300], p) > 0.7)

check("MMD 自身为 0", abs(cg_mmd(Z[, 1:100], Z[, 1:100])) < 1e-6)

## 共形：同分布时覆盖率应接近目标
set.seed(2)
Xc <- matrix(rnorm(30 * 2000), 30, 2000)
yc <- rep(c(0, 1), each = 1000)
s <- order(-abs(cg_tstat(Xc[, 1:1000], yc[1:1000])))[1:20]
bc <- cg_fit_logit(t(Xc[s, 1:1000]), yc[1:1000])
cf <- cg_conformal(cg_predict(t(Xc[s, 1001:1500]), bc), yc[1001:1500],
                   cg_predict(t(Xc[s, 1501:2000]), bc), yc[1501:2000], 0.1)
check("共形覆盖达标", cf$coverage >= 0.86)

cat(sprintf("\n%d passed, %d failed\n", ok, fail))
if (fail > 0) quit(status = 1)
