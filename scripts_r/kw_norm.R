r <- read.csv("data/processed/revision_metrics.csv")
m <- read.csv("data/processed/benchmark_R.csv")
d <- merge(m, r, by = "pair")
f <- function(x) summary(lm(d$drop ~ x))$r.squared
cat(sprintf("k05 norm %.3f | k1 norm %.3f | k2 norm %.3f | all norm %.3f\n",
            f(d$mmd_k05/(1+d$mmd_k05)), f(d$mmd_k1/(1+d$mmd_k1)),
            f(d$mmd_k2/(1+d$mmd_k2)), f(d$mmd_all/(1+d$mmd_all))))
