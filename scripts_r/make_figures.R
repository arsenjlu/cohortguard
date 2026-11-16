## Generate manuscript figures (ggplot2)
suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel); library(patchwork)
})

d <- read.csv(file.path("data", "processed", "benchmark_merged.csv"), stringsAsFactors = FALSE)

meta <- data.frame(
  pair = c("P1","P1r","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11","P12","P13","P14","P15"),
  label = c("DKD (A to B)","DKD (B to A)","LUAD","HCC","Breast","CRC","ccRCC",
            "IBD mucosa to blood","MS to RA","MS","SLE","SLE (small to large)",
            "UC colon","UC colon (cross-platform)","Asthma airway to blood","RA"),
  group = c("Mismatched","Mismatched","Matched","Matched","Matched","Matched","Matched",
            "Mismatched","Mismatched","Matched","Matched","Matched","Matched","Matched",
            "Mismatched","Matched"),
  stringsAsFactors = FALSE
)
d <- merge(d, meta, by = "pair", all.x = TRUE, sort = FALSE)
d$group[is.na(d$group)] <- "Matched"
d <- d[!is.na(d$drop) & !is.na(d$cv_auc), ]

pal <- c("Matched" = "#2C6FBB", "Mismatched" = "#C0392B")
theme_set(theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title = element_text(colour = "grey20"),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 10),
        legend.position = "top"))

d1 <- d[order(d$drop), ]
d1$label <- factor(d1$label, levels = d1$label)
med <- aggregate(drop ~ group, data = d1, median)
p1 <- ggplot(d1, aes(x = drop, y = label, colour = group)) +
  geom_vline(xintercept = 0, colour = "grey75", linetype = 2) +
  geom_segment(aes(x = 0, xend = drop, yend = label), colour = "grey88", linewidth = 0.6) +
  geom_point(size = 3.4) +
  geom_vline(data = med, aes(xintercept = drop, colour = group), linetype = 3, linewidth = 0.8) +
  scale_colour_manual(values = pal, name = NULL) +
  labs(x = "Drop in AUC (cross-validated minus external)", y = NULL,
       subtitle = "Dotted vertical lines: group medians. Positive means worse outside the training cohort.") +
  theme(panel.grid.major.y = element_line(colour = "grey93"))

mk <- function(var, xlab, ttl) {
  xv <- d[[var]]; yv <- d$drop
  fit <- lm(yv ~ xv)
  r2 <- summary(fit)$r.squared
  ggplot(d, aes(x = .data[[var]], y = drop, colour = group)) +
    geom_hline(yintercept = 0, colour = "grey85", linetype = 2) +
    geom_smooth(method = "lm", se = FALSE, colour = "grey30", linewidth = 0.8) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = pair), size = 3, show.legend = FALSE, seed = 1) +
    scale_colour_manual(values = pal, name = NULL) +
    labs(x = xlab, y = "Drop in AUC")
}
p2a <- mk("mmd_all_norm", "Kernel MMD, label-free (normalised)", "A. Label-free kernel distance")
p2b <- mk("mmd_norm", "Kernel MMD, model features (normalised)", "B. Feature-space kernel distance")
p2c <- mk("agree_delta", "Agreement drift", "C. Model-agreement drift")

p3 <- ggplot(d, aes(x = cov90, y = label, colour = group)) +
  geom_vline(xintercept = 0.9, colour = "grey40", linetype = 2) +
  geom_point(aes(size = abstain90), alpha = 0.9) +
  scale_size_continuous(range = c(2, 9), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
                        limits = c(0, 0.5),
                        labels = scales::percent, name = "Abstention") +
  scale_colour_manual(values = pal, name = NULL) +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(x = "Realised coverage at a nominal 90%", y = NULL,
       subtitle = "Dashed line: nominal level. Point size: fraction of abstentions.") +
  theme(panel.grid.major.y = element_line(colour = "grey93"))

p4 <- ggplot(d, aes(x = cv_auc, y = external_auc, colour = group)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey70", linetype = 2) +
  geom_hline(yintercept = 0.5, colour = "grey88", linetype = 3) +
  geom_point(size = 3.2) +
  geom_text_repel(aes(label = pair), size = 3, show.legend = FALSE, seed = 2) +
  scale_colour_manual(values = pal, name = NULL) +
  coord_cartesian(xlim = c(0.55, 1.02), ylim = c(0.25, 1.05)) +
  labs(x = "Cross-validated AUC (training cohort)", y = "External AUC (target cohort)",
       subtitle = "Dashed: identity. Dotted: chance level.")

ggsave(file.path("figures", "fig1_drop_by_contrast.png"), p1, width = 9, height = 5.8, dpi = 300, bg = "white")
ggsave(file.path("figures", "fig2_drift_vs_drop.png"), p2a | p2b | p2c, width = 16.5, height = 4.6, dpi = 300, bg = "white")
ggsave(file.path("figures", "fig3_conformal.png"), p3, width = 9, height = 5.6, dpi = 300, bg = "white")
ggsave(file.path("figures", "fig4_cv_vs_external.png"), p4, width = 6.8, height = 5.6, dpi = 300, bg = "white")
## 矢量版（投稿用）
ggsave(file.path("figures", "fig1_drop_by_contrast.pdf"), p1, width = 9, height = 5.8, device = cairo_pdf, bg = "white")
ggsave(file.path("figures", "fig2_drift_vs_drop.pdf"), p2a | p2b | p2c, width = 16.5, height = 4.6, device = cairo_pdf, bg = "white")
ggsave(file.path("figures", "fig3_conformal.pdf"), p3, width = 9, height = 5.6, device = cairo_pdf, bg = "white")
ggsave(file.path("figures", "fig4_cv_vs_external.pdf"), p4, width = 6.8, height = 5.6, device = cairo_pdf, bg = "white")
cat("figures written\n")
