#!/usr/bin/env Rscript
# Step 62: Response-letter panel (R3 #3) - sensitivity of reproducibility to the GR-TF encoding and covariates.
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(gridExtra)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
cl <- fread(file.path(base, "results", "calibration_linear.tsv"))
theme_set(theme_minimal(base_size = 9) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10.5),
  plot.subtitle = element_text(colour = "grey35", size = 7.5), panel.border = element_rect(fill = NA, colour = "grey85")))
cl[, scheme := factor(scheme, levels = c("ge1", "ge2", "ge3", "ge5", "log", "share"),
                      labels = c(">=1 read", ">=2 reads", ">=3 reads", ">=5 reads", "log dose\n(used)", "share of the\ncell's TF reads"))]
cl[, cov := factor(cov, levels = c("none", "depth", "depth+load"), labels = c("library only", "+ depth", "+ depth + GR-TF load / RFP"))]
fills <- c("#BBBBBB", "#56B4E9", "#1E7AB6")
pA <- ggplot(cl, aes(scheme, n_rep, fill = cov)) + geom_col(position = position_dodge(0.85), width = 0.8) +
  geom_point(aes(y = null_n_rep_mean, group = cov), position = position_dodge(0.85), shape = 4, size = 2, stroke = 0.8) +
  scale_fill_manual(values = fills, name = "Covariates") +
  labs(x = "GR-TF encoding", y = "TF -> gene effects significant in both\nexperiments (BH q < 0.05 each, same sign)",
       title = "A  Count-based encodings reproduce more TF effects than presence thresholds",
       subtitle = "Bars: real data; crosses: mean of 10 permutations of TF profiles within library x GR-TF-load strata") +
  theme(legend.position = "top")
pB <- ggplot(cl, aes(scheme, r_t, fill = cov)) + geom_col(position = position_dodge(0.85), width = 0.8) +
  geom_point(aes(y = null_r_t_mean, group = cov), position = position_dodge(0.85), shape = 4, size = 2, stroke = 0.8) +
  scale_fill_manual(values = fills, guide = "none") +
  labs(x = "GR-TF encoding", y = "Exp1 vs Exp2 correlation of\nTF x gene t-statistics",
       title = "B  Adjusting for depth and GR-TF load removes the correlation expected by chance",
       subtitle = "Crosses: permuted data. Without load adjustment, permuted data still correlate (spurious, load-driven signal)")
g <- arrangeGrob(pA, pB, ncol = 1, heights = c(1.1, 1))
ggsave(file.path(base, "figures", "paper", "Response_Fig_TF_encoding.png"), g, width = 8.5, height = 8, dpi = 300, bg = "white")
ggsave(file.path(base, "figures", "paper", "Response_Fig_TF_encoding.pdf"), g, width = 8.5, height = 8)
cat("ok\n")
