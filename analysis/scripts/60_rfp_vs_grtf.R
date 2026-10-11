#!/usr/bin/env Rscript
# Step 60: RFP reporter signal vs GR-TF junction signal per cell (QC cells), by vector type.
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(gridExtra)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
d <- fread(file.path(base, "results", "paper", "rfp_vs_grtf.tsv"))
theme_set(theme_minimal(base_size = 9) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10.5),
  plot.subtitle = element_text(colour = "grey35", size = 7.5), panel.border = element_rect(fill = NA, colour = "grey85"), strip.text = element_text(face = "bold")))
d[, rfp_cp := rfp / umi * 1e4][, tf_cp := tf / umi * 1e4]
pA <- ggplot(d, aes(tf_cp + 0.1, rfp_cp + 0.1)) + geom_point(size = 0.3, alpha = 0.3, colour = "#C2382F") + facet_wrap(~exp) +
  scale_x_log10() + scale_y_log10() + labs(x = "GR-TF junction reads per 10k UMIs (+0.1)", y = "RFP UMIs per 10k UMIs (+0.1)",
  title = "A  RFP and GR-TF signals agree cell by cell", subtitle = "Depth-adjusted Spearman rho = 0.71 (Exp1), 0.81 (Exp2)")
d[, nb := factor(pmin(n, 15), levels = 0:15, labels = c(0:14, "15+"))]
pB <- ggplot(d, aes(nb, rfp_cp + 0.1)) + geom_boxplot(outlier.size = 0.2, fill = "#F4B6A6", linewidth = 0.25) + facet_wrap(~exp) + scale_y_log10() +
  labs(x = "Number of GR-TFs detected in the cell", y = "RFP UMIs per 10k UMIs (+0.1)", title = "B  RFP rises with the number of constructs")
d[, grp := factor(grp, levels = c("no GR-TF", "GFP-vector TFs only", "RFP-vector TFs only", "both vector types"))]
pC <- ggplot(d[umi > median(umi), .SD, by = exp], aes(grp, rfp_cp + 0.1, fill = grp)) + geom_boxplot(outlier.size = 0.2, linewidth = 0.25) +
  facet_wrap(~exp) + scale_y_log10() + scale_fill_manual(values = c("grey70", "#519C3E", "#C2382F", "#7B3294"), guide = "none") +
  labs(x = NULL, y = "RFP UMIs per 10k UMIs (+0.1)", title = "C  Cells without GR-TF reads have ambient-level RFP",
       subtitle = "Cells deeper than the median of their experiment") + theme(axis.text.x = element_text(angle = 20, hjust = 1))
g <- arrangeGrob(pA, pB, pC, ncol = 1)
ggsave(file.path(base, "figures", "paper", "Response_Fig_RFP_vs_GRTF.png"), g, width = 8, height = 10, dpi = 300, bg = "white")
ggsave(file.path(base, "figures", "paper", "Response_Fig_RFP_vs_GRTF.pdf"), g, width = 8, height = 10)
cat("ok\n")
