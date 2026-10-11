#!/usr/bin/env Rscript
# Step 37: QC figure from raw CellRanger matrices (Reviewer 2 #9, Reviewer 3 minor 6, Reviewer 2 #2).
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(gridExtra)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
D <- file.path(base, "results", "qc_raw"); FG <- file.path(base, "figures")
theme_set(theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 11),
  plot.subtitle = element_text(colour = "grey35", size = 8.5), strip.text = element_text(face = "bold")))
q <- fread(file.path(D, "qc_plot_table.tsv"))
scol <- c(kept = "#4D4D4D", `low UMI/genes (MAD)` = "#D55E00", doublet = "#CC79A7", `removed by authors` = "#56B4E9")
q[, status := factor(status, levels = names(scol))]
pA <- ggplot(q, aes(umi, genes, colour = status)) + geom_point(size = 0.35, alpha = 0.5) + scale_x_log10() + scale_y_log10() +
  facet_wrap(~exp) + scale_colour_manual(values = scol, name = NULL) + guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(x = "UMI per cell (raw, CellRanger)", y = "Genes detected", title = "A  Cell quality from the raw CellRanger matrices",
       subtitle = "MAD flags computed per library; doublets called at the expected 10x multiplet rate (Scrublet-like simulation)")
m <- melt(q[, .(exp, status, `% mitochondrial` = pct_mito, `% plastid` = pct_plastid)], id.vars = c("exp", "status"))
pB <- ggplot(m, aes(value + 0.01, fill = exp)) + geom_histogram(bins = 60, position = "identity", alpha = 0.55) + scale_x_log10() +
  geom_vline(xintercept = 5, linetype = 2) + facet_wrap(~variable, scales = "free") + scale_fill_manual(values = c(Exp1 = "#0072B2", Exp2 = "#E69F00"), name = NULL) +
  labs(x = "% of UMIs (log scale, +0.01)", y = "Cells", title = "B  Organellar fractions are low; the 5% plastid filter removes few cells")
det <- fread(file.path(D, "Exp2_TF_detection_vs_depth.tsv"))
pC <- ggplot(det, aes(umi_median, frac_no_TF * 100)) + geom_line(colour = "#C2382F") + geom_point(colour = "#C2382F", size = 2) +
  geom_line(aes(y = (1 - frac_RFP) * 100), colour = "grey40", linetype = 2) + geom_point(aes(y = (1 - frac_RFP) * 100), colour = "grey40", size = 1.5) +
  scale_x_log10() + labs(x = "UMI per cell (decile median, Exp2)", y = "% cells without any read", title = "C  Detection depends on depth (Exp2)",
       subtitle = "Red: no GR-TF junction read; grey dashed: no RFP read")
pD <- ggplot(q[!is.na(n_TF_paper) & status != "removed by authors"], aes(status, n_TF_paper, fill = status)) + geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~exp) + scale_fill_manual(values = scol, guide = "none") +
  labs(x = NULL, y = "GR-TFs detected per cell", title = "D  Low-UMI cells carry almost all GR-TFs",
       subtitle = "High construct load, low endogenous transcriptome") + theme(axis.text.x = element_text(angle = 20, hjust = 1))
g <- arrangeGrob(pA, arrangeGrob(pB, pC, ncol = 2, widths = c(1.3, 1)), pD, ncol = 1, heights = c(1, 0.9, 0.85))
ggsave(file.path(FG, "QC_raw_cellranger.png"), g, width = 11, height = 12, dpi = 200, bg = "white")
cat("ok\n")
