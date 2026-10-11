#!/usr/bin/env Rscript
# Step 50: Supplemental figure - cell QC on raw CellRanger 7.1 counts and the SCTransform TF-imputation artefact.
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(gridExtra)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
FG <- file.path(base, "figures", "paper")
source(file.path(base, "scripts", "paper_style.R")); theme_set(theme_tpc(7))
q <- rbindlist(lapply(c("Exp1", "Exp2"), function(e) fread(file.path(base, "results", "qc", paste0(e, "_cell_qc.tsv.gz")))[, exp := e]))
q[, fate := fifelse(pass_qc, "kept", fifelse(low_umi | low_genes, "low UMI / genes", fifelse(high_mito | high_plastid, "mito or plastid > 5%", "doublet")))]
q[, fate := factor(fate, levels = c("kept", "low UMI / genes", "mito or plastid > 5%", "doublet"))]
fc <- c(kept = "grey60", `low UMI / genes` = "#C2382F", `mito or plastid > 5%` = "#519C3E", doublet = "#EA8535")
nk <- q[, .(lab = sprintf("%s kept / %s", format(sum(pass_qc), big.mark = ","), format(.N, big.mark = ","))), by = exp]
pA <- ggplot(q[order(fate)], aes(umi, genes, colour = fate)) + geom_point(size = 0.35, alpha = 0.6) + facet_wrap(~exp) +
  geom_text(data = nk, aes(x = 300, y = 12000, label = lab), inherit.aes = FALSE, hjust = 0, size = 6 * PT) +
  scale_x_log10() + scale_y_log10() + scale_colour_manual(values = fc, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  labs(x = "UMIs per cell", y = "Genes per cell", title = "A",
       subtitle = "Low UMI / genes: below median - 3 MAD per library; doublets: Scrublet-like score, expected rate 0.8% per 1,000 cells") + theme(legend.position = "bottom")
q[, bin := cut(umi, breaks = c(0, 1000, 2500, 5000, 10000, 20000, 40000, Inf), labels = c("<1k", "1-2.5k", "2.5-5k", "5-10k", "10-20k", "20-40k", ">40k"))]
d <- q[, .(no_tf = mean(n_tf == 0), med = as.numeric(median(n_tf)), n = .N), by = .(exp, bin)][n >= 20]
pB <- ggplot(d, aes(bin, no_tf * 100, fill = exp)) + geom_col(position = position_dodge(0.8), width = 0.75) +
  scale_fill_manual(values = c(Exp1 = "#0072B2", Exp2 = "#E69F00"), name = NULL) +
  labs(x = "UMIs per cell", y = "% cells with no GR-TF read", title = "B",
       subtitle = "All cells before QC; the no-TF fraction falls to ~5-15% in well-sequenced cells") + theme(legend.position = c(0.8, 0.8))
pt <- fread(file.path(base, "..", "reproducibility", "results", "qc_raw", "qc_plot_table.tsv"))[exp == "Exp2"]
m <- merge(q[exp == "Exp2", .(umi, genes, n_tf, fate)], pt[, .(umi, genes, n_TF_paper)], by = c("umi", "genes"))
m <- unique(m[fate %in% c("kept", "low UMI / genes") & !is.na(n_TF_paper)], by = c("umi", "genes"))
ml <- melt(m, measure.vars = c("n_tf", "n_TF_paper"))
ml[, variable := factor(variable, labels = c("raw counts", "SCTransform matrix"))]
pC <- ggplot(ml, aes(fate, value, fill = variable)) + geom_boxplot(outlier.size = 0.3, width = 0.7, linewidth = 0.3) +
  scale_fill_manual(values = c(`raw counts` = "grey80", `SCTransform matrix` = "#C2382F"), name = NULL) +
  labs(x = NULL, y = "GR-TFs detected per cell", title = "C",
       subtitle = sprintf("Exp2, %s cells. Low-UMI cells: median %d TFs in raw counts\nbut %d in the corrected matrix. TFs are now called on raw\ncounts in QC-passing cells only.",
                          format(nrow(m), big.mark = ","), as.integer(median(m[fate != "kept"]$n_tf)), as.integer(median(m[fate != "kept"]$n_TF_paper)))) +
  theme(legend.position = "bottom")
g <- arrangeGrob(pA, arrangeGrob(pB, pC, ncol = 2), ncol = 1, heights = c(1, 1))
save_tpc(g, FG, "SupFig_QC_raw_counts", 6.5, 6)
cat("matched", nrow(m), "\n")
