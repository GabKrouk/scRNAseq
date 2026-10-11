#!/usr/bin/env Rscript
# Step 22: box plots for genes whose ERF5 x CRF4 interaction agrees in bulk and single cells.
suppressPackageStartupMessages({library(ggplot2); library(data.table)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
D <- file.path(base, "results", "bulk_bridge"); FG <- file.path(base, "figures", "bulk_bridge"); dir.create(FG, showWarnings = FALSE)
theme_set(theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
  strip.text = element_text(face = "bold", size = 8.5), plot.title = element_text(face = "bold", size = 13),
  plot.subtitle = element_text(colour = "grey35"), panel.border = element_rect(fill = NA, colour = "grey85")))
# Fig 7 palette (Will's cbp1: EV, ERF5, CRF4, ERF5_CRF4)
cols <- c(neither = "#1E7AB6", ERF5 = "#EA8535", CRF4 = "#519C3E", both = "#C2382F")
sel <- fread(file.path(D, "agreeing_genes_ERF5xCRF4.tsv")); setnames(sel, 1, "gene")
bl <- fread(file.path(D, "agreeing_genes_bulk_values.tsv"))[, .(gene, exp, state, rel)]
sc <- fread(file.path(D, "agreeing_genes_sc_pseudobulk.tsv"))[, .(gene, exp, state, rel)]
d <- rbind(bl, sc)
d[, state := factor(state, levels = names(cols))]
d[, exp := factor(exp, levels = c("Bulk exp 1", "Bulk exp 2", "Single-cell Exp1", "Single-cell Exp2"))]
sel[, lab := sprintf("%s  [%s]   bulk p=%.2g, sc p=%.2g", gene, tier, bulk_p_inter, p_meta_inter)]
d <- merge(d, sel[, .(gene, lab, rank, p_meta_inter)], by = "gene")
ord <- sel[order(rank, p_meta_inter)]$lab
d[, lab := factor(lab, levels = ord)]
plotset <- function(dd, title, sub) ggplot(dd, aes(state, rel, fill = state)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.85, linewidth = 0.35) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 1), size = 0.9, alpha = 0.7) +
  facet_wrap(~ lab + exp, ncol = 4, scales = "free_y", labeller = function(l) list(paste0(sub("  .*", "", l$lab), " | ", l$exp))) +
  scale_fill_manual(values = cols, guide = "none") + scale_x_discrete(labels = c("EV / none", "ERF5", "CRF4", "ERF5+CRF4")) +
  labs(x = NULL, y = "Change vs 'neither'\n(bulk: log2 CPM per sample; single cell: adjusted log expr., 10 random cell pools)",
       title = title, subtitle = sub) + theme(strip.text = element_text(size = 8))
top <- d[lab %in% ord[1:8]]
g <- plotset(top, "ERF5 x CRF4: genes where the combination behaves the same in bulk DoubleTARGET and single cells",
             "Interaction same sign in bulk exp 1, bulk exp 2, single-cell Exp1 and Exp2. Rows 1-2: bulk p<0.05 & sc p<0.05; rows 3-8: bulk p<0.05 & sc p<0.1. Each panel has its own y-scale.")
ggsave(file.path(FG, "ERF5xCRF4_agreeing_genes_top8.png"), g, width = 12, height = 2.1 * 8 + 1.3, dpi = 200, bg = "white", limitsize = FALSE)
pdf(file.path(FG, "ERF5xCRF4_agreeing_genes_all.pdf"), width = 12, height = 2.1 * 9 + 1.3)
for (i in split(ord, ceiling(seq_along(ord) / 9))) print(plotset(d[lab %in% i], "ERF5 x CRF4: genes agreeing in bulk and single cells", "Tiers: strict; bulk p<0.05 & sc p<0.1; both p<0.1"))
invisible(dev.off()); cat("ok\n")
