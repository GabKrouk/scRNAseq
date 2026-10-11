#!/usr/bin/env Rscript
# Step 23: heatmap of the ERF5 x CRF4 genes agreeing in bulk and single cells (MinMax per gene per dataset).
suppressPackageStartupMessages({library(ggplot2); library(data.table)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
D <- file.path(base, "results", "bulk_bridge"); FG <- file.path(base, "figures", "bulk_bridge")
sel <- fread(file.path(D, "agreeing_genes_ERF5xCRF4.tsv")); setnames(sel, 1, "gene")
bl <- fread(file.path(D, "agreeing_genes_bulk_values.tsv"))[, .(m = mean(rel)), by = .(gene, exp, state)]
sc <- fread(file.path(D, "agreeing_genes_sc_pseudobulk.tsv"))[, .(m = mean(rel)), by = .(gene, exp, state)]
d <- rbind(bl, sc)
d[, mm := (m - min(m)) / (max(m) - min(m)), by = .(gene, exp)]
d[, state := factor(state, levels = c("neither", "ERF5", "CRF4", "both"), labels = c("--", "ERF5", "CRF4", "both"))]
d[, exp := factor(exp, levels = c("Bulk exp 1", "Bulk exp 2", "Single-cell Exp1", "Single-cell Exp2"))]
# order genes: combination-repressed block first, then combination-induced; clustered within block
sel[, grp := fifelse(meta_inter < 0, "Combination represses", "Combination induces")]
w <- dcast(d, gene ~ exp + state, value.var = "mm")
ord <- c()
for (g in c("Combination represses", "Combination induces")) {
  ww <- w[gene %in% sel[grp == g]$gene]
  ord <- c(ord, if (nrow(ww) > 2) ww$gene[hclust(dist(as.matrix(ww[, -1])), "average")$order] else ww$gene)
}
sel[, tier2 := fifelse(tier == "strict", "***", fifelse(tier == "bulk p<0.05, sc p<0.1", "**", "*"))]
d <- merge(d, sel[, .(gene, tier2, grp)], by = "gene")
d[, ylab := factor(paste(gene, tier2), levels = rev(paste(ord, sel[match(ord, gene)]$tier2)))]
d[, grp := factor(grp, levels = c("Combination represses", "Combination induces"))]
p <- ggplot(d, aes(state, ylab, fill = mm)) + geom_tile(colour = "white", linewidth = 0.4) + facet_grid(grp ~ exp, scales = "free_y", space = "free_y") +
  scale_fill_gradientn(colours = c("#2166AC", "#F7F7F7", "#B2182B"), name = "MinMax\nscaled") +
  labs(x = NULL, y = NULL, title = "ERF5 x CRF4: genes whose combinatorial response agrees in bulk DoubleTARGET and single cells",
       subtitle = "Same interaction sign in all four datasets.  *** bulk p<0.05 & sc p<0.05;  ** bulk p<0.05 & sc p<0.1;  * both p<0.1\nBulk: mean of replicates (log2 CPM); single cell: adjusted expression. Scaled per gene within each dataset.") +
  theme_minimal(base_size = 10) + theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"),
  plot.title = element_text(face = "bold"), plot.subtitle = element_text(colour = "grey35", size = 8.5), axis.text.y = element_text(family = "mono", size = 8))
ggsave(file.path(FG, "ERF5xCRF4_agreeing_genes_heatmap.png"), p, width = 10, height = 8, dpi = 220, bg = "white")
ggsave(file.path(FG, "ERF5xCRF4_agreeing_genes_heatmap.pdf"), p, width = 10, height = 8)
cat("ok\n")
