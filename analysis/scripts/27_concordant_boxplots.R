#!/usr/bin/env Rscript
# Step 27: Sup Fig box plots - combinatorial genes with dynamics shared by bulk DoubleTARGET and scTARGET-seq.
suppressPackageStartupMessages({library(ggplot2); library(data.table)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
D <- file.path(base, "results", "shared_dynamics"); FG <- file.path(base, "figures", "bulk_bridge")
theme_set(theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
  strip.text = element_text(face = "bold", size = 7), plot.title = element_text(face = "bold", size = 13),
  plot.subtitle = element_text(colour = "grey35", size = 9), panel.border = element_rect(fill = NA, colour = "grey88")))
cols <- c(neither = "#1E7AB6", ERF5 = "#EA8535", CRF4 = "#519C3E", both = "#C2382F")   # Fig 7 palette
d <- fread(file.path(D, "combinatorial_concordant_values.tsv"))
d[, state := factor(state, levels = names(cols))]
d[, exp := factor(exp, levels = c("Bulk exp 1", "Bulk exp 2", "Single-cell Exp1", "Single-cell Exp2"))]
gord <- unique(d[order(module_name, -score)][, .(gene, module_name, q)])
gord[, cls := sub("combination-", "comb.-", sub("^M[0-9]+: ", "", sub(" \\(n=.*", "", module_name)))]
gord[, lab := sprintf("%s - %s", gene, cls)]
d <- merge(d, gord[, .(gene, lab)], by = "gene"); d[, lab := factor(lab, levels = gord$lab)]
pp <- function(dd, title) ggplot(dd, aes(state, rel, fill = state)) + geom_hline(yintercept = 0, colour = "grey70") +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.9, linewidth = 0.3) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 1), size = 0.7, alpha = 0.6) +
  facet_wrap(~ lab + exp, ncol = 4, scales = "free_y", labeller = function(l) list(paste0(l$lab, "\n", l$exp))) +
  scale_fill_manual(values = cols, guide = "none") + scale_x_discrete(labels = c("EV", "ERF5", "CRF4", "ERF5+\nCRF4")) +
  labs(x = NULL, y = "Change vs EV / cells without ERF5 and CRF4", title = title,
       subtitle = "Bulk DoubleTARGET: log2 CPM per replicate. scTARGET-seq: expression adjusted for depth, load, cell type and other TFs, 10 random pools of cells per state.") +
  theme(axis.text.x = element_text(size = 7))
pages <- split(gord$lab, ceiling(seq_along(gord$lab) / 6))
pdf(file.path(FG, "SupFig_concordant_ERF5xCRF4_combinatorial_boxplots.pdf"), width = 12, height = 15)
for (i in seq_along(pages)) print(pp(d[lab %in% pages[[i]]], sprintf("ERF5 x CRF4 combinatorial genes shared by bulk and single cells (%d/%d)", i, length(pages))))
invisible(dev.off())
top <- gord[, head(.SD, 2), by = module_name]$lab
g <- pp(d[lab %in% top], "ERF5 x CRF4 combinatorial genes with the same dynamics in bulk DoubleTARGET and scTARGET-seq")
ggsave(file.path(FG, "SupFig_concordant_ERF5xCRF4_combinatorial_boxplots.png"), g, width = 12, height = 2.35 * length(top) + 1.2, dpi = 200, bg = "white")
cat(length(gord$lab), "genes;", length(pages), "pages\n")
