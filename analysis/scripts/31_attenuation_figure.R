#!/usr/bin/env Rscript
# Step 31: dose argument figure (Reviewer 3 #1): ERF5 attenuates CRF4-induced genes in single cells
# even though CRF4 dose is unchanged there, whereas in bulk the combined condition halves CRF4 dose.
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(gridExtra)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
D <- file.path(base, "results", "bulk_bridge"); FG <- file.path(base, "figures", "bulk_bridge")
theme_set(theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
  plot.title = element_text(face = "bold", size = 11.5), plot.subtitle = element_text(colour = "grey35", size = 8.5),
  strip.text = element_text(face = "bold", size = 8), panel.border = element_rect(fill = NA, colour = "grey88")))
cols <- c(neither = "#1E7AB6", ERF5 = "#EA8535", CRF4 = "#519C3E", both = "#C2382F")
labs4 <- c(neither = "EV", ERF5 = "ERF5", CRF4 = "CRF4", both = "ERF5+\nCRF4")
dd <- fread(file.path(D, "attenuation_examples_CRF4dose.tsv"))
dd[, state := factor(state, levels = c("CRF4", "both"))]
dd[, dataset := factor(dataset, levels = c("Bulk exp 1", "Bulk exp 2", "Single-cell Exp1", "Single-cell Exp2"))]
dd[, type := ifelse(grepl("Bulk", dataset), "log2 CPM (CRF4 locus)", "log1p junction reads per cell")]
pA <- ggplot(dd, aes(state, dose, fill = state)) + geom_violin(data = dd[grepl("Single", dataset)], scale = "width", width = 0.7, alpha = 0.6, colour = NA) +
  geom_boxplot(width = 0.35, outlier.shape = NA, alpha = 0.95) +
  geom_point(data = dd[grepl("Bulk", dataset)], size = 1.6) +
  facet_wrap(~dataset, nrow = 1, scales = "free_y") + scale_fill_manual(values = cols, guide = "none") +
  scale_x_discrete(labels = c("CRF4 only", "ERF5+CRF4")) +
  labs(x = NULL, y = "CRF4 dose", title = "A  CRF4 dose: halved in bulk DoubleTARGET, unchanged in single cells",
       subtitle = "Bulk: CRF4 locus expression per sample (~45% lower with ERF5). Single cell: CRF4 junction reads in CRF4+ cells with or without ERF5 (1.37 vs 1.35; 1.58 vs 1.54)")
r <- fread(file.path(D, "attenuation_all_CRF4_genes.tsv")); setnames(r, 1, "gene")
rt <- melt(r[, .(gene, t_Exp1, t_Exp2)], id.vars = "gene"); rt[, variable := sub("t_", "Single-cell ", variable)]
st <- rt[, .(lab = sprintf("%d%% of genes attenuated\n(sign test p = %.0e)", round(100 * mean(value < 0)), binom.test(sum(value < 0), .N)$p.value)), by = variable]
pB <- ggplot(rt, aes(value)) + geom_histogram(bins = 40, fill = "#C2382F", alpha = 0.8) + geom_vline(xintercept = 0, linetype = 2) +
  geom_text(data = st, aes(x = -Inf, y = Inf, label = lab), hjust = -0.05, vjust = 1.3, size = 3) + facet_wrap(~variable) +
  labs(x = "ERF5 effect in CRF4+ cells (both vs CRF4 only, t)", y = "Genes",
       title = "B  ERF5 attenuates the 283 CRF4-induced genes in single cells, at constant CRF4 dose")
v <- fread(file.path(D, "attenuation_examples_values.tsv"))
v[, state := factor(state, levels = names(cols))]
v[, exp := factor(exp, levels = c("Bulk exp 1", "Bulk exp 2", "Single-cell Exp1", "Single-cell Exp2"))]
pC <- ggplot(v, aes(state, rel, fill = state)) + geom_hline(yintercept = 0, colour = "grey75") +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.9, linewidth = 0.3) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 1), size = 0.7, alpha = 0.6) +
  facet_wrap(~ gene + exp, ncol = 4, scales = "free_y", labeller = function(l) list(paste(l$gene, "|", l$exp))) +
  scale_fill_manual(values = cols, guide = "none") + scale_x_discrete(labels = labs4) +
  labs(x = NULL, y = "Change vs EV / cells without ERF5 and CRF4",
       title = "C  Example genes: CRF4 induction lowered by ERF5 in bulk AND in single cells",
       subtitle = "Genes with significant attenuation (t < -1.96) in both single-cell experiments. Bulk: log2 CPM per replicate; single cell: adjusted expression, 10 random pools per state") +
  theme(axis.text.x = element_text(size = 7))
ng <- length(unique(v$gene))
g <- arrangeGrob(pA, pB, pC, ncol = 1, heights = c(1, 0.9, 0.55 * ng + 0.4))
ggsave(file.path(FG, "R3_dose_argument_ERF5_attenuation.png"), g, width = 11.5, height = 6.5 + 2.4 * ng, dpi = 200, bg = "white")
ggsave(file.path(FG, "R3_dose_argument_ERF5_attenuation.pdf"), g, width = 11.5, height = 6.5 + 2.4 * ng)
cat("ok\n")
