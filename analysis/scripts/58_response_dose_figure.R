#!/usr/bin/env Rscript
# Step 58: Response-letter figure for R3 #1 (bulk dose vs factorial models per FDR-defined gene class).
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(gridExtra)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
B <- file.path(base, "results", "bulk_bridge"); FG <- file.path(base, "figures", "paper")
theme_set(theme_minimal(base_size = 9) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10.5),
  plot.subtitle = element_text(colour = "grey35", size = 7.5), panel.border = element_rect(fill = NA, colour = "grey85"), strip.text = element_text(face = "bold")))
cols <- c(EV = "#1E7AB6", ERF5 = "#EA8535", CRF4 = "#519C3E", `ERF5+CRF4` = "#C2382F")
d <- fread(file.path(B, "bulk_TF_dose_per_sample.tsv"))
d[, cond := factor(c(N = "EV", E = "ERF5", C = "CRF4", EC = "ERF5+CRF4")[construct], levels = names(cols))]
m <- melt(d, id.vars = c("sample", "cond"), measure.vars = c("ERF5_log2CPM", "CRF4_log2CPM"))
m[, variable := sub("_log2CPM", "", variable)]
pA <- ggplot(m, aes(cond, value, fill = cond)) + geom_boxplot(outlier.shape = NA, width = 0.6) + geom_jitter(width = 0.1, size = 1) +
  facet_wrap(~variable) + scale_fill_manual(values = cols, guide = "none") +
  labs(x = NULL, y = "log2 CPM of the TF locus", title = "A  TF dose in bulk DoubleTARGET", subtitle = "Each TF is ~45% lower (-0.85 log2) in the combined condition") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
x <- merge(fread(file.path(B, "Fig6_gene_classes.tsv")), fread(file.path(B, "bulk_dose_models_per_gene.tsv"))[, .(gene = get(names(.SD)[1]), dose_better, fact_better)], by = "gene")
s <- x[, .(`explained by TF dose` = mean(dose_better), undecided = mean(!dose_better & !fact_better), `combination needed (factorial)` = mean(fact_better), n = .N), by = class]
sl <- melt(s, id.vars = c("class", "n")); sl[, lab := sprintf("%s (n=%d)", class, n)]
sl[, variable := factor(variable, levels = c("explained by TF dose", "undecided", "combination needed (factorial)"))]
pB <- ggplot(sl, aes(value * 100, lab, fill = variable)) + geom_col(width = 0.65) +
  scale_fill_manual(values = c("grey55", "grey88", "#C2382F"), name = NULL) +
  labs(x = "% of genes (AIC difference > 2)", y = NULL, title = "B  TF dose explains the weaker CRF4 response, not the combination-specific genes",
       subtitle = "Dose model: y ~ experiment + ERF5 level + CRF4 level. Factorial: y ~ experiment + ERF5 + CRF4 + ERF5:CRF4") + theme(legend.position = "bottom")
g <- arrangeGrob(pA, pB, ncol = 1, heights = c(1, 0.8))
ggsave(file.path(FG, "Response_R3-1_dose.png"), g, width = 8, height = 7, dpi = 250, bg = "white")
cat("ok\n")
