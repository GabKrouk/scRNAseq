#!/usr/bin/env Rscript
# Step 29: TF dose in bulk DoubleTARGET and whether it explains the responses.
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(gridExtra)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
D <- file.path(base, "results"); FG <- file.path(base, "figures", "bulk_bridge")
theme_set(theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 11.5),
  plot.subtitle = element_text(colour = "grey35", size = 8.5), strip.text = element_text(face = "bold")))
cols <- c(N = "#1E7AB6", E = "#EA8535", C = "#519C3E", EC = "#C2382F")
lab4 <- c(N = "EV", E = "ERF5", C = "CRF4", EC = "ERF5+CRF4")
d <- melt(fread(file.path(D, "bulk_bridge", "bulk_TF_dose_per_sample.tsv")), id.vars = c("sample", "construct"), variable.name = "TF", value.name = "v")
d[, TF := sub("_log2CPM", "", TF)]; d[, construct := factor(construct, levels = names(cols))]
pA <- ggplot(d, aes(construct, v, fill = construct)) + geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.9) +
  geom_point(position = position_jitter(width = 0.1, seed = 1), size = 1.3) + facet_wrap(~TF) +
  scale_fill_manual(values = cols, guide = "none") + scale_x_discrete(labels = lab4) +
  labs(x = NULL, y = "log2 CPM of the TF locus", title = "A  TF dose in DoubleTARGET",
       subtitle = "ERF5 and CRF4 are ~45% lower (-0.85 log2) in the combined condition")
r <- fread(file.path(D, "bulk_bridge", "bulk_dose_models_per_gene.tsv")); setnames(r, 1, "gene")
sd <- fread(file.path(D, "shared_dynamics", "concordant_genes_modules.tsv"))
x <- merge(sd[, .(gene, module_label)], r, by = "gene")
x[, verdict := fifelse(dose_better, "explained by TF dose", fifelse(fact_better, "combination needed (factorial)", "undecided"))]
x[, cls := fifelse(grepl("^combination", module_label), "Combinatorial genes (n=35)", fifelse(grepl("^CRF4-driven ind", module_label), "CRF4 induction, attenuated with ERF5 (n=283)", "Other single-TF genes (n=15)"))]
v <- x[, .N, by = .(cls, verdict)]
v[, verdict := factor(verdict, levels = c("explained by TF dose", "undecided", "combination needed (factorial)"))]
pB <- ggplot(v, aes(N, cls, fill = verdict)) + geom_col(position = "fill", width = 0.6) +
  scale_fill_manual(values = c("explained by TF dose" = "#999999", undecided = "#DDDDDD", "combination needed (factorial)" = "#C2382F"), name = NULL) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Share of genes (bulk, AIC difference > 2)", y = NULL, title = "B  Dose model vs factorial model, for the concordant genes",
       subtitle = "Dose model: y ~ ERF5 level + CRF4 level. Factorial: y ~ ERF5 + CRF4 + ERF5:CRF4") + theme(legend.position = "bottom")
rr <- r[q_any_response < 0.05]
pC <- ggplot(rr, aes(inter_fact, inter_after_dose)) + geom_abline(linetype = 2, colour = "grey50") +
  geom_point(size = 0.8, alpha = 0.5, colour = "#C2382F") +
  labs(x = "ERF5 x CRF4 interaction, factorial model", y = "Interaction after adjusting\nfor ERF5 and CRF4 levels",
       title = "C  The interaction is not absorbed by dose", subtitle = sprintf("%d responsive genes; r = %.2f", nrow(rr), cor(rr$inter_fact, rr$inter_after_dose)))
g <- arrangeGrob(arrangeGrob(pA, pC, ncol = 2, widths = c(1.2, 1)), pB, nrow = 2, heights = c(1, 0.75))
ggsave(file.path(FG, "bulk_TF_dose_models.png"), g, width = 11, height = 8, dpi = 200, bg = "white")
cat("ok\n")
