#!/usr/bin/env Rscript
# Step 25: figures for shared ERF5 x CRF4 dynamics across bulk and single-cell datasets.
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(gridExtra)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
D <- file.path(base, "results", "shared_dynamics"); FG <- file.path(base, "figures", "bulk_bridge")
source(file.path(base, "scripts", "paper_style.R")); theme_set(theme_tpc(7))
cols <- c(EV = "#1E7AB6", ERF5 = "#EA8535", CRF4 = "#519C3E", both = "#C2382F")   # Fig 7 palette
L <- fread(file.path(D, "module_profiles_long.tsv"))
L[, state := factor(state, levels = names(cols))]
L[, dataset := factor(dataset, levels = c("Bulk exp 1", "Bulk exp 2", "Single-cell Exp1", "Single-cell Exp2"))]
ordm <- L[, .N, by = module_name][order(-grepl("combination", module_name), module_name)]$module_name
L[, module_name := factor(module_name, levels = ordm)]
mean_ <- L[, .(m = mean(value), se = sd(value) / sqrt(.N)), by = .(module_name, dataset, state)]
pA <- ggplot(L, aes(state, value, group = gene)) + geom_hline(yintercept = 0, colour = "grey75") +
  geom_line(colour = "grey60", alpha = 0.25, linewidth = 0.3) +
  geom_ribbon(data = mean_, aes(x = as.integer(state), ymin = m - 1.96 * se, ymax = m + 1.96 * se), inherit.aes = FALSE, fill = "grey30", alpha = 0.18) +
  geom_line(data = mean_, aes(state, m, group = 1), inherit.aes = FALSE, linewidth = 1, colour = "black") +
  geom_point(data = mean_, aes(state, m, fill = state), inherit.aes = FALSE, shape = 21, size = 6 * PT, stroke = 0.4) +
  facet_grid(module_name ~ dataset, labeller = labeller(module_name = label_wrap_gen(22))) +
  scale_fill_manual(values = cols, guide = "none") + scale_x_discrete(labels = c("EV", "ERF5", "CRF4", "ERF5+\nCRF4")) +
  labs(x = NULL, y = "Scaled response (per gene, per dataset)",
       title = "A",
       subtitle = paste0(nrow(fread(file.path(base, "results", "shared_dynamics", "concordant_genes_modules.tsv"))), " genes whose response profile has the same shape in bulk exp 1, bulk exp 2, single-cell Exp1 and Exp2 (q<0.1)\nThin lines: genes; black: module mean +/- 95% CI")) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0), axis.text.x = element_text(size = 7))

R <- fread(file.path(D, "concordance_all_genes.tsv")); Nq <- fread(file.path(D, "null_quantiles.tsv"))
pB <- ggplot() + geom_density(data = Nq, aes(null_score, fill = "permuted genes"), alpha = 0.5, colour = NA) +
  geom_density(data = R, aes(score, fill = "observed"), alpha = 0.5, colour = NA) +
  scale_fill_manual(values = c(observed = "#C2382F", `permuted genes` = "grey55"), name = NULL) +
  geom_vline(xintercept = min(R[q < 0.1]$score), linetype = 2) +
  annotate("text", x = min(R[q < 0.1]$score), y = Inf, label = " q < 0.1", hjust = 0, vjust = 1.5, size = 6.5 * PT) +
  labs(x = "Cross-experiment concordance score (shape agreement x strength)", y = "Density",
       title = "B", subtitle = "Gene labels shuffled independently in each dataset (200x)") +
  theme(legend.position = c(0.8, 0.8))
g <- arrangeGrob(pA, pB, ncol = 1, heights = c(4.2, 1))
save_tpc(g, FG, "shared_dynamics_4experiments", 6.5, 8)
cat("ok\n")
