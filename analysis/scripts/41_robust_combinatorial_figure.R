#!/usr/bin/env Rscript
# Step 41: figure for robust TF x TF combinatorics on raw counts.
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(gridExtra)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
D <- file.path(base, "results", "combinatorial_push"); FG <- file.path(base, "figures"); dir.create(FG, showWarnings = FALSE)
theme_set(theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 11.5),
  plot.subtitle = element_text(colour = "grey35", size = 8.5), strip.text = element_text(face = "bold", size = 8),
  panel.border = element_rect(fill = NA, colour = "grey88")))
cols <- c(neither = "#1E7AB6", a = "#EA8535", b = "#519C3E", both = "#C2382F")
g <- fread(file.path(D, "grid_summary.tsv"))
gl <- melt(g, id.vars = "setting", measure.vars = c("strict_n", "medium_n", "relaxed_n"), variable.name = "tier", value.name = "observed")
gn <- melt(g, id.vars = "setting", measure.vars = c("strict_null", "medium_null", "relaxed_null"), variable.name = "tier", value.name = "null")
gl[, null := gn$null]; gl[, tier := factor(sub("_n", "", tier), levels = c("strict", "medium", "relaxed"))]
gl[, setting := factor(setting, levels = c("ge1", "ge1_loadmod", "ge2", "ge2_loadmod", "log", "log_loadmod"),
       labels = c(">=1 read", ">=1 read\n+ load-modulated", ">=2 reads", ">=2 reads\n+ load-modulated", "log dose", "log dose\n+ load-modulated"))]
pA <- ggplot(gl, aes(setting, observed + 1, fill = tier)) + geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_point(aes(y = null + 1, group = tier), position = position_dodge(0.8), shape = 4, size = 2.2) + scale_y_log10() +
  scale_fill_manual(values = c(strict = "#6A1B4D", medium = "#B5527F", relaxed = "#E8B4CB"), name = "Tier") +
  labs(x = NULL, y = "Reproducible TF x TF interactions + 1 (log)", title = "A  Search grid on raw counts (x = permutation null)",
       subtitle = "strict: q<0.05 & p<0.05 in each exp; medium: q<0.1 & p<0.1; relaxed: q<0.2 & p<0.2 (same sign in both experiments)") +
  theme(legend.position = "top", axis.text.x = element_text(size = 7.5))
p <- fread(file.path(D, "pair_replication.tsv"))[setting == "log_loadmod"][order(q)][1:20]
p[, pair := factor(pair, levels = rev(pair))]
pB <- ggplot(p, aes(r_exp1_exp2, pair, fill = q < 0.05)) + geom_col(width = 0.7) +
  geom_vline(xintercept = 0.060, linetype = 2) + annotate("text", x = 0.061, y = 1, label = "null 95%", hjust = 0, size = 2.8) +
  geom_text(aes(label = n_strict), hjust = -0.3, size = 2.7) +
  scale_fill_manual(values = c(`TRUE` = "#C2382F", `FALSE` = "grey65"), labels = c("q >= 0.05", "q < 0.05"), name = NULL) +
  labs(x = "Exp1 vs Exp2 correlation of the pair's interaction profile (all genes)", y = NULL,
       title = "B  TF pairs whose combinatorial profile replicates genome-wide", subtitle = "Number = strict gene-level interactions") +
  theme(legend.position = "bottom")
v <- fread(file.path(D, "examples_values.tsv"))
v[, state := factor(state, levels = names(cols))]
v[, lab := paste0(TFa, " x ", TFb, " -> ", target)]
v[, lab := factor(lab, levels = unique(lab))]
pC <- ggplot(v, aes(state, rel, fill = state)) + geom_hline(yintercept = 0, colour = "grey75") +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.9, linewidth = 0.3) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 1), size = 0.6, alpha = 0.6) +
  facet_wrap(~ lab + exp, ncol = 6, scales = "free_y", labeller = function(l) list(paste0(l$lab, "\n", sub("Single-cell ", "", l$exp)))) +
  scale_fill_manual(values = cols, guide = "none") + scale_x_discrete(labels = c("--", "a", "b", "a+b")) +
  labs(x = NULL, y = "Adjusted expression vs 'neither'", title = "C  Examples (interaction tested on TF dose; boxes show presence states for display)",
       subtitle = "10 random pools of cells per state; adjusted for depth, load (incl. per-TF load modulation), cell type and other TFs") +
  theme(strip.text = element_text(size = 6.5), axis.text.x = element_text(size = 7))
gg <- arrangeGrob(arrangeGrob(pA, pB, ncol = 2, widths = c(1.2, 1)), pC, ncol = 1, heights = c(1, 1.9))
ggsave(file.path(FG, "robust_combinatorics_raw.png"), gg, width = 14, height = 15, dpi = 180, bg = "white")
ggsave(file.path(FG, "robust_combinatorics_raw.pdf"), gg, width = 14, height = 15)
cat("ok\n")
