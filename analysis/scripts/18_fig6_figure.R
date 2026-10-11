#!/usr/bin/env Rscript
# Step 18: new Fig 6 - Fig 6A/6B DoubleTARGET-validated genes seen in Exp1, Exp2 and combined.
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(gridExtra); library(grid)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
D <- file.path(base, "results", "fig6"); FG <- file.path(base, "figures");
theme_set(theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
  plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(colour = "grey35", size = 9),
  strip.text = element_text(face = "bold")))
state_col <- c(neither = "#BBBBBB", a = "#56B4E9", b = "#E69F00", both = "#CC3311")
pair_col <- c(`ERF5-NAC4` = "#CC3311", `ERF5-CRF4` = "#0077BB")

# ---- A: AT1G65845, paper thresholds (ERF5 >= 1, CRF4 >= 3 reads), depth-normalised only
A <- fread(file.path(D, "fig6A_AT1G65845_depthnorm.tsv"))[CRF4_threshold == 3]
A[, state := factor(state, levels = c("neither", "a", "b", "both"))]
pA <- ggplot(A, aes(state, mean, fill = state)) + geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) + geom_hline(yintercept = 100, linetype = 2, colour = "grey40") +
  geom_text(aes(y = 4, label = paste0("n=", n)), size = 2.6, colour = "white") +
  facet_wrap(~exp) + scale_fill_manual(values = state_col, guide = "none") +
  scale_x_discrete(labels = c("--", "ERF5", "CRF4", "ERF5\n+CRF4")) +
  labs(x = NULL, y = "AT1G65845 expression\n(% of 'neither', mean +/- SE)", title = "A  AT1G65845 (Fig 6A gene), ERF5 x CRF4",
       subtitle = "Paper thresholds (ERF5 >= 1, CRF4 >= 3 reads); depth-normalised, no load correction")

# ---- B: heatmap of the 19 qPCR genes, min-max scaled per gene, per experiment
M <- fread(file.path(D, "fig6_state_means.tsv")); Tt <- fread(file.path(D, "fig6_interaction_tests.tsv"))
M[, mm := (mean - min(mean)) / (max(mean) - min(mean)), by = .(pair, target, exp)]
M[, state := factor(state, levels = c("neither", "a", "b", "both"), labels = c("--", "ERF5", "TF2", "both"))]
M[, exp := factor(exp, levels = c("Exp1", "Exp2", "Combined"))]
M[, pair := factor(pair, levels = c("ERF5-NAC4", "ERF5-CRF4"))]
ordg <- Tt[pair == "ERF5-NAC4"][order(-interaction_meta)]$target
M[, target := factor(target, levels = rev(ordg))]
Tt[, lab := sprintf("%s   %s %s", target, ifelse(same_sign, "same sign", "opposite "),
                    ifelse(p_interaction_meta < 0.001, "***", ifelse(p_interaction_meta < 0.01, "**", ifelse(p_interaction_meta < 0.05, "*", ""))))]
M <- merge(M, Tt[, .(pair, target = as.character(target), lab)], by = c("pair", "target"))
M[, pair := factor(pair, levels = c("ERF5-NAC4", "ERF5-CRF4"))]
M[, row := factor(paste(pair, target), levels = rev(c(paste("ERF5-NAC4", ordg), paste("ERF5-CRF4", ordg))))]
labs_map <- setNames(M$lab, M$row)
pB <- ggplot(M, aes(state, row, fill = mm)) + geom_tile(colour = "white", linewidth = 0.4) +
  facet_grid(pair ~ exp, scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) labs_map[x]) +
  scale_fill_gradientn(colours = c("#2166AC", "#F7F7F7", "#B2182B"), name = "MinMax\nscaled") +
  labs(x = NULL, y = NULL, title = "A  Fig 6B genes (qPCR-validated): four cell states, adjusted for load, depth, cell type and other TFs",
       subtitle = "TF2 = NAC4 (top) or CRF4 (bottom). Row label: whether the ERF5 x TF2 interaction has the same sign in Exp1 and Exp2; stars = combined interaction p") +
  theme(panel.grid = element_blank(), legend.position = "right", axis.text.y = element_text(family = "mono", size = 7.5))

# ---- C: interaction effect Exp1 vs Exp2
lim <- max(abs(c(Tt$interaction_Exp1, Tt$interaction_Exp2))) * 1.1
cnt <- Tt[, .(lab = sprintf("%s: same sign %d/%d", pair, sum(same_sign), .N)), by = pair]
pC <- ggplot(Tt, aes(interaction_Exp1, interaction_Exp2, colour = pair)) +
  annotate("rect", xmin = 0, xmax = lim, ymin = 0, ymax = lim, fill = "#FBE3D6", alpha = 0.5) +
  annotate("rect", xmin = -lim, xmax = 0, ymin = -lim, ymax = 0, fill = "#DCEAF6", alpha = 0.5) +
  geom_hline(yintercept = 0, colour = "grey60") + geom_vline(xintercept = 0, colour = "grey60") +
  geom_point(size = 2.4, alpha = 0.9) +
  geom_text(data = Tt[target == "AT1G65845"], aes(label = "AT1G65845"), size = 2.6, vjust = -1, show.legend = FALSE) +
  scale_colour_manual(values = pair_col, name = NULL, labels = setNames(cnt$lab, cnt$pair)) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(x = "Interaction (both - ERF5 - TF2), Exp1", y = "Interaction, Exp2", title = "B  Does the combinatorial effect replicate?",
       subtitle = paste0(paste(sapply(split(Tt, Tt$pair), function(x) sprintf("%s: %d/%d same sign (binomial p = %.2g)", x$pair[1], sum(x$same_sign), nrow(x), binom.test(sum(x$same_sign), nrow(x), alternative = "greater")$p.value)), collapse = "\n"), "\nRaw counts, QC cells: not replicated beyond chance")) +
  theme(legend.position = "bottom", legend.direction = "vertical")

# panel A (AT1G65845 bar plot) dropped at gab's request: heatmap + replication scatter only
g <- arrangeGrob(pB, arrangeGrob(pC, ncol = 1), ncol = 2, widths = c(2.3, 1))
ggsave(file.path(FG, "Fig6_revisited_Exp1_Exp2.png"), g, width = 16, height = 9.5, dpi = 220, bg = "white")
ggsave(file.path(FG, "Fig6_revisited_Exp1_Exp2.pdf"), g, width = 16, height = 9.5)
cat("ok\n")
