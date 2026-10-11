#!/usr/bin/env Rscript
# Step 15: figures for the four-state combinatorial analysis.
suppressPackageStartupMessages({library(ggplot2); library(data.table)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
D <- file.path(base, "results", "combinatorial"); FG <- file.path(base, "figures", "combinatorial")
dir.create(FG, showWarnings = FALSE, recursive = TRUE)
theme_set(theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
  strip.text = element_text(face = "bold", size = 8.5), plot.title = element_text(face = "bold", size = 13),
  plot.subtitle = element_text(colour = "grey35"), legend.position = "top"))
state_lab <- c(neither = "neither", a = "TFa only", b = "TFb only", both = "TFa + TFb")
state_col <- c(neither = "#BBBBBB", a = "#56B4E9", b = "#E69F00", both = "#CC3311")
exp_col <- c(Exp1 = "#0072B2", Exp2 = "#D55E00")
shape_col <- c(`AND-gate` = "#CC3311", reversal = "#7B3294", cancellation = "#009988", synergy = "#EE7733")

S <- fread(file.path(D, "state_means_for_plots.tsv"))
P <- fread(file.path(D, "proper_combinatorial_targets.tsv"))
A <- fread(file.path(D, "four_state_effects_all_consensus.tsv"))
S[, state := factor(state, levels = names(state_lab))]

# ---------- overview ----------
P[, tier := ifelse(core, "core: pair-specific & robust to number of TFs",
                    ifelse(pair_specific, "pair-specific only", "TF combinatorial with >=3 partners"))]
cnt <- P[, .N, by = .(shape, tier)]
cnt[, shape := factor(shape, levels = c("AND-gate", "reversal", "cancellation", "synergy"))]
cnt[, tier := factor(tier, levels = c("core: pair-specific & robust to number of TFs", "pair-specific only", "TF combinatorial with >=3 partners"))]
pA <- ggplot(cnt, aes(shape, N, fill = shape, alpha = tier)) + geom_col(width = 0.7) +
  scale_fill_manual(values = shape_col, guide = "none") + scale_alpha_manual(values = c(1, 0.55, 0.25), name = NULL) +
  geom_text(data = cnt[, .(N = sum(N)), by = shape], aes(shape, N, label = N), inherit.aes = FALSE, vjust = -0.4, size = 3.6) +
  labs(x = NULL, y = "TFa x TFb -> target cases", title = "Combinatorial cases by shape",
       subtitle = sprintf("%d reproducible interactions -> %d proper combinatorial (%d genes) -> %d core (%d genes)",
                          nrow(A), nrow(P), uniqueN(P$target), sum(P$core), uniqueN(P[core == TRUE]$target))) +
  theme(legend.position = "bottom", legend.direction = "vertical")

# median profile per shape (combined data, oriented so that 'both' is positive)
prof <- merge(S[exp == "Combined"], P[, .(pair, target, shape)], by = c("pair", "target"))
orient <- prof[state == "both", .(pair, target, sgn = sign(mean))]
prof <- merge(prof, orient, by = c("pair", "target"))
prof[, m := mean * sgn]
prof[, scl := max(abs(m)), by = .(pair, target)]; prof[, m := m / scl]
med <- prof[, .(m = median(m), lo = quantile(m, .25), hi = quantile(m, .75)), by = .(shape, state)]
med[, shape := factor(shape, levels = levels(cnt$shape))]
pB <- ggplot(med, aes(state, m, group = 1)) + geom_hline(yintercept = 0, colour = "grey70") +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = shape), alpha = 0.18) + geom_line(aes(colour = shape), linewidth = 1) +
  geom_point(aes(colour = shape), size = 2.4) + facet_wrap(~shape, nrow = 1) +
  scale_colour_manual(values = shape_col, guide = "none") + scale_fill_manual(values = shape_col, guide = "none") +
  scale_x_discrete(labels = c("--", "a", "b", "a+b")) +
  labs(x = NULL, y = "Scaled change vs neither\n(median, IQR)", title = "What each shape looks like (all cases, both experiments pooled)")
ggsave(file.path(FG, "overview_shapes.png"), pA, width = 6.4, height = 5, dpi = 220, bg = "white")
ggsave(file.path(FG, "overview_profiles.png"), pB, width = 9, height = 3, dpi = 220, bg = "white")

# pair matrix of proper combinatorial cases
pm <- P[core == TRUE, .N, by = .(TFa, TFb)]
pm2 <- rbind(pm, pm[, .(TFa = TFb, TFb = TFa, N)])
ord <- pm2[, .(t = sum(N)), by = TFa][order(-t)]$TFa
pm2[, TFa := factor(TFa, levels = ord)]; pm2[, TFb := factor(TFb, levels = rev(ord))]
pC <- ggplot(pm2, aes(TFa, TFb, fill = N)) + geom_tile(colour = "white") + geom_text(aes(label = N), size = 2.5) +
  scale_fill_gradient(low = "#FBE3D6", high = "#B2182B", trans = "sqrt", name = "cases") +
  labs(x = NULL, y = NULL, title = "Core combinatorial targets per TF pair") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), panel.grid = element_blank(), legend.position = "right")
ggsave(file.path(FG, "pair_matrix.png"), pC, width = 7.2, height = 6.2, dpi = 220, bg = "white")

# ---------- example panels ----------
panel <- function(keys, ncol = 4, title = NULL, subtitle = NULL) {
  d <- merge(S, keys, by = c("pair", "target"))
  d[, lab := factor(lab, levels = keys$lab)]
  comb <- d[exp == "Combined"]; reps <- d[exp != "Combined"]
  ggplot(mapping = aes(state, mean)) + geom_hline(yintercept = 0, colour = "grey75") +
    geom_col(data = comb, aes(fill = state), width = 0.72, alpha = 0.85) +
    geom_errorbar(data = comb, aes(ymin = mean - ci, ymax = mean + ci), width = 0.18, colour = "grey25") +
    geom_pointrange(data = reps, aes(ymin = mean - ci, ymax = mean + ci, colour = exp), shape = 21, fill = "white",
                    position = position_dodge(width = 0.55), size = 0.28, linewidth = 0.45, stroke = 0.9) +
    facet_wrap(~lab, ncol = ncol, scales = "free_y") +
    scale_fill_manual(values = state_col, labels = state_lab, name = NULL) +
    scale_colour_manual(values = exp_col, name = NULL) +
    scale_x_discrete(labels = c("--", "a", "b", "a+b")) +
    labs(x = NULL, y = "Change in log expression vs 'neither'\n(adjusted, mean +/- 95% CI)", title = title, subtitle = subtitle) +
    guides(fill = guide_legend(order = 1, nrow = 1), colour = guide_legend(order = 2))
}
mk <- function(x) x[, lab := sprintf("%s (a) x %s (b)\n-> %s  [%s%s]", TFa, TFb, target, shape, ifelse(core, ", core", ""))]
best <- P[core == TRUE][order(-score)][, head(.SD, 4), by = shape]
best[, shape := factor(shape, levels = levels(cnt$shape))]; setorder(best, shape, -score); mk(best)
g <- panel(best[, .(pair, target, lab)], ncol = 4, title = "Best reproducible combinatorial targets",
           subtitle = "Bars: Exp1 + Exp2 pooled.  Circles: each experiment separately.  a = first TF, b = second TF.  -- = cells with neither")
ggsave(file.path(FG, "gallery_best_by_shape.png"), g, width = 12, height = 3 * ceiling(nrow(best) / 4) + 1.2, dpi = 220, bg = "white")
ggsave(file.path(FG, "gallery_best_by_shape.pdf"), g, width = 12, height = 3 * ceiling(nrow(best) / 4) + 1.2)

# all proper cases, multi-page PDF (12 per page), ordered by shape then strength
allp <- copy(P); allp[, shape := factor(shape, levels = levels(cnt$shape))]; setorder(allp, shape, -core, -pair_specific, -score); mk(allp)
allp[, page := (seq_len(.N) - 1) %/% 12]
pdf(file.path(FG, "ALL_proper_combinatorial_targets.pdf"), width = 12, height = 9.5)
for (pg in unique(allp$page)) {
  k <- allp[page == pg]
  print(panel(k[, .(pair, target, lab)], ncol = 4,
              title = sprintf("Proper combinatorial targets: %s (page %d/%d)", paste(unique(k$shape), collapse = ", "), pg + 1, max(allp$page) + 1),
              subtitle = "Bars: Exp1 + Exp2 pooled. Circles: each experiment. Sorted by shape; core cases first, then strength."))
}
invisible(dev.off())
cat("done\n")
