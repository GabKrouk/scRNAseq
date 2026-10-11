#!/usr/bin/env Rscript
# Step 43: non-additive (ANOVA-type) TF x TF effects on raw counts.
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(gridExtra)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
D <- file.path(base, "results", "combinatorial_push"); FG <- file.path(base, "figures")
theme_set(theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 11.5),
  plot.subtitle = element_text(colour = "grey35", size = 8.5), strip.text = element_text(face = "bold", size = 8.5),
  panel.border = element_rect(fill = NA, colour = "grey88")))
SET <- Sys.getenv("SET", "strict")
d <- fread(file.path(D, paste0("nonadditivity_", SET, ".tsv")))
d[, type := ifelse(nonadd > 0, "synergy (more than the sum)", "antagonism (less than the sum)")]
d[, exp := paste("Single-cell", exp)]
lim <- quantile(abs(c(d$additive, d$dAB)), 0.99)
pA <- ggplot(d, aes(additive, dAB, colour = type)) + geom_abline(linetype = 2, colour = "grey40") +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(size = 1.6, alpha = 0.75) + facet_wrap(~exp) + coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  scale_colour_manual(values = c("synergy (more than the sum)" = "#C2382F", "antagonism (less than the sum)" = "#1E7AB6"), name = NULL) +
  labs(x = "Expected if additive: effect(TFa alone) + effect(TFb alone)", y = "Observed effect with both TFs",
       title = "A  The combined effect is not the sum of the single-TF effects",
       subtitle = sprintf("%d %s interactions (dose-based ANOVA-type interaction term, replicated in both experiments). Dashed line = additivity.", uniqueN(d[, paste(pair, target)]), ifelse(SET == "strict", "strict (q<0.05)", "strict + medium (q<0.1)"))) +
  theme(legend.position = "bottom")
w <- dcast(d, pair + target ~ exp, value.var = "nonadd"); setnames(w, c("pair", "target", "e1", "e2"))
pB <- ggplot(w, aes(e1, e2)) + geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(aes(colour = e1 > 0), size = 1.6, alpha = 0.8) + geom_abline(linetype = 2, colour = "grey50") +
  scale_colour_manual(values = c(`TRUE` = "#C2382F", `FALSE` = "#1E7AB6"), guide = "none") +
  labs(x = "Non-additivity, Exp1", y = "Non-additivity, Exp2", title = "B  Same deviation in both experiments",
       subtitle = sprintf("Same sign %d%%, r = %.2f", round(100 * mean(sign(w$e1) == sign(w$e2))), cor(w$e1, w$e2)))
l <- fread(file.path(D, "dose_landscapes_examples.tsv"))
l[, doseA := factor(doseA, levels = c("0", "low", "high"))]; l[, doseB := factor(doseB, levels = c("0", "low", "high"))]
l[, lab := paste0(sub("-.*", "", pair), " x ", sub(".*-", "", pair), " -> ", target)]
l[, exp := paste("Exp", sub("Exp", "", exp))]
keep <- unique(l$lab)[1:8]
pC <- ggplot(l[lab %in% keep], aes(doseA, doseB, fill = mean)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", mean)), size = 2.4) +
  facet_wrap(~ lab + exp, ncol = 8, labeller = function(x) list(paste0(x$lab, "\n", x$exp))) +
  scale_fill_gradient2(low = "#1E7AB6", mid = "white", high = "#C2382F", name = "change vs\nno TF") +
  labs(x = "Dose of TFa (junction reads: 0 / low / high)", y = "Dose of TFb", title = "C  Dose landscapes of example targets",
       subtitle = "Adjusted mean expression relative to cells with neither TF; the corner 'high/high' departs from the sum of the edges") +
  theme(strip.text = element_text(size = 6.5), axis.text = element_text(size = 7), panel.grid = element_blank())
g <- arrangeGrob(arrangeGrob(pA, pB, ncol = 2, widths = c(1.75, 1)), pC, ncol = 1, heights = c(1, 0.95))
ggsave(file.path(FG, paste0("nonadditive_TF_combinations", ifelse(SET == "strict", "", paste0("_", SET)), ".png")), g, width = 15, height = 11, dpi = 180, bg = "white")
ggsave(file.path(FG, paste0("nonadditive_TF_combinations", ifelse(SET == "strict", "", paste0("_", SET)), ".pdf")), g, width = 15, height = 11)
cat("ok\n")
