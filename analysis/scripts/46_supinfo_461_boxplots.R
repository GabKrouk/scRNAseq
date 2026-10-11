#!/usr/bin/env Rscript
# Step 46: Sup Info - box plots of the 461 reproducible non-additive TF x TF interactions,
# Exp1 and Exp2 side by side in each panel, normalised to cells carrying neither TF.
suppressPackageStartupMessages({library(ggplot2); library(data.table)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
D <- file.path(base, "results", "combinatorial_push"); FG <- file.path(base, "figures")
source(file.path(base, "scripts", "paper_style.R")); theme_set(theme_tpc(7) + theme(panel.grid.major.x = element_blank(), legend.position = "top", plot.title = element_text(size = 8, face = "bold"), strip.text = element_text(size = 6, face = "plain")))
cols <- c(neither = "#1E7AB6", a = "#EA8535", b = "#519C3E", both = "#C2382F")
v <- fread(file.path(D, "supinfo_461_values.tsv.gz")); s <- fread(file.path(D, "supinfo_461_list.tsv"))
s[, TFa := sub("-.*", "", pair)][, TFb := sub(".*-", "", pair)]
s[, lab := sprintf("%s x %s -> %s\n%s, q=%.2g", TFa, TFb, target, type, q_meta)]
setorder(s, pair, q_meta); s[, ord := .I]
v <- merge(v, s[, .(pair, target, lab, ord)], by = c("pair", "target"))
v[, state := factor(state, levels = names(cols))]
v[, exp := factor(exp, levels = c("Exp1", "Exp2"))]
v[, lab := factor(lab, levels = s$lab)]
per <- 12; pages <- split(s$lab, ceiling(seq_along(s$lab) / per))
MODE <- Sys.getenv("MODE", "combined")   # combined: one box per state (both experiments pooled); Exp1 = dots, Exp2 = crosses
pp <- function(dd, i) ggplot(dd, aes(state, rel, fill = state)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_boxplot(width = 0.65, outlier.shape = NA, linewidth = 0.25, alpha = 0.75) +
  geom_point(aes(shape = exp), position = position_jitter(width = 0.15, height = 0, seed = 1), size = 0.4, stroke = 0.3, alpha = 0.85) +
  scale_shape_manual(values = c(Exp1 = 16, Exp2 = 4), name = NULL, labels = c("Exp1 (dots)", "Exp2 (crosses)")) +
  facet_wrap(~lab, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = cols, labels = c("neither", "TFa only", "TFb only", "TFa + TFb"), name = NULL) +
  scale_x_discrete(labels = c("--", "a", "b", "a+b")) +
  labs(x = NULL, y = "Expression change vs cells with neither TF (adjusted log expression)",
       title = sprintf("Supplementary Figure 5 (page %d of %d)", i, length(pages)),
       subtitle = "Each panel: target gene in cells with neither TF, TFa only, TFb only or both (>=1 junction read); box = Exp1 + Exp2 combined; dots = Exp1, crosses = Exp2 (10 random pools of cells per state and experiment).\nExpression adjusted for depth, load (incl. per-TF load modulation), cell type and the 21 other TFs; normalised to cells carrying neither TF. Interaction tested on TF dose (q = meta-analysis FDR).")
pdf(file.path(FG, "SupFig5_461_nonadditive_interactions.pdf"), width = 6.5, height = 8, family = "ArialMT", useDingbats = FALSE)
for (i in seq_along(pages)) print(pp(v[lab %in% pages[[i]]], i))
invisible(dev.off())
png(file.path(FG, "SupFig5_461_page1_preview.png"), width = 6.5, height = 8, units = "in", res = 200, family = "Arial", bg = "white"); print(pp(v[lab %in% pages[[1]]], 1)); invisible(dev.off())
cat(length(pages), "pages\n")
