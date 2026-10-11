#!/usr/bin/env Rscript
# Step 59: Supplemental Figure 4 - linear models vs CART: cross-experiment prediction (former Fig 4E),
# stability of the CART root TF between experiments, and CART preference for TFs of intermediate prevalence (R1 M3).
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(gridExtra)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
R <- function(...) file.path(base, "results", ...); FG <- file.path(base, "figures", "paper")
col_exp <- c(Exp1 = "#0072B2", Exp2 = "#E69F00")
source(file.path(base, "scripts", "paper_style.R")); theme_set(theme_tpc(7))
tc <- fread(R("paper", "transfer_comparison.tsv"))
tc[, model := factor(model, levels = unique(model))]
tc[, direction := factor(direction, labels = c("learned in Exp1, tested in Exp2", "learned in Exp2, tested in Exp1"))]
pA <- ggplot(tc, aes(model, n_sig, fill = direction)) + geom_col(position = position_dodge(0.8), width = 0.75) +
  scale_fill_manual(values = unname(col_exp), name = NULL) +
  labs(x = NULL, y = "Genes predicted in the other\nexperiment (FDR 5%)", title = "A",
       subtitle = "Prediction of each gene from the TF profile of each cell, trained on one experiment and tested on the other; shuffled TF labels: 0 genes") +
  theme(legend.position = "top", axis.text.x = element_text(size = 7))
a <- fread(cmd = paste("gzip -dc", shQuote(R("cart", "cart_cells_faithful_Exp1.tsv.gz"))))[, .(gene, TF1)]
b <- fread(cmd = paste("gzip -dc", shQuote(R("cart", "cart_cells_faithful_Exp2.tsv.gz"))))[, .(gene, TF1)]
m <- merge(a, b, by = "gene", suffixes = c("_1", "_2"))[!is.na(TF1_1) & !is.na(TF1_2) & TF1_1 != "" & TF1_2 != ""]
f1 <- m[, .N / nrow(m), by = TF1_1]; f2 <- m[, .N / nrow(m), by = TF1_2]
chance <- sum(merge(f1, f2, by.x = "TF1_1", by.y = "TF1_2")[, V1.x * V1.y])
ag <- data.table(what = factor(c("observed", "expected by chance"), levels = c("observed", "expected by chance")), pct = 100 * c(mean(m$TF1_1 == m$TF1_2), chance))
pB <- ggplot(ag, aes(what, pct, fill = what)) + geom_col(width = 0.6) + geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.4, size = 6.5 * PT) +
  scale_fill_manual(values = c("#C2382F", "grey65"), guide = "none") + scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "% of genes", title = "B",
       subtitle = sprintf("Same first-split TF in Exp1 and Exp2\n(CART-CELLS settings, %s genes)", format(nrow(m), big.mark = ",")))
det <- fread(R("paper", "tf_detection.tsv"))
cnt <- rbind(a[, .(n = .N), by = .(tf = TF1)][, exp := "Exp1"], b[, .(n = .N), by = .(tf = TF1)][, exp := "Exp2"])
d <- merge(det, cnt, by = c("exp", "tf"), all.x = TRUE); d[is.na(n), n := 0]
rho <- d[, .(r = cor(abs(frac - 0.5), n, method = "spearman")), by = exp]
pC <- ggplot(d, aes(frac * 100, n, colour = exp)) + geom_point(size = 2) +
  geom_text(data = d[n > quantile(n, 0.85)], aes(label = tf), size = 6 * PT, vjust = -0.8, show.legend = FALSE) +
  facet_wrap(~exp, scales = "free_y") + scale_colour_manual(values = col_exp, guide = "none") +
  labs(x = "% of cells carrying the TF", y = "Genes for which the TF is the first CART split",
       title = "C",
       subtitle = sprintf("Spearman rho with distance of prevalence from 50%%: %.2f (Exp1), %.2f (Exp2)", rho[exp == "Exp1"]$r, rho[exp == "Exp2"]$r))
g <- arrangeGrob(pA, arrangeGrob(pB, pC, ncol = 2, widths = c(0.55, 1)), ncol = 1, heights = c(1, 1))
save_tpc(g, FG, "SupFig_linear_vs_CART", 6.5, 5.6)
cat("root agreement", round(mean(m$TF1_1 == m$TF1_2), 3), "chance", round(chance, 3), "\n"); print(rho)
