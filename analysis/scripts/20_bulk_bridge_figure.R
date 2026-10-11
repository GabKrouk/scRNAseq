#!/usr/bin/env Rscript
# Step 20: bulk DoubleTARGET reproducibility and single-cell vs bulk agreement.
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(gridExtra)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
D <- file.path(base, "results", "bulk_bridge")
theme_set(theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 12),
  plot.subtitle = element_text(colour = "grey35", size = 9), strip.text = element_text(face = "bold")))
b <- fread(file.path(D, "bulk_exp1_vs_exp2_effects.tsv"))
b[, effect := factor(effect, levels = c("CRF4 alone", "ERF5 alone", "ERF5+CRF4", "interaction"))]
st <- b[, .(rho = cor(E1, E2, method = "spearman"),
            top = {i <- order(-abs(E1))[1:200]; mean(sign(E1[i]) == sign(E2[i]))}), by = effect]
st[, lab := sprintf("rho = %.2f\ntop-200 genes same sign: %d%%", rho, round(100 * top))]
pA <- ggplot(b, aes(E1, E2)) + geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(size = 0.25, alpha = 0.25, colour = "#33415C") + geom_text(data = st, aes(x = -Inf, y = Inf, label = lab), hjust = -0.05, vjust = 1.2, size = 3) +
  facet_wrap(~effect, nrow = 1) + coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3)) +
  labs(x = "Bulk experiment 1 (log2 fold change vs EV)", y = "Bulk experiment 2",
       title = "A  Bulk DoubleTARGET (Fig 7 data): does each effect reproduce between the two bulk experiments?")
s <- fread(file.path(D, "sc_vs_bulk_effects.tsv"))
s[, effect := factor(effect, levels = c("CRF4", "ERF5", "inter"), labels = c("CRF4 alone", "ERF5 alone", "ERF5 x CRF4 interaction"))]
s[, singlecell := factor(singlecell, levels = c("Exp1", "Exp2", "meta"), labels = c("sc Exp1", "sc Exp2", "sc combined"))]
pB <- ggplot(s, aes(singlecell, sign_agreement_bulk_DEGs * 100, fill = singlecell)) + geom_col(width = 0.7) +
  geom_hline(yintercept = 50, linetype = 2) + facet_wrap(~effect) +
  geom_text(aes(label = sprintf("%d%%\np=%s", round(sign_agreement_bulk_DEGs * 100), formatC(binom_p, format = "g", digits = 1))), vjust = -0.2, size = 2.8) +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "grey30"), guide = "none") + coord_cartesian(ylim = c(0, 85)) +
  labs(x = NULL, y = "% of bulk-responsive genes (p<0.05)\nwith the same direction in single cells",
       title = "B  Single-cell vs bulk DoubleTARGET", subtitle = sprintf("Bulk genes tested: CRF4 %d, ERF5 %d, interaction %d", s[effect=="CRF4 alone"]$n_bulk_p05[1], s[effect=="ERF5 alone"]$n_bulk_p05[1], s[effect=="ERF5 x CRF4 interaction"]$n_bulk_p05[1]))
g <- arrangeGrob(pA, pB, nrow = 2, heights = c(1, 0.95))
ggsave(file.path(base, "figures", "bulk_vs_singlecell_ERF5_CRF4.png"), g, width = 13, height = 8.5, dpi = 200, bg = "white")
ggsave(file.path(base, "figures", "bulk_vs_singlecell_ERF5_CRF4.pdf"), g, width = 13, height = 8.5)
cat("ok\n")
