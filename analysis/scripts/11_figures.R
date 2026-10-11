#!/usr/bin/env Rscript
# Step 11: figures for the reproducibility report (ggplot2, Okabe-Ito colours).
suppressPackageStartupMessages({library(ggplot2); library(data.table)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
R <- function(...) file.path(base, "results", ...)
fig <- function(p, name, w = 7, h = 4.5) {
  ggsave(file.path(base, "figures", paste0(name, ".png")), p, width = w, height = h, dpi = 200, bg = "white")
  ggsave(file.path(base, "figures", paste0(name, ".pdf")), p, width = w, height = h)
}
col_exp <- c(Exp1 = "#0072B2", Exp2 = "#E69F00")
theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95")))

# ---- Fig 1: TF detection --------------------------------------------------
td <- fread(R("tf_detection_summary.tsv"))
td <- td[tf %in% fread(R("consensus_ge1", "per_TF_summary.tsv"))$tf]
ord <- td[exp == "Exp1"][order(frac_ge1)]$tf
td[, tf := factor(tf, levels = ord)]
p1a <- ggplot(td, aes(frac_ge1 * 100, tf, colour = exp)) + geom_point(size = 2.4) +
  scale_colour_manual(values = col_exp, name = NULL) +
  labs(x = "% of cells with >=1 GR-TF junction read", y = NULL, title = "A  TF detection per experiment")
nt <- fread(R("plot_ntf_per_cell.tsv"))
p1b <- ggplot(nt, aes(n_tfs, frac * 100, fill = exp)) + geom_col(position = position_dodge(0.8), width = 0.8) +
  scale_fill_manual(values = col_exp, name = NULL) +
  labs(x = "Number of the 23 GR-TFs detected in the cell", y = "% of cells", title = "B  TFs detected per cell")
fig(p1a, "fig1a_tf_detection", 5, 5.5); fig(p1b, "fig1b_tfs_per_cell", 6, 3.5)

# ---- Fig 2: calibration of TF calling / covariates -------------------------
cl <- fread(R("calibration_linear.tsv"))
cl[, scheme := factor(scheme, levels = c("ge1", "ge2", "ge3", "ge5", "log", "share"),
                      labels = c(">=1 read", ">=2 reads", ">=3 reads", ">=5 reads", "log dose", "share of TF reads"))]
cl[, cov := factor(cov, levels = c("none", "depth", "depth+load"), labels = c("no covariates", "+ depth", "+ depth + TF load/RFP"))]
p2a <- ggplot(cl, aes(scheme, r_t, fill = cov)) + geom_col(position = position_dodge(0.85), width = 0.8) +
  geom_point(aes(y = null_r_t_mean, group = cov), position = position_dodge(0.85), shape = 4, size = 2) +
  scale_fill_manual(values = c("#BBBBBB", "#56B4E9", "#009E73"), name = NULL) +
  labs(x = "GR-TF encoding", y = "Exp1 vs Exp2 correlation of\nTF x gene t-statistics",
       title = "Reproducibility of TF effects by TF-calling rule (x = permutation null)") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "top")
fig(p2a, "fig2_calibration_tf_calling", 7.5, 4.5)

# ---- Fig 3: per-TF reproducibility & consensus targets ---------------------
pt <- fread(R("consensus_ge1", "per_TF_summary.tsv"))
pt[, tf := factor(tf, levels = pt[order(n_consensus)]$tf)]
m3 <- melt(pt, id.vars = "tf", measure.vars = c("n_consensus_up", "n_consensus_down"))
m3[, variable := factor(variable, labels = c("induced", "repressed"))]
p3a <- ggplot(m3, aes(value, tf, fill = variable)) + geom_col() +
  scale_fill_manual(values = c(induced = "#D55E00", repressed = "#0072B2"), name = NULL) +
  labs(x = "Consensus targets (reproduced in both experiments)", y = NULL, title = "A  Reproducible targets per TF")
p3b <- ggplot(pt, aes(r_t_Exp1_Exp2, tf)) + geom_vline(xintercept = 0, colour = "grey60") +
  geom_point(size = 2.4, colour = "#009E73") +
  labs(x = "Genome-wide correlation of the TF's effects, Exp1 vs Exp2", y = NULL, title = "B  Per-TF agreement")
fig(p3a, "fig3a_consensus_targets_per_tf", 5.5, 5.5); fig(p3b, "fig3b_per_tf_agreement", 5, 5.5)

# ---- Fig 4: observed vs permutation null ----------------------------------
pn <- fread(R("consensus_ge1", "permutation_null.tsv"))
sm <- jsonlite_free <- NULL
obs <- c(main = 1972, inter = 1159)
sj <- readLines(R("consensus_ge1", "summary.json"))
get_n <- function(kind) { i <- grep(paste0('"', kind, '": \\{'), sj)[1]; j <- grep('"n_consensus"', sj); as.numeric(gsub("[^0-9]", "", sj[j[j > i][1]])) }
obs <- c(main = get_n("main"), inter = get_n("inter"))
d4 <- rbind(data.table(kind = "TF -> gene main effects", null = pn$main_n_consensus, obs = obs["main"]),
            data.table(kind = "TF x TF interactions", null = pn$inter_n_consensus, obs = obs["inter"]))
p4 <- ggplot(d4, aes(null + 1)) + geom_histogram(bins = 12, fill = "grey70") +
  geom_vline(aes(xintercept = obs + 1), colour = "#D55E00", linewidth = 1) + facet_wrap(~kind, scales = "free_y") +
  scale_x_log10() + labs(x = "Consensus edges + 1 (log scale)", y = "Permutations",
                         title = "Consensus edges: observed (orange) vs 50 TF-label permutations (grey)")
fig(p4, "fig4_permutation_null", 8, 3.5)

# ---- Fig 5: cross-experiment prediction -----------------------------------
tm <- fread(R("transfer_model_comparison.tsv"), skip = 3, header = FALSE)
setnames(tm, c("model", "r_1to2", "r_2to1", "n_1to2", "n_2to1"))
tm <- tm[!grepl("null", model)]
tm[, model := factor(model, levels = c("CART cart_cells_unadjusted", "CART stump", "CART tuned", "linear additive", "linear interactions"),
     labels = c("CART-CELLS settings,\nno covariate correction", "CART stump\n(1 split, corrected)",
                "CART depth 3, 10% leaf\n(corrected)", "Linear additive\n(corrected)", "Linear additive +\nTF x TF interactions"))]
m5 <- melt(tm, id.vars = "model", measure.vars = c("n_1to2", "n_2to1"))
m5[, variable := factor(variable, labels = c("trained Exp1 -> predicts Exp2", "trained Exp2 -> predicts Exp1"))]
p5 <- ggplot(m5, aes(model, value, fill = variable)) + geom_col(position = position_dodge(0.8), width = 0.75) +
  scale_fill_manual(values = c("#0072B2", "#E69F00"), name = NULL) +
  labs(x = NULL, y = "Genes predicted in the other experiment\n(BH q < 0.05)",
       title = "Do models learned in one experiment predict the other? (permuted labels: <= 35 genes)") +
  theme(legend.position = "top", axis.text.x = element_text(size = 8))
fig(p5, "fig5_cross_experiment_prediction", 8.5, 4.5)

fa <- fread(R("cart_faithful_TF1_agreement.tsv"))
m5b <- melt(fa, id.vars = c("min_log10p", "n_genes"), measure.vars = c("obs", "chance"))
m5b[, lab := paste0("-log10 p >= ", min_log10p, "\n(n = ", n_genes, ")")]
p5b <- ggplot(m5b, aes(lab, value * 100, fill = variable)) + geom_col(position = position_dodge(0.8), width = 0.75) +
  scale_fill_manual(values = c(obs = "#CC79A7", chance = "grey70"), labels = c("observed", "expected by chance"), name = NULL) +
  labs(x = "Genes whose TF1 split is significant in both experiments", y = "% genes with the same TF1",
       title = "CART-CELLS as published: is TF1 the same in Exp1 and Exp2?") + theme(legend.position = "top")
fig(p5b, "fig5b_cart_cells_tf1_agreement", 7, 4)

# ---- Fig 6: power vs cells ------------------------------------------------
pw <- fread(R("power_vs_cells.tsv"))
pw[, lab := ifelse(subsampled == "Exp2", "Exp2 down-sampled (Exp1 complete)", "Exp1 down-sampled (Exp2 complete)")]
p6 <- ggplot(pw, aes(n_cells, n_consensus, colour = lab)) + geom_point(size = 2) +
  stat_summary(fun = mean, geom = "line") + scale_colour_manual(values = c("#0072B2", "#E69F00"), name = NULL) +
  labs(x = "Cells in the down-sampled experiment", y = "Consensus TF -> gene edges",
       title = "Reproducible edges keep increasing with cell number") + theme(legend.position = "top")
fig(p6, "fig6_power_vs_cells", 6.5, 4.2)

# ---- Fig 7: TF x TF interaction matrix ------------------------------------
pp <- fread(R("consensus_ge1", "per_pair_summary.tsv"))
tfo <- levels(pt$tf)
mm <- rbind(pp[, .(a = tf_a, b = tf_b, n = n_consensus_interaction)], pp[, .(a = tf_b, b = tf_a, n = n_consensus_interaction)])
mm[, a := factor(a, levels = rev(tfo))]; mm[, b := factor(b, levels = rev(tfo))]
p7 <- ggplot(mm, aes(a, b, fill = n)) + geom_tile(colour = "white") +
  geom_text(data = mm[n >= 10], aes(label = n), size = 2.6) +
  scale_fill_gradient(low = "#F7F7F7", high = "#D55E00", name = "Consensus\ninteractions", trans = "sqrt") +
  labs(x = NULL, y = NULL, title = "Reproducible TF x TF interactions (both experiments, FDR 0.2%)") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), panel.grid = element_blank())
fig(p7, "fig7_interaction_matrix", 7.5, 6.5)

# ---- Fig 8: NAC4 & CRF4 effects Exp1 vs Exp2 --------------------------------
me <- fread(R("consensus_ge1", "main_effects_all_significant.tsv.gz"))
sel <- me[regulator %in% c("NAC4", "CRF4", "ERF5", "TGA1")]
sel[, regulator := factor(regulator, levels = c("NAC4", "CRF4", "TGA1", "ERF5"))]
p8 <- ggplot(sel, aes(beta_Exp1, beta_Exp2, colour = consensus)) + geom_hline(yintercept = 0, colour = "grey70") +
  geom_vline(xintercept = 0, colour = "grey70") + geom_point(size = 0.7, alpha = 0.6) +
  facet_wrap(~regulator, nrow = 1) + scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "#D55E00"), name = "Consensus") +
  labs(x = "Effect in Exp1 (log expression)", y = "Effect in Exp2",
       title = "TF effects on genes significant in either experiment") + coord_cartesian(xlim = c(-0.6, 0.6), ylim = c(-0.6, 0.6))
fig(p8, "fig8_effects_exp1_vs_exp2", 10, 3.4)
cat("figures written\n")
