#!/usr/bin/env Rscript
# Step 49: main-figure panels for the revised manuscript (raw counts, QC cells).
#   Fig 3  TF landscape          Fig 4  reproducibility          Fig 5  non-additive TF x TF interactions
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(gridExtra); library(jsonlite)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
R <- function(...) file.path(base, "results", ...)
FG <- file.path(base, "figures", "paper"); dir.create(FG, showWarnings = FALSE, recursive = TRUE)
source(file.path(base, "scripts", "paper_style.R"))
save <- function(g, name, w, h) save_tpc(g, FG, name, w, h)
col_exp <- c(Exp1 = "#0072B2", Exp2 = "#E69F00")
theme_set(theme_tpc(7))

# ================================================================ Fig 3: TF landscape
det <- fread(R("paper", "tf_detection.tsv"))
ord <- det[, .(m = mean(frac)), by = tf][order(m)]$tf
det[, tf := factor(tf, levels = ord)]
p3a <- ggplot(det, aes(frac * 100, tf, colour = exp, shape = exp)) + geom_point(size = 1.4) +
  scale_colour_manual(values = col_exp, name = NULL) + scale_shape_manual(values = c(Exp1 = 16, Exp2 = 17), name = NULL) +
  labs(x = "% of cells with >= 1 GR-TF junction read", y = NULL, title = "A") +
  theme(legend.position = "inside", legend.position.inside = c(0.8, 0.1), legend.background = element_rect(fill = "white", colour = NA))
nt <- fread(R("paper", "tfs_per_cell.tsv"))
med <- nt[, .(med = n_tfs[which(cumsum(frac) >= 0.5)[1]]), by = exp]
p3b <- ggplot(nt, aes(n_tfs, frac * 100, fill = exp)) + geom_col(position = position_dodge(0.85), width = 0.85) +
  geom_vline(data = med, aes(xintercept = med, colour = exp), linetype = 2, show.legend = FALSE) +
  scale_fill_manual(values = col_exp, name = NULL) + scale_colour_manual(values = col_exp) +
  labs(x = "Number of GR-TFs detected in the cell (of 23)", y = "% of cells", title = "B",
       subtitle = sprintf("Median %d (Exp1, %s cells) and %d (Exp2, %s cells); dashed lines = medians",
                          med[exp == "Exp1"]$med, format(sum(nt[exp == "Exp1"]$n_cells), big.mark = ","),
                          med[exp == "Exp2"]$med, format(sum(nt[exp == "Exp2"]$n_cells), big.mark = ","))) +
  theme(legend.position = "inside", legend.position.inside = c(0.85, 0.85))
co <- fread(R("paper", "tf_codetection.tsv"))
co[, a := factor(a, levels = rev(ord))][, b := factor(b, levels = rev(ord))]
co <- co[(exp == "Exp1" & as.integer(a) > as.integer(b)) | (exp == "Exp2" & as.integer(a) < as.integer(b)) | a == b]
co[a == b, log2_oe := NA]
p3c <- ggplot(co[a != b], aes(a, b, fill = log2_oe)) + geom_tile(colour = "white", linewidth = 0.2) +
  geom_tile(data = co[a == b & exp == "Exp1"], fill = "grey30") +
  scale_fill_gradient2(low = "#1E7AB6", mid = "white", high = "#C2382F", midpoint = 0, name = "log2\nobserved /\nexpected") +
  annotate("text", x = 21, y = 3, label = "Exp1", fontface = "bold", size = 7 * PT) +
  annotate("text", x = 3, y = 21, label = "Exp2", fontface = "bold", size = 7 * PT) +
  labs(x = NULL, y = NULL, title = "C") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6), axis.text.y = element_text(size = 6), panel.grid = element_blank(), legend.key.height = unit(14, "pt"))
save(arrangeGrob(arrangeGrob(p3a, p3b, ncol = 2, widths = c(0.8, 1)), p3c, ncol = 1, heights = c(0.85, 1.15)), "Fig3_TF_landscape", 6.5, 7.2)

# ================================================================ Fig 4: reproducibility
pt <- fread(R("consensus_log", "per_TF_summary.tsv"))
pt[, tf := factor(tf, levels = pt[order(n_consensus)]$tf)]
m4 <- melt(pt, id.vars = "tf", measure.vars = c("n_consensus_up", "n_consensus_down"))
m4[, variable := factor(variable, labels = c("induced", "repressed"))]
p4a <- ggplot(m4, aes(value, tf, fill = variable)) + geom_col(width = 0.75) +
  scale_fill_manual(values = c(induced = "#C2382F", repressed = "#1E7AB6"), name = NULL) +
  labs(x = "Targets reproduced in both experiments", y = NULL, title = "A",
       subtitle = sprintf("%s TF -> gene effects (FDR 5%%, p < 0.05 and same sign in each experiment)", format(sum(pt$n_consensus), big.mark = ","))) +
  theme(legend.position = "inside", legend.position.inside = c(0.75, 0.2))
me <- fread(R("consensus_log", "main_effects_all_significant.tsv.gz"))
sel <- me[regulator %in% c("NAC4", "CRF4", "COL5", "WRKY54")]
sel[, regulator := factor(regulator, levels = c("NAC4", "CRF4", "COL5", "WRKY54"))]
rr <- sel[, .(r = cor(beta_Exp1, beta_Exp2), n = sum(consensus)), by = regulator]
p4b <- ggplot(sel, aes(beta_Exp1, beta_Exp2)) + geom_hline(yintercept = 0, colour = "grey75") + geom_vline(xintercept = 0, colour = "grey75") +
  geom_point(aes(colour = consensus), size = 0.3, alpha = 0.6) + facet_wrap(~regulator, nrow = 1) +
  geom_text(data = rr, aes(x = -0.32, y = 0.4, label = sprintf("%d reproduced", n)), hjust = 0, size = 6 * PT) +
  scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "#C2382F"), labels = c("significant in one experiment", "reproduced"), name = NULL) +
  coord_cartesian(xlim = c(-0.35, 0.45), ylim = c(-0.35, 0.45)) +
  labs(x = "Effect in Exp1 (log expression per unit of log TF dose)", y = "Effect in Exp2", title = "D") + theme(legend.position = "bottom")
sj <- fromJSON(R("consensus_log", "summary.json")); fj <- fromJSON(R("final", "summary.json"))
pnm <- fread(R("consensus_log", "permutation_null.tsv")); pni <- fread(R("final", "permutation_null.tsv"))
nd <- rbind(data.table(kind = "TF -> gene\neffects", null = pnm$main_n_consensus, obs = sj$real$main$n_consensus),
            data.table(kind = "TF x TF\ninteractions (strict)", null = pni$inter_strict, obs = fj$real$inter_strict),
            data.table(kind = "TF x TF\ninteractions (medium)", null = pni$inter_medium, obs = fj$real$inter_medium))
ns <- nd[, .(obs = obs[1], mean = mean(null), max = max(null), nperm = .N), by = kind]
p4c <- ggplot(ns, aes(kind)) + geom_col(aes(y = obs + 1), fill = "#C2382F", width = 0.55) +
  geom_point(data = nd, aes(y = null + 1), position = position_jitter(width = 0.12, height = 0, seed = 2), colour = "grey35", size = 0.5) +
  geom_text(aes(y = obs + 1, label = format(obs, big.mark = ",")), vjust = -0.4, size = 6.5 * PT) +
  scale_y_log10(limits = c(1, 5000)) +
  labs(x = NULL, y = "Calls + 1 (log scale)", title = "B")
pw <- fread(R("power_vs_cells.tsv"))
pw[, lab := ifelse(subsampled == "Exp2", "Exp2 down-sampled, Exp1 complete", "Exp1 down-sampled, Exp2 complete")]
p4d <- ggplot(pw, aes(n_cells, n_consensus, colour = lab)) + geom_point(size = 1) + stat_summary(fun = mean, geom = "line", linewidth = 0.4) +
  scale_colour_manual(values = c("#0072B2", "#E69F00"), name = NULL) +
  labs(x = "Cells in the down-sampled experiment", y = "Reproduced TF -> gene effects", title = "C") +
  theme(legend.position = "inside", legend.position.inside = c(0.7, 0.18), legend.background = element_rect(fill = "white", colour = NA))
tc <- fread(R("paper", "transfer_comparison.tsv"))
tc[, model := factor(model, levels = unique(model))]
tc[, direction := factor(direction, labels = c("learned in Exp1, tested in Exp2", "learned in Exp2, tested in Exp1"))]
p4e <- ggplot(tc, aes(model, n_sig, fill = direction)) + geom_col(position = position_dodge(0.8), width = 0.75) +
  scale_fill_manual(values = unname(col_exp), name = NULL) +
  labs(x = NULL, y = "Genes predicted in the other\nexperiment (FDR 5%)", title = "E") +
  theme(legend.position = "top", axis.text.x = element_text(size = 7))
# panel E (linear vs CART transfer) moved to Supplemental Figure 4 (script 59)
g4 <- arrangeGrob(arrangeGrob(p4a, arrangeGrob(p4c, p4d, ncol = 1), ncol = 2, widths = c(0.9, 1)), p4b, ncol = 1, heights = c(1.25, 0.8))
save(g4, "Fig4_reproducibility", 6.5, 6.6)

# ================================================================ Fig 5: non-additive interactions
pp <- fread(R("final", "per_pair.tsv")); pp[, a := sub("-.*", "", pair)][, b := sub(".*-", "", pair)]
tfo <- levels(pt$tf)
mm <- rbind(pp[, .(a, b, n = interactions_medium, s = interactions_strict)], pp[, .(a = b, b = a, n = interactions_medium, s = interactions_strict)])
mm[, a := factor(a, levels = rev(tfo))][, b := factor(b, levels = rev(tfo))]
p5a <- ggplot(mm, aes(a, b, fill = n)) + geom_tile(colour = "white", linewidth = 0.2) +
  geom_text(data = mm[n >= 5], aes(label = n), size = 6 * PT) +
  scale_fill_gradient(low = "#F5F5F5", high = "#C2382F", name = "Reproduced\ninteractions", trans = "sqrt") +
  labs(x = NULL, y = NULL, title = "A") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6), axis.text.y = element_text(size = 6), panel.grid = element_blank(),
        legend.position = "bottom", legend.key.width = unit(14, "pt"), legend.key.height = unit(6, "pt"))
d <- fread(R("combinatorial_push", "nonadditivity_medium.tsv"))
d[, type := ifelse(nonadd > 0, "synergy", "antagonism")]
lim <- quantile(abs(c(d$additive, d$dAB)), 0.99)
p5b <- ggplot(d, aes(additive, dAB, colour = type)) + geom_abline(linetype = 2, colour = "grey40") +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(size = 0.4, alpha = 0.7) + facet_wrap(~exp) + coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  scale_colour_manual(values = c(synergy = "#C2382F", antagonism = "#1E7AB6"), name = NULL) +
  labs(x = "Sum of single-TF effects (additive expectation)", y = "Observed effect with both TFs",
       title = "C") + theme(legend.position = "bottom")
w <- dcast(d, pair + target ~ exp, value.var = "nonadd"); setnames(w, c("pair", "target", "e1", "e2"))
p5c <- ggplot(w, aes(e1, e2)) + geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(aes(colour = e1 > 0), size = 0.4, alpha = 0.75) + geom_abline(linetype = 2, colour = "grey50") +
  scale_colour_manual(values = c(`TRUE` = "#C2382F", `FALSE` = "#1E7AB6"), guide = "none") +
  labs(x = "Deviation from additivity, Exp1", y = "Deviation, Exp2", title = "D",
       subtitle = sprintf("Same sign %d%%, r = %.2f", round(100 * mean(sign(w$e1) == sign(w$e2))), cor(w$e1, w$e2)))
cols <- c(neither = "#1E7AB6", a = "#EA8535", b = "#519C3E", both = "#C2382F")
v <- fread(R("combinatorial_push", "supinfo_461_values.tsv.gz")); s <- fread(R("combinatorial_push", "supinfo_461_list.tsv"))
# Examples chosen for the clarity of the ANOVA-type interaction: same shape and sign in both experiments,
# ranked by the z-score of the non-additive deviation in the pooled data (results/paper/fig5d_candidates.tsv)
pick <- data.table(
  cls = rep(c("Emergent: only together", "Synergy: more than the sum", "Competition: cancel out", "Reversal: opposite together"), each = 2),
  pair = c("NAC4-RAV1", "HB6-TGA1", "GATA17-TGA4", "BEE2-NAC4", "CDF1-NAC4", "LBD38-TGA4", "NAP-VRN1", "HHO3-WRKY18"),
  target = c("AT4G28703", "AT2G20500", "AT1G69680", "AT3G06700", "AT3G04090", "AT1G19770", "AT1G20140", "AT3G54880"),
  sym = c("", "", "", "", "SIP1A", "PUP14", "SK4", ""))
pick[, a := sub("-.*", "", pair)][, b := sub(".*-", "", pair)]
pick[, lab := sprintf("%s\n%s x %s -> %s", cls, a, b, ifelse(sym == "", target, sym))]
vv <- merge(v, pick[, .(pair, target, lab, a, b)], by = c("pair", "target"))
vv[, xl := ifelse(state == "neither", "none", ifelse(state == "a", a, ifelse(state == "b", b, "both")))]
vv[, x := factor(paste(lab, state), levels = as.vector(t(outer(pick$lab, names(cols), paste))))]
vv[, state := factor(state, levels = names(cols))][, lab := factor(lab, levels = pick$lab)]
xlabs <- setNames(vv$xl, as.character(vv$x))
addv <- vv[, .(m = mean(rel)), by = .(lab, state, x)]
addv <- addv[, .(add = m[state == "a"] + m[state == "b"] - m[state == "neither"], x = x[state == "both"]), by = lab]
p5d <- ggplot(vv, aes(x, rel, fill = state)) + geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_boxplot(width = 0.65, outlier.shape = NA, linewidth = 0.25, alpha = 0.75) +
  geom_point(aes(shape = exp), position = position_jitter(width = 0.15, height = 0, seed = 1), size = 0.35, stroke = 0.3) +
  geom_errorbar(data = addv, aes(x = x, ymin = add, ymax = add), inherit.aes = FALSE, width = 0.8, linetype = "22", linewidth = 0.45) +
  scale_shape_manual(values = c(Exp1 = 16, Exp2 = 4), name = NULL, labels = c("Exp1 (dots)", "Exp2 (crosses)")) +
  facet_wrap(~lab, nrow = 2, scales = "free") +
  scale_fill_manual(values = cols, labels = c("neither TF", "TFa only", "TFb only", "both TFs"), name = NULL) +
  scale_x_discrete(labels = xlabs) +
  labs(x = NULL, y = "Change vs cells with neither TF (adjusted log expression)",
       title = "E") +
  theme(legend.position = "bottom", strip.text = element_text(size = 6, face = "plain"), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))
source(file.path(base, "scripts", "64_network_panel.R"))
p5net <- network_panel(base, "B")
g5 <- arrangeGrob(arrangeGrob(p5a, p5net, ncol = 2, widths = c(1, 1.05)), arrangeGrob(p5b, p5c, ncol = 2, widths = c(1.45, 1)), p5d,
                  ncol = 1, heights = c(1.5, 0.95, 1.45))
save(g5, "Fig5_nonadditive_interactions", 6.5, 8)
cat("paper figures written\n")
