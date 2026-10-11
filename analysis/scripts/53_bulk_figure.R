#!/usr/bin/env Rscript
# Step 53: Fig 6 (bulk ERF5 x CRF4 DoubleTARGET) rebuilt with FDR-controlled gene classes (R3 #7).
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(gridExtra); library(jsonlite)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
FG <- file.path(base, "figures", "paper"); B <- file.path(base, "results", "bulk_bridge")
W <- "/Users/gabkrouk/Library/CloudStorage/Dropbox/Science/Papers/Emanuel New TARGET paper/Figure Will"
source(file.path(base, "scripts", "paper_style.R")); theme_set(theme_tpc(7))
cols <- c(EV = "#1E7AB6", ERF5 = "#EA8535", CRF4 = "#519C3E", `ERF5+CRF4` = "#C2382F")
cnt <- read.table(file.path(W, "counts.txt"), header = TRUE, check.names = FALSE)
cnt <- cnt[, !(colnames(cnt) %in% c("E2_EC.2", "E2_E.1"))]
cpm <- t(t(cnt) / colSums(cnt) * 1e6); L <- log2(cpm[rowSums(cpm > 1) >= 3, ] + 1)
o <- fread(file.path(B, "bulk_condition_contrasts_per_gene.tsv")); setnames(o, 1, "gene")
cls_lab <- c(`CRF4 only` = "CRF4 alone", `shared single + combination` = "CRF4 and combination", `combination-specific` = "combination-specific")
o <- o[class %in% names(cls_lab)]; o[, cl := factor(cls_lab[class], levels = cls_lab)]
o[, dir := ifelse(class == "combination-specific", sign(log2FC_ERF5_CRF4), ifelse(class == "CRF4 only", sign(log2FC_CRF4), sign(log2FC_ERF5_CRF4 + log2FC_CRF4)))]
smp <- colnames(L); con <- sub(".*_(\\w+)\\..*", "\\1", smp); ex <- ifelse(startsWith(smp, "E1"), "bulk exp 1", "bulk exp 2")
lab <- c(N = "EV", E = "ERF5", C = "CRF4", EC = "ERF5+CRF4")
M <- L[o$gene, ]; for (e in unique(ex)) M[, ex == e] <- M[, ex == e] - rowMeans(M[, ex == e & con == "N", drop = FALSE])
Z <- t(scale(t(M)))
long <- data.table(gene = rep(rownames(Z), ncol(Z)), sample = rep(colnames(Z), each = nrow(Z)), z = as.vector(Z))
long[, cond := factor(lab[sub(".*_(\\w+)\\..*", "\\1", sample)], levels = lab)][, exp := ifelse(startsWith(sample, "E1"), "exp 1", "exp 2")]
long <- merge(long, o[, .(gene, cl, dir)], by = "gene")
ord <- o[order(cl, -dir, log2FC_ERF5_CRF4)]$gene
long[, gene := factor(gene, levels = rev(ord))]
long[, sample := factor(sample, levels = colnames(Z)[order(factor(lab[con], levels = lab), ex)])]
long[, grp := paste0(cl, ifelse(dir > 0, " (up)", " (down)"))]
n_lab <- o[, .N, by = .(cl, dir)][order(cl, -dir)]
pA <- ggplot(long, aes(sample, gene, fill = pmax(pmin(z, 2.5), -2.5))) + geom_raster() +
  facet_grid(cl ~ cond, scales = "free", space = "free") +
  scale_fill_gradient2(low = "#1E7AB6", mid = "white", high = "#C2382F", name = "z-score") +
  labs(x = NULL, y = NULL, title = "A",
       subtitle = sprintf("Genes significant vs EV at FDR 5%% (moderated t, experiment as covariate). CRF4 alone: %d; CRF4 and combination: %d;\ncombination-specific (combination FDR 5%%, neither TF alone p < 0.05): %d. ERF5 alone: 2 genes.",
                          sum(o$cl == "CRF4 alone"), sum(o$cl == "CRF4 and combination"), sum(o$cl == "combination-specific"))) +
  theme(axis.text = element_blank(), panel.grid = element_blank(), strip.text.y = element_text(angle = 0, hjust = 0), panel.spacing = unit(1, "pt"))
fc <- melt(o[, .(gene, cl, ERF5 = log2FC_ERF5, CRF4 = log2FC_CRF4, `ERF5+CRF4` = log2FC_ERF5_CRF4)], id.vars = c("gene", "cl"))
fc[, variable := factor(variable, levels = c("ERF5", "CRF4", "ERF5+CRF4"))]
pB <- ggplot(fc, aes(variable, value, fill = variable)) + geom_hline(yintercept = 0, colour = "grey60") +
  geom_violin(scale = "width", linewidth = 0.2, alpha = 0.8) + geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.25) +
  facet_wrap(~cl, nrow = 1) + scale_fill_manual(values = cols[-1], guide = "none") + coord_cartesian(ylim = c(-2.5, 2.5)) +
  labs(x = NULL, y = "log2 fold change vs EV", title = "B") + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ho <- fromJSON(file.path(B, "bulk_combination_specific_heldout.json"))
hd <- rbindlist(lapply(names(ho), function(k) data.table(set = k, frac = ho[[k]]$same_sign_frac, n = ho[[k]]$n, p = ho[[k]]$binom_p)))
hd[, cls := ifelse(grepl("^comb", set), "combination-specific", "CRF4 alone")][, dir := sub(".* ", "", set)]
hd[, dir := sub("E1->E2", "found in exp 1,\ntested in exp 2", sub("E2->E1", "found in exp 2,\ntested in exp 1", dir))]
pC <- ggplot(hd, aes(dir, frac * 100, fill = cls)) + geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_hline(yintercept = 50, linetype = 2) + geom_text(aes(label = sprintf("n=%d\np=%.0e", n, p)), position = position_dodge(0.8), vjust = -0.2, size = 6 * PT) +
  scale_fill_manual(values = c(`combination-specific` = "#C2382F", `CRF4 alone` = "#519C3E"), name = NULL) + coord_cartesian(ylim = c(0, 115)) +
  labs(x = NULL, y = "% same direction in the other experiment", title = "C",
       subtitle = "Genes found in one experiment (p < 0.01)\ntested in the other; dashed line = chance") + theme(legend.position = "bottom")
go <- fread(file.path(base, "results", "go", "go_enrichment_all.tsv"))[grepl("^bulk", set) & q_BH < 0.05]
go <- go[order(q_BH)][, .SD[!duplicated(overlap)][1:min(.N, 4)], by = set]
go[, set := sub("^bulk ", "", set)][, set := factor(set, levels = c("CRF4 alone up", "CRF4 alone down", "CRF4 and combination up", "CRF4 and combination down", "combination-specific up", "combination-specific down"))]
go[, term := factor(term, levels = rev(unique(term[order(set, q_BH)])))]
pD <- ggplot(go, aes(set, term, size = overlap, colour = -log10(q_BH))) + geom_point() +
  scale_colour_gradient(low = "#F4B6A6", high = "#C2382F", name = "-log10 q") + scale_size(range = c(0.8, 3.2), name = "genes") +
  labs(x = NULL, y = NULL, title = "D",
       subtitle = "Universe: genes expressed in the bulk experiments; top non-redundant terms per class") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.text.y = element_text(size = 6.5))
g <- arrangeGrob(pA, arrangeGrob(pB, pC, ncol = 2, widths = c(1.5, 1)), pD, ncol = 1, heights = c(1.5, 1, 0.95))
save_tpc(g, FG, "Fig6_ERF5_CRF4_bulk_FDR", 6.5, 8)
fwrite(o[, .(gene, class = cl, log2FC_ERF5, q_ERF5, log2FC_CRF4, q_CRF4, log2FC_ERF5_CRF4, q_ERF5_CRF4, p_EC_vs_C, p_EC_vs_E)], file.path(B, "Fig6_gene_classes.tsv"), sep = "\t")
cat("ok\n")
