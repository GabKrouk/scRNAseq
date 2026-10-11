#!/usr/bin/env Rscript
# Step 56: Supplemental figure - cell identity: composition, GR-TF detection per cell type, robustness of effects.
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(gridExtra); library(jsonlite)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
R <- function(...) file.path(base, "results", ...); FG <- file.path(base, "figures", "paper")
source(file.path(base, "scripts", "paper_style.R")); theme_set(theme_tpc(7))
tb <- rbindlist(lapply(c("Exp1", "Exp2"), function(e) { d <- fread(R(paste0(e, "_TF_by_celltype.tsv"))); setnames(d, 1, "row"); d[, exp := e] }), fill = TRUE)
m <- melt(tb, id.vars = c("row", "exp"), variable.name = "cell_type")
m[, cell_type := gsub("_", " ", cell_type)]
nc <- m[row == "n_cells"]; nc[, frac := value / sum(value), by = exp]
nc[, cell_type := factor(cell_type, levels = nc[exp == "Exp2"][order(value)]$cell_type)]
pA <- ggplot(nc, aes(frac * 100, cell_type, fill = exp)) + geom_col(position = position_dodge(0.8), width = 0.75) +
  scale_fill_manual(values = c(Exp1 = "#0072B2", Exp2 = "#E69F00"), name = NULL) +
  labs(x = "% of QC-passing cells", y = NULL, title = "A",
       subtitle = "Protoplasts assigned when one marker score dominates;\nmost cells remain unassigned, so scores enter the models as covariates") + theme(legend.position = c(0.75, 0.2))
tf <- m[!(row %in% c("n_cells", "mean_n_TFs")) & exp == "Exp2"]
tf <- merge(tf, nc[exp == "Exp2", .(cell_type, n = value)], by = "cell_type")
tf <- tf[n >= 30]; tf[, cell_type := factor(cell_type, levels = levels(nc$cell_type))]
pB <- ggplot(tf, aes(cell_type, row, fill = value)) + geom_tile(colour = "white") + geom_text(aes(label = round(value)), size = 6 * PT) +
  scale_fill_gradient(low = "white", high = "#C2382F", name = "% cells\nwith TF") +
  labs(x = NULL, y = NULL, title = "B", subtitle = "Cell types with >= 30 cells") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), panel.grid = element_blank(), axis.text.y = element_text(size = 6.5))
cr <- fromJSON(R("paper", "celltype_robustness.json"))
d <- data.table(kind = rep(c("TF -> gene effects", "TF x TF interactions"), each = 3),
                what = rep(c("with cell-type\ncovariates", "without", "found\nboth ways"), 2),
                n = c(cr$main_with, cr$main_without, cr$main_overlap, cr$inter_with, cr$inter_without, cr$inter_overlap))
d[, what := factor(what, levels = c("with cell-type\ncovariates", "without", "found\nboth ways"))]
pC <- ggplot(d, aes(what, n, fill = what)) + geom_col(width = 0.7) + geom_text(aes(label = format(n, big.mark = ",")), vjust = -0.3, size = 6 * PT) +
  facet_wrap(~kind, scales = "free_y") + scale_fill_manual(values = c("#C2382F", "grey65", "#1E7AB6"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Reproduced calls", title = "C",
       subtitle = sprintf("Correlation of estimates with vs without cell-type covariates: r = %.2f (main), %.2f (interactions)", cr$main_beta_r, cr$inter_beta_r)) +
  theme(axis.text.x = element_text(size = 5.5))
g <- arrangeGrob(arrangeGrob(pA, pC, ncol = 1), pB, ncol = 2, widths = c(1, 1.05))
save_tpc(g, FG, "SupFig_celltype", 6.5, 5.2)
cat("ok\n")
