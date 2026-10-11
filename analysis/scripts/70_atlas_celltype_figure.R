#!/usr/bin/env Rscript
# Step 70: figure - the paper's network connections in the Shahan 2022 atlas, cell type by cell type, and with
# less cell-type information in the adjustment (step 69). gab's style: letters only, Arial, 6.5 in wide.
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(gridExtra)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
R <- function(...) file.path(base, "results", ...)
FG <- file.path(base, "figures", "atlas"); dir.create(FG, showWarnings = FALSE, recursive = TRUE)
source(file.path(base, "scripts", "paper_style.R"))
theme_set(theme_tpc(7))
scol <- c(recovered = "#519C3E", opposite = "grey25", `not significant` = "grey85", `not testable` = "white")
cols <- c(neither = "#1E7AB6", A = "#EA8535", B = "#519C3E", both = "#C2382F")
ct_order <- c("Columella", "Lateral Root Cap", "Atrichoblast", "Trichoblast", "Cortex", "Endodermis", "Quiescent Center",
              "Xylem Pole Pericycle", "Phloem Pole Pericycle", "Procambium", "Protoxylem", "Metaxylem", "Protophloem",
              "Metaphloem & Companion Cell")
ct_lab <- c("Columella", "Lateral root cap", "Atrichoblast", "Trichoblast", "Cortex", "Endodermis", "QC", "XPP", "PPP",
            "Procambium", "Protoxylem", "Metaxylem", "Protophloem", "MP & CC")
var_lab <- c(full = "Sample x cell type x stage", no_celltype = "Sample x stage\n(no cell type)", sample_only = "Sample only\n(no cell type, no stage)")
S <- fread(R("atlas", "celltype_and_variant_summary.tsv"))
C <- S[variant %in% ct_order][, ct := factor(variant, levels = rev(ct_order), labels = rev(ct_lab))]
V <- S[variant %in% names(var_lab)][, v := factor(variant, levels = rev(names(var_lab)), labels = rev(var_lab))]

# ---------------------------------------------------------------- A: combined effect (both TFs vs neither), per cell type
A <- melt(C[, .(ct, observed = rho_both, chance = rho_both_chance, sig = p_rho_both < 0.05, testable)], id.vars = c("ct", "sig", "testable"))
pA <- ggplot(A, aes(value, ct)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.25) +
  geom_line(aes(group = ct), colour = "grey70", linewidth = 0.35) +
  geom_point(data = A[variable == "chance"], shape = 21, fill = "white", colour = "grey40", size = 1.3, stroke = 0.35) +
  geom_point(data = A[variable == "observed"], aes(fill = sig), shape = 21, colour = "grey15", size = 1.6, stroke = 0.3) +
  geom_text(data = A[variable == "observed"], aes(x = 0.72, label = testable), size = 5.5 * PT, colour = "grey35", hjust = 1) +
  scale_fill_manual(values = c(`TRUE` = "#C2382F", `FALSE` = "white"), guide = "none") +
  scale_x_continuous(limits = c(-0.12, 0.72), breaks = c(0, 0.2, 0.4, 0.6)) +
  labs(x = "Agreement with scTARGET-seq, both TFs vs neither\n(Spearman rho; open circle: chance)", y = NULL, title = "A")

# ---------------------------------------------------------------- B, C: gates recovered (right) and opposite (left)
gbar <- function(D, yv, letter, xlab = "Gates (connections, BH q < 0.05)") {
  L <- rbind(D[, .(y = get(yv), side = "recovered", n = recovered, ch = recovered_chance)],
             D[, .(y = get(yv), side = "opposite", n = -opposite, ch = -opposite_chance)])
  lim <- max(abs(c(L$n, L$ch))) * 1.08
  ggplot(L, aes(n, y, fill = side)) + geom_col(width = 0.68) +
    geom_point(aes(x = ch), shape = 124, size = 2.4, colour = "black") +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
    scale_fill_manual(values = scol, name = NULL, breaks = c("recovered", "opposite"),
                      labels = c("same direction as scTARGET-seq", "opposite direction")) +
    scale_x_continuous(labels = abs, limits = c(-lim, lim)) +
    labs(x = xlab, y = NULL, title = letter) + theme(legend.position = "bottom", legend.direction = "vertical", panel.grid.major.y = element_blank())
}
pB <- gbar(C, "ct", "B") + theme(legend.position = "none")
pC <- gbar(V, "v", "C")

# ---------------------------------------------------------------- D: connections recovered in >= 1 cell type, across cell types
T <- fread(R("atlas", "connections_by_variant_and_celltype.tsv.gz"))
T[, key := paste0(sub("-", " x ", pair), " -> ", target)]
keys <- T[is_celltype == TRUE & status == "recovered", .(k = uniqueN(variant)), by = key][order(-k, key)]
M <- CJ(key = keys$key, variant = c("full", ct_order))
M <- merge(M, T[variant %in% c("full", ct_order), .(key, variant, status)], by = c("key", "variant"), all.x = TRUE)
M[is.na(status), status := "not testable"]
M[, status := factor(status, levels = names(scol))]
M[, col := factor(variant, levels = c("full", ct_order), labels = c("All cells", ct_lab))]
M[, row := factor(key, levels = rev(keys$key))]
pD <- ggplot(M, aes(col, row, fill = status)) + geom_tile(colour = "grey92", linewidth = 0.2) +
  geom_vline(xintercept = 1.5, colour = "grey30", linewidth = 0.35) +
  scale_fill_manual(values = scol, name = NULL, drop = FALSE) +
  labs(x = NULL, y = NULL, title = "D") + scale_x_discrete(position = "top") +
  theme(axis.text.x.top = element_text(angle = 60, hjust = 0, vjust = 0, size = 6), axis.text.y = element_text(size = 6),
        panel.grid = element_blank(), panel.border = element_blank(), legend.position = "bottom", legend.direction = "horizontal",
        legend.key.size = unit(6, "pt"), legend.text = element_text(size = 6), plot.margin = margin(3, 16, 3, 3)) +
  guides(fill = guide_legend(nrow = 2))

# ---------------------------------------------------------------- E: endodermis examples, scTARGET vs atlas endodermis
pick <- data.table(pair = c("COL5-HB6", "LBD38-NAC4", "HAT22-TGA4", "HHO3-WRKY18"),
                   target = c("AT3G53110", "AT1G27730", "AT4G22590", "AT2G30600"))
P5 <- fread(R("paper", "fig5d_candidates.tsv"))[, gate := ifelse(shape == "cancellation", "competition", shape)]
pick <- merge(pick, P5[, .(pair, target, gate)], by = c("pair", "target"), sort = FALSE)
pick[, lab := sprintf("%s: %s x %s -> %s", tools::toTitleCase(gate), sub("-.*", "", pair), sub(".*-", "", pair), target)]
sc <- rbindlist(lapply(c("Exp1", "Exp2"), function(e) P5[, .(pair, target, exp = e, A = get(paste0("dA_", e)), B = get(paste0("dB_", e)), both = get(paste0("dAB_", e)))]))
sc <- melt(sc, id.vars = c("pair", "target", "exp"), variable.name = "state", value.name = "est")
sc <- rbind(sc, unique(sc[, .(pair, target, exp, state = "neither", est = 0)]))[, `:=`(se = NA_real_, data = "scTARGET-seq")]
en <- T[variant == "Endodermis"]
en <- rbind(en[, .(pair, target, state = "A", est = dA, se = se_dA)], en[, .(pair, target, state = "B", est = dB, se = se_dB)],
            en[, .(pair, target, state = "both", est = dAB, se = se_dAB)], unique(en[, .(pair, target, state = "neither", est = 0, se = 0)]))
en[, `:=`(exp = "Atlas", data = "Atlas, endodermis")]
ex <- merge(rbind(sc, en, fill = TRUE), pick[, .(pair, target, lab)], by = c("pair", "target"))
ex[, state := factor(state, levels = names(cols))][, lab := factor(lab, levels = pick$lab)][, data := factor(data, levels = c("scTARGET-seq", "Atlas, endodermis"))]
pick[, top := sprintf("%s: %s x %s\n-> %s (scTARGET-seq)", tools::toTitleCase(gate), sub("-.*", "", pair), sub(".*-", "", pair), target)]
pick[, bot := sprintf("%s x %s -> %s\n(atlas, endodermis)", sub("-.*", "", pair), sub(".*-", "", pair), target)]
ex <- merge(ex, pick[, .(pair, target, top, bot)], by = c("pair", "target"))
ex[, panel := factor(ifelse(data == "scTARGET-seq", top, bot), levels = c(pick$top, pick$bot))]
addv <- ex[, .(add = mean(est[state == "A"]) + mean(est[state == "B"])), by = .(panel)]
pE <- ggplot(ex, aes(state, est, colour = state)) + geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.25) +
  geom_errorbar(data = addv, aes(x = "both", ymin = add, ymax = add), inherit.aes = FALSE, width = 0.8, linetype = "22", linewidth = 0.4) +
  geom_point(data = ex[data == "scTARGET-seq"], aes(shape = exp), size = 1.1, stroke = 0.45, position = position_dodge(width = 0.45)) +
  geom_pointrange(data = ex[data != "scTARGET-seq"], aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se), size = 0.18, linewidth = 0.4) +
  facet_wrap(~panel, nrow = 2, scales = "free_y") +
  scale_colour_manual(values = cols, guide = "none") + scale_shape_manual(values = c(Exp1 = 16, Exp2 = 4), name = NULL) +
  scale_x_discrete(labels = c(neither = "none", A = "TFa", B = "TFb", both = "both")) +
  labs(x = NULL, y = "Change vs cells with neither TF\n(adjusted log expression)", title = "E") +
  theme(legend.position = "bottom", legend.margin = margin(-4, 0, 0, 0), strip.text = element_text(size = 6, face = "plain"),
        axis.text.x = element_text(size = 6))

left <- arrangeGrob(pA, pB, pC, ncol = 1, heights = c(2.05, 2.05, 1.6))
g <- arrangeGrob(arrangeGrob(left, pD, ncol = 2, widths = c(1, 1.05)), pE, ncol = 1, heights = c(5.45, 2.5))
save_tpc(g, FG, "Atlas_connections_by_celltype", 6.5, 7.95)
cat("cell-type figure written\n")
