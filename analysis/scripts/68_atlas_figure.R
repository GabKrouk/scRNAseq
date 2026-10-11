#!/usr/bin/env Rscript
# Step 68: figure - the network connections of the paper (Fig 5B) tested with endogenous TFs in the
# Shahan et al. (2022) root atlas (steps 65-67). gab's style: letters only, Arial, 6.5 in wide.
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(gridExtra); library(jsonlite)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
R <- function(...) file.path(base, "results", ...)
FG <- file.path(base, "figures", "atlas"); dir.create(FG, showWarnings = FALSE, recursive = TRUE)
source(file.path(base, "scripts", "paper_style.R"))
theme_set(theme_tpc(7))
gcol <- c(emergent = "#7B3294", synergy = "#C2382F", competition = "#1E7AB6", reversal = "#EA8535")
scol <- c(recovered = "#519C3E", opposite = "grey25", `not significant` = "grey82")
cols <- c(neither = "#1E7AB6", A = "#EA8535", B = "#519C3E", both = "#C2382F")
T <- fread(R("atlas", "network_connections_atlas.tsv"))
T[, status := factor(status, levels = names(scol))][, gate_paper := factor(gate_paper, levels = names(gcol))]

# ---------------------------------------------------------------- A: endogenous detection of the 23 TFs
d <- fread(R("atlas", "tf_detection_by_celltype.tsv")); setnames(d, 1, "celltype")
ct_order <- c("Columella", "Lateral Root Cap", "Atrichoblast", "Trichoblast", "Cortex", "Endodermis", "Quiescent Center",
              "Xylem Pole Pericycle", "Phloem Pole Pericycle", "Procambium", "Protoxylem", "Metaxylem", "Protophloem",
              "Metaphloem & Companion Cell", "All cells")
ct_lab <- c("Columella", "Lateral root cap", "Atrichoblast", "Trichoblast", "Cortex", "Endodermis", "QC", "XPP", "PPP",
            "Procambium", "Protoxylem", "Metaxylem", "Protophloem", "MP & CC", "All cells")
dl <- melt(d, id.vars = "celltype", variable.name = "tf", value.name = "frac")
tf_ord <- dl[celltype == "All cells"][order(frac)]$tf
dl[, tf := factor(tf, levels = tf_ord)][, celltype := factor(celltype, levels = ct_order, labels = ct_lab)]
pA <- ggplot(dl, aes(celltype, tf, fill = 100 * frac)) + geom_tile(colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 14.5, colour = "grey30", linewidth = 0.3) +
  scale_fill_gradientn(colours = c("white", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"), trans = "sqrt",
                       breaks = c(1, 10, 30, 60), name = "Cells with\nTF detected (%)") +
  labs(x = NULL, y = NULL, title = "A") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 5.5), axis.text.y = element_text(size = 5.5),
        panel.grid = element_blank(), legend.position = "right", legend.key.height = unit(14, "pt"))

# ---------------------------------------------------------------- B: network, gate -> target edges coloured by atlas result
N <- fread(R("paper", "network_nodes.tsv")); Ed <- fread(R("paper", "network_edges.tsv"))
Ed <- merge(Ed, N[, .(from = node, x, y, ftype = type)], by = "from"); Ed <- merge(Ed, N[, .(to = node, xend = x, yend = y)], by = "to")
Ed[ftype != "TF", `:=`(pair = sub("\\|.*", "", sub("^G:", "", from)), target = sub("^T:", "", to))]
Ed <- merge(Ed, T[, .(pair, target, status)], by = c("pair", "target"), all.x = TRUE)
Ed[ftype != "TF" & is.na(status), status := "not significant"]
tfn <- N[type == "TF"]; gt <- N[!(type %in% c("TF", "target"))]; tg <- N[type == "target"]
gt[, type := factor(type, levels = names(gcol))]
eg <- Ed[ftype != "TF"][order(match(status, c("not significant", "opposite", "recovered")))]
pB <- ggplot() +
  geom_segment(data = Ed[ftype == "TF"], aes(x, y, xend = xend, yend = yend), colour = "grey80", linewidth = 0.1, alpha = 0.6) +
  geom_segment(data = eg, aes(x, y, xend = xend, yend = yend, colour = status, linewidth = status)) +
  geom_point(data = tg, aes(x, y), size = 0.15, colour = "grey55") +
  geom_point(data = gt, aes(x, y, fill = type), size = 1.1, shape = 23, colour = "white", stroke = 0.15) +
  geom_point(data = tfn, aes(x, y, size = degree * 0.5), shape = 21, fill = "#F2C14E", colour = "grey20", stroke = 0.3, show.legend = FALSE) +
  geom_text(data = tfn[, .(lx = x * 1.1, ly = y * 1.1, label, h = ifelse(x > 0.05, 0, ifelse(x < -0.05, 1, 0.5)))],
            aes(lx, ly, label = label, hjust = h), size = 5.5 * PT) +
  scale_colour_manual(values = scol, name = "Gate in the atlas", drop = FALSE) +
  scale_linewidth_manual(values = c(recovered = 0.45, opposite = 0.3, `not significant` = 0.12), guide = "none") +
  scale_fill_manual(values = gcol, name = "Gate (paper)") + scale_size_continuous(range = c(0.6, 2.6)) +
  guides(colour = guide_legend(override.aes = list(linewidth = 0.8), order = 1), fill = guide_legend(override.aes = list(size = 2.5), order = 2)) +
  coord_equal(clip = "off") + scale_x_continuous(expand = expansion(mult = 0.18)) + scale_y_continuous(expand = expansion(mult = 0.08)) +
  labs(title = "B") + theme_void(base_size = 7) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0), plot.title.position = "plot",
        legend.position = "right", legend.text = element_text(size = 5.5), legend.title = element_text(size = 6),
        legend.key.size = unit(7, "pt"), legend.spacing.y = unit(2, "pt"), plot.margin = margin(3, 3, 3, 3))

# ---------------------------------------------------------------- C: effects of TFa alone, TFb alone and both, paper vs atlas
eff <- rbind(T[, .(panel = "TFa alone", paper = paper_dA, atlas = dA)], T[, .(panel = "TFb alone", paper = paper_dB, atlas = dB)],
             T[, .(panel = "Both TFs", paper = paper_dAB, atlas = dAB)])
eff[, panel := factor(panel, levels = c("TFa alone", "TFb alone", "Both TFs"))]
rh <- eff[, .(rho = cor(paper, atlas, method = "spearman"), same = mean(sign(paper) == sign(atlas))), by = panel]
pcol <- c(`TFa alone` = "#EA8535", `TFb alone` = "#519C3E", `Both TFs` = "#C2382F")
pC <- ggplot(eff, aes(paper, atlas, colour = panel)) + geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.25) + geom_point(size = 0.45, alpha = 0.6, stroke = 0) +
  geom_text(data = rh, aes(x = -Inf, y = Inf, label = sprintf("rho = %.2f\nsame sign %d%%", rho, round(100 * same))),
            hjust = -0.08, vjust = 1.15, size = 5.5 * PT, colour = "grey15", lineheight = 0.9) +
  facet_wrap(~panel, nrow = 1, scales = "free") + scale_colour_manual(values = pcol, guide = "none") +
  labs(x = "Effect in scTARGET-seq (vs neither TF, Exp1 + Exp2)", y = "Effect in the atlas\n(endogenous TFs)", title = "C")

# ---------------------------------------------------------------- D: the gate (non-additivity), paper vs atlas
pD <- ggplot(T[order(match(status, c("not significant", "opposite", "recovered")))], aes(paper_nonadd, nonadd, colour = status)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.25) + geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.25) +
  geom_point(size = 0.6, stroke = 0) + scale_colour_manual(values = scol, guide = "none") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.2, size = 5.5 * PT, colour = "grey15",
           label = sprintf("rho = %.2f", cor(T$paper_nonadd, T$nonadd, method = "spearman"))) +
  labs(x = "Gate in scTARGET-seq\n(deviation from additivity)", y = "Gate in the atlas", title = "D")

# ---------------------------------------------------------------- E: outcome per gate type
cnt <- T[, .N, by = .(gate_paper, status)]
pE <- ggplot(cnt, aes(N, gate_paper, fill = status)) + geom_col(width = 0.7, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = scol, name = NULL) + scale_y_discrete(limits = rev(names(gcol))) +
  labs(x = "Connections", y = NULL, title = "E") + theme(legend.position = "bottom", legend.direction = "vertical", panel.grid.major.y = element_blank())

# ---------------------------------------------------------------- F: atlas gates follow the atlas single-TF effects
pF <- ggplot(T[order(match(status, c("not significant", "opposite", "recovered")))], aes(dA + dB, nonadd, colour = status)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.25) + geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.25) +
  geom_point(size = 0.6, stroke = 0) + scale_colour_manual(values = scol, guide = "none") +
  annotate("text", x = Inf, y = Inf, hjust = 1.08, vjust = 1.2, size = 5.5 * PT, colour = "grey15",
           label = sprintf("rho = %.2f", cor(T$dA + T$dB, T$nonadd, method = "spearman"))) +
  labs(x = "Sum of single-TF effects\nin the atlas", y = "Gate in the atlas", title = "F")

# ---------------------------------------------------------------- G: examples, paper vs atlas
pick <- data.table(pair = c("COL5-HB6", "GATA17-TGA4", "CRF4-HB6", "HB6-NAC4"),
                   target = c("AT3G53110", "AT1G69680", "AT1G49140", "AT1G04270"))
pick[, `:=`(a = sub("-.*", "", pair), b = sub(".*-", "", pair))]
pick <- merge(pick, T[, .(pair, target, gate_paper)], by = c("pair", "target"), sort = FALSE)
pick[, lab := sprintf("%s: %s x %s -> %s", tools::toTitleCase(as.character(gate_paper)), a, b, target)]
P5 <- fread(R("paper", "fig5d_candidates.tsv"))
sc <- rbindlist(lapply(c("Exp1", "Exp2"), function(e) P5[, .(pair, target, exp = e, A = get(paste0("dA_", e)), B = get(paste0("dB_", e)), both = get(paste0("dAB_", e)))]))
sc <- melt(sc, id.vars = c("pair", "target", "exp"), variable.name = "state", value.name = "est")
sc <- rbind(sc, unique(sc[, .(pair, target, exp, state = "neither", est = 0)]))
sc[, se := NA_real_][, data := "scTARGET-seq"]
at <- fread(R("atlas", "network_connections_states.tsv"))[, `:=`(exp = "Atlas", data = "Atlas (endogenous)")]
ex <- merge(rbind(sc, at, fill = TRUE), pick[, .(pair, target, lab, a, b)], by = c("pair", "target"))
ex[, state := factor(state, levels = names(cols))][, lab := factor(lab, levels = pick$lab)]
ex[, data := factor(data, levels = c("scTARGET-seq", "Atlas (endogenous)"))]
addv <- ex[, .(add = mean(est[state == "A"]) + mean(est[state == "B"])), by = .(lab, data)]
pG <- ggplot(ex, aes(state, est, colour = state)) + geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.25) +
  geom_errorbar(data = addv, aes(x = "both", ymin = add, ymax = add), inherit.aes = FALSE, width = 0.8, linetype = "22", linewidth = 0.4) +
  geom_point(data = ex[data == "scTARGET-seq"], aes(shape = exp), size = 1.1, stroke = 0.45, position = position_dodge(width = 0.45)) +
  geom_pointrange(data = ex[data != "scTARGET-seq"], aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se), size = 0.18, linewidth = 0.4) +
  facet_grid(data ~ lab, scales = "free_y", labeller = labeller(lab = label_wrap_gen(22))) +
  scale_colour_manual(values = cols, guide = "none") + scale_shape_manual(values = c(Exp1 = 16, Exp2 = 4), name = NULL) +
  scale_x_discrete(labels = c(neither = "none", A = "TFa", B = "TFb", both = "both")) +
  labs(x = NULL, y = "Change vs cells with neither TF\n(adjusted log expression)", title = "G") +
  theme(legend.position = "bottom", legend.margin = margin(-4, 0, 0, 0), strip.text.x = element_text(size = 5.5, face = "plain"),
        strip.text.y = element_text(size = 5.5), axis.text.x = element_text(size = 5.5))

g <- arrangeGrob(arrangeGrob(pA, pB, ncol = 2, widths = c(1, 1.15)), pC,
                 arrangeGrob(pD, pE, pF, ncol = 3, widths = c(1, 1, 1)), pG,
                 ncol = 1, heights = c(2.35, 1.45, 1.75, 2.4))
save_tpc(g, FG, "Atlas_network_connections", 6.5, 7.95)
cat("atlas figure written\n")
