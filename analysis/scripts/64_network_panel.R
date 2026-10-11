#!/usr/bin/env Rscript
# Step 64: network panel (TFa, TFb -> gate -> target) of the 461 reproduced non-additive interactions; used in Fig 5.
suppressPackageStartupMessages({library(ggplot2); library(data.table)})
base <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
PT <- 1 / 2.845276
network_panel <- function(base, title = "B  Network of the 461 non-additive interactions") {
  N <- fread(file.path(base, "results", "paper", "network_nodes.tsv")); Ed <- fread(file.path(base, "results", "paper", "network_edges.tsv"))
  gcol <- c(emergent = "#7B3294", synergy = "#C2382F", competition = "#1E7AB6", reversal = "#EA8535")
  Ed <- merge(Ed, N[, .(from = node, x, y, tf = type)], by = "from"); Ed <- merge(Ed, N[, .(to = node, xend = x, yend = y, gt = type)], by = "to")
  Ed[, kind := ifelse(tf == "TF", "tf", gt)][tf != "TF", kind := N$type[match(from, N$node)]]
  tfn <- N[type == "TF"]; gt <- N[!(type %in% c("TF", "target"))]; tg <- N[type == "target"]
  gt[, type := factor(type, levels = names(gcol))]
  ggplot() +
    geom_segment(data = Ed[tf == "TF"], aes(x, y, xend = xend, yend = yend), colour = "grey70", linewidth = 0.1, alpha = 0.6) +
    geom_segment(data = Ed[tf != "TF"], aes(x, y, xend = xend, yend = yend, colour = kind), linewidth = 0.12, alpha = 0.55, show.legend = FALSE) +
    geom_point(data = tg, aes(x, y), size = 0.2, colour = "grey45") +
    geom_point(data = gt, aes(x, y, fill = type, size = n_targets), shape = 23, colour = "white", stroke = 0.2) +
    geom_point(data = tfn, aes(x, y, size = degree * 0.5), shape = 21, fill = "#F2C14E", colour = "grey20", stroke = 0.3, show.legend = FALSE) +
    geom_text(data = tfn[, .(lx = x * 1.09, ly = y * 1.09, label, h = ifelse(x > 0.05, 0, ifelse(x < -0.05, 1, 0.5)))], aes(lx, ly, label = label, hjust = h), size = 6 * PT) +
    geom_label(data = gt[order(-n_targets)][1:2], aes(x, y, label = sub("\\|.*", "", sub("-", " x ", label))), size = 6 * PT,
               label.padding = unit(0.06, "lines"), fill = alpha("white", 0.8), colour = "grey15", vjust = 1.6) +
    scale_colour_manual(values = gcol) + scale_fill_manual(values = gcol, name = "Gate (TFa x TFb)") +
    scale_size_continuous(range = c(0.4, 3.2), name = "Targets per gate", breaks = c(1, 5, 20)) +
    guides(fill = guide_legend(override.aes = list(size = 3)), size = guide_legend(override.aes = list(shape = 23, fill = "grey55", colour = "white"))) + coord_equal(clip = "off") + scale_x_continuous(expand = expansion(mult = 0.16)) +
    labs(title = title, subtitle = "TFa, TFb (yellow) -> gate (diamond: one per TF pair and interaction shape) -> target genes (grey).\nLabels: the two gates with the most targets") +
    theme_void(base_size = 7) + theme(plot.title = element_text(face = "bold", size = 12, hjust = 0), plot.title.position = "plot", plot.subtitle = element_blank(),
                                      legend.position = "bottom", legend.box = "vertical", legend.text = element_text(size = 6), legend.title = element_text(size = 6.5), legend.key.size = unit(7, "pt"), legend.margin = margin(0, 0, 0, 0), legend.spacing.y = unit(1, "pt"))
}
if (sys.nframe() == 0) {
  p <- network_panel(base)
  ggsave(file.path(base, "figures", "paper", "Fig5_network_test.png"), p, width = 6.5, height = 7, dpi = 250, bg = "white")
  cat("ok\n")
}
