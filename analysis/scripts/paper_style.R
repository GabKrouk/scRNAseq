# Shared style for the revised paper figures (The Plant Cell): Arial, 7-8 pt text (>= 6 pt), bold 12 pt panel
# letters with no panel titles, two-column width 6.5 in, height < 8 in, PDF (vector) + TIFF/PNG at 300 dpi, RGB.
suppressPackageStartupMessages({library(ggplot2); library(grid)})
TPC_W <- 6.5
PT <- 1 / 2.845276            # ggplot text size (mm) per point
theme_tpc <- function(base = 7) {
  theme_minimal(base_size = base) + theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0, margin = margin(0, 0, 2, 0)), plot.title.position = "plot",
    plot.subtitle = element_blank(), plot.caption = element_blank(),
    panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA, colour = "grey80", linewidth = 0.3),
    axis.text = element_text(size = base - 1, colour = "grey15"), axis.title = element_text(size = base),
    strip.text = element_text(face = "bold", size = base - 0.5, margin = margin(1, 0, 1, 0)),
    legend.text = element_text(size = base - 1), legend.title = element_text(size = base - 0.5),
    legend.key.size = unit(8, "pt"), legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"),
    plot.margin = margin(3, 4, 3, 3))
}
save_tpc <- function(g, dir, name, w = TPC_W, h, png = TRUE) {
  stopifnot(w <= 6.5 + 1e-9, h <= 8 + 1e-9)
  pdf(file.path(dir, paste0(name, ".pdf")), width = w, height = h, family = "ArialMT", useDingbats = FALSE); grid.newpage(); grid.draw(g); invisible(dev.off())
  if (png) { png(file.path(dir, paste0(name, ".png")), width = w, height = h, units = "in", res = 300, family = "Arial", bg = "white"); grid.newpage(); grid.draw(g); invisible(dev.off()) }
}
