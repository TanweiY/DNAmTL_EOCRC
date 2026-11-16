source("functions.R", local = TRUE)
require_packages(c("readxl", "ggplot2", "ggpubr", "officer", "rvg"))
# ================================================================
# DNAmTL: Time-dependent AUC - eight-panel figure (A4 portrait)
# Rows 1-2: Overall survival; rows 3-4: Progression-free survival
# Each row: EOCRC (left), LOCRC (right)
# Single common model legend; panels A-H
# Exports a vector PDF and an editable PowerPoint slide
# ================================================================

library(readxl)
library(ggplot2)
library(ggpubr)
library(grid)
library(officer)
library(rvg)

# Run 13_time_dependent_AUC.R first to create both cutoff-specific outputs.
old_dir <- dirname(result_file("auc_cutoff55", "placeholder"))
new_dir <- dirname(result_file("auc_cutoff50", "placeholder"))
out_dir <- dirname(result_file("placeholder"))
a4_template <- file.path("data", "A4_portrait_blank_template.pptx")

# --- 2. Common settings ---
model_colors <- c(ref = "steelblue3", mtl = "#e5383b")
model_shapes <- c(ref = 1, mtl = 2)
model_labels <- c(ref = "Age + Sex + stage",
                  mtl = "Age + Sex + stage + DNAmTL")

# Figure-specific plot maker: the legend is off in all eight panels.
make_auc <- function(path, title, left_col = TRUE, bottom_row = FALSE) {
  dat <- as.data.frame(read_excel(path))
  needed <- c("times", "AUC", "model")
  if (!all(needed %in% names(dat))) {
    stop("Missing columns in ", path, ": ",
         paste(setdiff(needed, names(dat)), collapse = ", "))
  }
  dat$times <- as.numeric(dat$times)
  dat$AUC <- as.numeric(dat$AUC)
  dat <- dat[dat$times >= 1 & dat$times <= 10, , drop = FALSE]
  dat$model <- factor(as.character(dat$model), levels = c("ref", "mtl"))
  if (any(is.na(dat$model))) {
    stop("Unexpected model values in ", path,
         ". Expected: ref and mtl.")
  }

  ggplot(dat, aes(x = times, y = AUC, colour = model, shape = model)) +
    geom_line(linewidth = 0.45) +
    geom_point(size = 1.7, stroke = 0.65, fill = "white") +
    scale_colour_manual(values = model_colors, breaks = c("ref", "mtl"),
                        labels = model_labels, name = "Model") +
    scale_shape_manual(values = model_shapes, breaks = c("ref", "mtl"),
                       labels = model_labels, name = "Model") +
    scale_x_continuous(breaks = 1:10) +
    scale_y_continuous(breaks = seq(0.60, 1.00, 0.05)) +
    coord_cartesian(xlim = c(1, 10), ylim = c(0.60, 1.00)) +
    labs(title = title,
         x = if (bottom_row) "Follow-up time (years)" else NULL,
         y = if (left_col) "Time-dependent AUC" else NULL) +
    guides(colour = guide_legend(override.aes = list(
      shape = unname(model_shapes), linewidth = 0.6)),
      shape = "none") +
    theme_minimal(base_family = "sans", base_size = 10) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 11, hjust = 0),
      axis.title = element_text(size = 10.5),
      axis.text = element_text(size = 9),
      axis.title.y = element_text(margin = margin(r = 4)),
      axis.title.x = element_text(margin = margin(t = 4)),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.45),
      panel.grid.minor = element_line(colour = "grey93", linewidth = 0.2),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.25),
      plot.margin = margin(5, 6, 5, 5)
    )
}

# --- 3. Eight plots in the requested order ---
p1 <- make_auc(file.path(old_dir, "AUCeocrc_os.xlsx"), "EOCRC (<55)", TRUE, FALSE)
p2 <- make_auc(file.path(old_dir, "AUClocrc_os.xlsx"), "LOCRC (>=55)", FALSE, FALSE)
p3 <- make_auc(file.path(new_dir, "AUCeocrc_os.xlsx"), "EOCRC (<50)", TRUE, FALSE)
p4 <- make_auc(file.path(new_dir, "AUClocrc_os.xlsx"), "LOCRC (>=50)", FALSE, FALSE)
p5 <- make_auc(file.path(old_dir, "AUCeocrc_pfs.xlsx"), "EOCRC (<55)", TRUE, FALSE)
p6 <- make_auc(file.path(old_dir, "AUClocrc_pfs.xlsx"), "LOCRC (>=55)", FALSE, FALSE)
p7 <- make_auc(file.path(new_dir, "AUCeocrc_pfs.xlsx"), "EOCRC (<50)", TRUE, TRUE)
p8 <- make_auc(file.path(new_dir, "AUClocrc_pfs.xlsx"), "LOCRC (>=50)", FALSE, TRUE)

# --- 4. Extract one shared legend from a temporary version of p1 ---
legend_source <- p1 + theme(
  legend.position = "top",
  legend.text = element_text(size = 10, family = "sans"),
  legend.title = element_text(size = 10, family = "sans"),
  legend.key.width = unit(0.8, "cm"),
  legend.background = element_blank()
)
shared_legend <- ggpubr::get_legend(legend_source)

# --- 5. Combine A-D for OS, E-H for PFS ---
# ggarrange() adds panel letters in the plot layout, separate from plot titles.
os_grid <- ggpubr::ggarrange(
  p1, p2, p3, p4,
  ncol = 2, nrow = 2,
  labels = c("A", "B", "C", "D"),
  font.label = list(size = 12, face = "bold", family = "sans"),
  label.x = 0.005, label.y = 1.00,
  common.legend = FALSE
)
pfs_grid <- ggpubr::ggarrange(
  p5, p6, p7, p8,
  ncol = 2, nrow = 2,
  labels = c("E", "F", "G", "H"),
  font.label = list(size = 12, face = "bold", family = "sans"),
  label.x = 0.005, label.y = 1.00,
  common.legend = FALSE
)

# --- 6. Section headers and one legend ---
os_header <- ggpubr::text_grob("Overall survival",
                               face = "bold", size = 13, family = "sans")
pfs_header <- ggpubr::text_grob("Progression-free survival",
                                face = "bold", size = 13, family = "sans")

sectioned <- ggpubr::ggarrange(
  os_header, os_grid, pfs_header, pfs_grid,
  ncol = 1, nrow = 4,
  heights = c(0.07, 1, 0.07, 1),
  align = "v"
)
final_fig <- ggpubr::ggarrange(
  shared_legend, sectioned,
  ncol = 1, nrow = 2,
  heights = c(0.08, 1),
  align = "v"
)

print(final_fig)

# --- 7. Export PDF at true A4 portrait size (210 x 297 mm) ---
pdf(file.path(out_dir, "DNAmTL_AUC_8panel_A4.pdf"),
    width = 210 / 25.4, height = 297 / 25.4,
    family = "Helvetica", onefile = TRUE)
print(final_fig)
dev.off()

# --- 8. Editable PowerPoint ---
# An optional blank A4 template preserves true A4 slide size. Without it, use
# officer's standard slide size and fit the A4 figure without distortion.
ppt <- if (file.exists(a4_template)) officer::read_pptx(a4_template) else officer::read_pptx()
if (file.exists(a4_template)) {
  size <- officer::slide_size(ppt)
  stopifnot(abs(size$width - 210 / 25.4) < 0.02, abs(size$height - 297 / 25.4) < 0.02)
}
while (length(ppt) > 0) ppt <- officer::remove_slide(ppt, index = 1)
layout <- subset(officer::layout_summary(ppt), layout == "Blank")
if (!nrow(layout)) stop("PowerPoint template must contain a Blank layout.")
ppt <- officer::add_slide(ppt, layout = "Blank", master = layout$master[1])
size <- officer::slide_size(ppt)
margin <- 0.12
height <- min(size$height - 2 * margin, (size$width - 2 * margin) * 297 / 210)
width <- height * 210 / 297
ppt <- officer::ph_with(ppt,
  value = rvg::dml(code = print(final_fig), editable = TRUE, fonts = list(sans = "Arial")),
  location = officer::ph_location(left = (size$width - width) / 2,
                                  top = (size$height - height) / 2, width = width, height = height))
print(ppt, target = file.path(out_dir, "DNAmTL_AUC_8panel_editable.pptx"))
save_session_info("14_AUC_plots")
