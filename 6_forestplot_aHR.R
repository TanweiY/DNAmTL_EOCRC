source("functions.R", local = TRUE)
require_packages(c("readxl", "forestplot", "officer", "rvg"))
# DNAmTL prognosis forest plot: 55-year versus 50-year cutoffs
# Uses forestplot; portable sans-serif typography, larger point boxes, italic P,
# dashed HR=1 line, header rules and manual EOCRC/LOCRC colour legend.

library(readxl)
library(forestplot)
library(grid)

# --- 1. Read workbook -----------------------------------------------------
df <- as.data.frame(read_excel(
  if (!is.null(getOption("dnamtl.forest_file"))) getOption("dnamtl.forest_file") else input_file("forestdf.xlsx")
))

df$Group <- trimws(as.character(df$Group))
for (v in c("N", "HR", "LL", "UL", "P for interaction")) {
  df[[v]] <- suppressWarnings(as.numeric(df[[v]]))
}

# Heading rows contain no HR estimates.
outcomes <- c("Overall mortality", "Disease progression or death")
row_type <- ifelse(
  df$Group %in% outcomes, "outcome",
  ifelse(grepl("^(50|55)[-–]year cut-off$", df$Group), "cutoff", "estimate")
)
is_estimate <- row_type == "estimate"

# Check confidence intervals, without changing the source data.
stopifnot(all(is.finite(df$HR[is_estimate])),
          all(is.finite(df$LL[is_estimate])),
          all(is.finite(df$UL[is_estimate])),
          all(df$LL[is_estimate] <= df$HR[is_estimate] &
                df$HR[is_estimate] <= df$UL[is_estimate]))

# --- 2. Text columns ------------------------------------------------------
group_label <- ifelse(is_estimate, paste0("   ", df$Group), df$Group)
n_label <- ifelse(is.na(df$N), "", format(df$N, big.mark=",", trim=TRUE))
hr_label <- ifelse(
  is_estimate,
  sprintf("%.2f (%.2f–%.2f)", df$HR, df$LL, df$UL),
  ""
)
p_value <- df[["P for interaction"]]
p_label <- ifelse(
  is.na(p_value), "",
  ifelse(p_value < 0.0001, "<0.0001", sprintf("%.2f", p_value))
)

# A list of columns allows mathematical notation in only the P header.
# Unlike a character matrix, this preserves italic(P).
fp_table <- list(
  as.list(c("Group", group_label)),
  as.list(c("N", n_label)),
  as.list(c("aHR (95% CI)", hr_label)),
  c(list(expression(italic(P) ~ "value for interaction")), as.list(p_label))
)

fp_mean <- c(NA_real_, df$HR)
fp_lower <- c(NA_real_, df$LL)
fp_upper <- c(NA_real_, df$UL)
summary_rows <- c(TRUE, row_type != "estimate")

# --- 3. Colours and forest plot styling ----------------------------------
eocrc_col <- "#E41A1C"
locrc_col <- "#377EB8"
row_col <- c(
  "#444444", # header placeholder
  ifelse(df$Group == "EOCRC", eocrc_col,
         ifelse(df$Group == "LOCRC", locrc_col, "#444444"))
)

shape_gp <- fpShapesGp(
  box = lapply(row_col, function(z) gpar(col = z, fill = z)),
  lines = lapply(row_col, function(z) gpar(col = z, lwd = 1.6)),
  vertices = lapply(row_col, function(z) gpar(col = z, lwd = 1.5)),
  # Style the vertical grid line at HR=1 (rather than the default zero line).
  grid = gpar(col = "black", lwd = 1, lty = "dashed")
)

# Black, 1-point horizontal rule immediately under the column headings,
# plus one separator between the two outcome sections.
rule_rows <- unique(c(2L, which(row_type == "outcome") + 1L))
header_rules <- setNames(
  lapply(rule_rows, function(i) gpar(col = "black", lwd = 1, lty = 1)),
  as.character(rule_rows)
)

# A grid line at HR=1 can span the plotted graph region, unlike the
# usual zero line, which may stop above the first estimate.
ref_grid <- 1
attr(ref_grid, "gp") <- gpar(col = "black", lwd = 1, lty = "dashed")

p <- forestplot(
  labeltext = fp_table,
  mean = fp_mean,
  lower = fp_lower,
  upper = fp_upper,
  is.summary = summary_rows,
  graph.pos = 3,
  align = c("l", "r", "r", "c"),
  zero = NA_real_,   # prevent double-drawing HR=1
  grid = ref_grid,   # black dashed line at HR=1
  xlog = TRUE,
  clip = c(0.35, 1.25),
  xticks = c(0.4, 0.6, 0.8, 1.0, 1.2),
  xticks.digits = 1,
  graphwidth = unit(69, "mm"),
  lineheight = unit(9.0, "mm"),
  colgap = unit(5, "mm"),
  boxsize = 0.38, # larger squares; try 0.42 if more emphasis is desired
  ci.vertices = TRUE,
  ci.vertices.height = 0.15,
  hrzl_lines = header_rules,
  xlab = "Adjusted HR per 1 SD increase in DNAmTL",
  shapes_gp = shape_gp,
  col = fpColors(
    box = c(eocrc_col, locrc_col),
    lines = c(eocrc_col, locrc_col),
    zero = "black",
    axes = "black"
  ),
  txt_gp = fpTxtGp(
    label = gpar(fontsize = 12.5, fontfamily = "sans"),
    summary = gpar(fontsize = 12.5, fontface = "bold", fontfamily = "sans"),
    ticks = gpar(fontsize = 12.5, fontfamily = "sans"),
    xlab = gpar(fontsize = 14, fontfamily = "sans"),
    legend = gpar(fontsize = 12.5, fontfamily = "sans")
  ),
  mar = unit(c(5, 5, 4, 5), "mm"),
  new_page = FALSE
)

# --- 4. Draw a legend separately: EOCRC and LOCRC are different rows, ----
# not multiple CI series in the forestplot package.
draw_final_plot <- function(p) {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = 2, ncol = 1,
    heights = unit.c(unit(13, "mm"), unit(1, "null"))
  )))

  pushViewport(viewport(layout.pos.row = 1))
  # Small coloured squares, matching forest-plot markers.
  grid.rect(x = unit(0.40, "npc"), y = unit(0.5, "npc"),
            width = unit(4.0, "mm"), height = unit(4.0, "mm"),
            gp = gpar(col = eocrc_col, fill = eocrc_col))
  grid.text("EOCRC", x = unit(0.40, "npc") + unit(4, "mm"),
            y = unit(0.5, "npc"), just = "left",
            gp = gpar(fontfamily = "sans", fontsize = 12.5))
  grid.rect(x = unit(0.57, "npc"), y = unit(0.5, "npc"),
            width = unit(4.0, "mm"), height = unit(4.0, "mm"),
            gp = gpar(col = locrc_col, fill = locrc_col))
  grid.text("LOCRC", x = unit(0.57, "npc") + unit(4, "mm"),
            y = unit(0.5, "npc"), just = "left",
            gp = gpar(fontfamily = "sans", fontsize = 12.5))
  popViewport()

  pushViewport(viewport(layout.pos.row = 2))
  print(p, new_page = FALSE)
  popViewport(2)
}

# --- 5. Export ------------------------------------------------------------
# Use the portable sans-serif family. Cairo produces a vector PDF.
grDevices::cairo_pdf(
  filename = result_file("DNAmTL_forestplot_final.pdf"),
  width = 11.69, height = 8.27, family = "sans"
)
draw_final_plot(p)
dev.off()

png(
  filename = result_file("DNAmTL_forestplot_final.png"),
  width = 3508, height = 2480, res = 300,
  type = "cairo", bg = "white"
)
draw_final_plot(p)
dev.off()

# Also show in RStudio
#if (interactive()) draw_final_plot(p)


# --- 6. Save as editable PowerPoint ---------------------------------------
# forestplot objects are grid graphics, NOT ggplot objects.
# Use rvg::dml(code = ...) instead of dml(ggobj = ...).
library(officer)
library(rvg)

ppt <- read_pptx()
ppt <- add_slide(ppt, layout = "Blank", master = "Office Theme")

# Maintain A4-landscape aspect ratio on the PowerPoint slide.
slide_dims <- slide_size(ppt)
aspect_ratio <- 11.69 / 8.27
margin_in <- 0.30
plot_width <- min(
  slide_dims$width - 2 * margin_in,
  (slide_dims$height - 2 * margin_in) * aspect_ratio
)
plot_height <- plot_width / aspect_ratio

ppt <- ph_with(
  ppt,
  value = dml(
    code = draw_final_plot(p),
    editable = TRUE,
    fonts = list(sans = "sans")
  ),
  location = ph_location(
    left = (slide_dims$width - plot_width) / 2,
    top = (slide_dims$height - plot_height) / 2,
    width = plot_width,
    height = plot_height
  )
)

print(ppt, target = result_file("DNAmTL_forestplot_editable.pptx"))

save_session_info("06_forest")
