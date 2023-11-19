####╔════  ═══╗####
####💠 Theme 💠####
####╚════  ═══╝####

cli_h2("┗ [CONFIG] Setting project's plots & tables' themes")

#---------------#
####🔺Colors ####
#---------------#

bg_color_light <- "white"
primary_color_light <- "black"
secondary_color_light <- "#0d6efd"

bg_color_dark <- "#222"
primary_color_dark <- "white"
secondary_color_dark <- "#20c997"

strip_color <- "#adb5bd"

color_text <- "black"
color_text_bi <- "#a8aeb4"

colors_cond <- c("#02ccae", "#9965db") # N, IH (factor)
colors_fold <- c("#f02b00", "#01944b") # Downreg, Upreg
colors_effect <- c("#b625bb", "#0264cc") # Bad, Good

my_palettes_d <- list(colors_cond, scales::viridis_pal()(3), scales::viridis_pal()(5), scales::viridis_pal()(10), scales::viridis_pal()(15))

sunburst_pcr_colors <- list(
  # Stages
  P4 = "#95bfed", P8 = "#6498d0", P12 = "#4681bf", P21 = "#316cab", P70 = "#1e5795",
  # Layers
  EGL = "#b398d4", EGLo = "#ab8dcf", EGLi = "#9978c2", ML = "#8b66b7", MLPC = "#7d57ac", PC = "#6d479d", IGL = "#653b98", WM = "#532a86", Whole = "#321359",
  # Expression
  Downregulated = "#f02b00", Upregulated = "#01944b",
  # ND pathway families
  `Proliferation and Repair` = "#eb4f2d", `Migration and Response` = "#83c108", `Cellular Differentiation` = "#f5a402", `Cell Communication` = "#2cc969",
  # ND pathways
  `Fate Mapping` = "#f02b00", Survival = "#faa491", `Myelin Sheath` = "#f59002", `Neurite Growth` = "#ffda8f", Guidance = "#81c700", Motility = "#d1eb52", Membrane = "#01944b", Soluble = "#78e3a1",
  # OS pathways
  `Apoptotic Pathways` = "#FFE792", `Autophagy and Mitophagy` = "#E78AC3", `Cell Death and Protection` = "#66C2A4", 
  `Inflammation Pathways` = "#FC6A4A", `Antioxidant Response` = "#A5D853", `Redox Enzymes` = "#3690C0",  `ROS Production` = "#FEB24C"
)

#----------------------#
####🔺ggplot themes ####
#----------------------#

library(ggplot2)

options(
  ggplot2.discrete.colour = my_palettes_d,
  ggplot2.discrete.fill = my_palettes_d,
  ggplot2.continuous.colour = \() scale_color_viridis_c(),
  ggplot2.continuous.fill = \() scale_fill_viridis_c(),
  ggplot2.binned.colour = \() scale_color_viridis_b(),
  ggplot2.binned.fill = \() scale_fill_viridis_b()
)

base_theme_mar <- ggplot2::theme_minimal() +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "transparent", colour = NA),
    ## Panel
    panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
    panel.border = ggplot2::element_rect(fill = "transparent"),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    ## Titles
    plot.title = ggtext::element_markdown(size = 18, face = "bold"),
    plot.subtitle = ggtext::element_markdown(size = 15, face = "italic"),
    ## Legend
    legend.title = ggplot2::element_text(size = 16, face = "bold"),
    legend.text = ggplot2::element_text(size = 15),
    legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
    legend.key = ggplot2::element_rect(fill = "transparent", colour = NA),
    ## Facets
    strip.text = ggplot2::element_text(size = 16, face = "bold"),
    ## Axes
    axis.title.x = ggtext::element_markdown(size = 15, face = "bold", hjust = 0.5),
    axis.title.y = ggtext::element_markdown(size = 15, face = "bold", hjust = 0.5),
    axis.text.x = ggplot2::element_text(size = 13),
    axis.text.y = ggplot2::element_text(size = 13),
  )

light_addon_mar <- ggplot2::theme(
  text = ggplot2::element_text(color = primary_color_light),
  ## Panel
  panel.border = ggplot2::element_rect(colour = primary_color_light),
  ## Titles
  plot.title = ggtext::element_markdown(colour = primary_color_light),
  plot.subtitle = ggtext::element_markdown(colour = primary_color_light),
  ## Legend
  legend.title = ggplot2::element_text(colour = primary_color_light),
  ## Facets
  strip.background = ggplot2::element_rect(fill = strip_color),
  strip.text = ggplot2::element_text(colour = primary_color_light),
  ## Axes
  axis.text = ggplot2::element_text(colour = primary_color_light),
  axis.title.x = ggtext::element_markdown(colour = primary_color_light),
  axis.title.y = ggtext::element_markdown(colour = primary_color_light)
)

theme_light_mar <- base_theme_mar + light_addon_mar

dark_addon_mar <- ggplot2::theme(
  text = ggplot2::element_text(color = primary_color_dark),
  ## Panel
  panel.border = ggplot2::element_rect(colour = primary_color_dark),
  ## Titles
  plot.title = ggtext::element_markdown(colour = primary_color_dark),
  plot.subtitle = ggtext::element_markdown(colour = primary_color_dark),
  ## Legend
  legend.title = ggplot2::element_text(colour = primary_color_dark),
  ## Facets
  strip.background = ggplot2::element_rect(fill = strip_color),
  strip.text = ggplot2::element_text(colour = primary_color_dark),
  ## Axes
  axis.text = ggplot2::element_text(colour = primary_color_dark),
  axis.title.x = ggtext::element_markdown(colour = primary_color_dark),
  axis.title.y = ggtext::element_markdown(colour = primary_color_dark)
)

theme_dark_mar <- base_theme_mar + dark_addon_mar

ggplot2::theme_set(theme_light_mar)

#------------------#
####🔺gt themes ####
#------------------#

library(gt)

format_gt <- function(gt_tbl) {
  
  gt_tbl <- gt::fmt(
    gt_tbl,
    columns = select(gt_tbl[["_data"]], matches("p.val|^pr|pr\\(.*\\)|^p$")) |> colnames(),
    fns = \(x) purrr::map_chr(x, \(v) ifelse(!is.na(v) && utils::type.convert(v, as.is = TRUE) |> is.numeric(), label_pval(as.numeric(v)), v))
  )
  
  gt_tbl <- gt::fmt_number(
    gt_tbl,
    columns = select(gt_tbl[["_data"]], where(\(v) is.numeric(v))) |> colnames(),
    decimals = 3, drop_trailing_zeros = TRUE # n_signif = 2
  )
  
  # gt_tbl <- gt::fmt_markdown(
  #   gt_tbl,
  #   columns = select(gt_tbl[["_data"]], where(\(v) is.character(v)) & !matches("p.val|^pr|pr\\(.*\\)|^p$")) |> colnames()
  # )
  
  return(gt_tbl)
}

gt_style_light <- function(gt_tbl) {
  gt_tbl <- (gt_tbl 
    |> format_gt()
    |> gt::tab_style(
     style = list(
       cell_fill(color = bg_color_light, alpha = 1),
       cell_text(color = secondary_color_light, weight = "bold"),
       cell_borders(sides = c("top", "bottom"), color = secondary_color_light, style = "solid", weight = px(2))
     ),
     locations = list(cells_title(), cells_column_labels())
    ) 
    |> gt::tab_style(
     style = list(
       cell_fill(color = bg_color_light, alpha = 1),
       cell_text(color = primary_color_light)
     ),
     locations = list(cells_stub(), cells_body(), cells_row_groups(), cells_footnotes(), cells_source_notes())
    )
    |> gt::tab_style(
     style = list(cell_text(weight = "bold")),
     locations = list(cells_row_groups())
    )
    |> gt::tab_options(container.width = pct(100), table.width = pct(100))
  )
  
  return(gt_tbl)
}

#------------------------#
####🔺reactable theme ####
#------------------------#

options(
  reactable.theme = reactableTheme(
    headerStyle = list(
      "&:hover[aria-sort]" = list(background = "hsl(0, 0%, 96%)"),
      "&[aria-sort='ascending'], &[aria-sort='descending']" = list(background = "hsl(0, 0%, 96%)"),
      borderColor = "#555"
    )
  )
)