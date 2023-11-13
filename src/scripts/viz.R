####╔═════      ═════╗####
####💠Visualizations💠####
####╚═════      ═════╝####

cli_h2("┗ [SCRIPTS] Loading vizualisation functions")

#----------------------------------#
####🔺Correlation / composition ####
#----------------------------------#

corr_matrix_plot <- function(df, vars, title = "") {
  return(df
    |> mutate(
      across(where(is.character), factor),
      across(where(is.factor), label_encoding)
    )
    |> correlation(select = vars, include_factors = TRUE, redundant = TRUE, method = "auto")
    |> rename(R = matches("^r$|^rho$"))
    |> mutate(across(matches("Parameter[1-2]"), \(x) factor(x, levels = vars)))
    |> ggplot(aes(x = Parameter1, y = Parameter2))
      + geom_tile(aes(fill = R), colour = "white", linewidth = 1.2, stat = "identity")
      + geom_text(aes(label = round(R, 2), colour = abs(R) > 0.5), size = rel(4.5))
      + scale_color_manual(values = c("black", "white"))
      + scale_fill_gradient2(na.value = "white", breaks = seq(-1, 1, 0.2), limits = c(-1, 1))
      + scale_x_discrete(position = "top")
      + scale_y_discrete(limits = rev)
      + guides(fill = guide_colourbar(title = "R", barheight = rel(17), title.hjust = 0.15), colour = "none")
      + labs(title = title)
      + theme(
        plot.title = element_markdown(hjust = 0.5)
        , axis.title.x = element_blank()
        , axis.title.y = element_blank()
        , axis.text.x = element_text(face = "bold", angle = 30, hjust = 0, size = 8)
        , axis.text.y = element_text(face = "bold", angle = 45, hjust = 1, size = 8)
      )
  )
}

compositional_plot <- function(dat, responses, prefix) {
  
  pattern <- paste0("^(", paste0(prefix, collapse = "|"), ")_(", ".*)$")
  pattern_prop <- paste0("^Prop_(", paste0(prefix, collapse = "|"), ")_(", ".*)$")
  responses <- str_subset(responses, paste0("^(", paste0(prefix, collapse = "|"), ")", "_[^Tot]"))
  
  dat <- dat |> select(Sample, Mouse, Stage, Condition, any_of(responses))
  
  dat_long <- (left_join(
      dat |>
        pivot_longer(cols = any_of(responses), names_pattern = pattern, names_to = c("Type", "Layer"), values_to = "Value") |>
        filter(Layer %ni% c("Tot", "Total")) |> 
        mutate(Type = factor(Type, levels = prefix)),
      dat |>
        group_by(Stage, Sample) |>
        mutate(map_dfc(prefix, \(pre) across(matches(str_glue("^{pre}_[^Tot]")), \(x) x / sum(c_across(matches(str_glue("^{pre}_[^Tot]"))), na.rm = T), .names = "Prop_{.col}"))) |>
        ungroup() |>
        pivot_longer(cols = matches("Prop_"), names_pattern = pattern_prop, names_to = c("Type", "Layer"), values_to = "Percentage") |> 
        select(Stage, Sample, Mouse, Condition, Type, Layer, Percentage)
    )
    |> group_by(Stage, Layer, Type, Condition) 
    |> summarize(
        Value = mean(Value, na.rm = TRUE),
        Percentage = mean(Percentage, na.rm = TRUE)
      ) 
    |> ungroup()
  )
  
  plot <- ggplot(dat_long, aes(x = Condition, y = Value, fill = Layer)) +
    geom_col() +
    geom_text(aes(label = scales::percent(Percentage, accuracy	= 0.1)), position = position_stack(vjust = 0.5), color = "white") +
    labs(x = " ", y = " ")
  
  if (length(unique(dat_long$Stage)) == 1) plot <- plot + facet_wrap(vars(Type), scales = "free_y", ncol = length(unique(dat_long$Type)))
  if (length(unique(dat_long$Stage)) > 1) plot <- plot + facet_grid(vars(Type), vars(Stage), scales = "free_y")
  
  return(plot)
}

#------------------#
####🔺Timelines ####
#------------------#

make_fold_timeline_plot <- function(
    dat, facet_rows = "Pathway", trans = "identity", 
    color_by = NULL, colors = colors_effect, size_boost = 1
) {
  
  origin <- do.call(trans, list(1))
  
  dat <- (
    dat
    |> mutate(Fold_trans = do.call(trans, list(Fold)))
    |> mutate(Fold_Amp = ifelse(
      max(Fold_trans, na.rm = TRUE) - min(Fold_trans, na.rm = TRUE) != 0, 
      max(Fold_trans, na.rm = TRUE) - min(Fold_trans, na.rm = TRUE), 
      mean(Fold_trans, na.rm = TRUE)) * 0.1,
      .by = all_of(c(facet_rows, "Stage"))
    )
  )
  
  timeline <- (
    ggplot(dat)
    + { if(is.null(color_by)) aes(x = Gene, color = Fold >= 1) else aes(x = Gene, color = .data[[color_by]]) }
    + geom_linerange(aes(ymax = Fold_trans), ymin = origin, size = 2 + (size_boost * 0.5))
    + geom_hline(yintercept = origin, size = 0.3, linetype = "dotted")
    + geom_text(aes(
        label = str_c(round(Fold, 2), stars.pval(p.value) |> str_replace(fixed("."), ""), sep = " "), 
        y = ifelse(Fold_trans > origin, Fold_trans + Fold_Amp, Fold_trans - Fold_Amp),
        hjust = ifelse(Fold > 1, 0, 1)
      ),
      vjust = 0.5, angle = 0, size = 2 + (size_boost * 0.25), check_overlap = TRUE
    )
    + scale_color_manual(" ", values = colors)
    + scale_y_continuous(breaks = c(0,1,2,3), expand = expansion(mult = 1.01 * (1 + (size_boost/100))))
    + scale_x_discrete(expand = expansion(add = 1 * size_boost), limits = \(x) rev(x))
    + labs(
      x = "",
      y = ifelse(trans != "identity", str_glue("Fold Change *({trans} scale)*"), "Fold Change")
    )
    + coord_flip()
    + facet_grid(
      vars(.data[[facet_rows]]), vars(Stage), 
      scales = "free_y", space = "free_y", labeller = label_wrap_gen(width = 12, multi_line = TRUE)
    )
    + { if (!is.null(color_by)) guides(color = guide_legend(title = color_by)) }
    + theme(
      legend.position = ifelse(is.null(color_by), "none", "bottom")
      , axis.text.x = element_blank()
      , axis.title.x = element_markdown(size = 9)
      , axis.text.y = element_text(size = 7)
      , strip.text = element_text(size = 5 * size_boost)
      , plot.title = element_markdown(size = 9, face = "plain", vjust = 1, hjust = 0.5)
    )
  )
  
  return(timeline)
}

#-----------------#
####🔺Circlize ####
#-----------------#

make_circlize_plot <- function(dat, resolution = 15000, save_to = NULL) {
  
  stage_list <- read_excel(configs$data$PCR$data_dict, sheet = "Stage")$Name
  layer_list <- read_excel(configs$data$PCR$data_dict, sheet = "Layer")$Name |> discard(\(x) x == "Whole")
  pathway_families <- read_excel(configs$data$PCR$data_dict, sheet = "Pathway_family")$Name
  
  crossed_layer_order <- (
    crossing(Stage = stage_list, Layer = layer_list)
    |> mutate(across(c(Stage, Layer), \(x) to_factor(x, configs$data$PCR$data_dict)))
    |> arrange(Stage, Layer)
    |> unite(c(Stage, Layer), col = "Layer", sep = "_", remove = FALSE)
    |> pull(Layer)
  )
  
  circlize_data <- (
    dat
    |> filter(p.value <= .05)
    |> select(Stage, Layer, Gene, Pathway, Pathway_family, Fold, p.value, Expression)
    |> unite(c(Stage, Layer), col = "Layer", sep = "_", remove = FALSE)
    |> mutate(Layer = factor(Layer, levels = crossed_layer_order))
    |> arrange(Layer)
    |> droplevels()
  )
  
  df_chord <- (
    circlize_data 
    |> distinct(Layer, Pathway, Gene)
    |> summarize(N_reg = n(), .by = c(Layer, Pathway))
    |> mutate(Layer = as.character(Layer))
  )
  
  mat_chord <- (
    xtabs(
      N_reg ~ Layer + Pathway, 
      data = df_chord |> 
        mutate(
          Layer = factor(Layer, levels = crossed_layer_order),
          Pathway = to_factor(Pathway, configs$data$PCR$data_dict)
        ) |> 
        arrange(Layer),
      drop.unused.levels = TRUE
    )
  )
  
  df_layer <- (
    tibble(Sector = rownames(mat_chord))
    |> mutate(N_Connection = rowSums(mat_chord)) 
    |> filter(N_Connection > 0)
  )
  
  df_pathway <- (
    tibble(Sector = colnames(mat_chord))
    |> mutate(N_Connection = colSums(mat_chord)) 
    |> filter(N_Connection > 0)
  )
  
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  df_extra <- (
    bind_rows(df_layer, df_pathway)
    |> mutate(
      Color = col_vector[1:n()],
      Layer = ifelse(str_detect(Sector, "P[0-9]+_"), str_remove(Sector, "P[0-9]+_"), NA),
      Stage = str_extract(Sector, "P[0-9]+"),
    )
    |> select(Sector, Layer, Stage, Color, N_Connection)
  )
  
  get_genes_barplot <- function(index, theta) {
    return(
      circlize_data
      |> filter(Layer == index)
      |> distinct(Gene, Fold, p.value, Expression)
      |> mutate(
        Orientation = theta,
        Label = ifelse(
          Orientation < 90 | Orientation > 270,
          as.character(str_glue("{Gene} ({round(Fold, 2)}) {gtools::stars.pval(p.value)}")),
          as.character(str_glue("{gtools::stars.pval(p.value)} ({round(Fold, 2)}) {Gene}"))
        ),
        Fold = ifelse(
          Fold >= 1,
          scales::rescale(Fold, to = c(0, 1), from = c(1, max(circlize_data$Fold))),
          scales::rescale(Fold, to = c(-1, 0), from = c(min(circlize_data$Fold), 1))
        ),
        Color = ifelse(str_detect(Expression, regulation_type$UPREG), colors_fold[2], colors_fold[1])
      )
      |> arrange(Gene)
      |> select(Fold, Label, Color)
    )
  }
  
  # Drawing the plot
  
  if (!is.null(save_to)) png(save_to, width = resolution, height = resolution)
  
  circos.clear()
  
  circos.par(cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE, start.degree = 180)
  
  colors <- c(
    "#03937e","#02ccae","#86e3d5","#b9f0e7",
    "#0264cc","#3680cf","#66a0de","#95bce6","#bed1e6",
    "#5204b3","#7835cc","#9965db","#ac85de","#d1bbed",
    "#76187a", "#b625bb","#e37ae6","#eeb1f0", 
    "#8c013b", "#b8024e","#f589b6",
    "#f02b00","#faa491","#f59002","#ffda8f","#81c700", "#d1eb52", "#01944b", "#78e3a1"
  )
  
  ## Chord diagram
  
  h0 <- resolution / 220
  h1 <- resolution / 180
  h2 <- resolution / 90
  
  chord_plot <- circlize::chordDiagram(
    as.matrix(mat_chord),
    grid.col = colors,
    directional = 1,
    diffHeight = mm_h(h0*1.1),
    target.prop.height = mm_h(h0),
    direction.type = c("diffHeight", "arrows"),
    link.arr.type = "big.arrow",
    # big.gap = 10,
    small.gap = 0.15,
    transparency = 0.20,
    annotationTrack = "grid",
    preAllocateTracks = list(
      list(
        track.height = mm_h(h1)
      ), # Track 1 (or 3 ?)
      list(
        track.height = mm_h(h2),
        track.margin = c(mm_h(h2), 0)
      ), # Track 2
      list(
        track.height = mm_h(h1)
      ) # Track 3 (or 1 ?)
    ),
    annotationTrackHeight = mm_h(h1), # Track 4
    scale = TRUE
    
  )
  
  ## Track 4 : Layer/Stage legend
  
  track_chord <- 4
  cex_chord <- resolution / 1350
  
  circos.track(
    track.index = track_chord,
    panel.fun = function(x, y) {
      
      theta <- (circlize(mean(CELL_META$xlim), 1.3)[1, 1] %% 360)[[1]]
      facing <- ifelse(theta > 220 && theta < 320,  c("outside", "bending.outside"), c("inside", "bending.inside"))
      
      circos.text(
        CELL_META$xcenter,
        CELL_META$ycenter,
        str_remove(CELL_META$sector.index, "P[0-9]+_"),
        col = "black",
        font = 1,
        cex = cex_chord,
        adj = c(0.5, 0.5),
        facing = facing
      )
    },
    bg.border = NA
  )
  
  ## Track 3 : Stage & Pathway family legend
  
  track_families <- 3
  cex_families <- resolution / 1200
  
  colors.stages <- c("#03937e", "#0264cc", "#5204b3", "#76187a", "#8c013b")
  
  for (stage in stage_list) {
    highlight.sector(
      sector.index = df_extra |> filter(Stage == stage) |> pull(Sector), # purrr::keep(df_extra$Sector, \(x) str_detect(x, "P4"))
      track.index = track_families,
      col = colors.stages[which(stage_list == stage)[[1]]],
      text = stage,
      cex = cex_families,
      font = 2,
      text.col = "white",
      facing = "bending.inside",
      niceFacing = FALSE
    )
  }
  
  colors.pf <- c("#eb4f2d","#83c108", "#f5a402", "#2cc969")
  
  for (pf in pathway_families) {
    
    circlize_data |> filter(Pathway_family == pf) |> distinct(Pathway) |> pull()
    
    highlight.sector(
      sector.index = circlize_data |> filter(Pathway_family == pf) |> distinct(Pathway) |> pull(),
      track.index = track_families,
      col = colors.pf[which(pathway_families == pf)[[1]]],
      text = pf,
      cex = cex_families,
      font = 2,
      text.col = "white",
      facing = "bending.outside",
      niceFacing = TRUE
    )
  }
  
  ## Track 2 : barplot
  
  track_barplot <- 2
  cex_barplot <- resolution / 2000
  
  circos.track(
    track.index = track_barplot,
    panel.fun = \(x, y) {
      
      if (str_detect(CELL_META$sector.index, "P[0-9]+_")) {
        
        circos.segments(
          x0 = 0,
          x1 = 1,
          y0 = 0,
          y1 = 0,
          lty = "solid"
        )
        
        draw.sector(
          get.cell.meta.data("cell.start.degree", sector.index = CELL_META$sector.index),
          get.cell.meta.data("cell.end.degree", sector.index = CELL_META$sector.index),
          rou1 = get.cell.meta.data("cell.top.radius", track.index = track_barplot),
          rou2 = get.cell.meta.data("cell.top.radius", track.index = track_barplot + 1) + mm_h(2),
        )
        
        theta <- (circlize(mean(CELL_META$xlim), 1.3)[1, 1] %% 360)[[1]]
        genes_barplot <- get_genes_barplot(CELL_META$sector.index, theta)
        # xpos <- seq(from = 0.1, to = 0.9, length.out = nrow(genes_barplot))
        # xamp <- abs(CELL_META$cell.xlim[2] - CELL_META$cell.xlim[1])
        xpos <- seq(from = CELL_META$cell.xlim[1] + 0.05 , to = CELL_META$cell.xlim[2] - 0.05 , length.out = nrow(genes_barplot))
        facing <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
        adjust <- ifelse(theta < 90 || theta > 270, c(0, 0.5), c(1, 0.5)) 
        
        
        if (nrow(genes_barplot) > 0) {
          circos.barplot(
            pos = xpos,
            value = genes_barplot$Fold,
            bar_width = 0.025,
            col = genes_barplot$Color,
            border = genes_barplot$Color
          )
          
          circos.text(
            x = xpos,
            y = get.cell.meta.data("cell.top.radius", track.index = track_barplot) + 0.2,
            labels = genes_barplot$Label,
            col = genes_barplot$Color,
            font = 2,
            cex = cex_barplot,
            adj = c(adjust),
            facing = c(facing)
          )
        }
      }
    },
    bg.border = NA
  )
  
  ## Track 1 : Empty (gene labels)
  
  if (!is.null(save_to)) dev.off()
}

#-----------------#
####🔺Boxplots ####
#-----------------#

## Generating a boxplot for an individual gene, showing the main effect of a predictor (using the model fitted to this gene's data as input)
make_signif_boxplot <- function(
    mod, xaxis = "Condition", facet = NULL, cluster = "Mouse", add_cluster_averages = TRUE, subtitle = NULL, caption = NULL, 
    invert_DCq = TRUE, scale = "link", adjust = "none", method = "pairwise", resp_name = NULL
) {
  
  get_n_units <- function(df) {
    if(!is.null(cluster) && cluster %in% colnames(df)) return(length(unique(df[[cluster]])))
    else return(dplyr::tally(df))
  }
  
  dat <- insight::get_data(mod)
  
  if (!is.null(cluster) && cluster %ni% colnames(dat)) {
    cluster <- NULL
    add_cluster_averages <- FALSE
  }
  
  if (!is.null(cluster) 
      && cluster %in% colnames(dat) 
      && dat |> group_by(across(any_of(c(xaxis, facet, cluster)))) |> count() |> filter(n > 1) |> nrow() == 0
  ) {
    cluster <- NULL
    add_cluster_averages <- FALSE
  }
  
  resp <- insight::find_response(mod)
  if (is.null(resp_name)) resp_name <- get_response_name(resp, "IHC")
  
  is_DCq <- resp %in% c("DCq", "DCt", "dct", "dcq")
  if (is_DCq && invert_DCq) {
    dat <- dat |> mutate(DCq = -1 * DCq)
    resp_name <- "-1 * DCq"
  }
  
  ## Making sure the variables of interest are contrasts for emmeans
  dat <- dat |> mutate(across(c(any_of(c(xaxis, facet)) & where(\(c) !is.factor(c))), as.factor))
  
  extra_dat <- dat |> group_by(across(any_of(c(xaxis, facet)))) |> summarize(N = str_glue("N = {get_n_units(pick(everything()))}")) |> ungroup()
  
  max <- max(dat[[resp]])
  min <- min(dat[[resp]])
  amp <- abs(max - min)
  
  if(adjust == "none") correction <- "(uncorrected)"
  else correction <- str_glue("({adjust} corrected)")
  
  # -----------[ Contrasts ]----------- #
  
  specs <- paste0(" ~ ", xaxis)
  if(!is.null(facet)) specs <- paste0(specs, " | ", facet)
  specs <- as.formula(specs)
  
  emms <- emmeans::emmeans(mod, specs = specs, type = "response", data = insight::get_data(mod))
  if (tolower(scale) %in% c("response", "resp")) emm <- regrid(emm, transform = "response")
  
  contrasts <- emmeans::contrast(emms, method = method, adjust = adjust, infer = TRUE) |> 
    as_tibble() |> 
    rename(Contrast = contrast) |> 
    tidyr::extract(col = Contrast, into = c("X1", "X2"), regex = "(.*) [- | /] (.*)", remove = FALSE)
  
  p_data_contrasts <- (
    contrasts
    |> group_by(across(any_of(c(facet))))
    |> mutate(
      x1 = match(X1, levels(dat[[xaxis]])),
      x2 = match(X2, levels(dat[[xaxis]])),
      p.signif = label_pval(p.value)
    ) 
    |> arrange(x.diff := abs(x2 - x1))
    |> mutate(
      step = 1:n(),
      pos.x = (x2 + x1) * 0.5,
      pos.y = max + step * 0.1 * (max - min)
    ) 
    |> ungroup()
    |> filter(p.signif <= alpha)
  )
  
  # -----------[ Plot ]----------- #
  
  plot <- (
    ggplot(dat, aes(x = .data[[xaxis]], y = .data[[resp]], color = .data[[xaxis]], fill = .data[[xaxis]]))
    + geom_boxplot(outlier.alpha = 0, size = 1.1, fill = NA)
    + stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)), width = 0.75, size = 1.1, linetype = "dotted")
    + { if (!is.null(cluster)) geom_jitter(
     data = \(x) x |> group_by(across(any_of(c(xaxis, facet)))) |> group_modify(\(d, g) slice_sample(d, n = min(nrow(d), 50))) |> ungroup(), 
     size = 1.5, width = 0.1, alpha = 0.3
    )
     else geom_jitter(
       data = \(x) x |> group_by(across(any_of(c(xaxis, facet)))) |> group_modify(\(d, g) slice_sample(d, n = min(nrow(d), 50))) |> ungroup(), 
       mapping = aes(fill = .data[[xaxis]]), shape = 23, color = color_text, size = 3, width = 0.1, alpha = 0.9
     )
    }
    + {if (add_cluster_averages) stat_summary(
     aes(group = .data[[cluster]], fill = .data[[xaxis]]), geom = "point", fun = mean, 
     size = ifelse(is.null(facet), 4, 3), shape = 23, color = color_text, alpha = 0.9, position = position_dodge(0.2)
    )}
    + geom_errorbarh(
     data = p_data_contrasts, aes(xmin = x1, xmax = x2, y = pos.y), inherit.aes = FALSE, 
     color = "black", height = 0.03 * amp, linewidth = 0.5
    )
    + geom_text(
     data = p_data_contrasts, aes(x = pos.x, y = pos.y, label = p.signif), inherit.aes = FALSE,
     size = 5, color = "black", fontface = "bold", vjust = 0, hjust = 0.5, position = position_nudge(y = 0.02 * amp)
    )
    + geom_label(
     aes(y = min - 0.05 * amp, fontface = "bold", label = N, color = .data[[xaxis]]),
     data = extra_dat, fill = NA, size = 5, alpha = 0.7
    )
    + theme(
     legend.position = "none", 
     panel.grid.major = element_blank(), 
     panel.grid.minor = element_blank(), 
     axis.title.x = element_blank(),
     plot.subtitle = ggtext::element_markdown(hjust = 0.5, face = "plain"),
     plot.caption = element_text(hjust = 0.5, face = "plain", size = 9)
    )
    + labs(y = resp_name)
    + {if(!is.null(subtitle)) labs(subtitle = subtitle)}
    + {if(!is.null(caption)) labs(caption = caption)}
    + {if (!is.null(facet)) facet_wrap( ~ .data[[facet]])}
    + {if (add_cluster_averages) labs(caption = str_glue("Small round points are individual measurements\n Diamonds represent {cluster}-averages"))}
  )
  
  return(plot)
}

## Generating a boxplot for an individual gene, showing interaction effects between two predictors (using the model fitted to this gene's data as input)
make_signif_boxplot_inter <- function(
    mod, pred1 = "Condition", pred2, cluster = NULL, add_cluster_averages = FALSE, invert_DCq = TRUE, stage = NULL,
    scale = "link", adjust = "none", resp_name = NULL
) {
  
  get_n_units <- function(df) {
    if(!is.null(cluster) && cluster %in% colnames(df)) return(length(unique(df[[cluster]])))
    else return(dplyr::tally(df))
  }
  
  dat <- insight::get_data(mod)
  resp <- insight::find_response(mod)
  
  is_DCq <- resp %in% c("DCq", "DCt", "dct", "dcq")
  
  if (is_DCq && invert_DCq) dat <- dat |> mutate({{ resp }} := -1 * .data[[resp]]) # -1 * DCq for better legibility
  
  if (is.null(resp_name)) {
    if (is_DCq) resp_name <- ifelse(invert_DCq, str_glue("-1 * {resp}"), str_glue("{resp}"))
    else resp_name <- get_response_name(resp, "IHC")
  }
  
  ## Making sure the variables of interest are contrasts for emmeans
  dat <- dat |> mutate(across(c(any_of(c(pred1, pred2)) & where(\(c) !is.factor(c))), as.factor))
  
  extra_dat <- dat |> group_by(across(any_of(c(pred1, pred2)))) |> summarize(N = str_glue("N = {get_n_units(pick(everything()))}")) |> ungroup()
  
  max <- max(dat[[resp]])
  min <- min(dat[[resp]])
  amp <- abs(max - min)
  
  # -----------[ Contrasts ]----------- #
  
  specs <- paste0(" ~ ", pred1)
  if(!is.null(pred2)) specs <- paste0(specs, " | ", pred2)
  specs <- as.formula(specs)
  
  emmeans <- emmeans::emmeans(mod, specs = specs, type = "response", data = insight::get_data(mod))
  if (tolower(scale) %in% c("response", "resp")) emmeans <- regrid(emmeans, transform = "response")
  
  contrasts <- emmeans::contrast(emmeans, method = "pairwise", adjust = adjust, infer = TRUE) |> 
    as.data.frame() |> 
    rename(Contrast = contrast) |> 
    tidyr::extract(col = Contrast, into = c("X1", "X2"), regex = "(.*) [- | /] (.*)", remove = FALSE)
  
  p_data_contrasts <- contrasts |>
    group_by(across(any_of(c(pred2)))) |>
    mutate(
      x1 = (match(.data[[pred2]], levels(dat[[pred2]])) - 1) * length(unique(dat[[pred1]])) + match(X1, levels(dat[[pred1]])),
      x2 = (match(.data[[pred2]], levels(dat[[pred2]])) - 1) * length(unique(dat[[pred1]])) + match(X2, levels(dat[[pred1]])),
      p.signif = label_pval(p.value)
    ) |>
    arrange(x.diff := abs(x2 - x1)) |>
    mutate(
      step = 1:n(),
      pos.x = (x2 + x1) * 0.5,
      pos.y = max + step * 0.1 * (max - min)
    ) |>
    ungroup()
  
  contrasts_interactions <- emmeans::contrast(emmeans, interaction = c("pairwise"), by = NULL, adjust = "none", infer = TRUE) |> 
    as.data.frame() |> 
    tidyr::extract(col = paste0(pred1, "_pairwise"), into = c("X1", "X2"), regex = "(.*) [- | /] (.*)", remove = FALSE) |> 
    tidyr::extract(col = paste0(pred2, "_pairwise"), into = c("F1", "F2"), regex = "(.*) [- | /] (.*)", remove = FALSE)
  
  p_data_interactions <- contrasts_interactions |>
    mutate(
      x1 = 0.5 * ((match(F1, levels(dat[[pred2]])) - 1) * length(unique(dat[[pred1]])) + match(X1, levels(dat[[pred1]])) +
                    (match(F1, levels(dat[[pred2]])) - 1) * length(unique(dat[[pred1]])) + match(X2, levels(dat[[pred1]]))),
      x2 = 0.5 * ((match(F2, levels(dat[[pred2]])) - 1) * length(unique(dat[[pred1]])) + match(X1, levels(dat[[pred1]])) +
                    (match(F2, levels(dat[[pred2]])) - 1) * length(unique(dat[[pred1]])) + match(X2, levels(dat[[pred1]]))),
      p.signif = label_pval(p.value)
    ) |>
    arrange(x.diff := abs(x2 - x1)) |>
    mutate(
      step = 1:n() + choose(length(unique(dat[[pred1]])), 2),
      pos.x = (x2 + x1) * 0.5,
      pos.y = max + step * 0.1 * (max - min)
    )
  
  # -----------[ Plot ]----------- #
  
  plot <- (
    ggplot(dat, aes(x = interaction(.data[[pred1]], .data[[pred2]], sep = "_"), y = .data[[resp]], color = .data[[pred1]]))
    + geom_boxplot(outlier.alpha = 0, size = 1.1, fill = NA)
    + stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)), width = 0.75, size = 1.1, linetype = "dotted")
    + { 
     if (!is.null(cluster)) geom_jitter(
       data = \(x) x |> group_by(across(any_of(c(pred1, pred2)))) |> group_modify(\(d, g) slice_sample(d, n = min(nrow(d), 50))) |> ungroup(), 
       size = 1.5, width = 0.1, alpha = 0.3
     )
     else geom_jitter(
       data = \(x) x |> group_by(across(any_of(c(pred1, pred2)))) |> group_modify(\(d, g) slice_sample(d, n = min(nrow(d), 50))) |> ungroup(), 
       mapping = aes(fill = .data[[pred1]]), shape = 23, color = color_text, size = 3, width = 0.1, alpha = 0.9
     )
    }
    + {
     if (add_cluster_averages) stat_summary(
       aes(group = .data[[cluster]], fill = .data[[pred1]]), geom = "point", fun = mean, 
       size = 3, shape = 23, color = color_text, alpha = 0.9, position = position_dodge(0.2)
     )
    }
    + geom_errorbarh(
     data = p_data_contrasts, aes(xmin = paste(X1, .data[[pred2]], sep = "_"), xmax = paste(X2, .data[[pred2]], sep = "_"), y = pos.y), inherit.aes = FALSE,
     color = "black", height = 0.02 * amp, size = 0.5
    )
    + geom_text(
     data = p_data_contrasts, aes(x = pos.x, y = pos.y, label = p.signif), inherit.aes = FALSE,
     size = 5, color = "black", fontface = "bold", vjust = 0, hjust = 0.5, position = position_nudge(y = 0.02 * amp)
    )
    + geom_label(
     aes(y = min - 0.05 * amp, fontface = "bold", label = N, color = .data[[pred1]]), 
     data = extra_dat, fill = NA, size = 5, alpha = 0.7
    )
    ## Interactions
    + geom_errorbarh(
     data = p_data_interactions, aes(xmin = x1, xmax = x2, y = pos.y), inherit.aes = FALSE,
     color = "black", height = 0.02 * amp, size = 0.5
    )
    + geom_text(
     data = p_data_interactions, aes(x = pos.x, y = pos.y, label = p.signif), inherit.aes = FALSE,
     size = 5, color = "black", fontface = "bold", vjust = 0, hjust = 0.5, position = position_nudge(y = 0.02 * amp)
    )
    + theme(
     legend.position = "none",
     plot.subtitle = ggtext::element_markdown(hjust = 0.5, face = "plain")
    )
    + labs(y = resp_name, x = str_c(pred1, " by ", pred2))
    + {if(!is.null(stage)) labs(subtitle = str_glue("{stage}"))}
    + {if (add_cluster_averages) labs(caption = str_glue("Small round points are individual measurements\n Diamonds represent {cluster}-averages"))}
    + scale_x_discrete(labels = \(l) str_replace(l, "_", "\n"))
  )
  
  return(plot)
}

#--------------------------#
####🔺Model diagnostics ####
#--------------------------#

make_acf_plot <- function(mod) {
  forecast::ggAcf(residuals(mod, type = "response", retype = "normalized"), color = "#1b6ca8") + 
    geom_point() + 
    labs(
      title = "Autocorrelation of residuals",
      subtitle = "Data (lines) should be inside the dashed area"
    ) + 
    see::theme_lucid() +
    theme(
      plot.title = element_markdown(size = 15, margin = margin(0, 0, 0, 1)),
      plot.subtitle = element_markdown(size = 12, margin = margin(1, 0, 0, 0))
    )
}

ppc_plots <- function(mod, simulations, term = "Condition", type = "fixed", is_count = NULL, max_cols_per_plot = 3) {
  
  Y <- insight::get_response(mod)
  n_unique <- n_distinct(insight::get_data(mod)[[term]])
  
  if(is.null(is_count)) is_count <- ifelse(insight::get_family(mod)$family |> str_detect("binom|poiss"), TRUE, FALSE)
  
  # ppc_fun <- ifelse(is_count, bayesplot::ppc_bars, bayesplot::ppc_dens_overlay)
  ppc_fun_grouped <- ifelse(is_count, bayesplot::ppc_bars_grouped, bayesplot::ppc_dens_overlay_grouped)
  ppc_fun_pred_grouped <- bayesplot::ppc_intervals_grouped
  
  if(type %in% c("fixed", "fe")) {
    .term <- insight::get_predictors(mod)[[term]]
    
    # ppc_global <- ppc_fun(y = Y, yrep = simulations) 
    if(is_count) ppc_root <- bayesplot::ppc_rootogram(Y, simulations, style = "suspended")
    
    ppc_grouped <- ppc_fun_grouped(Y, simulations, group = .term) + 
      facet_wrap(~ group, ncol = min(max_cols_per_plot, n_unique), scales = "free")
    ppc_pred_grouped <- ppc_fun_pred_grouped(Y, simulations, group = .term, prob_outer = 0.95) + 
      facet_wrap(~ group, ncol = min(max_cols_per_plot, n_unique), scales = "free")
  }
  else if(type %in% c("random", "re")) {
    .term <- insight::get_random(mod)[[term]]
    
    ppc_grouped <- ppc_fun_grouped(Y, simulations, group = .term) + 
      facet_wrap(~ group, ncol = min(max_cols_per_plot, n_unique), scales = "free")
    ppc_pred_grouped <- ppc_fun_pred_grouped(Y, simulations, group = .term, prob_outer = 0.95) + 
      facet_wrap(~ group, ncol = min(max_cols_per_plot, n_unique), scales = "free")
  }
  
  return(
    if(type %in% c("fixed", "fe")) {
      if(is_count) { 
        (ppc_root / ppc_grouped / ppc_pred_grouped) + plot_layout(guides = 'collect', ncol = 1, nrow = 3) +
          # plot_annotation(title = "Simulation-based Posterior Predictive Checks", subtitle = str_glue("For [{term}]")) & 
          theme(legend.position = 'right', axis.title.x = element_blank())
      } else {
        (ppc_grouped / (ppc_pred_grouped + theme(axis.title.x = element_blank()))) + plot_layout(ncol = 1, nrow = 2) + 
          # plot_annotation(title = "Simulation-based Posterior Predictive Checks", subtitle = str_glue("For [{term}]")) & 
          theme(legend.position = 'right')
      }
    }
    else list(ppc_grouped, ppc_pred_grouped)
  )
}

ppc_stat_plots <- function(mod, simulations, term = "Condition", type = "fixed", stats = c("min", "max", "mean", "sd"), n_cols = 2, max_cols_per_plot = 5) {
  
  n_unique <- n_distinct(insight::get_data(mod)[[term]])
  
  if(type %in% c("fixed", "fe")) .term <- insight::get_predictors(mod)[[term]]
  else if(type %in% c("random", "re")) .term <- insight::get_random(mod)[[term]]
  
  return(
    patchwork::wrap_plots(
      purrr::map(
        stats, 
        \(.x) bayesplot::ppc_stat_grouped(
          insight::get_response(mod), 
          simulations, group = .term, stat = .x,
          facet_args = list(ncol = min(max_cols_per_plot, n_unique))
        ) + scale_x_continuous(labels = \(l) signif(l, digits = 2))
      ), 
      ncol = n_cols, guides = 'auto'
    ) + 
      # plot_annotation(title = "Simulation-based Predictive Checks (on statistics)", subtitle = str_glue("For [{term}]")) & 
      theme(legend.position = 'bottom', axis.text.x = element_text(size = rel(1.5), angle = 30, hjust = 1))
  )
}

#-----------------#
####🔺Sunburst ####
#-----------------#

make_suburst_plot <- function(dat, layers, tooltips = NULL, colors = NULL, plot_options = list()) {
  
  if (!is.null(colors)) colors <- set_names(colors, \(x) x |> str_replace_all(" ", "\n") |> str_replace_all("and\n", "and "))
  
  .make_sunburst_data <- function(dat, layers_part) {
    
    if (all(layers == layers_part))
      sun_dat <- summarize(dat, across(any_of(c(layers_part, tooltips)), first), .by = all_of(layers_part))
    else 
      sun_dat <- summarize(dat, across(any_of(layers_part), first), .by = all_of(layers_part))
    
    return(
      bind_cols(
        transmute(sun_dat, ids = pick(any_of(layers_part)) |> pmap_chr(\(...) str_c(..., sep = ">>"))),
        transmute(
          sun_dat, 
          parents = pick(any_of(layers_part)) |> 
            select(-last_col()) |> 
            pmap_chr(\(...) str_c(..., sep = ">>")) |> 
            bind(x, if (is_empty(x)) "" else x)
        ),
        transmute(sun_dat, labels = pick(any_of(layers_part)) |> select(last_col()) |> pull(1))
      )
      |> bind(
        x, 
        if (all(layers == layers_part) & !is.null(tooltips)) 
          bind_cols(
            x, 
            transmute(
              sun_dat, 
              hovertext = pmap_chr(
                pick(all_of(tooltips)), 
                \(...) str_c(
                  names(c(...)), 
                  map_if(list(...), is.numeric, \(x) round(x, 3)) |> map_if(is.factor, as.character), 
                  sep = ": ", 
                  collapse = "\n"
                )
              )
            )
          )
        else x
      )
      |> mutate(labels = labels |> str_replace_all(" ", "\n") |> str_replace_all("and\n", "and "))
    )
  }
  
  params <- map_dfr(1:length(layers), \(x) .make_sunburst_data(dat, layers[1:x])) |> 
    mutate(bg_color = map_chr(labels, \(x) colors[[x]] %||% NA))
  
  rlang::exec(
    plotly::plot_ly,
    !!!select(params, -bg_color),
    type = "sunburst",
    branchvalues = "total",
    hoverinfo = "text",
    sort = FALSE,
    marker = list(colors = ~params$bg_color),
    !!!plot_options
  )
}


#----------------#
####🔺Heatmap ####
#----------------#

make_heatmap <- function(data, xaxis, facet) {
  
  max_upreg <- data |> filter(Fold >= 1) |> pull(Fold) |> max()
  
  return(
    ggplot(data, aes(x = .data[[xaxis]], y = Gene))
    + scale_fill_gradient(name = regulation_type$UPREG, low = "seagreen2", high = "seagreen4", limits = c(1, max_upreg), trans = scales::log10_trans(), labels = \(x) round(x, 2))
    + geom_tile(data = \(d) filter(d, Fold >= 1), aes(fill = Fold), colour = "white")
    + new_scale("fill")
    + scale_fill_gradient(name = regulation_type$DOWNREG, low = "firebrick4", high = "firebrick2", limits = c(0, 1), labels = \(x) round(x, 2))
    + geom_tile(data = \(d) filter(d, Fold < 1), aes(fill = Fold), colour = "white")
    + geom_text(aes(label = round(Fold, 2)), size = 3, colour = "white", fontface = "bold", check_overlap = TRUE)
    + facet_grid(cols = vars(.data[[facet]]), scales = "free_x", space = "free_x")
    + theme_light_mar
    + theme(
      legend.title = element_markdown(face = "bold", vjust = 0.80),
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(hjust = 0),
      axis.ticks = element_line(linewidth = 0.4)
    )
    + labs(x = xaxis, y = "")
  )
}


#---------------#
####🔺Tables ####
#---------------#

make_data_dict_reactable <- function(data, data_dict_path) {
  data_dict <- read_excel(data_dict_path, sheet = 1) |> 
    filter(Name %in% colnames(data))
  
  get_var_details <- function(idx) {
    sheet <- pluck(data_dict[idx, ], 1)
    
    if (sheet %in% excel_sheets(data_dict_path))
      reactable(read_excel(data_dict_path, sheet = sheet), fullWidth = TRUE, compact = TRUE)
  }
  
  reactable(
    data_dict
    , defaultColDef = colDef(sortNALast = TRUE, vAlign = "center", headerVAlign = "center", align = "center")
    , columns = list(Description = colDef(html = TRUE), Label = colDef(html = TRUE))
    , details = \(idx) get_var_details(idx)
    , outlined = TRUE
    , striped = TRUE
    , highlight = TRUE
    , compact = TRUE
    , fullWidth = TRUE
    , defaultPageSize = 20
  )
}


make_diag_reactable <- function(dat) {
    
  make_diagnostic_plots <- function(mod) {
    
    diag_list <- c("normality", "homogeneity", "qq", "pp_check")
    if (!is.null(find_random(mod))) diag_list <- c("normality", "homogeneity", "qq", "reqq", "pp_check")
    
    return(
      map(diag_list, \(x) plot(suppressMessages({check_model(mod, check = x, panel = FALSE)}))[[1]])
      |> wrap_plots(ncol = 2)
      |> girafe(ggobj = _)
    )
  }
  
  diagnostic_plot_list <- map(dat$Mod, \(mod) make_diagnostic_plots(mod), .progress = "Generating diagnostic plots:")
  
  reactable(
    dat
    , defaultColDef = colDef(sortNALast = TRUE, vAlign = "center", headerVAlign = "center", align = "center")
    , columns = list(
      Mod = colDef(
        details = \(idx) diagnostic_plot_list[[idx]],
        cell = \(x, idx) "Click for diagnostics"
      ),
      Fold = colDef(format = colFormat(digits = 3))
    )
    , outlined = TRUE
    , striped = TRUE
    , highlight = TRUE
    , compact = TRUE
    , fullWidth = TRUE
    , defaultPageSize = 5
  )
}