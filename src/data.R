####╔═════   ═════╗####
####💠Loading Data💠####
####╚═════   ═════╝####

## Factor levels
layer_list <- c("EGL", "EGLo", "EGLi", "ML", "MLPC", "PC", "IGL", "WM", "Tot")
stage_list <- c("P4", "P8", "P12", "P21", "P70")

## ⚠️ Only used for the Diff panel
pathway_families <- c("Proliferation and Repair", "Migration and Response", "Cellular Differentiation", "Cell Communication")
pathways <- c("Fate Mapping", "Survival", "Myelin Sheath", "Neurite Growth", "Guidance", "Motility", "Membrane", "Soluble")

#-------------------------#
####🔺Helper functions ####
#-------------------------#

## Helper functions to associate a regulation status (down or up-regulated) to each gene
regulation_type_enum <- function() {
  list(
    NOT_REG = "Not Regulated", 
    MAYBE_UPREG = "Maybe Upregulated", 
    UPREG = "Upregulated",
    MAYBE_DOWNREG = "Maybe Downregulated",
    DOWNREG = "Downregulated"
  )
}
regulation_type <- regulation_type_enum()

get_regulation_type <- function(fold_change, p_value) {
  case_when(
    p_value <= alpha ~ ifelse(
      fold_change < threshold.reg,
      regulation_type$DOWNREG,
      regulation_type$UPREG
    ),
    p_value <= trend ~ ifelse(
      fold_change < threshold.reg,
      regulation_type$MAYBE_DOWNREG,
      regulation_type$MAYBE_UPREG
    ),
    is.na(p_value) | is.na(fold_change) ~ NA_character_,
    .default = regulation_type$NOT_REG
  )
}

#-------------------------------#
####🔺Data loading functions ####
#-------------------------------#

#### ┗ PCR ------

## Loading and shaping the PCR data for the Differentiation panel
load_Diffs <- function(files, max_cq = 33) {
  return(
    map_dfr(
      map(files, \(file_name) rmatch(config, file_name)), 
      \(file) load_Diff(file, max_cq)
    )
  )
}

load_Diff <- function(path, max_cq = 33) {
  stage <- str_extract(path_ext_remove(path), "P(\\d{1,2}|Ad)")
  return(
    map_dfr(excel_sheets(path) |> set_names(), \(x) read_excel(path = path, sheet = x), .id = "Layer")
    |> tidyr::extract(Mouse, into = c("Litter", NA, "Condition"), regex = "^(\\w{2})(\\d{1,2})(H|N)$", convert = TRUE, remove = FALSE)
    |> normalize_litter_names()
    |> drop_na(Mean_Cq, DCq)
    |> filter(Mean_Cq <= max_cq)
    |> mutate(Stage = stage)
    |> mutate(
      Condition = factor(Condition, levels = c("N", "H")),
      Stage = factor(Stage, levels = stage_list),
      Layer = factor(Layer, levels = layer_list)
    )
    #### Adding additional information
    |> add_sex_and_weight()
    # |> add_gene_family()
    |> add_layer_families()
    |> add_pathways(sheet = "Diff")
    |> add_fold_change()
    #### Selecting / arranging data
    |> select(Mouse, Litter, Stage, Layer, Layer_Family, Pathway, Pathway_family, Gene, Condition, Sex, Weight, Litter, Mean_Cq, DCq, Fold)
    |> mutate(
      Layer_Family = factor(Layer_Family, levels = layer_list),
      Pathway = factor(Pathway, levels = pathways),
      Pathway_family = factor(Pathway_family, levels = pathway_families)
    )
    |> arrange(across(any_of(c("Experiment", "Stage", "Layer", "Gene", "Mouse"))))
    |> mutate(Id = row_number(), .before = 1)
  )
}

## Loading and shaping the PCR data for the Oxidative Stress panel (for the whole cerebellum)
load_OS_Tot <- function(path = config$data$PCR$OS_Tot_path, layer = "Tot", max_cq = 33, hk_genes = c("Hsp90", "Ppia"), update_with_cal = TRUE) {
  return(
    map_dfr(
      excel_sheets(path) |> set_names(),
      \(sheet) read_excel(path = path, sheet = sheet), 
      .id = "Stage"
    )
    |> tidyr::extract(Mouse, into = c("Litter", NA, "Condition"), regex = "^(\\w{2})(\\d{1,2})(H|N)$", convert = TRUE, remove = FALSE)
    |> normalize_litter_names()
    |> drop_na(Mean_Cq)
    |> filter(Mean_Cq <= max_cq)
    |> mutate(Layer = layer) # Adding "Tot" as default layer
    |> mutate(
      Condition = factor(Condition, levels = c("N", "H")),
      Stage = factor(Stage, levels = stage_list),
      Layer = factor(Layer, levels = layer_list),
      Plate = factor(str_c(Plate, Experiment, sep = " "))
    )
    #### Replacing Cq values with calibration data
    |> bind(
      x, 
      if (update_with_cal) 
        rows_update(
          x, 
          read_excel(config$data$calibration_path) |> 
            filter(Gene %in% unique(x$Gene)) |> 
            left_join(x |> select(Mouse, Gene, Stage), by = c("Mouse", "Gene")) |> 
            rename(Mean_Cq = Cq) |> 
            filter(Stage %in% c("P21", "P70")), 
          by = c("Mouse", "Gene")) 
      else x
    )
    #### Adding additional information
    |> add_dcq(hk_genes)
    |> add_sex_and_weight()
    |> add_fold_change()
    |> add_pathways(sheet = "OS")
    #### Selecting / arranging data
    |> mutate(Experiment = factor(Experiment))
    |> select(Experiment, Plate, Mouse, Litter, Sex, Weight, Stage, Layer, Gene, Pathway, Effect, Figure, Condition, Cq = Mean_Cq, DCq, Fold)
    |> arrange(Stage, Gene, Mouse)
    |> mutate(Id = row_number(), .before = 1)
  )
}

## Loading and shaping the PCR data for the Oxidative Stress panel (Purkinje Cells only)
load_OS_PC <- function(path = config$data$PCR$OS_PC_path, layer, max_cq = 33) {
  return(
    map_dfr(
      excel_sheets(path) |> set_names(),
      \(sheet) read_excel(path = path, sheet = sheet), 
      .id = "Stage"
    )
    |> tidyr::extract(Mouse, into = c("Litter", NA, "Condition"), regex = "^(\\w{2})(\\d{1,2})(H|N)$", convert = TRUE, remove = FALSE)
    |> normalize_litter_names()
    |> drop_na(Mean_Cq, DCq)
    |> filter(Mean_Cq <= max_cq)
    |> mutate(Layer = layer) # Adding "Tot" as default layer
    |> mutate(
      Condition = factor(Condition, levels = c("N", "H")),
      Layer = factor(Layer, levels = layer_list),
      Stage = factor(Stage, levels = stage_list)
    )
    #### Adding additional information
    |> add_sex_and_weight()
    |> add_fold_change()
    |> add_pathways(sheet = "OS")
    #### Selecting / arranging data
    |> select(Mouse, Litter, Sex, Stage, Layer, Gene, Pathway, Effect, Figure, Condition, DCq, Fold)
    |> arrange(Stage, Gene, Mouse)
    |> mutate(Id = row_number(), .before = 1)
  )
}

## Loading and shaping the PCR data for the Vascularization panel
load_Vasc <- function(path = config$data$PCR$Vasc_path, layer = "Tot", max_cq = 33) {
  return(
    map_dfr(
      excel_sheets(path) |> set_names(),
      \(sheet) read_excel(path = path, sheet = sheet), 
      .id = "Stage"
    )
    |> filter(!Outlier)
    |> tidyr::extract(Mouse, into = c("Litter", NA, "Condition"), regex = "^(\\w{2})(\\d{1,2})(H|N)$", convert = TRUE, remove = FALSE)
    |> normalize_litter_names()
    |> drop_na(Mean_Cq, DCq)
    |> filter(Mean_Cq <= max_cq)
    |> mutate(Layer = layer) # Adding "Tot" as default layer
    |> mutate(
      Condition = factor(Condition, levels = c("N", "H")),
      Stage = factor(Stage, levels = stage_list),
      Layer = factor(Layer, levels = layer_list),
      Experiment = factor(Experiment)
    )
    |> filter(Experiment != "B")
    #### Adding additional information
    |> add_sex_and_weight()
    |> add_fold_change()
    |> add_pathways(sheet = "Vasc", cols = "Pathway")
    #### Selecting / arranging data
    |> select(Mouse, Litter, Experiment, Stage, Sex, Layer, Gene, Pathway, Condition, Mean_Cq, DCq, Fold)
    |> arrange(Stage, Gene, Mouse)
    |> mutate(Id = row_number(), .before = 1)
  )
}

#### ┗ IHC ------

## Loading and shaping the Caspase IHC data
load_Casp <- function(path = config$data$IHC$casp_path, stage = NULL) {
  if (is.null(stage) || str_to_sentence(stage) == "All") {
    return(
      map_dfr(
        excel_sheets(path) |> set_names(), 
        \(x) load_Casp(path = path, stage = x), 
        .id = "Stage"
      )
    )
  } 
  else {
    return(
      read_xlsx(path, sheet = stage)
      |> tidyr::extract(
        Mouse,
        into = c("Bloodline", "MouseID", "Condition"),
        regex = "^(\\w{2})(\\d{1})([hH]+|[nN]+)$",
        convert = TRUE, remove = TRUE
      )
      |> tidyr::unite("Sample", c(Bloodline, MouseID, Condition, Slice), sep = "", remove = FALSE)
      |> tidyr::unite("Mouse", c(Bloodline, MouseID, Condition), sep = "", remove = FALSE)
      |> mutate(
        across(matches("^A_|^C_"), \(x) x / 1e6),
        Dens_EGL = N_EGL / A_EGL,
        Dens_ML_PCL = N_ML_PCL / A_ML_PCL ,
        Dens_IGL = N_IGL / A_IGL,
        Dens_WM = N_WM / A_WM,
        Dens_Tot = N_Tot / A_Tot,
        Prop_C_EGL = C_EGL / A_EGL,
        Prop_C_ML_PCL = C_ML_PCL / A_ML_PCL,
        Prop_C_IGL = C_IGL / A_IGL,
        Prop_C_WM = C_WM / A_WM,
        Prop_C_Tot = C_Tot / A_Tot,
        Stage = factor(stage, levels = stage_list),
        Condition = factor(Condition, levels = c("N", "H"), labels = c("N", "IH")),
        Z = factor(str_to_sentence(Z))
      )
      |> arrange(Condition, Sample)
      |> select(Sample, Mouse, Stage, Condition, Z, matches("_EGL"), matches("_ML_PCL"), matches("_IGL"), matches("_WM"), matches("_Tot"))
    )
  }
}

## Loading and shaping the Thickness IHC data
load_Thickness <- function(path = config$data$IHC$thickness_path, stage = "All") {
  # Loading data for all stages
  if (is.null(stage) || str_to_sentence(stage) == "All") {
    return(
      map_df(
        excel_sheets(path) |> set_names(), 
        \(x) load_Thickness(path = path, stage = x)
      )
    )
  }
  # Loading data for a specific stage
  else {
    return(
      read_xlsx(path, sheet = stage)
      |> filter(!Outlier)
      |> group_by(Stage, Mouse, Slice)
      |> mutate(MeasureID = seq(1:n()))
      |> ungroup()
      |> mutate(
        Slice = str_c(Mouse, Slice),
        Sample = str_c(Slice, MeasureID),
        across(c(Mouse, Slice), \(.x) as.factor(.x)),
        Condition = factor(Condition, levels = c("N", "H"), labels = c("N", "IH"))
      )
      |> arrange(Condition, Mouse)
      |> select(Sample, Stage, Slice, MeasureID, Mouse, Condition, everything(), -Outlier)
    )
  }
}


#-------------------------------------#
####🔺Extra data loading functions ####
#-------------------------------------#

## Load housekeeping genes data
load_HK <- function(path = config$data$HK_path, hk_genes = c("Hsp90", "Ppia")) {
  return(
    read_excel(path)
    |> bind(x, if (!is.null(hk_genes)) filter(x, Gene %in% hk_genes) else filter(x, Gene != "Ywhaz")) # Ignoring Ywhaz HK calibration data by default (very few samples)
    |> rename(Cq_hk = Cq)
    |> group_by(Mouse)
    |> summarize(
      Cq_hk_avg = mean(Cq_hk, na.rm = TRUE),
      Experiment = unique(Experiment)
    )
    |> tidyr::extract(Mouse, into = c("Litter", NA, NA), regex = "^(\\w{2})(\\d{1,2})(H|N)$", convert = TRUE, remove = FALSE)
    |> normalize_litter_names()
  )
}

## Compute DCq values from HK Cq values
add_dcq <- function(data, hk_genes = NULL) {
  return(
    load_HK(hk_genes = hk_genes)
    |> right_join(data, by = intersect(colnames(data), c("Experiment", "Mouse", "Litter")))
    |> validate_hk_join()
    |> mutate(DCq = Mean_Cq - Cq_hk_avg)
    |> check_dcq_diff()
  )
}

## Computing fold change values
add_fold_change <- function(data) {
  return(data 
         |> mutate(Condition = factor(Condition, levels = c("N", "H")))
         |> group_by(across(any_of(c("Stage", "Layer", "Gene", "Condition"))))
         |> mutate(avg.DCq = mean(DCq, na.rm = TRUE), avg.N = mean(DCq, na.rm = TRUE))
         |> ungroup()
         |> mutate(avg.N = ifelse(Condition == "N", avg.N, NA))
         |> group_by(across(any_of(c("Stage", "Layer", "Gene"))))
         |> arrange(across(any_of(c("Stage", "Layer", "Gene", "Condition"))))
         |> fill(avg.N, .direction = "down")
         |> ungroup()
         |> mutate(dDCq = avg.DCq - avg.N, two = 2**(-dDCq))
         |> group_by(across(any_of(c("Stage", "Layer", "Gene", "Condition")))) # Experiment
         |> mutate(Fold = mean(two, na.rm = TRUE))
         |> ungroup()
  )
}

## Computing expression type based on fold change sign and its statistical significance
add_expression <- function(data, p.name = "p.val") {
  return(data 
         |> group_by(across(any_of(c("Stage", "Layer", "Gene"))))
         |> mutate(Expression = ifelse(Condition == "H", get_regulation_type(Fold, .data[[p.name]]), NA))
         |> fill(Expression, .direction = "updown")
         |> ungroup()
  )
}

## Adding sex and weight information for each mouse
add_sex_and_weight <- function(data) {
  if (file.exists(config$data$animal_data_path)) {
    return(
      read_xlsx(config$data$animal_data_path)
      |> select(
        Mouse = `ID souris`, Animalerie, Experiment = `Manip HI`,
        Litter = Portée, Stage = Stade, Sex = Sexe, Condition,
        Date_HI = `Date HI`, Weight = `Poids au jour du sacrifice`, 
        Date_Extract = `Date prélèvement`
      )
      |> mutate(
        across(where(is.character), \(c) str_replace_all(c, "\\?", NA_character_) |> utils::type.convert(as.is = TRUE)),
        Condition = factor(Condition, levels = c("N", "H")),
        Stage = factor(Stage, levels = stage_list),
        Experiment = factor(Experiment),
        Sex = toupper(Sex)
      )
      |> filter(Experiment != "Sarah", !is.na(Stage))
      |> select(Mouse, Stage, Condition, Sex, Weight)
      |> right_join(data, by = intersect(colnames(data), c("Mouse", "Stage", "Condition")))
      |> mutate(Sex = factor(Sex))
    )
  } else {
    cat(log.warn("[DATA] No `mouse_info` file found, skipping.\n"))
    return(data)
  }
}

## Adding gene family information
add_gene_family <- function(data) {
  if(file.exists(config$data$gene_families_path)) {
    return(read_excel(path = config$data$gene_families_path)
           |> select(Gene = gene, Gene_Family = family, Gene_Role = role) 
           |> right_join(data, by = "Gene")
    )
  } else {
    cat(log.warn("\n[DATA] No `family` file found, skipping.\n"))
    return(data)
  }
}

## Adding pathway information
add_pathways <- function(data, path = config$data$gene_data_path, sheet, cols = c("Pathway", "Pathway_family", "Gene", "Effect", "Figure")) {
  if(file.exists(path)) {
    return(
      read_excel(path = path, sheet = sheet)
      |> select(any_of(c("Gene", cols)))
      |> mutate(Gene = str_extract(Gene, "\\w+"))
      |> right_join(data, by = "Gene")
    )
  } else {
    cat(log.warn("\n[DATA] No `pathways` file found, skipping.\n"))
    return(data)
  }
}

## Adding layer family information
add_layer_families <- function(data) {
  return(
    bind_rows(
      crossing(Layer_Family = "EGL", Layer = c("EGL", "EGLi", "EGLo")),
      tibble(Layer_Family = "ML", Layer = "ML"),
      crossing(Layer_Family = "PC", Layer = c("PC", "MLPC")),
      tibble(Layer_Family = "IGL", Layer = "IGL"),
      tibble(Layer_Family = "WM", Layer = "WM")
    )
    |> right_join(data, by = "Layer")
  )
}

## Normalizing litter names: some litters had different codes for males and females -> we're setting everything to one code for both
normalize_litter_names <- function(data) {
  return(
    mutate(data, 
           Litter = case_when(
             Litter == "SX" ~ "SO",
             Litter == "SQ" ~ "SY",
             Litter == "EB" ~ "EA",
             TRUE ~ Litter
           )
    )
  )
}

#--------------------------------#
####🔺Data checking functions ####
#--------------------------------#

check_na <- function(data, col = NULL) {
  nrow_na <- data |> filter(is.na({{col}})) |> nrow()
  if (nrow_na > 0) cli::cli_inform(c("\n", "i" = glue("[INFO] There are {nrow_na} rows where {col} is NA."), "\n"))
  return(data)
}


check_dcq_diff <- function(data, threshold = 0.01) {
  if ("DCq_manual" %in% colnames(data)) {
    nrow_diff <- mutate(data, DCq_diff = abs(DCq_manual - DCq)) |> filter(abs(DCq_diff) >= threshold) |> nrow()
    cli::cli_inform(c("\n", "i" = glue("[INFO] There are {nrow_diff} rows where DCq differs from DCq_manual by more than {threshold}."), "\n"))
  }
  return(data)
}


validate_hk_join <- function(data) {
  nrow_unmatched <- filter(data, is.na(Cq_hk_avg)) |> nrow()
  if (nrow_unmatched > 0)
    cli::cli_alert_warning(c("\n", "!" = glue("[WARN] There are {nrow_unmatched} data points that could not be matched to any HK gene."), "\n"))
  return(data)
}
