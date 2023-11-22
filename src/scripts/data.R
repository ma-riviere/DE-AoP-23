####╔═════   ═════╗####
####💠Loading Data💠####
####╚═════   ═════╝####

cli_h2("┗ [SCRIPTS] Loading data ingestion functions")

#-------------------------#
####🔺Helper functions ####
#-------------------------#

## Turn a variable into a factor, ordered based on the level order defined within the relevant sheet of the provided data dictionary
to_factor <- \(var, dict_path) factor(var, read_excel(dict_path, sheet = deparse(substitute(var)))$Name)

## Helper functions to associate a regulation status (down or up-regulated) to each gene

regulation_type <- list(
  NOT_REG = "Not Regulated", 
  MAYBE_UPREG = "Maybe Upregulated", 
  UPREG = "Upregulated",
  MAYBE_DOWNREG = "Maybe Downregulated",
  DOWNREG = "Downregulated"
)

get_regulation_type <- function(fold, p_value) {
  case_when(
    p_value <= .05 & fold < 1 ~ regulation_type$DOWNREG,
    p_value <= .05 & fold > 1 ~ regulation_type$UPREG,
    is.na(p_value) | is.na(fold) ~ NA_character_,
    .default = regulation_type$NOT_REG
  )
}

## Computing fold change between N & IH conditions, for a given Gene, in a given Layer, at a given Stage
# add_fold_change <- function(data) {
#   cli_alert_info("[DATA] Computing Fold change values")
#   return(
#     data
#     |> select(Stage, Layer, Gene, Condition, DCq)
#     |> pivot_wider(id_cols = c(Stage, Layer, Gene), names_from = Condition, values_from = DCq, values_fn = \(x) mean(x, na.rm = TRUE))
#     |> summarize(Fold = 2**(-1 * (IH - N)), .by = c(Stage, Layer, Gene))
#     |> left_join(data, y = _, join_by(Stage, Layer, Gene))
#     |> mutate(Fold = if_else(Condition == "N", 1, Fold))
#   )
# }

#-------------------------------#
####🔺Data loading functions ####
#-------------------------------#

## Loading the supplementary data (animal_data, gene_data, & layer_families)
load_supplementary_data <- function() {
  
  res <- list()
  
  res$animal_data <- read_excel(configs$data$animal_data) |> 
    mutate(across(c(Stage, Condition), \(x) to_factor(x, configs$data$PCR$data_dict)))
  
  res$gene_data <- map(
    set_names(excel_sheets(configs$data$gene_data)), 
    \(sheet) read_excel(configs$data$gene_data, sheet)
  )
  
  res$layer_families <- (
    bind_rows(
      crossing(Layer_family = "EGL", Layer = c("EGL", "EGLi", "EGLo")),
      tibble(Layer_family = "ML", Layer = "ML"),
      crossing(Layer_family = "PC", Layer = c("PC", "MLPC")),
      tibble(Layer_family = "IGL", Layer = "IGL"),
      tibble(Layer_family = "WM", Layer = "WM")
    )
  )
  
  return(res)
}

#### ┗ PCR ------

## Loading the PCR data
# @target: which data to load (i.e. OS or ND)
# @reprocess: if TRUE, re-process the raw data, else simply load it from the processed data
# @refit: if TRUE, refit the models and extract their predictions (careful, refitting takes upward to 10 minutes)
# @max_cq_clean: maximum Cq value allowed for the clean data
# @model: the model to fit to the data. Only required if reprocess = TRUE
load_pcr_data <- function(target, reprocess = FALSE, max_cq_clean = 33, refit = FALSE, model) {
  
  possible_targets <- configs$data$PCR |> 
    names() |> 
    keep(\(x) str_detect(x, "_raw$|_processed$")) |> 
    str_split_i("_", 1) |> 
    unique()
  
  if (missing(target) || target %ni% possible_targets) 
    cli_abort("[DATA] Incorrect {.var target} provided: possible values are {.pkg {possible_targets}}")
  
  if (reprocess) {
    
    path <- configs$data$PCR[[str_glue("{target}_raw")]]
    
    res <- list()
    
    ## Raw data
    res$raw <- (
      map(set_names(excel_sheets(path)), \(x) read_excel(path, sheet = x)) 
      |> list_rbind(names_to = "Stage")
      |> mutate(across(c(Condition, Stage, Layer), \(x) to_factor(x, configs$data$PCR$data_dict)))
      |> arrange(Stage, Layer, Gene, Mouse)
      |> select(Mouse, any_of("Experiment"), Stage, Layer, Gene, Condition, Mean_Cq, DCq)
      |> mutate(ID = row_number(), .before = 1)
    )
    
    ## Processed data (i.e. without unusable data points/groups)
    res$clean <- (
      res$raw
      ## Filtering useless data points
      |> drop_na(Mean_Cq, DCq)
      |> filter(Mean_Cq <= max_cq_clean)
      |> filter(n() >= 3, .by = c(Stage, Layer, Gene, Condition)) # Checking that there are at least 3 values per Condition
      |> filter(n_distinct(Condition) == 2, .by = c(Stage, Layer, Gene)) # Checking that each Gene has been measured in both Conditions
      ## Arranging the data
      |> select(ID, Mouse, any_of("Experiment"), Stage, Layer, Gene, Condition, Mean_Cq, DCq)
      |> mutate(across(c(Condition, Stage, Layer), \(x) to_factor(x, configs$data$PCR$data_dict)))
      |> arrange(Stage, Layer, Gene, Mouse)
    )
    
    if (refit) {
      
      if (missing(model)) 
        cli_abort("[DATA] Please provide a model to fit the data on.")
      
      compute_fold_change <- function(mod) {
        return(
          get_data(mod) 
          |> select(Condition, DCq) 
          |> pivot_wider(names_from = Condition, values_from = DCq, values_fn = \(x) mean(x, na.rm = TRUE)) 
          |> summarize(Fold = 2**(-1 * (IH - N))) 
          |> pull(Fold) 
        )
      }
      
      ## Fitting the provided model to each Gene, for each Layer and Stage
      res$models <- (
        res$clean
        |> group_split(Stage, Layer, Gene)
        |> map_dfr(
          \(d) suppressMessages({summarize(d, Mod = list(model(pick(everything()))), .by = c(Stage, Layer, Gene))}), 
          .progress = "Fitting models:"
        )
        |> filter(!has_na_coefs(Mod))
        |> mutate(Fold = map_dbl(Mod, compute_fold_change))
        |> select(Stage, Layer, Gene, Fold, Mod)
      )
      
      get_emmeans_data <- function(mod) {
        return(
          emmeans(mod, specs = "Condition", type = "response")
          |> contrast(method = "pairwise", adjust = "none", infer = TRUE)
          |> as.data.frame()
          |> pivot_wider(names_from = contrast, values_from = estimate)
          |> select(last_col(), LCB = lower.CL, UCB = upper.CL, p.value)
          |> mutate(across(where(is.character), \(x) na_if(x, "NaN")))
        )
      }
      
      ## Extracting model predictions
      res$predictions <- (
        res$models
        |> group_split(Stage, Layer, Gene)
        |> map_dfr(\(d) mutate(d, get_emmeans_data(Mod[[1]])), .progress = "Extracting model predictions:")
        |> filter(!is.na(p.value))
        |> mutate(Expression = get_regulation_type(Fold, p.value))
        |> select(Stage, Layer, Gene, Fold, Expression, matches("-|/"), LCB, UCB, p.value)
      )
    }
    
    return(res)
    
  } else {
    readRDS(configs$data$PCR[[str_glue("{target}_processed")]])
  }
}

#### ┗ IHC ------

## Loading and shaping Calbindin IHC data
load_calb_data <- function(path = configs$data$IHC$calb) {
  
  res <- list()
  
  res$raw <- (
    map(set_names(excel_sheets(path)), \(sheet) read_excel(path, sheet))
    |> list_rbind(names_to = "Stage")
    |> arrange(Sample)
    |> mutate(
      across(c(Stage, Condition), \(x) to_factor(x, configs$data$IHC$data_dict)),
      ID = row_number()
    )
    |> select(ID, Sample, Mouse, Stage, Condition, N_CC, Vol_PC)
  )
  
  res$clean <- (
    res$raw
    |> mutate(Vol_PC_per_cell = Vol_PC / N_CC)
    |> select(ID, Sample, Mouse, Stage, Condition, N_CC, Vol_PC_per_cell)
  )
  
  return(res)
}

## Loading and shaping the Caspase IHC data
load_casp_data <- function(path = configs$data$IHC$casp) {
  
  res <- list()
  
  res$raw <- (
    map(set_names(excel_sheets(path)), \(sheet) read_excel(path, sheet))
    |> list_rbind(names_to = "Stage")
    |> arrange(Sample)
    |> mutate(
      across(c(Stage, Condition, Z), \(x) to_factor(x, configs$data$IHC$data_dict)),
      ID = row_number()
    )
    |> select(ID, Sample, Mouse, Stage, Condition, Z, matches("_Tot"), matches("_EGL"), matches("_ML_PCL"), matches("_IGL"), matches("_WM"))
  )
  
  res$clean <- (
    res$raw
    |> mutate(Dens_Tot = N_Tot / A_Tot)
    |> select(ID, Sample, Mouse, Stage, Condition, Z, matches("_Tot"), matches("_EGL"), matches("_ML_PCL"), matches("_IGL"), matches("_WM"))
  )
  
  return(res)
}

#----------------------------------#
####🔺Data validation functions ####
#----------------------------------#

check_na <- function(data, col = NULL) {
  nrow_na <- data |> filter(is.na({{col}})) |> nrow()
  if (nrow_na > 0) cli::cli_inform(c("\n", "i" = str_glue("[INFO] There are {nrow_na} rows where {col} is NA."), "\n"))
  return(data)
}
