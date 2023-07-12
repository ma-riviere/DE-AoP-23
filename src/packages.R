####╔═══════     ══════╗####
####💠Project Packages💠####
####╚═══════     ══════╝####

project_pkgs <- c(
  ### Base packages
  "renv", 
  "here", 
  "config", 
  "rlang", 
  "fs", 
  "crayon", 
  "usethis", 
  "glue", 
  "magrittr", 
  "pipebind",
  
  ### Data wrangling
  "tibble",
  "janitor",
  "readxl",
  "stringr",
  "purrr",
  "tidyr",
  "dplyr",

  ### Model fitting
  "car",
  "glmmTMB@1.1.5",
  "lme4",
  "afex",
  "optimx",
  
  ### Model analysis
  "broom",
  "insight",
  "datawizard",
  "performance",
  "qqplotr",           # Required by performance
  "correlation",
  "psych",             # For categorical correlations
  "parameters",
  "DHARMa",
  "emmeans",
  
  ### Visualizations
  "ggplot2",
  "ggtext",
  "patchwork",
  "see",
  "ggdist",
  "bayesplot",
  "ggnewscale",
  "RColorBrewer",
  "circlize",

  ### Reporting
  "gt",
  "gtExtras",
  "knitr",
  "rmarkdown",
  
  ### Misc
  "styler", 
  "miniUI", 
  "gtools"
)