####╔═════     ═════╗####
####💠Stats helpers💠####
####╚═════     ═════╝####

cli_h2("┗ [SCRIPTS] Loading stats functions")

#------------------------------------#
####🔺Summarizing data or a model ####
#------------------------------------#

distribution_summary <- function(data, dvs, between = "Condition") {
  data |> select(all_of(between), all_of(dvs)) |> 
    pivot_longer(all_of(dvs), names_to = "DV", values_to = "Value") |> 
    group_by(across(any_of(between)))  |> 
    group_map(
      \(d, g) datawizard::describe_distribution(group_by(d, DV), verbose = FALSE) |>
        mutate(
          Variance = SD^2,
          CoV = ifelse(SD / Mean > 1e4, NA_real_, SD / Mean),
          Variable = str_remove(.group, fixed("DV="))
        ) |> 
        add_column(g, .after = 1) |> 
        select("Variable", all_of(between), "Mean", "SD", "Variance", "CoV", "IQR", "Min", "Max", "Skewness", "Kurtosis", "n")
    ) |> 
    reduce(full_join, by = c("Variable", between, "Mean", "SD", "Variance", "CoV", "IQR", "Min", "Max", "Skewness", "Kurtosis", "n")) |> 
    arrange(Variable, across(any_of(between)))
}

get_model_based_outliers <- function(data, mod, mod_dharma, responses) {
  
  outliers <- get_data(mod) |> 
    rownames_to_column("ID") |> 
    filter(ID %in% DHARMa::outliers(mod_dharma)) |> 
    utils::type.convert(as.is = TRUE) 
  
  if (nrow(outliers) > 0) outliers <- semi_join(data, y = outliers) |> select(-setdiff(responses, find_response(mod)))
  
  return(outliers)
}

#--------------------------------------------#
####🔺Extracting information from a model ####
#--------------------------------------------#

# Should we exponentiate the coefficients of a model (based on its link function)
should_exp <- \(mod) insight::get_family(mod)$link %in% c("log", "logit")

## Check which (if any) models have any NA as fixed effect coefficients (which signals that the model fitting failed silently)
has_na_coefs <- function(mods) {
  map_lgl(
    mods,
    \(mod) as.data.frame(summary(mod)$coefficients$cond) |> 
      mutate(across(where(is.character), \(x) na_if(x, "NaN"))) |> 
      lapply(\(x) any(is.na(x))) |> 
      flatten_lgl() |> 
      any()
  )
}

find_formula_formatted <- function(mod) {
  return(
    c(
      find_formula(mod)$conditional,
      map(
        find_formula(mod)$random |> unclass(),
        \(x) paste0("(", reduce(deparse(x), \(y) paste(y)) |> str_replace_all("~", ""), ")")
      )
    ) |> keep(\(z) z != "()") |> paste(collapse = " + ")
  )
}

find_coefficient_count_mod <- function(mod) {
  return(find_coefficient_count(
    get_data(mod),
    find_predictors(mod)$conditional,
    find_interactions(mod)$conditional
  ))
}

find_coefficient_count <- function(data, predictors, interactions = NULL) {
  
  coef_count <- 0
  coef_count <- map(predictors, \(pred) length(unique(data[[pred]])) - 1) |> reduce(\(x, y) sum(x, y, na.rm = TRUE), .init = 0)
  
  if (!is.null(interactions)) {
    coef_count <- coef_count + map(
      str_split(interactions, ":"), 
      \(int_preds) {
        map(int_preds, \(int_pred) length(unique(data[[int_pred]])) - 1) |> reduce(\(x, y) prod(x, y, na.rm = TRUE), .init = 1)
      }
    ) |> reduce(\(x, y) sum(x, y, na.rm = TRUE), .init = 0)
  }
  
  return(coef_count + 1) # Add one for the intercept
}

#---------------------------#
####🔺Variable selection ####
#---------------------------#

## Likelihood Ratio Test helper function: 
### - Attempts to refit the provided model without the predictor of interest (i.e. the reduced model), and without using REML 
### - Runs an LRT between the original model (without REML) and the reduced model
LRT <- function(mod, pred = "Condition", print_eq = FALSE) {
  data <- get_data(mod)
  resp <- find_response(mod)
  link <- link_function(mod)
  preds <- find_predictors(mod)$conditional
  inters <- find_interactions(mod)$conditional
  
  preds_reduced <- str_subset(preds, pred, negate = TRUE)
  inters_reduced <- str_subset(inters, pred, negate = TRUE)
  
  reduced_formula <- str_glue("{resp} ~ . -{paste0(stringr::str_subset(preds, pred))}")
  if (!is.null(inters)) reduced_formula <- str_c(reduced_formula, str_glue("-{paste0(stringr::str_subset(inters, pred))}"))
  
  if (toupper(find_algorithm(mod)$algorithm) == "REML") 
    cat(crayon::yellow(str_glue("\n{crayon::bold('\n[LRT]')} Full model was fit with REML --> Refitting with ML.\n")))
  
  mod_full <- tryCatch(
    update(mod, REML = FALSE),
    
    warning = \(w) {
      cat(crayon::red(str_glue("{crayon::bold('\n[LRT]')} Model convergence warning (full model):\n\n")))
      print(w)
      cat(crayon::yellow(str_glue("{crayon::bold('[LRT]')} Retrying with better starting values ...\n\n")))
      
      beta_start_full <- c(link(mean(insight::get_response(mod))), rep(0, find_coefficient_count(data, preds, inters) - 1))
      
      return(
        tryCatch(
          update(mod, REML = FALSE, start = list(beta = beta_start_full)),
          warning = \(w) {
            update(mod, REML = FALSE, start = list(beta = beta_start_full), control = glmmTMBControl())
          },
          error = \(e) {
            update(mod, REML = FALSE, start = list(beta = beta_start_full), control = glmmTMBControl())
          }
        )
      )
    },
    error = \(e) {
      cat(crayon::red(str_glue("{crayon::bold('\n[LRT]')} Model convergence error (full model):\n\n")))
      print(e)
      cat(crayon::yellow(str_glue("{crayon::bold('[LRT]')} Retrying with better starting values ...\n\n")))
      
      beta_start_full <- c(link(mean(insight::get_response(mod))), rep(0, find_coefficient_count(data, preds, inters) - 1))
      
      return(
        tryCatch(
          update(mod, REML = FALSE, start = list(beta = beta_start_full)),
          warning = \(w) {
            update(mod, REML = FALSE, start = list(beta = beta_start_full), control = glmmTMBControl())
          },
          error = \(e) {
            update(mod, REML = FALSE, start = list(beta = beta_start_full), control = glmmTMBControl())
          }
        )
      )
    }
  )
  
  mod_reduced <- tryCatch(
    update(mod, formula. = reduced_formula, REML = FALSE, start = NULL), 
    
    warning = \(w) {
      cat(crayon::red(str_glue("{crayon::bold('\n[LRT]')} Model convergence warning (reduced model):\n\n")))
      print(w)
      cat(crayon::yellow(str_glue("{crayon::bold('[LRT]')} Retrying with better starting values ...\n\n")))
      
      beta_start_reduced <- c(link(mean(insight::get_response(mod))), rep(0, find_coefficient_count(data, preds_reduced, inters_reduced) - 1))
      
      return(
        tryCatch(
          update(mod, formula. = reduced_formula, REML = FALSE, start = list(beta = beta_start_reduced)),
          warning = \(w) {
            update(mod, formula. = reduced_formula, REML = FALSE, start = list(beta = beta_start_reduced), control = glmmTMBControl())
          },
          error = \(e) {
            update(mod, formula. = reduced_formula, REML = FALSE, start = list(beta = beta_start_reduced), control = glmmTMBControl())
          }
        )
      )
    },
    error = \(e) {
      cat(crayon::red(str_glue("{crayon::bold('\n[LRT]')} Model convergence error (reduced model):\n\n")))
      print(e)
      cat(crayon::yellow(str_glue("{crayon::bold('[LRT]')} Retrying with better starting values ...\n\n")))
      
      beta_start_reduced <- c(link(mean(insight::get_response(mod))), rep(0, find_coefficient_count(data, preds_reduced, inters_reduced) - 1))
      
      return(
        tryCatch(
          update(mod, formula. = reduced_formula, REML = FALSE, start = list(beta = beta_start_reduced)),
          warning = \(w) {
            update(mod, formula. = reduced_formula, REML = FALSE, start = list(beta = beta_start_reduced), control = glmmTMBControl())
          },
          error = \(e) {
            update(mod, formula. = reduced_formula, REML = FALSE, start = list(beta = beta_start_reduced), control = glmmTMBControl())
          }
        )
      )
    }
  )
  
  formula_full <- find_formula_formatted(mod_full)
  formula_reduced <- find_formula_formatted(mod_reduced)
  
  res <- stats::anova(mod_full, mod_reduced, test = "LRT") |> as.data.frame() |> rownames_to_column("Model") |> janitor::clean_names()
  # res <- performance::test_lrt(mod_full, mod_reduced, estimator = "ML") # estimator = "OLS" <=> anova(test = "LRT")
  
  if (print_eq) {
    cat(crayon::blue(str_glue("{crayon::bold('\n[LRT]')} Full formula: {formula_full}\n")))
    cat(crayon::blue(str_glue("{crayon::bold('\n[LRT]')} Reduced formula: {formula_reduced}\n")))
    
    cat(paste0(
      "\n\n$\\mathcal{X}_", str_glue("{res$chi_df[2]}"), "^2", str_glue(" = {round(res$chisq[2], 3)}; "),
      scales::pvalue(res$pr_chisq[2], add_p = T, prefix = c("p < ", "p = ", "p > ")), "$"
    ))
  }
  
  return(res)
}