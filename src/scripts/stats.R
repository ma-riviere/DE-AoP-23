####╔═════     ═════╗####
####💠Stats helpers💠####
####╚═════     ═════╝####

cli_h2("┗ [SCRIPTS] Loading stats functions")

#------------------------------------#
####🔺Summarizing data or a model ####
#------------------------------------#

distribution_summary <- function(data, dvs, between = "Condition") {
    data |>
        select(all_of(between), all_of(dvs)) |>
        pivot_longer(all_of(dvs), names_to = "DV", values_to = "Value") |>
        group_by(across(any_of(between))) |>
        group_map(
            \(d, g) {
                datawizard::describe_distribution(group_by(d, DV), verbose = FALSE) |>
                    mutate(
                        Variance = SD^2,
                        CoV = ifelse(SD / Mean > 1e4, NA_real_, SD / Mean),
                        Variable = str_remove(.group, fixed("DV="))
                    ) |>
                    add_column(g, .after = 1) |>
                    select(
                        "Variable",
                        all_of(between),
                        "Mean",
                        "SD",
                        "Variance",
                        "CoV",
                        "IQR",
                        "Min",
                        "Max",
                        "Skewness",
                        "Kurtosis",
                        "n"
                    )
            }
        ) |>
        reduce(
            full_join,
            by = c(
                "Variable",
                between,
                "Mean",
                "SD",
                "Variance",
                "CoV",
                "IQR",
                "Min",
                "Max",
                "Skewness",
                "Kurtosis",
                "n"
            )
        ) |>
        arrange(Variable, across(any_of(between)))
}

get_model_based_outliers <- function(data, mod, mod_dharma, responses) {
    outliers <- get_data(mod) |>
        rownames_to_column("ID") |>
        filter(ID %in% DHARMa::outliers(mod_dharma)) |>
        utils::type.convert(as.is = TRUE)

    if (nrow(outliers) > 0) {
        outliers <- semi_join(data, y = outliers) |> select(-setdiff(responses, find_response(mod)))
    }

    return(outliers)
}

#--------------------------------------------#
####🔺Extracting information from a model ####
#--------------------------------------------#

# Should we exponentiate the coefficients of a model (based on its link function)
should_exp <- \(mod) insight::get_family(mod)$link %in% c("log", "logit")

## Recompute Nakagawa's R2 for a Generalized Poisson model, and overwrite the value that
## `parameters(include_info = TRUE)` shows in its table footer.
## glmmTMB's genpois variance is mu * phi^2, where `sigma()` returns the index of dispersion phi^2.
## insight's `.variance_distributional()` groups `genpois` with the negative-binomial families and
## uses (1/mu + 1/phi) for the observation-level variance instead. On an under-dispersed fit that
## inflates it ~24x and collapses the R2 (0.089 -> 0.004 for the Calbindin N_CC model).
fix_genpois_r2 <- function(params, mod) {
    if (!identical(insight::get_family(mod)$family, "genpois")) {
        cli::cli_abort("{.fn fix_genpois_r2} only applies to {.val genpois} models.")
    }

    variances <- insight::get_variance(mod)
    # Nakagawa's lognormal approximation, log(1 + V(mu) / mu^2), with V(mu) / mu^2 = phi^2 / mu
    var_distribution <- mean(log1p(sigma(mod) / fitted(mod)))
    total_variance <- variances$var.fixed + variances$var.random + var_distribution

    r2 <- attr(params, "r2")
    r2$R2_conditional[[1]] <- (variances$var.fixed + variances$var.random) / total_variance
    r2$R2_marginal[[1]] <- variances$var.fixed / total_variance
    attr(params, "r2") <- r2

    return(params)
}

## Check which (if any) models have any NA as fixed effect coefficients (which signals that the model fitting failed silently)
has_na_coefs <- function(mods) {
    map_lgl(
        mods,
        \(mod) {
            as.data.frame(summary(mod)$coefficients$cond) |>
                mutate(across(where(is.character), \(x) na_if(x, "NaN"))) |>
                lapply(\(x) anyNA(x)) |>
                flatten_lgl() |>
                any()
        }
    )
}
