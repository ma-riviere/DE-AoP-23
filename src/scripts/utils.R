####╔═══    ══╗####
####💠 Utils 💠####
####╚═══    ══╝####

cli_h2("┗ [SCRIPTS] Loading helper functions")

#--------------#
####🔺Pipes ####
#--------------#

"%ni%" <- Negate("%in%")

"%s+%" <- \(lhs, rhs) paste0(lhs, rhs)

"%ne%" <- \(lhs, rhs) {
    if (is.null(lhs) || rlang::is_empty(lhs) || (length(lhs) == 1 && !nzchar(lhs))) return(rhs) else return(lhs)
}


#-------------------------#
####🔺Formatting utils ####
#-------------------------#

label_pval <- function(p) {
    str_c(
        scales::label_pvalue()(p) |> str_remove(">") |> str_trim(),
        gtools::stars.pval(p) |> str_replace(fixed("."), ""),
        sep = " "
    )
}

get_response_name <- function(var, from, col = "Label") {
    res <- read_xlsx(configs$data[[from]]$data_dict, sheet = 1) |> filter(Name == var) |> pull(col)

    return(res %ne% var)
}

get_var_level_name <- function(var, level, from, col = "Description") {
    res <- read_excel(configs$data[[from]]$data_dict, sheet = var) |> filter(Name == level) |> pull(col)

    return(res %ne% level)
}

get_model_family <- function(mod) {
    family <- get_family(mod)$family |> str_to_sentence()
    link <- get_family(mod)$link

    model_tag <- str_glue("{family} ('{link}')")

    cov_struct <- str_match(get_call(mod)$formula |> toString(), "\\s(\\w{2,3})\\(.*\\)")[[2]]
    if (!is.null(cov_struct) && !is.na(cov_struct) && nzchar(cov_struct)) {
        model_tag <- str_glue("{model_tag} + {toupper(cov_struct)}")
    }

    return(model_tag)
}

get_model_tag <- function(mod) {
    resp <- insight::find_response(mod)
    return(str_glue("{resp} - {get_model_family(mod)}"))
}

### From: https://michaelbarrowman.co.uk/post/getting-a-variable-name-in-a-pipeline/
get_var_name <- function(x) {
    lhs <- get_lhs()
    if (is.null(lhs)) {
        lhs <- rlang::ensym(x)
    }
    return(rlang::as_name(lhs))
}

get_lhs <- function() {
    calls <- sys.calls()

    #pull out the function or operator (e.g. the `%>%`)
    call_firsts <- lapply(calls, `[[`, 1)

    #check which ones are equal to the pipe
    pipe_calls <- vapply(call_firsts, identical, logical(1), quote(`%>%`))

    #if we have no pipes, then get_lhs() was called incorrectly
    if (!any(pipe_calls)) {
        NULL
    } else {
        #Get the most recent pipe, lowest on the
        pipe_calls <- which(pipe_calls)
        pipe_calls <- pipe_calls[length(pipe_calls)]

        #Get the second element of the pipe call
        this_call <- calls[[c(pipe_calls, 2)]]

        #We need to dig down into the call to find the original
        while (is.call(this_call) && identical(this_call[[1]], quote(`%>%`))) {
            this_call <- this_call[[2]]
        }
        this_call
    }
}

#--------------------#
####🔺Stats utils ####
#--------------------#

label_encoding <- function(var) {
    vals <- unique(var)
    car::Recode(
        var,
        glue::glue_collapse(str_glue("'{vals}' = {as.vector(seq.int(1, length(vals)))}"), sep = "; ") |> as_string()
    ) |>
        as.character() |>
        as.numeric()
}

#-------------#
####🔺Misc ####
#-------------#

## Get element by name from a list:
rmatch <- function(x, name) {
    pos <- match(name, names(x))
    if (!is.na(pos)) {
        return(x[[pos]])
    }
    for (el in x) {
        if (inherits(el, "list")) {
            out <- Recall(el, name)
            if (!is.null(out)) return(out)
        }
    }
}

# Returns the source lines of a function (comments included), to display it in a qmd chunk through the knitr
# `code` option: `#| code: !expr get_function_code("fn_name")`. Requires `options(keep.source = TRUE)`.
get_function_code <- function(fn_name) {
    fn_source <- attr(get(fn_name), "srcref")

    if (is.null(fn_source)) {
        cli_abort("[UTILS] No source reference for {.fn {fn_name}}: source it with {.code keep.source = TRUE}")
    }

    fn_source <- as.character(fn_source)
    fn_source[1] <- paste(fn_name, "<-", fn_source[1])

    return(fn_source)
}
