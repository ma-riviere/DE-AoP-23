####╔═════    ══════╗####
####💠Project Config💠####
####╚═════    ══════╝####

cli_h2("┗ [CONFIG] Loading external {.file _configs.yml} file")

configs <- load_yml(here("_configs.yml"))

options(
    scipen = 999L,
    digits = 4L,
    mc.cores = max(1L, parallel::detectCores(logical = TRUE)),
    glmmTMB.cores = max(1L, parallel::detectCores(logical = FALSE)),
    na.action = "na.omit",
    contrasts = c("contr.sum", "contr.poly"),
    seed = 256,
    keep.source = TRUE # Keeps function source references, so `get_function_code()` can display them in the qmd files
)

set.seed(getOption("seed"))

#--------------------------#
####🔺Knitr & RMarkdown ####
#--------------------------#

cli_h2("┗ [CONFIG] Setting knitr options")

knitr::opts_chunk$set(
    warning = FALSE,
    message = FALSE
    # , fig.align = "center"
    # , fig.retina = 2
    # , dpi = 300
    # , dev = 'svg'
    # , dev.args = list(bg = "transparent")
)

#------------------------#
####🔺Package options ####
#------------------------#

cli_h2("┗ [CONFIG] Setting packages options")

emmeans::emm_options(
    lmer.df = "kenward-roger",
    opt.digits = 4,
    back.bias.adj = FALSE
)

#----------------#
####🔺Masking ####
#----------------#

get <- base::get
