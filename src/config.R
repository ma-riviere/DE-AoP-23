####╔═════    ══════╗####
####💠Project Config💠####
####╚═════    ══════╝####

log.main("[CONFIG] Configuring project ...")

config <- config::get(file = "_config.yml")

options(
  scipen = 999L, 
  digits = 4L,
  mc.cores = max(1L, parallel::detectCores(logical = TRUE)),
  glmmTMB.cores = max(1L, parallel::detectCores(logical = FALSE)),
  na.action = "na.omit",
  contrasts = c("contr.sum", "contr.poly"),
  seed = 256
)

set.seed(getOption("seed"))

#--------------------------#
####🔺Knitr & RMarkdown ####
#--------------------------#

knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  fig.align = "center",
  # fig.retina = 2,
  dpi = 250,
  dev = 'svg',
  dev.args = list(bg = "transparent")
)

dpi_save_png <- knitr::opts_chunk$get("dpi")

#------------------------#
####🔺Package options ####
#------------------------#

emmeans::emm_options(
  lmer.df = "kenward-roger",
  opt.digits = 4,
  back.bias.adj = FALSE 
)

#----------------#
####🔺Masking ####
#----------------#

get <- base::get

#------------------#
####🔺Constants ####
#------------------#

alpha <- 0.05
trend <- 0.1
threshold.reg <- 1