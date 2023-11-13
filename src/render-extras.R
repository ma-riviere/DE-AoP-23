####╔═════      ═════╗####
####💠 Extra-render 💠####
####╚═════      ═════╝####

suppressPackageStartupMessages({library(here)})
suppressPackageStartupMessages({library(fs)})

#----------------------------------#
####🔺Generating data README.md ####
#----------------------------------#

knitr::knit(
  text = "# Contents

```{r echo = FALSE}
read.csv(here::here('data', '_data_file_desc.csv')) |> insight::print_md()
```
",
output = here("data", "README.md")
)

#-----------------------------------#
####🔺Generating data .zip files ####
#-----------------------------------#

configs <- yaml::read_yaml(here("_configs.yml"), eval.expr = TRUE)

output_dir <- Sys.getenv("QUARTO_PROJECT_OUTPUT_DIR")

dir_create(here(output_dir, "res", "data"))

zip(
  here("res", "data", "supplementary.zip"), 
  configs$data |> purrr::keep_at(\(x) stringr::str_detect(x, "PCR|IHC", TRUE)) |> unlist(), 
  extras = '-j'
)
file_copy(here("res", "data", "supplementary.zip"), here(output_dir, "res", "data", "supplementary.zip"), overwrite = TRUE)

zip(here("res", "data", "IHC.zip"), unlist(configs$data$IHC), extras = '-j')
file_copy(here("res", "data", "IHC.zip"), here(output_dir, "res", "data", "IHC.zip"), overwrite = TRUE)

zip(here("res", "data", "PCR.zip"), unlist(configs$data$PCR), extras = '-j')
file_copy(here("res", "data", "PCR.zip"), here(output_dir, "res", "data", "PCR.zip"), overwrite = TRUE)
