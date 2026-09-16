####╔═══════     ═══════╗####
####💠 Package Manager 💠####
####╚═══════     ═══════╝####

cli_h2("┗ [PACKAGES] Loading external {.file _dependencies.yml} file")

dependencies <- load_yml(here("_dependencies.yml"))

options(
    renv.config.snapshot.inference = FALSE,
    renv.config.snapshot.validate = FALSE,
    renv.config.synchronized.check = FALSE,
    renv.config.updates.parallel = max(1L, parallel::detectCores(logical = TRUE) - 1L),
    Ncpus = max(1L, parallel::detectCores(logical = TRUE) - 1L),
    install.opts = "--no-lock",
    repos = unlist(dependencies$repos),
    pkgType = ifelse(Sys.info()[["sysname"]] == "Windows", "both", "source"),
    verbose = FALSE
)

Sys.setenv(MAKEFLAGS = paste0("-j", getOption("Ncpus")))

#----------------------------#
####🔺Installing packages ####
#----------------------------#

install_project_packages <- function() {
    cli_h2("[PACKAGES] Installing project packages")

    ## To make sure we have the most up-to-date version of the yml
    dependencies <- load_yml(here("_dependencies.yml"))

    install(dependencies$packages, prompt = FALSE)

    cli_h2("[PACKAGES] Snapshotting project packages")

    snapshot(prompt = FALSE, type = "all")

    load_project_packages()
}

#---------------------------#
####🔺Restoring packages ####
#---------------------------#

restore_project_packages <- function() {
    cli_h2("┗ [PACKAGES] Restoring project packages")

    renv::restore(prompt = FALSE)

    load_project_packages()
}

#-------------------------#
####🔺Loading packages ####
#-------------------------#

load_project_packages <- function() {
    cli_h2("┗ [PACKAGES] Loading project packages")

    suppressWarnings({
        dependencies$packages |>
            sapply(\(pkg) load_package_if_installed(get_pkg_name(pkg))) |>
            invisible()
    })
}

get_pkg_name <- function(pkg) {
    pkg_name <- pkg
    if (grepl("/", pkg_name, fixed = TRUE)) {
        pkg_path <- strsplit(pkg_name, "/", fixed = TRUE)[[1]]
        pkg_name <- pkg_path[length(pkg_path)]
    }
    if (grepl("@", pkg_name, fixed = TRUE)) {
        pkg_path <- strsplit(pkg_name, "@", fixed = TRUE)[[1]]
        pkg_name <- pkg_path[1]
    }
    return(pkg_name)
}

get_pkg_version <- function(pkg) {
    if (grepl("@", pkg, fixed = TRUE)) {
        pkg_path <- strsplit(pkg, split = "@", fixed = TRUE)[[1]]
        return(pkg_path[length(pkg_path)])
    }
    return(NA_character_)
}
