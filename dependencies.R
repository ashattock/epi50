###########################################################
# DEPENDENCIES
#
# Deal with all package dependencies in one place.
#
###########################################################

# ---- R version check ----

# # R versions for which this project has been tested and is stable
# stable_versions = "4.3.0"
# 
# # R versions for which this project is stable (as a string)
# stable_str = paste(stable_versions, collapse = ", ")
# 
# # Get details of R version currently running
# version_info = R.Version()
# 
# # Construct version number from list details
# version_num = paste0(version_info$major, ".",  version_info$minor)
# 
# # Throw an error if this R version is unsuitable
# if (!version_num %in% stable_versions)
#   stop("This software is stable with R version(s): ", stable_str,
#        " (currently running ", version_num, ")")

# Clear global environment
rm(list = ls())

# ---- Source files ----

# Scripts that should bot be sourced
no_src = c("launch.R", "dependencies.R")

# All R files, and those to source
all_files = list.files(pattern = ".+\\.R$")
src_files = setdiff(all_files, no_src)

# Source each of these files
for (file in src_files)
  source(file)

# ---- Define packages ----

# Complete list of all R packages required for this project
packages = c("tidyverse",      # Includes ggplot2, dplyr, tidyr (www.tidyverse.org/packages/)
             "data.table",     # Next generation dataframes
             # "gsubfn",         # Output multiple variables from functions
             "magrittr",       # Additional pipe operators, such as %<>%
             # "rlist",          # List-related operation functions
             "wrapr",          # Convenience functions (eg qc)
             "stats",          # Statistical calculations and random number generation
             "stats4",         # MLE algorithm
             "matrixStats",    # Matrix row and column operations
             "tgp",            # Latin hypercube sampler
             # "pspline",        # Spline fitting functions
             "wpp2019",        # Demogrphic data
             "yaml",           # Data loading functionality
             # "xlsx",           # Data loading functionality (or use readxl?)
             "lubridate",      # Data formatting functionality
             "naniar",         # Data formatting functionality
             "ggpubr",         # Plotting functionality
             "scales",         # Plotting functionality
             "pals",           # Colour palettes
             "colorspace")     # Colour palettes

# Other packages that are likely needed...
# DBI
# DT
# RSQLite
# rnaturalearthdata
# bigrquery
# dbplyr
# knitr
# rmarkdown
# testthat

# ---- Install and/or load R packages with pacman ----

message("* Installing required packages")

# Check whether pacman itself has been installed
pacman_installed = "pacman" %in% rownames(installed.packages())

# If not, install it
if (!pacman_installed) 
  install.packages("pacman")

# Load pacman
library(pacman) 

# Load all required packages, installing them if required
pacman::p_load(char = packages)

# ---- Redefine or unmask particular functions ----

# Unmask certain commonly used functions
select  = dplyr::select
filter  = dplyr::filter
rename  = dplyr::rename
recode  = dplyr::recode
count   = dplyr::count
union   = dplyr::union
predict = stats::predict

# ---- Clean up ----

# Tidy up console after package loading
if (interactive()) clf()  # Close figures
if (interactive()) clc()  # Clear console

