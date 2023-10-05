
# Mainly doing this for %>% pipe operator
#' @rawNamespace import(dplyr, except = c(last, first, between))
NULL

# Assignment pipe operator
#' @importFrom magrittr "%<>%"
NULL

# All data.table functionality
#
# NOTE: data.table is generally careful to minimize the scope for namespace
# conflicts (i.e., functions with the same name as in other packages);
# a more conservative approach using @importFrom should be careful to
# import any needed data.table special symbols as well, e.g., if you
# run DT[ , .N, by='grp'] in your package, you'll need to add
# @importFrom data.table .N to prevent the NOTE from R CMD check.
# See ?data.table::`special-symbols` for the list of such symbols
# data.table defines; see the 'Importing data.table' vignette for more
# advice (vignette('datatable-importing', 'data.table')).
#
#' @import data.table
NULL

# Many helpful tidy functions including pivot_longer/wider, expand_grid, and unite
#' @import tidyr
NULL

# Many helpful string manipulation functions
#' @import stringr
NULL

# Helpful date functions
#
# NOTE: Has numerous conflicts with datatable
#' @importFrom lubridate parse_date_time floor_date ymd hours
NULL

# For shiny plots
#' @import highcharter
#' @import shiny
NULL

# For static plots
#' @import ggplot2
#' @import scales
#' @import pals
NULL

# Based on devtools::check() recommendations
#' @importFrom utils data head tail write.csv
NULL
