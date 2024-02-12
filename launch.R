###########################################################
# LAUNCH
#
# Main launch function for WHO EPI50 analysis.
#
# Note that this software requires an internet connection.
#
# Author: A.J.Shattock
# Contact: shattocka@who.int, andrewjames.shattock@unibas.ch
###########################################################

# Set working directory to sourced file
if (interactive()) setwd(getSrcDirectory(function() {}))
setwd("C:\\Users\\helen\\Documents\\GitHub\\epi50-vaccine-impact")


# Load all required packages and functions
source("dependencies.R")

message("Running EPI50 pipeline")

# Set options (see options.R)
o = set_options(do_step = c(6))

# Step 1) Prepare all inputs (only needs to be done once)
run_prepare()  # See prepare.R

# Step 2) Interface with external polio and measles models
run_external()  # See external.R

# Step 3) Estimate impact for static models using GBD
run_static()  # See static.R

# Step 4) Impute missing countries for VIMC-modelled pathogens
run_impute()  # See impute.R

# Step 5) Fit and select impact-FVP functions
run_impact()  # See impact.R

# Step 6) Apply impact functions to historical coverage
run_history()  # See history.R

# Step 7) Re-fit time series regression models to infer predictors
run_regression()  # See impute.R (second call)

# Step 8) Generate uncertainty draws
run_uncertainty()  # See uncertainty.R

# Step 9) Produce results
run_results()  # See results.R

# Finish up
message("* Finished!")

