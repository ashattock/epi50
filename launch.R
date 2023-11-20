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

# Load all required packages and functions
source("dependencies.R")

message("Running EPI50 pipeline")

# Set options (see options.R)
o = set_options(do_step = 4)

# Step 1) Prepare all inputs (only needs to be done once)
run_prepare()  # See prepare.R

# Step 2) Estimate impact for non-modelled pathogens using GBD
run_non_modelled()  # See non_modelled.R

# Step 3) Impute missing countries for VIMC-modelled pathogens
run_impute()  # See impute.R

# Step 4) Fit and select impact-FVP functions
run_impact()  # See impact.R

# Step 5) Apply impact functions to historical coverage
run_history()  # See history.R

# Step 6) Generate uncertainty draws
run_uncertainty()  # See uncertainty.R

# Step 7) Produce results
run_results()  # See results.R

# Finish up
message("* Finished!")

