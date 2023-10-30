###########################################################
# LAUNCH
#
# Main launch function for WHO EPI50 analysis.
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

# Step 2) Calculate and impute relative risk
run_relative_risk()  # See relative_risk.R

# Step 3) Calculate impact per FVP
run_impact()  # See impact_nonlinear.R

# Step 4) Extrapolate impact for missing countries
run_extrapolate()  # See extrapolation.R

# # Step 4) Calculate DALYs
# run_dalys()  # See dalys.R
# 
# # Step 5) Generate uncertainty draws
# run_uncertainty()  # See uncertainty.R
# 
# # Step 6) Produce results
# run_results()  # See results.R

# Finish up
message("* Finished!")

