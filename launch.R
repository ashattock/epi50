###########################################################
# LAUNCH
#
# Main launch function for WHO EPI50 analysis.
#
# Project lead: A.J.Shattock
###########################################################

# Set working directory to sourced file
if (interactive()) setwd(getSrcDirectory(function() {}))

# Load all required packages and functions
source("dependencies.R")

message("Running EPI50 pipeline")

# Set options (see options.R)
o = set_options(do_step = 0)

# Step 0) Prepare all inputs (only needs to be done once)
run_prepare()  # See prepare.R

# Step 1) Calculate and impute relative risk
run_relative_risk()  # See relative_risk.R

# Step 2) Calculate impact factors
run_impact_factors()  # See impact_factors.R
# run_impact_nonlinear()  # See impact_nonlinear.R

# Step 3) Calculate DALYs
run_dalys()  # See dalys.R

# Step 4) Generate uncertainty draws
run_uncertainty()  # See uncertainty.R

# Step 5) Produce results
run_results()  # See results.R

# Finish up
message("* Finished!")

