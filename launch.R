###########################################################
# LAUNCH
#
# Main launch function for WHO EPI50 analysis.
#
###########################################################

# Set working directory to sourced file
if (interactive()) setwd(getSrcDirectory(function() {}))

# Load all required packages and functions
source("dependencies.R")

message("Running EPI50 pipeline")

# Set options (see options.R)
o = set_options(do_step = 0)

# Other possible pre-steps: 
#  - Generate database
#  - Population projection 

# Additional functionality:
# prep_sia()
# explore_coverage()
# explore_dalys()
# explore_nonlinear()

# Step 0) Prepare input (only needs to be done once)
run_prepare()  # See prepare.R

# Step 1) Calculate and impute relative risk
run_relative_risk()  # See relative_risk.R

# Step 2) Impact factors
run_impact_factors()  # See impact_factors.R

# Step 3) Generate uncertainty draws
run_uncertainty()  # See uncertainty.R

# Step 4) Produce results
run_results()  # See results.R

# Finish up
message("* Finished!")