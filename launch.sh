#!/bin/bash

############################################################
# LAUNCH
#
# Launch EPI50 pipeline from command line.
#
# Command line usage:
#   sh launch.sh
#
############################################################

# Load R
module purge
ml R/4.3.0-foss-2021a

# Call main launch script
Rscript launch.R

