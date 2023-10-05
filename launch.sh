#!/bin/bash

############################################################
# LAUNCH
#
# Launch IA2030 pipeline from command line.
#
# Command line usage:
#   sh launch.sh
#
############################################################

# Load R
module purge
ml R/4.1.0-foss-2018b

# Call main launch script
Rscript launch.R

