###########################################################
# SET DIRECTORIES
#
# Set and get directories in one place in the name of consistency
# and ease. Creates any directories that do not currently exist.
#
# OUTPUTS:
#	- A list of relevant directories (within o$pth) which can
#   be referenced throughout in the pipeline
#
###########################################################

# ---------------------------------------------------------
# Define paths for project inputs and outputs
# ---------------------------------------------------------
set_dirs = function(o) {
  
  # Initiate file path lists
  pth = out = list()
  
  # We've already moved to main directory
  pth$main = getwd()
  
  # Path to cluster log files and data cache
  out$code  = file.path(pth$main, "R")
  out$cache = file.path(pth$main, "cache")
  
  # Path for temporary files (needed for uploading to database)
  out$temp = getwd() # file.path(pth$main, "temp")
  
  # ---- Input and configuration files ----
  
  # Parent path of all input files
  pth$config = file.path(pth$main, "config")
  
  # Paths to specific configuration files
  pth$credentials = file.path(pth$config, "gcs_credentials.json")
  pth$parameters  = file.path(pth$config, "parameters.json")
  
  # ---- Output directories ----
  
  # Parent path of all output files
  pth_output = file.path(pth$main, "output")
  
  # Path to test run files
  out$testing = file.path(pth_output, "0_testing")
  
  # Path to relative risk and impact factor files
  out$relative_risk  = file.path(pth_output, "1_relative_risk",  o$analysis_name)
  out$impact_factors = file.path(pth_output, "2_impact_factors", o$analysis_name)
  
  # Path to relative risk and impact factor files
  out$uncertainty = file.path(pth_output, "3_uncertainty", o$analysis_name)
  
  # Path to figures and other output resusts
  out$results     = file.path(pth_output, "4_results", o$analysis_name)
  out$figures     = file.path(pth_output, "5_figures", o$analysis_name)
  out$diagnostics = file.path(out$figures, "diagnostics")
  
  # ---- Create directory structure ----
  
  # Make all output directories
  make_out_dirs(out)
  
  # Append paths to o list
  o = append_dirs(o, pth, out)
  
  return(o)
}

# ---------------------------------------------------------
# Make all output directories if they do not already exist
# ---------------------------------------------------------
make_out_dirs = function(out) {
  
  # Extract all path names in list
  pth_names = names(out)
  
  # Loop through these path names
  for (pth_name in pth_names) {
    this_pth = out[[pth_name]]
    
    # If it does not already exist, create it
    if (!dir.exists(this_pth) & !grepl("\\*", this_pth))
      dir.create(this_pth, recursive = TRUE)
  }
}

# ---------------------------------------------------------
# Concatenate separators and append directories to o list
# ---------------------------------------------------------
append_dirs = function(o, pth, out) {
  
  # Extract all path names in list
  pth_names = names(out)
  
  # Loop through these path names
  for (pth_name in pth_names) {
    this_pth = out[[pth_name]]
    
    # We use * to denote a partial path
    if (grepl("\\*", this_pth)) {
      out[[pth_name]] = substr(this_pth, 1, nchar(this_pth) - 1)
      
    } else {  # Otherwise add a file separator to end of output paths
      out[[pth_name]] = paste0(this_pth, .Platform$file.sep)
    }
  }
  
  # Concatenate lists
  o$pth = c(pth, out)
  
  return(o)
}

