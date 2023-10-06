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
  
  # Initiate list with reference to main code directory
  #
  # NOTE: We should already be in this code directory
  pth = list(code = getwd())
  
  # ---- Input and configuration files ----
  
  # Parent path of all input files
  pth$input  = file.path(pth$code, "input")
  pth$config = file.path(pth$code, "config")
  
  # Path to cached data tables
  pth$cache  = file.path(pth$code, "cache")
  # pth$tables = file.path(pth$code, "tables")
  
  # Paths to specific configuration files
  # pth$parameters  = file.path(pth$config, "parameters.json")
  
  # ---- Output directories ----
  
  # Parent path of all output files
  pth_output = file.path(pth$code, "output")
  
  # Path to test run files
  pth$testing = file.path(pth_output, "0_testing")
  
  # Path to relative risk and impact factor files
  pth$relative_risk  = file.path(pth_output, "1_relative_risk",  o$analysis_name)
  pth$impact_factors = file.path(pth_output, "2_impact_factors", o$analysis_name)
  
  # Path to relative risk and impact factor files
  pth$uncertainty = file.path(pth_output, "3_uncertainty", o$analysis_name)
  
  # Path to figures and other output resusts
  pth$results     = file.path(pth_output, "4_results", o$analysis_name)
  pth$figures     = file.path(pth_output, "5_figures", o$analysis_name)
  pth$diagnostics = file.path(pth$figures, "diagnostics")
  
  # ---- Create directory structure ----
  
  # Make all output directories
  make_out_dirs(pth)
  
  # Append paths to o list
  o = append_dirs(o, pth)
  
  return(o)
}

# ---------------------------------------------------------
# Make all output directories if they do not already exist
# ---------------------------------------------------------
make_out_dirs = function(pth) {
  
  # Extract all path names in list
  pth_names = names(pth)
  
  # Loop through these path names
  for (pth_name in pth_names) {
    this_pth = pth[[pth_name]]
    
    # If it does not already exist, create it
    if (!dir.exists(this_pth) & !grepl("\\*", this_pth))
      dir.create(this_pth, recursive = TRUE)
  }
}

# ---------------------------------------------------------
# Concatenate separators and append directories to o list
# ---------------------------------------------------------
append_dirs = function(o, pth) {
  
  # Extract all path names in list
  pth_names = names(pth)
  
  # Loop through these path names
  for (pth_name in pth_names) {
    this_pth = pth[[pth_name]]
    
    # We use * to denote a partial path
    if (grepl("\\*", this_pth)) {
      pth[[pth_name]] = substr(this_pth, 1, nchar(this_pth) - 1)
      
    } else {  # Otherwise add a file separator to end of output paths
      pth[[pth_name]] = paste0(this_pth, .Platform$file.sep)
    }
  }
  
  # Append to o list
  o$pth = pth
  
  return(o)
}

