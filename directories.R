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
prepare_dirs = function(o) {
  
  # Initiate list with reference to main code directory
  #
  # NOTE: We should already be in this code directory
  pth = list(code = getwd())
  
  # ---- Input and configuration directories ----
  
  # Parent path of all input files
  pth$input  = file.path(pth$code, "input")
  pth$config = file.path(pth$code, "config")
  pth$data   = file.path(pth$code, "data")
  
  # ---- Output directories ----
  
  # Parent path of all output files
  pth_output = file.path(pth$code, "output")
  
  # Path to cached data tables
  pth$tables = file.path(pth_output, "0_tables")
  
  # Path to non-modelled assumptions and results
  pth$non_modelled = file.path(pth_output, "1_non_modelled")
  
  # Path to imputatation and impact calculation files
  pth$impute = file.path(pth_output, "2_impute")
  pth$impact = file.path(pth_output, "3_impact")
  
  # Path to relative risk and impact factor files
  # pth$uncertainty = file.path(pth_output, "3_uncertainty")
  
  # Path to figures and other output resusts
  pth$results = file.path(pth_output, "4_results")
  pth$figures = file.path(pth_output, "5_figures")
  
  # ---- Create directory structure ----
  
  # Append paths to o list
  o = set_dirs(o, pth)
  
  return(o)
}

# ---------------------------------------------------------
# Make all output directories and append to o list
# ---------------------------------------------------------
set_dirs = function(o, pth) {
  
  # Iterate through dirs
  for (pth_name in names(pth)) {
    this_pth = pth[[pth_name]]
    
    # If it does not already exist, create it
    if (!dir.exists(this_pth)) {
      dir.create(this_pth, recursive = TRUE)
    }
    
    # Add a file separator to end of dir path
    pth[[pth_name]] = paste0(this_pth, file_sep())
  }
  
  # Append to o list
  o$pth = pth
  
  return(o)
}

