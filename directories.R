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
  
  # Parent path of all external model output files
  pth$extern = file.path(pth$code, "extern")
  
  # Path to log files
  pth$log = file.path(pth$code, "log")
  
  # ---- Output directories ----
  
  # Parent path of all output files
  pth$output = file.path(pth$code, "output")
  
  # Path to cached data tables
  pth$tables = file.path(pth$output, "0_tables")
  
  # Path to non-modelled assumptions and results
  pth$non_modelled   = file.path(pth$output, "1_non_modelled")
  pth$non_modelled_d = file.path(pth$non_modelled, "disease")
  pth$non_modelled_v = file.path(pth$non_modelled, "vaccine")
  pth$non_modelled_t = file.path(pth$non_modelled, "vaccine_type")
  pth$non_modelled_w = file.path(pth$non_modelled, "weights")
  
  # Path to imputatation and impact calculation files
  pth$impute      = file.path(pth$output, "2_impute")
  pth$impact      = file.path(pth$output, "3_impact")
  pth$impact_runs = file.path(pth$impact, "runs")
  
  # Path to relative risk and impact factor files
  # pth$uncertainty = file.path(pth$output, "3_uncertainty")
  
  # Path to figures and other output resusts
  pth$results = file.path(pth$output, "4_results")
  pth$figures = file.path(pth$output, "5_figures")
  
  # Append paths to o list
  o = set_dirs(o, pth)
  
  # ---- Data URLs ----
  
  # stem = "https://ghdx.healthdata.org/sites/default/files/record-attached-files/"
  # 
  # o$data = list(
  #   sdi = c("IHME_GBD_2019_SDI_1990_2019_Y2020M10D15.CSV", 
  #           "IHME_GBD_2019_SDI_1970_1989_Y2020M10D15.CSV"), 
  #   haqi = "IHME_GBD_2019_HAQ_1990_2019_DATA.zip")
    
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

