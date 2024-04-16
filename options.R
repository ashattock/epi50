###########################################################
# OPTIONS
#
# Set key options for all things model related. The output
# of this function, o (a list), lives in the global 
# environment, so can be referenced throughout the pipeline.
#
###########################################################

# ---------------------------------------------------------
# Set model options and assumptions
# ---------------------------------------------------------
set_options = function(run_module = NA) {

  message("* Setting options")

  # Several global R settings to make life easier
  default_R_options()  # See auxiliary.R

  # Initiate options list
  o = list(run_module = run_module)

  # Prepare output directory system
  o = prepare_dirs(o)  # See directories.R
  
  # Convert config yaml files to datatables
  o = parse_general_options(o)
  
  # ---- Parallelisation settings ----
  
  # Detect number of cores available to this user
  o$n_cores = detectCores()
  
  # Use multiple cores to speed up several processes
  #
  # NOTE: Not available on windows - flags ignored
  o$parallel = list(
    interp  = TRUE,
    impute  = TRUE, 
    impact  = FALSE,  # NOTE: Having issues with shared memory
    history = TRUE)
  
  # Apply parallelisation depending on operating system
  o = set_parallel(o)
  
  # ---- Loading flags ----
  
  # Flag to force (re)download Gapminder data from github
  o$force_download_gapminder = FALSE
  
  # ---- Simulation flags ----
  
  # Directly simulate Dynamice model
  o$simulate_dynamice = FALSE

  # ---- Plotting flags ----

  # Turn figures on or off
  o$plot_inputs     = TRUE
  o$plot_static     = TRUE
  o$plot_imputation = TRUE
  o$plot_impact     = TRUE
  o$plot_history    = TRUE
  
  # Main results table
  o$results_table = TRUE

  # ---- Plotting settings ----

  # Saved figure size
  o$save_width  = 14
  o$save_height = 10

  # Units of figures sizes
  o$save_units = "in"

  # Plotting resolution (in dpi)
  o$save_resolution = 300

  # Image format for saving figure
  #
  # NOTE: Use a character vector to save with multiple formats at once
  o$figure_format     = "png"
  o$manuscript_format = "svg"  # Lancet requires: "eps" or "svg"

  return(o)
}

# ---------------------------------------------------------
# Apply (or trivalise) parallelisation depending on operating system
# ---------------------------------------------------------
set_parallel = function(o) {
  
  # Auto turn off parallelisation on Windows
  if (.Platform$OS.type == "windows")
    o$parallel[] = FALSE
  
  return(o)
}

# ---------------------------------------------------------
# Load and parse general options from config/general yaml
# ---------------------------------------------------------
parse_general_options = function(o) {
  
  # Function to parse general inputs
  parse_fn = function(exp) {
    
    # Attempt to evaluate expression
    value = tryCatch(
      eval_str(exp), 
      error = function(e) exp)
    
    # Do not interpret strings as functions
    if (is.function(value)) 
      value = exp
    
    return(value)
  }
  
  # Load general options
  general = o$pth$config %>%
    paste0("general.yaml") %>%
    read_yaml() %>%
    lapply(parse_fn)
  
  # Append to global options list
  o = c(o, general)
  
  return(o)
}

