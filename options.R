###########################################################
# OPTIONS
#
# Set key options for all things model related. The output 
# of this function, o (a list), lives in the global environment, 
# so can be referenced throughout the pipeline.
#
###########################################################

# ---------------------------------------------------------
# Set model options and assumptions
# Called by: launch.R
# ---------------------------------------------------------
set_options = function(do_step = NA, quiet = FALSE) {
  
  if (!quiet) message("* Setting options")
  
  # Reset R's most annoying default options
  options(stringsAsFactors = FALSE, scipen = 999, dplyr.summarise.inform = FALSE)
  
  # Initiate options list
  o = list(do_step = do_step)
  
  # ---- Analysis settings ----
  
  # Name of analysis to run
  o$analysis_name = "v01"
  
  # Create output directory system
  o = set_dirs(o)  # See directories.R
  
  # ---- Non-linear impact assumptions ----
  
  o$per_person = 1
  
  # A very good fit is required to go non-linear
  #
  # NOTE: We also require a better AICc than the simple linear model
  o$r2_threshold0 = 0.98  # Above this guarentees a linear model
  o$r2_threshold1 = 0.95  # Non-linear models must be above this
  
  # Multiply impact when fitting for more consistent FVP-impact scales
  o$impact_scaler = 1000
  
  o$eval_x_scale = 2
  
  # Metric with which to select best fitting model
  o$model_metric = "r2"  # OPTIONS: "aicc" or "r2"
  
  # ---- Time settings ----
  
  # Years to analyse
  o$analysis_years = 1974 : 2024  # Vaccine deployed across these dates
  
  # Year and age ranges stored in coverage database
  o$data_years = 1980 : 2024  # Vaccine effect calculated across these dates
  o$data_ages  = 0 : 95
  
  # ---- Uncertainty settings ----
  
  # Flag for reproducible uncertainty draws (consistent randomly sampled seeds)
  o$uncertainty_reproducible = TRUE  # TODO: Implement this
  
  # Number of draws to sample
  o$n_draws = 200
  
  # Parameter bounds for fitting beta distribution to GBD disease vaccine efficacy
  o$par_lower = log(1)
  o$par_upper = log(10)
  
  # Statistical summary to use for 'best estimate' projection
  #
  # OPTIONS:
  #  "median" := Median of uncertainty simulations (stochastic and parameter uncertainty)
  #    "mean" := Mean of uncertainty simulations (stochastic and parameter uncertainty)
  # o$best_estimate_simulation = "mean"
  
  # Quantiles for credibility intervals
  # o$quantiles = c(0.025, 0.975)
  
  # ---- Results flags ----
  
  # Turn results generation on or off
  o$results_markdown = TRUE  # Full markdown results document
  o$results_upload   = FALSE  # Upload reference results to database
  
  # ---- Plotting flags ----
  
  # Turn figures on or off
  o$plot_diagnostics = TRUE  # All diagnostic figures
  
  # ---- Plotting settings ----
  
  # # Lower bound of age groups for plotting - bounded above by maximum age
  # o$plot_ages = c(0, 18, 60)  # Captures 3 age groups as per ECDC request
  
  # Colour palette for SIA data exploration plots
  o$palette_sia = "pals::kovesi.rainbow"
  
  # Colour packages and palettes (see colour_scheme in auxiliary.R)
  o$palette_disease = "pals::kovesi.rainbow"  # ~15 values needed
  o$palette_country = "pals::kovesi.rainbow"  # ~190 values needed
  o$palette_region  = "brewer::paired"  # 6 values
  o$palette_economy = "brewer::dark2"  # 4 values
  o$palette_gavi    = "brewer::greys"  # 2 values (yes or no)
  
  # # Define some nice properties for baseline metric plots
  # o$baseline_name   = "Baseline scenario"
  # o$baseline_colour = "grey50"  # Light grey
  # 
  # # Grey colour for current date dashed line
  # o$data_colour = "#555555"  # Dark grey
  # o$dash_colour = "#808080"  # Even darker grey
  
  # Font sizes: title, axis, tick, strip, legend, key
  o$font_size = c(34, 28, 16, 24, 20, 18)
  
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
  o$figure_format = "png" # Classic options: "png", "pdf", or "svg"
  
  # Display analysis details
  if (!quiet) message(" - Analysis name: ", o$analysis_name)
  
  return(o)
}

