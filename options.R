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
# ---------------------------------------------------------
set_options = function(do_step = NA) {

  message("* Setting options")

  # Several global R settings to make life easier
  default_R_options()  # See auxiliary.R

  # Initiate options list
  o = list(do_step = do_step)

  # Prepare output directory system
  o = prepare_dirs(o)  # See directories.R
  
  # ---- Time settings ----
  
  # Years to analyse
  o$analysis_years = 1974 : 2024  # Vaccine deployed across these dates
  
  # Year and age ranges stored in coverage database
  o$data_years = 1980 : 2024
  o$data_ages  = 0 : 95
  
  # ---- Data settings ----
  
  # Smooth non-modelled coverage to better align with GBD burden estimates
  #
  # NOTE: Set to NULL to turn smoothing off
  o$gbd_coverage_smoother = "kernel"  # OPTIONS: "kernel" or "spline"
  
  # Bound coverage values below x%
  o$max_coverage = 0.99

  # ---- Non-linear impact assumptions ----
  
  # Use multiple cores to speed up impact function fitting
  #
  # NOTE: Set to 1 to turn off (and fit sequentially)
  o$n_cores = 1 # 28  # Set to NA to autodetect cores
  
  # Multiply impact when fitting for more consistent FVP-impact scales
  o$impact_scaler = 1000

  o$eval_x_scale = 2
  
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

  # ---- Plotting flags ----

  # Turn figures on or off
  o$plot_inputs       = TRUE
  o$plot_non_modelled = TRUE
  o$plot_imputation   = TRUE
  o$plot_impact       = TRUE
  o$plot_uncertainty  = FALSE
  o$plot_history      = TRUE  # Primary result

  # ---- Plotting settings ----

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

  return(o)
}

