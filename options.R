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
  o$years = 1974 : 2024
  
  # Ages modelled
  o$ages = 0 : 95
  
  # ---- Coverage settings ----
  
  # Bound coverage values below x%
  o$max_coverage = 0.99
  
  # Year that high-income countries switch to acellular pertussis vaccine
  o$wholecell_acellular_switch = 1995
  
  # Smooth non-modelled coverage to better align with GBD burden estimates
  #
  # NOTE: Set to NULL to turn smoothing off
  o$gbd_coverage_smoother = "kernel"  # OPTIONS: "kernel" or "spline"

  # ---- Impact function settings ----
  
  # Use multiple cores to speed up impact function fitting
  o$parallel = FALSE
  
  # Multiply impact when fitting for more consistent FVP-impact scales
  o$impact_scaler = 1000

  # Defalt x scale for evaluating impact functions
  o$eval_x_scale = 2  # Not a critical value - often overwritten with actual FVPs
  
  # Number of initial years to average over for back-projecting impact ratio
  o$init_impact_years = 1 # 3
  
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
  o$plot_inputs       = FALSE
  o$plot_non_modelled = FALSE
  o$plot_imputation   = FALSE
  o$plot_impact       = TRUE
  o$plot_uncertainty  = FALSE
  o$plot_history      = TRUE  # Primary result

  # ---- Plotting settings ----
  
  # Colour packages and palettes (see colour_scheme in auxiliary.R)
  o$palette = list(
    disease = "pals::kovesi.rainbow",  # 11 values needed
    region  = "brewer::paired",        # 6 values needed
    income  = "brewer::dark2")         # 4 values needed

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

