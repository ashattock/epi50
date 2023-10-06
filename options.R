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
  o$analysis_years = 2000 : 2030  # Vaccine deployed across these dates
  o$future_years   = 2021 : 2030  # Cohort years
  
  # Year and age ranges stored in coverage database
  o$data_years = 2000 : 2039  # Vaccine effect calculated across these dates
  o$data_ages  = 0 : 95
  
  # ---- Database settings ----
  
  # Version number of google database
  o$db_version = "v2"
  
  # Force a fresh database pull even if cache loading available
  o$force_db_pull = FALSE
  
  # Only load from cache if pulled within the last n hours
  o$cache_hour_limit = 168 * 52  # 168 hours = 1 week
  
  # # ECDC data links
  # o$ecdc_api = 
  #   list(cases = "https://opendata.ecdc.europa.eu/covid19/nationalcasedeath_eueea_daily_ei/csv", 
  #        hosp  = "https://opendata.ecdc.europa.eu/covid19/hospitalicuadmissionrates/csv")
  # 
  # # ETH effective reproduction number estimates
  # o$eth_api = "https://raw.githubusercontent.com/covid-19-Re/dailyRe-Data/master/<country>-estimates.csv"
  # 
  # # API endpoint for national-level Oxford Stringency Index data
  # o$osi_api = "https://covidtrackerapi.bsg.ox.ac.uk/api/v2/stringency/date-range/"
  # 
  # # Data dictionary: ECDC hospital & ICU indicators
  # o$data_dict$ecdc = c(hospital_beds = "Daily hospital occupancy", 
  #                      icu_beds      = "Daily ICU occupancy",
  #                      hospital_admissions = "Weekly new hospital admissions per 100k", 
  #                      icu_admissions      = "Weekly new ICU admissions per 100k")
  
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
  
  # ---- Prepare output ----
  
  # Append helpful properties
  # o = append_shortcuts(o)
  
  # Display analysis details
  if (!quiet) message(" - Analysis name: ", o$analysis_name)
  
  return(o)
}

# ---------------------------------------------------------
# Override options set in my_options file
# Called by: set_options()
# ---------------------------------------------------------
append_shortcuts = function(o) {
  
  # Disease source (VIMC or GBD)
  for (i in unique(disease_table$source))
    o$disease[[i]] = disease_table[source == i, disease]
  
  return(o)
}

