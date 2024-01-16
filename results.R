###########################################################
# RESULTS
#
# Call plotting functions as defined by o$plot_xxx flags (see
# options.R). All plotting functions themselves live in 
# plotting.R.
#
###########################################################

# ---------------------------------------------------------
# Call plots as defined by o$plot_xxx flags
# ---------------------------------------------------------
run_results = function() {
  
  # Only continue if specified by do_step
  if (!is.element(8, o$do_step)) return()
  
  message("* Producing results")
  
  # ---- Input data plots ----

  # Check plotting flag
  if (o$plot_inputs) {

    # Methodology pathogen-country-scope figure
    plot_scope()

    # Total number of FVP over time by source
    plot_total_fvps()

    # Plot coverage density by disease
    plot_coverage()

    # Coverage data density by age
    plot_coverage_age_density()

    # GBD death estimates
    plot_gbd_estimates()
  }

  # ---- Static model plots ----

  # Check plotting flag
  if (o$plot_static) {

    # Plot vaccine efficacy profiles for static model pathogens
    plot_vaccine_efficacy()

    # Effective coverage with waning immunity for static model pathogens
    plot_effective_coverage()

    # Deaths and DALYs averted for static model pathogens
    plot_static()
  }
  
  # ---- Imputation plots ----
  
  # Check plotting flag
  if (o$plot_imputation) {
    
    # Plot predictor and response relationships
    plot_covariates()

    # Plot imputation quality of fit
    plot_impute_quality()
    
    # Plot train-predict countries
    plot_impute_countries()
  }
  
  # ---- Impact function plots ----
  
  # Check plotting flag
  if (o$plot_impact) {
    
    # Exploratory plots of data used to fit impact functions
    plot_impact_data()
    
    # Plot all-time impact per FVPs
    plot_impact_fvps(scope = "all_time")
    
    # Plot impact vs coverage by vaccine, income, and decade 
    # plot_impact_coverage()
    
    # Plot function selection statistics
    plot_model_selection()
    
    # Plot impact function evaluation
    plot_model_fits()
  }
  
  # ---- Historical results ----
  
  # Check plotting flag
  if (o$plot_history) {
    
    # Plot inital impact ratios used to back project
    plot_impact_fvps(scope = "initial")

    # Main results plot - historical impact over time
    plot_historical_impact()
    
    # Another key plot - change in child mortality rates over time
    plot_child_mortality()
  }
  
  # ---- Uncertainty plots ----
  
  # Check plotting flag
  if (o$plot_uncertainty) {
    
    # Plot annual totals to check alignment of means
    plot_annual_total()
    
    # Plot uncertainty draws for all diseases
    plot_uncertainty_draws()
    
    # Plot fit parameters around vaccine efficacy (GBD diseases only)
    plot_gbd_uncertainty_dist()
  }
}

