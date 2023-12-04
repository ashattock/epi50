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
  if (!is.element(7, o$do_step)) return()
  
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
    
    # SDI - HAQi relationship
    plot_sdi_haqi()
  }
  
  # ---- Non-modelled plots ----
  
  # Check plotting flag
  if (o$plot_non_modelled) {
    
    # Plot vaccine efficacy profiles for non-modelled pathogens
    plot_vaccine_efficacy()
    
    # Effective coverage with waning immunity for non-modelled pathogens
    plot_effective_coverage()
    
    # Deaths and DALYs averted for non-modelled pathogens
    plot_non_modelled()
  }
  
  # ---- Imputation plots ----
  
  # Check plotting flag
  if (o$plot_imputation) {
    
    # Plot predictor and response relationships
    # plot_target()
    plot_covariates()

    # Plot imputation outcomes
    plot_impute_fit()
    plot_impute_countries()
  }
  
  # ---- Impact function plots ----
  
  # Check plotting flag
  if (o$plot_impact) {
    
    # # Exploratory plots of data used to fit impact functions
    plot_impact_data()

    # Plot function selection statistics
    plot_model_selection()
    
    # Plot impact function evaluation
    plot_model_fits()
  }
  
  # ---- Historical results ----
  
  # Check plotting flag
  if (o$plot_history) {
    
    # Main results plot - impact over time
    plot_history()
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

