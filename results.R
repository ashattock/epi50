###########################################################
# RESULTS
#
# xxxxxxx
#
###########################################################

# ---------------------------------------------------------
# xxxxxxxx
# ---------------------------------------------------------
run_results = function() {
  
  # Only continue if specified by do_step
  if (!is.element(6, o$do_step)) return()
  
  message("* Producing results")
  
  # ---- Input data plots ----
  
  # Check plotting flag
  if (o$plot_inputs) {
    
    # Coverage data density by age
    plot_coverage_age_density()
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
  
  # ---- Uncertainty plots ----
  
  # Check plotting flag
  if (o$plot_uncertainty) {
    
    # Plot parameters of fitted beta distribution to vaccine efficacy (GBD diseases only)
    fig_name = "Uncertainty distributions - GBD diseases"
    plot_gbd_uncertainty_dist(fig_name)  # See plotting.R
    
    # Plot uncertainty draws for all diseases
    fig_name = "Uncertainty draws - All diseases"
    plot_draws(fig_name)  # See plotting.R
    
    # Plot annual totals to check alignment of means
    fig_name = "Uncertainty bounds - Annual total"
    plot_annual_total(fig_name)  # See plotting.R
  }
}

