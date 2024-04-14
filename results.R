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
  
  # Only continue if specified by run_module
  if (!is.element(8, o$run_module)) return()
  
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
  }
  
  # ---- External model plots ----
  
  # Check plotting flag
  if (o$plot_external) {
    
    # Plot outcomes from each external model
    plot_external_models()
    
    # Missing coverage data by country
    plot_missing_data()
  }
  
  # ---- Static model plots ----
  
  # Check plotting flag
  if (o$plot_static) {
    
    # Global Burden of Disease death estimates by age
    plot_gbd_estimates()

    # Proportion of GBD burden we have coverage data for
    plot_gbd_missing()

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
    
    # Repeat for deaths and DALYs
    for (metric in "deaths") { # o$metrics) {
      
      # # Plot model choice by region
      # plot_model_choice(metric)

      # Plot predicted vs. observed for all countries
      plot_impute_quality(metric)
      
      # Plot predicted vs observed for each country
      plot_impute_perform(metric)
      
      # # Plot fit to data in train-predict countries
      # plot_impute_fit(metric)
      # 
      # # Plot validation
      # plot_validation(metric)
    }
  }
  
  # ---- Impact function plots ----
  
  # Check plotting flag
  if (o$plot_impact) {
    
    # Repeat for deaths and DALYs
    for (metric in o$metrics) {
      
      # Exploratory plots of data used to fit impact functions
      plot_impact_data(metric)

      # Plot all-time impact per FVPs
      # plot_impact_fvps(metric, scope = "all_time")

      # Plot function selection statistics
      plot_model_selection(metric)
      
      # Plot impact function evaluation
      plot_model_fits(metric)
      
      # Plot impact vs coverage by vaccine, income, and decade 
      # plot_impact_coverage(metric)
    }
  }
  
  # ---- Historical results ----
  
  # Check plotting flag
  if (o$plot_history) {
    
    # Create full results table with bounds
    all_results_table()
    
    # Inital impact ratios used to back project
    for (metric in o$metrics)
      plot_impact_fvps(metric, scope = "initial")

    # Main results plot - historical impact over time
    plot_historical_impact()

    # Equivalent plot for each region and income
    plot_historical_impact(by = "region")
    plot_historical_impact(by = "income")

    # Non-cumulative, pathogen specific results
    for (metric in o$metrics)
      plot_temporal_impact(metric)
    
    # Infant mortality rates over time with and without vaccination
    plot_infant_mortality()
    
    # xxx
    plot_mortality_region()
    
    # Regional differences in child mortality changes
    plot_mortality_change()
    
    # Measles deaths in context of all cause deaths
    plot_measles_in_context()
    
    # Plot absolute and relative probability of death in 2024
    plot_prob_death_age()
    
    # Plot absolute and relative probability of death in 2024
    plot_survival_increase()
    
    # Plot comparison of EPI50 outcomes vs VIMC outcomes
    plot_vimc_comparison()
  }
}

# ---------------------------------------------------------
# Full table of disease-specific results with bounds
# ---------------------------------------------------------
all_results_table = function() {
  
  message("  - Creating main results table")
  
  # Initiate results list
  results_list = list()
  
  # Function to format all numbers
  fmt = function(x) thou_sep(round(x, -3))
  
  # Iterate through key metrics
  for (metric in o$metrics) {
    
    message("   ~ ", metric)
    
    # Load results and summarise bounds
    results_list[[metric]] = 
      read_rds("history", "all_samples", metric) %>%
      # Summarising uncertainty...
      summarise_uncertainty(cumulative = TRUE) %>%
      left_join(y  = table("d_v_a"), 
                by = "d_v_a_id") %>%
      left_join(y  = table("disease_name"), 
                by = "disease") %>%
      # Results for each disease and each region...
      lazy_dt() %>%
      filter(year == max(year)) %>%
      group_by(disease_name) %>%
      summarise(x  = sum(impact), 
                lb = sum(lower), 
                ub = sum(upper)) %>%
      ungroup() %>%
      # Format strings...
      mutate(result = paste0(
        fmt(x), "  [", fmt(lb), " - ", fmt(ub), "]")) %>%
      mutate(metric = metric) %>%
      select(disease = disease_name, metric, result) %>%
      as.data.table()
  }
  
  # Pivot metrics wider for pretty table
  results_dt = rbindlist(results_list) %>%
    pivot_wider(names_from  = metric, 
                values_from = result)
  
  # File path to save to 
  save_path = paste0(o$pth$figures, "manuscript")
  save_file = file.path(save_path, "Table 1.csv")
  
  # Save in output folder
  fwrite(results_dt, file = save_file)
}

