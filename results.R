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
    plot_smooth_fvps()

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
    
    # Plot predicted vs observed for all countries
    plot_impute_quality("deaths")
    
    # Plot predicted vs observed for each country
    plot_impute_perform("deaths")
    
    # xxx
    plot_model_choice("deaths")
  }
  
  # ---- Impact function plots ----
  
  # Check plotting flag
  if (o$plot_impact) {
    
    # Repeat for deaths and DALYs
    for (metric in o$metrics) {
      
      # Exploratory plots of data used to fit impact functions
      plot_impact_data(metric)

      # Plot function selection statistics
      plot_model_selection(metric)
      
      # Plot impact function evaluation
      plot_model_fits(metric)
    }
  }
  
  # ---- Historical results ----
  
  # Check plotting flag
  if (o$plot_history) {

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
    
    # Inital impact ratios used to back project
    for (metric in o$metrics)
      plot_impact_fvps(metric, scope = "initial")
    
    # Plot comparison of EPI50 outcomes vs VIMC outcomes
    plot_vimc_comparison()
  }
  
  # ---- Main results table ----
  
  # Create full results table with bounds
  if (o$results_table)
    all_results_table()
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

    # Load results and summarise uncertainty
    results_dt = read_rds("history", "all_samples", metric) %>%
      summarise_uncertainty(cumulative = TRUE) %>%
      filter(year == max(year)) %>%
      # Append disease details...
      left_join(y  = table("d_v_a"), 
                by = "d_v_a_id") %>%
      left_join(y  = table("disease_name"), 
                by = "disease") %>%
      # Append region...
      append_region_name() %>%
      select(disease = disease_name, region, 
             country, impact, lower, upper)
    
    # Results for each disease
    global_dt = results_dt %>%
      group_by(disease) %>%
      summarise(x  = sum(impact), 
                lb = sum(lower), 
                ub = sum(upper)) %>%
      ungroup() %>%
      mutate(region = "Gloabl", 
             .before = 1) %>%
      as.data.table()
    
    # Results by region
    regional_dt = results_dt %>%
      group_by(region, disease) %>%
      summarise(x  = sum(impact), 
                lb = sum(lower), 
                ub = sum(upper)) %>%
      ungroup() %>%
      as.data.table()
    
    # Format strings
    format_dt = rbind(global_dt, regional_dt) %>%
      mutate(result = paste0(
        fmt(x), "  [", fmt(lb), " - ", fmt(ub), "]")) %>%
      mutate(metric = metric) %>%
      select(region, disease, metric, result)
    
    # Store result
    results_list[[metric]] = format_dt
  }
  
  # Pivot metrics wider for pretty table
  results_dt = rbindlist(results_list) %>%
    pivot_wider(names_from  = metric, 
                values_from = result)
  
  # Save result
  fwrite(results_dt, file = paste0(o$pth$figures, "Table 1.csv"))
}

