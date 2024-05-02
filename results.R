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
    
    # Coverage data density by age
    plot_coverage_age_density()
  }
  
  # ---- Static model plots ----
  
  # Check plotting flag
  if (o$plot_static) {
    
    # Global Burden of Disease death estimates by age
    plot_gbd_estimates()
    
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
    
    # Plot for deaths only
    for (metric in o$metrics[1]) {
      
      # Plot predicted vs observed for all countries
      plot_impute_quality(metric)
      
      # Plot predicted vs observed for each country
      plot_impute_perform(metric)
    }
  }
  
  # ---- Impact function plots ----
  
  # Check plotting flag
  if (o$plot_impact) {
    
    # Plot for deaths only
    for (metric in o$metrics[1]) {
      
      # Plot impact function evaluation
      plot_model_fits(metric)
      
      # Plot function selection statistics
      plot_model_selection(metric)
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
    
    # Regional differences in child mortality changes
    plot_mortality_change()
    
    # Plot absolute and relative probability of death in 2024
    plot_survival_increase()
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
  
  # Convert to pathogen names for consistency with paper
  pathogen_dict = c(
    "Meningitis A" = "Neisseria meningitidis A", 
    "Pneumococcal" = "Streptococcus pneumoniae")
  
  # Define age groups to produce results for
  age_groups = c(max(o$ages), 5)
  
  # Initiate results list
  results_list = list()
  
  # Function to format all numbers
  fmt = function(x) thou_sep(round(x, -3))
  
  # Iterate through key metrics
  for (metric in o$metrics) {
    
    # Load results for this metric
    samples_dt = read_rds("history", "all_samples", metric)
    
    # Repeat for each age group
    for (age_group in age_groups) {
      
      message("   ~ ", metric, ": under ", age_group)
      
      # Proportion of impact in this age group
      age_dt = table("impact_age_multiplier") %>%
        filter(age <= age_group) %>%
        group_by(d_v_a_id) %>%
        summarise(scaler = sum(scaler)) %>%
        ungroup() %>%
        as.data.table()
      
      # Subset ages and summarise uncertainty
      results_dt = samples_dt %>%
        # Apply age effect...
        left_join(y  = age_dt,
                  by = "d_v_a_id") %>%
        mutate(impact = impact * scaler) %>%
        select(-scaler) %>%
        # Summarise uncertainty in final year...
        summarise_uncertainty(cumulative = TRUE) %>%
        filter(year == max(year)) %>%
        # Append disease details...
        left_join(y  = table("d_v_a"),
                  by = "d_v_a_id") %>%
        left_join(y  = table("disease_name"),
                  by = "disease") %>%
        mutate(disease_name = recode(disease_name, !!!pathogen_dict)) %>%
        # Append region...
        append_region_name() %>%
        select(disease = disease_name, region,
               country, impact, lower, upper) %>%
        arrange(region, disease)
      
      # Results for each disease
      global_dt = results_dt %>%
        group_by(disease) %>%
        summarise(x  = sum(impact), 
                  lb = sum(lower), 
                  ub = sum(upper)) %>%
        ungroup() %>%
        mutate(region = "Global", 
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
      
      # Total by region
      total_dt = rbind(global_dt, regional_dt) %>%
        group_by(region) %>%
        summarise(x  = sum(x), 
                  lb = sum(lb), 
                  ub = sum(ub)) %>%
        ungroup() %>%
        mutate(disease = "~Total~", 
               .after = region) %>%
        as.data.table()
      
      # Format strings
      format_dt = global_dt %>%
        rbind(regional_dt) %>%
        rbind(total_dt) %>%
        mutate(result = paste0(
          fmt(x), "\n[", fmt(lb), " - ", fmt(ub), "]")) %>%
        mutate(metric = metric, 
               age    = paste("Under", age_group)) %>%
        select(region, disease, metric, age, result) %>%
        arrange(region, disease)
      
      # Store result
      results_list[[paste1(metric, age_group)]] = format_dt
    }
  }
  
  # Pivot metrics wider for pretty table
  results_dt = rbindlist(results_list) %>%
    append_metric_name() %>%
    mutate(group = paste0(metric_impact, ": ", age)) %>%
    select(region, disease, group, result) %>%
    pivot_wider(names_from  = group, 
                values_from = result) %>%
    as.data.table()
  
  # Save result
  fwrite(results_dt, file = paste0(o$pth$figures, "Table 1.csv"))
}

