###########################################################
# NON-MODELLED
#
# Estimate impact (deaths and DALYs averted) for non-modelled
# pathogens using Global Burden of Disease estimates.
#
###########################################################

# ---------------------------------------------------------
# Parent function for calculating impact for non-modelled pathogens
# ---------------------------------------------------------
run_non_modelled = function() {
  
  # Only continue if specified by do_step
  if (!is.element(2, o$do_step)) return()
  
  message("* Estimating vaccine impact for non-modelled pathogens")
  
  # All diseases of interest (everything non-modelled)
  diseases = table("disease") %>%
    filter(source == "gbd") %>%
    pull(disease)
  
  # Iterate through these diseases
  for (disease in diseases) {

    message(" - ", disease)

    # Effective coverage considering waning immunity and boosters
    effective_coverage(disease)

    # Deaths averted considering effective coverage and GBD disease burden
    deaths_averted(disease)

    # Use deaths averted to calculate DALYs ----
    # dalys_averted(disease)
  }
  
  # Compile all results
  outputs = qc(effective_coverage, deaths_averted) #, dalys_averted)
  compile_outputs(outputs, diseases)
  
  # ---- Plot results ----
  
  # Effective coverage with waning immunity for non-modelled pathogens
  plot_effective_coverage()
  
  # Deaths and DALYs averted for non-modelled pathogens
  plot_non_modelled()
}

# ---------------------------------------------------------
# Effective coverage considering waning immunity and boosters
# ---------------------------------------------------------
effective_coverage = function(disease) {
  
  # CALCULATION PROCESS: 
  # 1) Number of immune people over time - considering waning immunity
  # 2) Number of people covered with different vaccine / schedules
  # 3) Weight immunity by primary & booster dosing
  
  # TODO: Why are so many values hitting this cap?...
  
  effective_capped = 0.95
  
  # ---- Set up ----
  
  # Vaccines targeting this disease
  vaccines = table("d_v_a") %>% 
    filter(disease == !!disease) %>%
    pull(vaccine)
  
  # Vaccine immunity efficacy profiles
  profile_dt = table("vaccine_efficacy_profiles") %>%
    pivot_wider(id_cols = time, 
                names_from  = vaccine,
                values_from = profile) %>%
    select(time, all_of(vaccines)) %>%
    as.data.table()
  
  # Vaccine coverage
  coverage_dt = table("coverage") %>%
    left_join(y  = table("v_a"), 
              by = "v_a_id") %>%
    filter(vaccine %in% vaccines) %>%
    select(country, v_a_id, vaccine, year, age, fvps)
  
  # Full factorial country-year-age grid (we join results to this)
  full_dt = expand_grid(
    country = all_countries(),
    year    = o$data_years, 
    age     = o$data_ages) %>%
    as.data.table()
  
  # ---- Waning immunity per vaccine / shedule ----
  
  # Initiate list for effective immunity
  immunity_list = list()
  
  # Iterate through these vaccines
  for (vaccine in vaccines) {
    
    message("  > ", vaccine)
    
    # Immunity effiacy profile for this vaccine
    profile = profile_dt[[vaccine]]
    
    # Whether vaccine is primary series or booster dose
    schedule = ifelse(
      test = grepl("_BX$", vaccine), 
      yes  = "booster", 
      no   = "primary")
    
    # For each entry, apply waning immunity over time
    immunity_list[[vaccine]] = coverage_dt %>%
      filter(vaccine == !!vaccine) %>%
      dtapply(waning_immunity, profile) %>%
      rbindlist() %>%
      # Summarise where multiple sources of immunity...
      group_by(country, year, age) %>%
      summarise(covered   = sum(covered), 
                effective = sum(effective)) %>%
      ungroup() %>%
      # Join with full factorial grid...
      full_join(y  = full_dt, 
                by = names(full_dt)) %>%
      replace_na(list(covered   = 0, 
                      effective = 0)) %>%
      # Append vaccine details...
      mutate(schedule = schedule) %>%
      arrange(country, year, age) %>%
      as.data.table()
  }
  
  # Bind results of each vaccine / schedule
  immunity_dt = rbindlist(immunity_list)
  
  # ---- Weight between primary & booster ----
  
  # Flag for whether vaccine has booster
  has_booster = "booster" %in% immunity_dt$schedule
  
  # If booster, apply primary-booster weighting
  if (has_booster == TRUE)
    weighted_dt = weight_booster(immunity_dt)
  
  # If no booster, nothing to do here
  if (has_booster == FALSE) {
    weighted_dt = immunity_dt %>%
      select(country, year, age, effective)
  }
  
  # ---- Overall effective coverage ----
  
  # TODO: Can we avoid this conversion back to coverage?
  
  # Convert effective FVPs to effective coverage
  effective_dt = weighted_dt %>%
    left_join(y  = table("wpp_pop"),
              by = c("country", "year", "age")) %>%
    mutate(coverage = pmin(effective / pop, effective_capped)) %>%
    replace_na(list(coverage = 0)) %>%
    # Append disease details and tidy up...
    mutate(disease = disease) %>%
    select(country, disease, year, age, coverage) %>%
    arrange(country, disease, year, age) %>%
    as.data.table()
  
  # Save this result to file
  save_rds(effective_dt, "non_modelled", "effective_coverage", disease)
}

# ---------------------------------------------------------
# People effectively covered over time considering waning immunity
# ---------------------------------------------------------
waning_immunity = function(data, profile) {
  
  # Indicies for years and ages
  year_idx = match(data$year, o$data_years) : length(o$data_years)
  age_idx  = match(data$age,  o$data_ages)  : length(o$data_ages)
  
  # Index upto only the smallest of these two vectors
  vec_idx = 1 : min(length(year_idx), length(age_idx))
  
  # Effective FVPs: initial efficacy and immunity decay
  effective_fvps = data$fvps * profile
  
  # These form the only non-trivial entries
  waning_immunity_dt = data.table(
    country   = data$country, 
    year      = o$data_years[year_idx[vec_idx]],
    age       = o$data_ages[age_idx[vec_idx]], 
    covered   = data$fvps, 
    effective = effective_fvps[vec_idx])
  
  return(waning_immunity_dt)
}

# ---------------------------------------------------------
# Weight between primary & booster for total effective FVPs
# ---------------------------------------------------------
weight_booster = function(immunity_dt) {
  
  # Weighting towards booster efficacy (and away from primary efficacy)
  weighting_dt = immunity_dt %>%
    select(country, year, age, covered, schedule) %>%
    pivot_wider(names_from  = schedule, 
                values_from = covered) %>%
    # Proportion of primary cases that have booster (capped at 100%)...
    mutate(primary = pmax(primary, booster), 
           weight  = pmin(booster / primary, 1)) %>%
    filter(!is.nan(weight)) %>%
    select(country, year, age, weight) %>%
    as.data.table()
  
  # Effective number of FVP after booster considerations
  weighted_dt = immunity_dt %>%
    select(country, year, age, effective, schedule) %>%
    pivot_wider(names_from  = schedule, 
                values_from = effective) %>%
    # Append booster weighting details...
    left_join(y  = weighting_dt, 
              by = c("country", "year", "age")) %>%
    replace_na(list(weight = 0)) %>%
    # Weight efficacies based on number covered...
    mutate(effective = primary * (1 - weight) + booster * weight) %>%
    select(country, year, age, effective) %>%
    as.data.table()
  
  return(weighted_dt)
}

# ---------------------------------------------------------
# Deaths averted considering effective coverage and GBD disease burden
# ---------------------------------------------------------
deaths_averted = function(disease) {
  
  # Load effective coverage for this disease from file
  effective_dt = read_rds("non_modelled", "effective_coverage", disease)
  
  # Load disease deaths, append coverage, and estimate deaths averted
  averted_dt = effective_dt %>%
    inner_join(y  = table("gbd_estimates"),
               by = c("disease", "country", "year", "age")) %>%
    # Use this to estimate deaths without a vaccine...
    mutate(deaths_without = deaths_disease / (1 - coverage), 
           deaths_averted = deaths_without - deaths_disease) %>%
    select(country, disease, year, age, deaths_disease, deaths_averted)
  
  # Save this result to file
  save_rds(averted_dt, "non_modelled", "deaths_averted", disease)
}

# ---------------------------------------------------------
# Compile and save all non-modelled outputs 
# ---------------------------------------------------------
compile_outputs = function(outputs, diseases) {
  
  message(" - Compiling results")

  # Function to read rds files
  load_fn = function(file)
    read_rds("non_modelled", file)
  
  # Loop through all outputs
  for (output in outputs) {
    
    # Load all files related to this output
    output_dt = paste1(output, diseases) %>%
      lapply(load_fn) %>%
      rbindlist() %>%
      arrange(country, disease, year, age)
    
    # Save compiled results to file
    save_rds(output_dt, "non_modelled", output)
  }
}

