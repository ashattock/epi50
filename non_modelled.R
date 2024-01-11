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
  if (!is.element(3, o$do_step)) return()
  
  message("* Estimating vaccine impact for non-modelled pathogens")
  
  # All diseases of interest (everything non-modelled)
  diseases = table("disease") %>%
    filter(source == "gbd") %>%
    pull(disease)
  
  # Iterate through these diseases
  for (disease in diseases) {
    
    # Effective coverage considering waning immunity and boosters
    effective_coverage(disease)
    
    # Deaths averted considering effective coverage and GBD disease burden
    deaths_averted(disease)
    
    # Use deaths averted to calculate DALYs
    # dalys_averted(disease)  # See dalys.R
  }
  
  # Compile all results
  compile_outputs(diseases)
  
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
  
  # ---- Set up ----
  
  # Vaccines targeting this disease
  vaccine_dt = table("d_v_a") %>% 
    filter(disease == !!disease) %>%
    mutate(type = str_remove(vaccine, "([0-9]+|_BX)$")) %>% 
    select(vaccine, type)
  
  # Shorthand for vaccines for this disease
  vaccines = vaccine_dt$vaccine
  
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
    year    = o$years, 
    age     = o$ages) %>%
    as.data.table()
  
  # ---- Waning immunity per vaccine / shedule ----
  
  # Initiate list for effective immunity
  immunity_list = list()
  
  # Iterate through these vaccines
  for (vaccine in vaccines) {
    
    message(" - ", disease, ", ", vaccine)
    
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
      mutate(vaccine  = vaccine, 
             schedule = schedule) %>%
      left_join(y  = vaccine_dt, 
                by = "vaccine") %>%
      arrange(country, year, age) %>%
      as.data.table()
  }
  
  # Bind results of each vaccine / schedule
  immunity_dt = rbindlist(immunity_list)
  
  # And save the result in it's own right
  save_rds(immunity_dt, "non_modelled_d", "immunity", disease)
  
  # ---- Overall effective coverage ----
  
  message(" - Calculating effective coverage")
  
  # TODO: Can we avoid this conversion back to coverage?
  
  # Inititate list for effective coverage
  effective_list = list()
  
  # Iterate through vacicne types
  for (vaccine_type in unique(vaccine_dt$type)) {
    
    # Convert effective FVPs to effective coverage
    effective_v_dt = immunity_dt %>%
      filter(type == vaccine_type) %>%
      select(-vaccine, -type) %>%
      # Weight between primary & booster...
      weight_booster(vaccine_type) %>%
      # Convert to effective coverage...
      left_join(y  = table("wpp_pop"),
                by = c("country", "year", "age")) %>%
      mutate(coverage = pmin(effective / pop, o$max_coverage)) %>%
      replace_na(list(coverage = 0)) %>%
      # Append vaccine and disease details and tidy up...
      mutate(disease = disease, 
             type    = vaccine_type) %>%
      select(country, disease, type, 
             year, age, effective, pop, coverage) %>%
      arrange(country, disease, type, year, age) %>%
      as.data.table()
    
    # Store result for this vaccine type in list
    effective_list[[vaccine_type]] = effective_v_dt
    
    # And save the result in it's own right
    save_rds(effective_v_dt, "non_modelled_t", "effective_coverage", vaccine_type)
  }
  
  # Combine vaccine types (considered additive)
  effective_d_dt = rbindlist(effective_list) %>%
    group_by(country, disease, year, age) %>%
    summarise(effective = sum(effective), 
              pop       = mean(pop)) %>%
    ungroup() %>%
    mutate(coverage = pmin(effective / pop, o$max_coverage)) %>%
    replace_na(list(coverage = 0)) %>%
    select(-effective, -pop) %>%
    as.data.table()
  
  # Save this result to file
  save_rds(effective_d_dt, "non_modelled_d", "effective_coverage", disease)
}

# ---------------------------------------------------------
# People effectively covered over time considering waning immunity
# ---------------------------------------------------------
waning_immunity = function(data, profile) {
  
  # Indicies for years and ages
  year_idx = match(data$year, o$years) : length(o$years)
  age_idx  = match(data$age,  o$ages)  : length(o$ages)
  
  # Index upto only the smallest of these two vectors
  vec_idx = 1 : min(length(year_idx), length(age_idx))
  
  # Effective FVPs: initial efficacy and immunity decay
  effective_fvps = data$fvps * profile
  
  # These form the only non-trivial entries
  waning_immunity_dt = data.table(
    country   = data$country, 
    year      = o$years[year_idx[vec_idx]],
    age       = o$ages[age_idx[vec_idx]], 
    covered   = data$fvps, 
    effective = effective_fvps[vec_idx])
  
  return(waning_immunity_dt)
}

# ---------------------------------------------------------
# Weight between primary & booster for total effective FVPs
# ---------------------------------------------------------
weight_booster = function(immunity_dt, type) {
  
  # Flag for whether vaccine has booster
  has_booster = "booster" %in% immunity_dt$schedule
  
  # Number covered by schedule (primary and booster)
  covered_dt = immunity_dt %>%
    select(country, year, age, covered, schedule) %>%
    pivot_wider(names_from  = schedule, 
                values_from = covered) %>%
    # Append trivial column if no booster doses...
    {if (!has_booster) mutate(., booster = 0) else .} %>%
    as.data.table()
  
  # Proportion of primary cases that have booster (capped at 100%)
  weighting_dt = covered_dt %>%
    mutate(primary = pmax(primary, booster), 
           weight  = pmin(booster / primary, 1)) %>%
    filter(!is.nan(weight)) %>%
    select(country, year, age, weight) %>%
    as.data.table()
  
  # Save this weighting to file
  save_rds(weighting_dt, "non_modelled_w", "booster_weight", type)
  
  # TODO: I'm now questioning this, I think it should be:
  #       effective = primary * (1 - weight) + booster
  #       The idea is those being boosted should not be double counted
  
  # Effective number of FVP after booster considerations
  weighted_dt = immunity_dt %>%
    select(country, year, age, effective, schedule) %>%
    pivot_wider(names_from  = schedule, 
                values_from = effective) %>%
    # Append trivial column if no booster doses...
    {if (!has_booster) mutate(., booster = 0) else .} %>%
    # Append booster weighting details...
    left_join(y  = weighting_dt, 
              by = c("country", "year", "age")) %>%
    replace_na(list(weight = 0)) %>%
    # Weight efficacies based on number covered...
    mutate(effective = primary * (1 - weight) + booster) %>% # * weight) %>%  # ?? See comment above ??
    select(country, year, age, effective) %>%
    as.data.table()
  
  return(weighted_dt)
}

# ---------------------------------------------------------
# Deaths averted considering effective coverage and GBD disease burden
# ---------------------------------------------------------
deaths_averted = function(disease) {
  
  # ---- Deaths averted for this disease ----
  
  # Load effective coverage for this disease from file
  effective_dt = read_rds("non_modelled_d", "effective_coverage", disease)
  
  # Load disease deaths, append coverage, and estimate deaths averted
  averted_disease = effective_dt %>%
    inner_join(y  = table("gbd_estimates"),
               by = c("disease", "country", "year", "age")) %>%
    # Estimate deaths without a vaccine and deaths averted...
    mutate(deaths_without = deaths_disease / (1 - coverage), 
           deaths_averted = deaths_without - deaths_disease) %>%
    select(country, disease, year, age, deaths_disease, deaths_averted) %>%
    as.data.table()
  
  # Save this result to file
  save_rds(averted_disease, "non_modelled_d", "deaths_averted", disease)
  
  # ---- Attribute impact to vaccine schedule ----
  
  # Relative weighting of effective coverage for each vaccine
  weight_dt = read_rds("non_modelled_d", "immunity", disease) %>%
    left_join(y  = table("d_v_a"),
              by = "vaccine") %>%
    group_by(country, year, age) %>%
    mutate(weight = effective / sum(effective)) %>%
    ungroup() %>%
    filter(!is.nan(weight)) %>%
    select(d_v_a_id, vaccine, country, year, age, weight) %>%
    as.data.table()
  
  # We'll split impact across schedules by previously calculated weighting
  averted_vaccine = averted_disease %>%
    # Attribute impact to each shedule...
    inner_join(y  = weight_dt, 
               by = c("country", "year", "age")) %>%
    # mutate(deaths_averted = deaths_averted * weight) %>%
    # Summarise over age ...
    group_by(country, d_v_a_id, year) %>%
    summarise(impact = sum(deaths_averted * weight)) %>%
    ungroup() %>%
    # Tidy up...
    select(country, d_v_a_id, year, impact) %>%
    arrange(country, d_v_a_id, year) %>%
    as.data.table()
  
  # Number of 'FVP' by vaccine - yeah, confusing concept in the context of boosters
  fvps_dt = table("coverage") %>%
    left_join(y  = table("v_a"), 
              by = "v_a_id") %>%
    left_join(y  = table("d_v_a"), 
              by = c("vaccine", "activity")) %>%
    filter(disease == !!disease) %>%
    # Summarise over age ...
    group_by(country, d_v_a_id, year) %>%
    summarise(fvps = sum(fvps)) %>%
    ungroup() %>%
    as.data.table()
  
  # Full vaccine-specific results 
  result_list = averted_vaccine %>%
    inner_join(y  = fvps_dt, 
               by = c("country", "d_v_a_id", "year")) %>%
    select(country, d_v_a_id, year, fvps, impact) %>%
    split(.$d_v_a_id)
  
  # Function for saving vaccine-specific results
  save_result_fn = function(x, name)
    save_rds(x, "non_modelled_v", "deaths_averted", name)
  
  # Apply save function to each vaccine
  napply(result_list, save_result_fn)
}

# ---------------------------------------------------------
# Compile and save all non-modelled outputs 
# ---------------------------------------------------------
compile_outputs = function(diseases) {
  
  message(" - Compiling results")
  
  # IDs targeting these disease
  d_v_a_id = table("d_v_a") %>% 
    filter(disease %in% diseases) %>%
    pull(d_v_a_id)
  
  # Vaccine types targeting these disease
  vaccine_types = table("d_v_a") %>% 
    filter(disease %in% diseases) %>%
    mutate(type = str_remove(vaccine, "([0-9]+|_BX)$")) %>% 
    pull(type) %>% 
    unique()
  
  # IDs for loading disease/vaccine results
  ids = list(
    disease = diseases,
    vaccine = d_v_a_id, 
    type    = vaccine_types)
  
  # All outputs to compile
  outputs = list(
    disease = qc(deaths_averted, effective_coverage), #, dalys_averted),
    vaccine = qc(deaths_averted), #, dalys_averted),
    type    = qc(effective_coverage))
  
  for (group in names(outputs)) {
    
    # Sub directory to load disease/vaccine results from
    load_dir = paste1("non_modelled", str_sub(group, 1, 1))
    
    # Function to read rds files
    load_fn = function(file) 
      read_rds(load_dir, file)
    
    # Loop through all outputs
    for (output in outputs[[group]]) {
      
      # Load all files related to this output
      output_dt = paste1(output, ids[[group]]) %>%
        lapply(load_fn) %>%
        rbindlist()
      
      # Save compiled results to file
      save_rds(output_dt, "non_modelled", output, group)
    }
  }
}

