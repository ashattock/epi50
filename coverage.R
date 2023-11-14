###########################################################
# COVERAGE
#
# All coverage related functionality in one place.
#
###########################################################

# ---------------------------------------------------------
# Parent function for preparing vaccine coverage data
# ---------------------------------------------------------
prepare_coverage = function() {
  
  message(" - Coverage data")
  
  # Extract coverage for VIMC pathogens
  vimc_dt = coverage_vimc()
  
  # However not every country is covered by VIMC for these pathogens
  vimc_countries_dt = vimc_dt %>%
    left_join(y  = table("v_a"),
              by = "v_a_id") %>%
    select(country, vaccine, year, source) %>%
    arrange(vaccine, country, year) %>%
    unique()
  
  # For everything remaining, extract coverage from WIISE database
  wiise_dt = coverage_wiise(vimc_countries_dt)
  
  # Finally, incorporate SIA data (from WIISE)
  sia_dt = coverage_sia(vimc_countries_dt)  # See sia.R
  
  # Combine sources
  rbind(vimc_dt, wiise_dt, sia_dt) %>%
    filter(fvps > 0) %>%  # Remove trivial values
    arrange(country, v_a_id, year, age) %>%
    save_table("coverage")
  
  # ---- Data visualisation plots ----
  
  # Plot total number of FVP over time
  plot_total_fvps()
  
  # Coverage data density by age
  plot_coverage_age_density()
}

# ---------------------------------------------------------
# Extract coverage from VIMC outputs
# ---------------------------------------------------------
coverage_vimc = function() {
  
  message("  > VIMC coverage")
  
  # Extract VIMC vaccine coverage data
  vimc_dt = fread(paste0(o$pth$input, "vimc_coverage.csv")) %>%
    select(country, disease, vaccine, activity = activity_type, 
           gender, year, age, fvps_adjusted, cohort_size) %>%
    filter(year %in% o$analysis_years) %>%
    # Combine gender where necessary...
    mutate(gender = ifelse(gender == "Both", "b", "x")) %>%
    group_by(country, disease, vaccine, activity, gender, year, age) %>%
    summarise(fvps     = sum(fvps_adjusted),
              cohort   = sum(cohort_size),
              coverage = fvps / cohort) %>%
    ungroup() %>%
    # Append v_a ID...
    inner_join(y  = table("v_a"), 
               by = c("vaccine", "activity")) %>%
    select(country, v_a_id, year, age, fvps, cohort, coverage) %>%
    arrange(country, v_a_id, year, age) %>%
    mutate(source = "vimc") %>%
    as.data.table()
  
  return(vimc_dt)
}

# ---------------------------------------------------------
# Extract coverage from WIISE database
# ---------------------------------------------------------
coverage_wiise = function(vimc_countries_dt) {
  
  message("  > WIISE coverage")
  
  # ---- Load data ----
  
  # File path for already-downloaded WIISE coverage data
  raw_file = paste0(o$pth$data, "wiise_raw.rds")
  
  # If file has already been downloaded, read it now
  if (file.exists(raw_file)) {
    raw_dt = read_rds(raw_file)
    
  } else {  # Otherwise we'll need to download
    
    # Non-VIMC coverage taken from WIISE database
    raw_url = "https://whowiise.blob.core.windows.net/upload/coverage--2021.xlsx"
    raw_dt  = read_url_xls(raw_url, sheet = 1) 
    
    # Save as an RDS file for easy future loading
    save_rds(raw_dt, raw_file)
  }
  
  # ---- Wrangle WIISE data ----
  
  # Parse 'interventions' into EPI50 vaccines
  data_dt = raw_dt %>%
    # Convert to lower case...
    setnames(names(.), tolower(names(.))) %>% 
    mutate_if(is.character, tolower) %>%
    # Reduce columns...
    select(country = code, intervention = antigen, 
           year, coverage, source = coverage_category) %>% 
    # Remove any unknown countries...
    mutate(country = toupper(country)) %>%
    filter(country %in% table("country")$country, 
           year    %in% o$analysis_years) %>%
    # Convert coverage to proportion...
    mutate(coverage = coverage / 100) %>%
    filter(coverage > 0) %>%
    # Use WUENIC data as primary source...
    mutate(wuenic   = ifelse(source == "wuenic", coverage, NA), 
           coverage = ifelse(source != "wuenic", coverage, NA)) %>%
    # Compare against average of all other sources...
    group_by(country, intervention, year) %>%
    summarise(wuenic = mean(wuenic,   na.rm = TRUE),
              other  = mean(coverage, na.rm = TRUE)) %>%
    ungroup() %>%
    # Salvage coverage from non-WUENIC sources...
    mutate(wuenic = ifelse(is.nan(wuenic), other, wuenic)) %>%
    select(country, intervention, year, coverage = wuenic) %>%
    filter(coverage > 0) %>%
    # Interpret 'intervention'...
    left_join(y  = table("vaccine_dict"), 
              by = "intervention", 
              relationship = "many-to-many") %>%
    filter(!is.na(vaccine)) %>%
    # Tidy up...
    select(country, vaccine, year, coverage) %>%
    arrange(country, vaccine, year) %>%
    as.data.table()
  
  # Function for parsing and expanding to single age bins
  expand_age_fn = function(x)
    y = expand_grid(x, expand_age = eval_str(x$age))
  
  # Expanded datatable of ages per vaccine
  age_dt = table("vaccine") %>%
    select(vaccine, age) %>%
    dtapply(expand_age_fn) %>%
    rbindlist() %>%
    select(vaccine, age = expand_age)
  
  # Routine activities (or 'all' for GBD pathogens)
  v_a_dt = table("v_a") %>%
    filter(activity %in% c("routine", "all"))
  
  # Append age and calculate FVPs
  wiise_dt = data_dt %>%
    # Remove countries and years already covered by VIMC...
    left_join(y  = vimc_countries_dt, 
              by = c("country", "vaccine", "year")) %>%
    filter(is.na(source)) %>%
    select(-source) %>%
    # Append ages...
    left_join(y  = age_dt, 
              by = "vaccine", 
              relationship = "many-to-many") %>%
    # Calculate fully vaccinated people...
    left_join(y  = table("wpp_pop"), 
              by = c("country", "year", "age")) %>%
    rename(cohort = pop) %>%
    mutate(fvps = coverage * cohort) %>%
    # Summarise in case of multiple interventions...
    group_by(country, vaccine, year, age) %>%
    summarise(cohort = mean(cohort), 
              fvps   = sum(fvps)) %>%
    ungroup() %>%
    mutate(coverage = pmin(fvps / cohort, 1)) %>%
    # Apply v_a ID...
    left_join(y  = v_a_dt, 
              by = c("vaccine")) %>%
    # Tidy up...
    select(country, v_a_id, year, age, fvps, cohort, coverage) %>%
    arrange(country, v_a_id, year, age) %>%
    mutate(source = "wiise") %>%
    as.data.table()
  
  return(wiise_dt)
}

# ---------------------------------------------------------
# Calculate effective coverage for each cohort
# ---------------------------------------------------------
total_coverage = function(coverage_dt, d_v_a) {
  
  # TODO: Allow each d_v_a to be 'targeted' or 'non-targeted'
  
  # Create full combination table
  #
  # NOTE: Final result is sparse => most age-year values will be zero
  full_dt = expand_grid(
    country = unique(coverage_dt$country),
    year    = o$data_years, 
    age     = o$data_ages) %>%
    as.data.table()
  
  # Extract waning immunity profile for this d_v
  profile = table("vaccine_efficacy_profiles") %>%
    filter(disease == d_v_a$disease, 
           vaccine == d_v_a$vaccine) %>%
    pull(profile)
  
  # Normalise profile between 0 and 1
  norm_profile = profile / max(profile)
  
  if (length(norm_profile) == 0)
    stop("No waning immunity profile detected")
  
  # Function to extract total coverage for each data point
  total_coverage_fn = function(i) {
    
    # Data relating to this row of coverage_dt
    data = coverage_dt[i, ]
    
    # Indicies for years and ages
    year_idx = match(data$year, o$data_years) : length(o$data_years)
    age_idx  = match(data$age,  o$data_ages)  : length(o$data_ages)
    
    # Index upto only the smallest of these two vectors
    vec_idx = 1 : min(length(year_idx), length(age_idx))
    
    # Represent immunity decay with waning coverage
    waning_coverage = data$coverage * norm_profile
    
    # These form the only non-trivial entries
    total_dt = data.table(
      country = data$country, 
      year    = o$data_years[year_idx[vec_idx]],
      age     = o$data_ages[age_idx[vec_idx]], 
      value   = waning_coverage[vec_idx])
    
    return(total_dt)
  }
  
  # Apply total coverage function to each row of coverage datatable
  total_dt = seq_row(coverage_dt) %>%
    lapply(total_coverage_fn) %>%
    rbindlist() %>%
    # First summarise for cumulative total coverage...
    group_by(country, year, age) %>%
    # summarise(total_coverage = 1 - prod(1 - value)) %>%  # Assumes non-targeted vaccination
    summarise(total_coverage = min(sum(value), 1)) %>%   # Assumes targeted vaccination
    ungroup() %>%
    # Then join with full grid...
    full_join(y  = full_dt, 
              by = names(full_dt)) %>%
    replace_na(list(total_coverage = 0)) %>%
    arrange(country, year, age) %>%
    as.data.table()
  
  return(total_dt)
}

