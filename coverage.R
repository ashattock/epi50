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
    select(vaccine, country, year, source) %>%
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
  coverage_file = paste0(o$pth$data, "wiise_coverage.rds")
  
  # If file has already been downloaded, read it now
  if (file.exists(coverage_file)) {
    data_dt = readRDS(coverage_file)
    
  } else {  # Otherwise we'll need to download
    
    # Non-VIMC coverage taken from WIISE database
    data_url = "https://whowiise.blob.core.windows.net/upload/coverage--2021.xlsx"
    data_dt  = read_url_xls(data_url, sheet = 1) 
    
    # Save as an RDS file for easy future loading
    save_rds(data_dt, coverage_file)
  }
  
  # Load WIISE-related vaccine details 
  wiise_info = table("wiise_vaccine")
  
  # ---- Extract coverage ----
  
  # Initiate list to store results
  wiise_list = list()
  
  # Iterate through WIISE vaccines
  for (i in seq_row(wiise_info)) {
    v = wiise_info[i, ]
    
    # Countries and years already covered by VIMC
    vimc_countries = vimc_countries_dt %>%
      filter(vaccine == v$vaccine) %>%
      select(-vaccine)
    
    # Filter data by coverage_category type
    wiise_coverage = data_dt %>%
      setnames(names(.), tolower(names(.))) %>% 
      rename(country  = code, 
             wiise_id = antigen) %>%
      # Select only vaccine and years of interest...
      filter(coverage_category == v$coverage_category, 
             wiise_id %in% v$wiise_id,
             year     %in% o$analysis_years, 
             country  %in% table("country")$country) %>% 
      select(country, year, coverage) %>%
      # Remove countries and years already covered by VIMC...
      left_join(y  = vimc_countries, 
                by = c("country", "year")) %>%
      filter(is.na(source)) %>%
      select(-source) %>%
      # Format coverage...
      replace_na(list(coverage = 0)) %>%
      mutate(coverage = coverage / 100) %>%
      select(country, year, coverage) %>%
      arrange(country, year)
    
    # In some cases coverage is repeated for multiple ages
    wiise_list[[i]] = wiise_info %>% 
      # Expand for each age considered...
      slice(i) %>% 
      expand_grid(expand_age = eval_str(age)) %>%
      select(vaccine, activity, gender, age = expand_age) %>%
      # Repeat coverage values for each age...
      expand_grid(wiise_coverage) %>%
      # Apply v_a ID...
      left_join(y  = table("v_a"), 
                by = c("vaccine", "activity")) %>%
      select(country, v_a_id, gender, year, age, coverage) %>%
      as.data.table()
  }
  
  # ---- Calculate FVPs ----
  
  # TODO: Check whether HPV coverage is actually modelled this way...
  
  # Calculate FVPs from coverage 
  wiise_dt = rbindlist(wiise_list) %>%
    # Combine gender where necessary...
    group_by(country, v_a_id, year, age) %>%
    summarise(coverage = mean(coverage)) %>%
    ungroup() %>%
    # Calculate number of fully vaccinated people...
    inner_join(y  = table("wpp_pop"), 
               by = c("country", "year", "age")) %>%
    rename(cohort = pop) %>%
    mutate(fvps = coverage * cohort) %>%
    # Final formatting...
    select(country, v_a_id, year, age, fvps, cohort, coverage) %>%
    arrange(country, v_a_id, year, age) %>%
    mutate(source = "wiise") %>%
    as.data.table()
  
  return(wiise_dt)
}

# ---------------------------------------------------------
# Calculate total lifetime coverage for each cohort
# ---------------------------------------------------------
total_coverage = function(coverage_dt) {
  
  # TODO: Allow each d_v_a to be 'targeted' or 'non-targeted'
  # TODO: Set a cap on BCG effect at age 15
  # TODO: May be use an efficacy decay for each pathogen
  
  browser()
  
  # Create full combination table
  #
  # NOTE: Final result is sparse => most age-year values will be zero
  full_dt = expand_grid(
    country = unique(coverage_dt$country),
    year    = o$data_years, 
    age     = o$data_ages) %>%
    as.data.table()
  
  # Function to extract total coverage for each data point
  total_coverage_fn = function(i) {
    
    # Data relating to this row of coverage_dt
    data = coverage_dt[i, ]
    
    # Indicies for years and ages
    year_idx = match(data$year, o$data_years) : length(o$data_years)
    age_idx  = match(data$age,  o$data_ages)  : length(o$data_ages)
    
    # Index upto only the smallest of these two vectors
    vec_idx = 1 : min(length(year_idx), length(age_idx))
    
    # These form the only non-trivial entries
    total_dt = data.table(
      country = data$country, 
      year    = o$data_years[year_idx[vec_idx]],
      age     = o$data_ages[age_idx[vec_idx]], 
      value   = data$coverage)
    
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

