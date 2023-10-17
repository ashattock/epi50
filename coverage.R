###########################################################
# COVERAGE
#
# xxxxxxxxxxxx
#
###########################################################

# ---------------------------------------------------------
# Parent function for preparing vaccine coverage data
# ---------------------------------------------------------
prepare_coverage = function() {
  
  # TODO: Incorporate SIA data
  # prepare_sia()  # See sia.R
  
  message(" - Coverage data")
  
  # Extract coverage for VIMC pathogens
  vimc_dt = coverage_vimc()
  
  # However not every country is covered by VIMC for these pathogens
  vimc_countries = vimc_dt %>%
    left_join(y  = table("v_a"), 
              by = "v_a_id") %>%
    select(vaccine, country) %>%
    arrange(vaccine, country) %>%
    unique() %>%
    split(.$vaccine) %>%
    lapply(function(x) x$country)
  
  # For everything remaining, extract coverage from WIISE database
  wiise_dt = coverage_wiise(vimc_countries)
  
  # Combine sources
  rbind(vimc_dt, wiise_dt) %>%
    filter(fvps > 0) %>%  # Remove trivial values
    arrange(country, v_a_id, year, age) %>%
    save_table("coverage")
}

# ---------------------------------------------------------
# Extract coverage from VIMC outputs
# ---------------------------------------------------------
coverage_vimc = function() {
  
  # Extract VIMC vaccine coverage data
  vimc_dt = fread(paste0(o$pth$input, "vimc_coverage.csv")) %>%
    select(country, disease, vaccine, activity = activity_type, 
           gender, year, age, fvps_adjusted, cohort_size) %>%
    # Combine gender where necessary...
    mutate(gender = ifelse(gender == "Both", "b", "x")) %>%
    group_by(country, disease, vaccine, activity, gender, year, age) %>%
    summarise(fvps     = sum(fvps_adjusted),
              coverage = fvps / sum(cohort_size)) %>%
    ungroup() %>%
    # Append v_a ID...
    left_join(y  = table("v_a"), 
              by = c("vaccine", "activity")) %>%
    select(country, v_a_id, year, age, fvps, coverage) %>%
    arrange(country, v_a_id, year, age) %>%
    as.data.table()
  
  return(vimc_dt)
}

# ---------------------------------------------------------
# Extract coverage from WIISE database
# ---------------------------------------------------------
coverage_wiise = function(vimc_countries) {
  
  # ---- Load data ----
  
  temp_file = "temp/wiise_coverage.rds"
  if (file.exists(temp_file)) {
    data_dt = readRDS(temp_file)
  } else {
    
    # Non-VIMC coverage taken from WIISE database
    data_url = "https://whowiise.blob.core.windows.net/upload/coverage--2021.xlsx"
    data_dt  = read_url_xls(data_url, sheet = 1) 
    
    saveRDS(data_dt, temp_file)
  }
  
  # Load WIISE-related vaccine details 
  wiise_info = table("wiise_vaccine")
  
  # ---- Extract coverage ----
  
  # Initiate list to store results
  wiise_list = list()
  
  # Iterate through WIISE vaccines
  for (i in seq_len(nrow(wiise_info))) {
    v = wiise_info[i, ]
    
    # Countries to consider (ignoring VIMC countries)
    wiise_countries = setdiff(
      x = table("country")$country, 
      y = vimc_countries[[v$vaccine]])
    
    # Filter data by coverage_category type
    wiise_coverage = data_dt %>%
      setnames(names(.), tolower(names(.))) %>% 
      rename(country  = code, 
             wiise_id = antigen) %>%
      # Select only vaccine and countries of interest...
      filter(coverage_category == v$coverage_category, 
             wiise_id %in% v$wiise_id,
             country  %in% wiise_countries) %>%  
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
  
  # Total population by country, year, and age
  pop_dt = table("wpp_input") %>%
    group_by(country, year, age) %>%
    summarise(cohort_size = sum(nx)) %>%
    ungroup() %>%
    as.data.table()
  
  # Calculate FVPs from coverage 
  wiise_dt = rbindlist(wiise_list) %>%
    # Combine gender where necessary...
    group_by(country, v_a_id, year, age) %>%
    summarise(coverage = mean(coverage)) %>%
    ungroup() %>%
    # Calculate number of fully vaccinated people...
    inner_join(y  = pop_dt, 
               by = c("country", "year", "age")) %>%
    mutate(fvps = coverage * cohort_size) %>%
    # Final formatting...
    select(country, v_a_id, year, age, fvps, coverage) %>%
    arrange(country, v_a_id, year, age) %>%
    as.data.table()
  
  return(wiise_dt)
}

# ---------------------------------------------------------
# Calculate total lifetime coverage for each cohort
# ---------------------------------------------------------
total_coverage = function(coverage_dt) {
  
  # TODO: Allow each d_v_a to be 'targeted' or 'non-targeted'
  
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
  
  # Coverage data values to work through  
  total_idx = seq_len(nrow(coverage_dt))
  
  # Apply total coverage function to each row
  total_dt = lapply(total_idx, total_coverage_fn) %>%
    rbindlist() %>%
    # First summarise for cumulative total coverage...
    group_by(country, year, age) %>%
    # summarise(value = 1 - prod(1 - value)) %>%  # Assumes non-targeted vaccination
    summarise(value = min(sum(value), 1)) %>%   # Assumes targeted vaccination
    ungroup() %>%
    # Then join with full grid...
    full_join(y  = full_dt, 
              by = names(full_dt)) %>%
    replace_na(list(value = 0)) %>%
    arrange(country, year, age) %>%
    # Append d_v_a info...
    mutate(d_v_a_id = unique(coverage_dt$v_a_id), 
           .after = 1) %>%
    as.data.table()
  
  # TODO: Set a cap on BCG effect at age 15
  
  return(total_dt)
}

# ---------------------------------------------------------
# Calculate FVPs using coverage and demographic data
# ---------------------------------------------------------
cov2fvp = function(coverage_dt) {
  
  browser() # Needs updating...
  
  # Load demographic data
  wpp_input = table("wpp_input")
  
  # Total number of people per country (both genders combined)
  #
  # NOTE: nx := number of people
  both_dt = wpp_input[, .(nx = sum(nx)), .(country, age, year)]
  both_dt[, gender := 3]
  
  # Combine so we have both genders seperate and combined
  pop_dt = rbind(wpp_input, both_dt, fill = T)
  
  # Join with coverage details
  fvp_dt = merge(coverage_dt, pop_dt[, .(country, age, gender, year, nx)])
  
  # Then just a simple calculation for FVPs
  fvp_dt[, fvps := coverage * nx]
  fvp_dt[, nx := NULL]
  
  return(fvp_dt)
}

