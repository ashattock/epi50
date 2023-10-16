###########################################################
# COVERAGE
#
# xxxxxxxxxxxx
#
###########################################################

# ---------------------------------------------------------
# xxxxxxxxx
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
    arrange(country, v_a_id, sex, year, age) %>%
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
    # Combine sex where necessary...
    mutate(sex = ifelse(gender == "Both", "b", "x")) %>%
    group_by(country, disease, vaccine, activity, sex, year, age) %>%
    summarise(fvps     = sum(fvps_adjusted),
              coverage = fvps / sum(cohort_size)) %>%
    ungroup() %>%
    mutate(sex = "b") %>%  # As both genders nbow combined
    # Append v_a ID...
    left_join(y  = table("v_a"), 
              by = c("vaccine", "activity")) %>%
    select(country, v_a_id, sex, year, age, fvps, coverage) %>%
    arrange(country, v_a_id, year, age) %>%
    as.data.table()
  
  return(vimc_dt)
}

# ---------------------------------------------------------
# Extract coverage from WIISE database
# ---------------------------------------------------------
coverage_wiise = function(vimc_countries) {
  
  # TODO: Combine genders for HPV ??
  
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
      select(vaccine, activity, sex, age = expand_age) %>%
      # Repeat coverage values for each age...
      expand_grid(wiise_coverage) %>%
      # Final formatting...
      left_join(y  = table("v_a"), 
                by = c("vaccine", "activity")) %>%
      select(country, v_a_id, sex, year, age, coverage) %>%
      as.data.table()
  }
  
  # ---- Calculate FVPs ----
  
  wiise_dt = table("wpp_input") %>%
    # Total number of people...
    group_by(country, year, age) %>%
    summarise(cohort_size = sum(nx)) %>%
    ungroup() %>%
    # Calculate number of fully vaccinated people...
    inner_join(y  = rbindlist(wiise_list), 
               by = c("country", "year", "age")) %>%
    mutate(fvps = coverage * cohort_size) %>%
    # Final formatting...
    select(country, v_a_id, sex, year, age, fvps, coverage) %>%
    arrange(country, v_a_id, sex, year, age) %>%
    as.data.table()
  
  return(wiise_dt)
}

