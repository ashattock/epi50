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
  
  # Extract coverage for non-VIMC pathogens
  wiise_dt = coverage_wiise(vimc_dt)
  
  browser() # v_at_table is now v_a_table
  
  # Combine sources
  coverage = table("wpp_input") %>%
    # Total number of people...
    group_by(country, year, age) %>%
    summarise(cohort_size = sum(nx)) %>%
    ungroup() %>%
    # Calculate FVPs for non-VIMC diseases...
    inner_join(y  = wiise_dt, 
               by = c("country", "year", "age")) %>%
    mutate(fvps = coverage * cohort_size) %>%
    select(-cohort_size) %>% 
    # Combine VIMC diseases and format...
    bind_rows(vimc_dt) %>%
    filter(fvps > 0) %>%  # Remove trivial values
    left_join(y  = table("v_a"), 
              by = c("vaccine", "activity")) %>%
    select(country, v_at_id, year, age, sex_id, fvps, coverage) %>%
    arrange(country, v_at_id, year, age, sex_id) %>%
    setDT()
  
  # coverage[, age := as.integer(age)]
  # coverage[, sex_id := as.integer(sex_id)]
  
  # Upload table to database
  upload_object(coverage, "coverage")
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
    arrange(country, disease, vaccine, activity, year, age) %>%
    as.data.table()
  
  browser()
  
  left_join(y  = table("v_a"), 
            by = c("vaccine", "activity")) %>%
  
  return(vimc_dt)
}

# ---------------------------------------------------------
# Extract coverage from WIISE database
# ---------------------------------------------------------
coverage_wiise = function(vimc_dt) {
  
  # TODO: Combine genders for HPV ??
  
  temp_file = "temp/wiise_coverage.rds"
  if (file.exists(temp_file)) {
    data_dt = readRDS(temp_file)
  } else {
    
    # Non-VIMC coverage taken from WIISE database
    data_url = "https://whowiise.blob.core.windows.net/upload/coverage--2021.xlsx"
    data_dt  = read_url_xls(data_url, sheet = 1) 
    
    saveRDS(data_dt, temp_file)
  }
  
  all_countries = table("country")$country
  
  wiise_table = table("wiise")
  
  wiise_list = list()
  
  for (i in seq_len(nrow(wiise_table))) {
    j = wiise_table[i, ]
    
    # Some countries may already be covered for this vaccine by VIMC
    vimc_countries = vimc_dt %>%
      filter(vaccine == j$vaccine) %>%
      pull(country) %>%
      unique()
    
    # Filter data by coverage_category type
    wiise_list[[i]] = data_dt %>%
      setnames(names(.), tolower(names(.))) %>% 
      rename(country  = code, 
             wiise_id = antigen) %>%
      # Select only vaccine and countries of interest...
      filter(coverage_category == j$coverage_category, 
             wiise_id %in% j$wiise_id,
             country  %in% all_countries, 
             !country %in% vimc_countries) %>%  # Ignoring VIMC countries
      mutate(vaccine = j$vaccine) %>%
      # Format coverage...
      replace_na(list(coverage = 0)) %>%
      mutate(coverage = coverage / 100) %>%
      # Append additional vaccination details...
      expand_grid(activity = j$activity, 
                  sex      = j$sex, 
                  age      = eval_str(j$age)) %>%
      select(country, vaccine, activity, 
             sex, year, age, coverage) %>%
      as.data.table()
  }
  
  wiise_dt = rbindlist(wiise_list) %>%
    arrange(country, vaccine, year, age)
  
  return(wiise_dt)
}

