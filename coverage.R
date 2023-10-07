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
  # prepare_sia()
  
  message(" - Coverage data")
  
  stop("Under development...")
  
  # Load VIMC coverage data
  vimc_dt = coverage_vimc()
  
  # Non-VIMC coverage taken from WIISE database
  data_url = "https://whowiise.blob.core.windows.net/upload/coverage--2021.xlsx"
  data_dt  = read_url_xls(data_url, sheet = 1) 
  
  browser()
  
  # Extract coverage for each non-VIMC vaccine
  non_vimc_dt = 1 : nrow(wiise_vaccine_table) %>%
    lapply(coverage_wiise, data_dt, vimc_dt) %>%
    rbindlist()
  
  browser() # v_at_table is now v_a_table
  
  # Combine sources
  coverage = table("wpp_input") %>%
    # Calculate FVPs for non-VIMC diseases...
    group_by(country, year, age) %>%
    summarise(cohort_size = sum(nx)) %>%
    ungroup() %>%
    inner_join(y  = non_vimc_dt, 
               by = c("country", "year", "age")) %>%
    mutate(fvps = coverage * cohort_size) %>%
    select(-cohort_size) %>% 
    # Combine VIMC diseases and format...
    bind_rows(vimc_dt) %>%
    filter(fvps > 0) %>%  # Remove trivial values
    left_join(y  = v_at_table, 
              by = c("vaccine", "activity_type")) %>%
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
  
  # Dictionary for recoding gender to sex_id
  gender_dict = setNames(1 : 3, c("Male", "Female", "Both"))
  
  # Path to VIMC output datatable
  vimc_path = system.file("extdata", "vimc_coverage.csv", package = "vieIA2030")
  
  # NOTE: This isn't quite right. Need to consider cases that have 'both' gender
  #       values AND values for 'male' and 'female'. See get_temporal_fvps().
  
  # Extract coverage
  vimc_dt = fread(vimc_path) %>%
    mutate(sex_id = recode(gender, !!!gender_dict)) %>%
    group_by(country, disease, vaccine, activity_type, sex_id, year, age) %>%
    summarise(fvps     = sum(fvps_adjusted),
              coverage = fvps / sum(cohort_size)) %>%
    ungroup() %>%
    select(-disease) %>%
    setDT()
  
  return(vimc_dt)
}

# ---------------------------------------------------------
# Extract coverage from WIISE for a given vaccine 
# ---------------------------------------------------------
coverage_wiise = function(v_idx, data_dt, vimc_dt) {
  
  browser()
  
  # Details of vaccine at index v_idx
  v = wiise_vaccine_table[v_idx, ]
  
  # Some countries may already be covered for this vaccine by VIMC
  vimc_country = vimc_dt %>%
    filter(vaccine == v$vaccine) %>%
    pull(country) %>%
    unique()
  
  # Filter data by coverage_category type
  coverage_dt = data_dt %>%
    setnames(names(.), tolower(names(.))) %>% 
    rename(country  = code, 
           wiise_id = antigen) %>%
    # Select only vaccine and countries of interest...
    filter(coverage_category == v$coverage_category, 
           wiise_id %in% v$wiise_id,
           country  %in% country_table$country, 
           !country %in% vimc_country) %>%  # Ignoring VIMC countries
    mutate(vaccine = v$vaccine) %>%
    # Format coverage...
    replace_na(list(coverage = 0)) %>%
    mutate(coverage = coverage/ 100) %>%
    # Append additional vaccination details...
    expand_grid(activity_type = v$activity, 
                sex_id = v$sex, 
                age    = eval_str(v$age)) %>%
    select(country, vaccine, activity_type, 
           sex_id, year, age, coverage) %>%
    arrange(country, vaccine, year, age) %>%
    setDT()
  
  return(coverage_dt)
}

