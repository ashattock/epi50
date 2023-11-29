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
  
  # For other countries and years, extract coverage from WIISE database
  wiise_dt = coverage_wiise(vimc_countries_dt)
  
  # Incorporate non-routine SIA data (from WIISE)
  sia_dt = coverage_sia(vimc_countries_dt)  # See sia.R
  
  # Combine sources and deal with pertussis special case
  source_dt = rbind(vimc_dt, wiise_dt, sia_dt) %>%
    wholecell_acellular_switch()
  
  # Sanity check that no zero entries remain
  if (any(source_dt$fvps <= 1e-10))
    stop("Trival coverage entries identified")
  
  # Summarise, assuming partially targeted SIAs
  coverage_dt = source_dt %>%
    group_by(country, v_a_id, year, age) %>%
    summarise(fvps     = max(fvps),  # Essentially a placeholder until next calculation
              cohort   = mean(cohort), 
              coverage = 1 - prod(1 - coverage)) %>%  # Key assumption
    ungroup() %>%
    # Use combined coverage - unless FVPs already eclipses 100% (unlikely)
    mutate(fvps = pmax(cohort * coverage, fvps)) %>%
    as.data.table()
  
  # Save both datatables to file
  save_table(source_dt,   "coverage_source")
  save_table(coverage_dt, "coverage")
  
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
    filter(country %in% all_countries(), 
           year    %in% o$analysis_years) %>%
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
  raw_file = paste0(o$pth$input, "wiise_coverage.csv")
  
  # If file has already been downloaded, read it now
  if (file.exists(raw_file)) {
    raw_dt = fread(raw_file)
    
  } else {  # Otherwise we'll need to download
    
    # Non-VIMC coverage taken from WIISE database
    raw_url = "https://whowiise.blob.core.windows.net/upload/coverage--2021.xlsx"
    raw_dt  = read_url_xls(raw_url, sheet = 1) 
    
    # Save csv file locally for easy re-loading
    write_delim(raw_dt, raw_file, delim = ",")
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
    filter(country %in% all_countries(), 
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
    # Remove countries and years already covered by VIMC...
    left_join(y  = vimc_countries_dt, 
              by = c("country", "vaccine", "year")) %>%
    filter(is.na(source)) %>%
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
    # Apply v_a ID...
    left_join(y  = v_a_dt, 
              by = c("vaccine")) %>%
    select(-activity) %>%
    # Append ages...
    left_join(y  = age_dt, 
              by = "vaccine", 
              relationship = "many-to-many") %>%
    # Calculate fully vaccinated people...
    left_join(y  = table("wpp_pop"), 
              by = c("country", "year", "age")) %>%
    mutate(doses = coverage * pop) %>%
    calculate_fvps() %>%
    # Tidy up...
    arrange(country, v_a_id, year, age) %>%
    mutate(source = "wiise") %>%
    as.data.table()
  
  return(wiise_dt)
}

# ---------------------------------------------------------
# Distribute across pertussis vaccine types by country and year
# ---------------------------------------------------------
wholecell_acellular_switch = function(coverage_dt) {
  
  # NOTE: A first attempt to defining when countries switched - to be improved
  
  # Details of who switched to acellular pertussis and when
  #
  # TODO: Do a more thorough job of this
  switch = list(
    year    = 1995,
    economy = c("hic"), 
    region  = c("EURO", "PAHO", "WPRO"))
  
  # Construct datatable of switch details 
  switch_dt = table("country") %>%
    filter(economy %in% switch$economy, 
           region  %in% switch$region) %>%
    mutate(switch_year = switch$year) %>%
    select(country, switch_year)
  
  # IDs of both wholecell and acellular pertussis vaccines
  id = list(
    wp = table("v_a")[vaccine == "wPer3", v_a_id], 
    ap = table("v_a")[vaccine == "aPer3", v_a_id])
    
  # Only a subset of that defined should be acelluar
  acellular_dt = coverage_dt %>%
    filter(v_a_id == id$ap) %>%
    left_join(y  = switch_dt, 
              by = "country") %>%
    replace_na(list(switch_year = Inf)) %>%
    filter(year > switch_year) %>%
    select(-switch_year)
  
  # Everything else should be wholecell
  wholecell_dt = acellular_dt %>%
    select(country, year, age, source) %>%
    mutate(remove = TRUE) %>%
    full_join(y  = coverage_dt, 
              by = c("country", "year", "age", "source")) %>%
    filter(v_a_id %in% unlist(id), 
           is.na(remove)) %>%
    select(-remove) %>%
    # Covert to wholecell...
    mutate(v_a_id = id$wp) %>%
    group_by(country, v_a_id, year, age, source) %>%
    summarise(fvps   = sum(fvps), 
              cohort = mean(cohort)) %>%
    ungroup() %>%
    # Recalculate coverage...
    mutate(coverage = pmin(fvps / cohort, 1)) %>%
    select(all_of(names(coverage_dt))) %>%
    as.data.table()

  # Recombine all data
  switched_dt = coverage_dt %>%
    filter(!v_a_id %in% unlist(id)) %>%
    rbind(wholecell_dt) %>%
    rbind(acellular_dt) %>%
    arrange(country, v_a_id, year, age)
  
  # Sanity check that we haven't altered total number of FVPs
  diff = sum(coverage_dt$fvps) - sum(switched_dt$fvps)
  if (abs(diff) > 1e-6)
    stop("FVPs have been lost/gained through wholecell-acellular switch")
    
  return(switched_dt)
}

# ---------------------------------------------------------
# Effect of multiple booster doses for DTP
# ---------------------------------------------------------
calculate_fvps = function(coverage_dt) {
  
  # NOTES: 
  #  - Using mean for pop as all values should all equal
  #  - Coverage bounded by 100%, but retain true FVPs
  
  # For primary schedule, assume all new FVPs
  primary_dt = coverage_dt %>% 
    filter(!grepl("_BX$", vaccine)) %>%
    group_by(country, v_a_id, year, age) %>%
    summarise(fvps   = sum(doses),  # Using sum
              cohort = mean(pop)) %>%
    ungroup() %>%
    as.data.table()
  
  # For boosters, consecutive doses are for the same person
  booster_dt = coverage_dt %>% 
    filter(grepl("_BX$", vaccine)) %>%
    group_by(country, v_a_id, year, age) %>%
    summarise(fvps   = max(doses),  # Using max
              cohort = mean(pop)) %>%
    ungroup() %>%
    as.data.table()
  
  # Re-bind everything together and calculate coverage
  fvps_dt = rbind(primary_dt, booster_dt) %>%
    mutate(coverage = pmin(fvps / cohort, 1))
  
  return(fvps_dt)
}

