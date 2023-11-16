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
  
  # ---- Set up ----
  
  # All diseases of interest (everything non-modelled)
  diseases = table("disease") %>%
    filter(source == "gbd") %>%
    pull(disease)
    
  # Vaccine and shedule info for these diseases
  d_v_dt = table("d_v_a") %>%
    filter(disease %in% diseases) %>%
    mutate(schedule = ifelse(
      test = grepl("_BX$", vaccine), 
      yes  = "booster", 
      no   = "primary")) %>%
    select(disease, vaccine, schedule)
  
  # Vaccine efficacy profiles
  efficacy_dt = table("vaccine_efficacy_profiles") %>%
    pivot_wider(id_cols = time, 
                names_from  = vaccine,
                values_from = profile) %>%
    select(time, all_of(d_v_dt$vaccine)) %>%
    as.data.table()
  
  # Vaccine coverage
  coverage_dt = table("coverage") %>%
    left_join(y  = table("v_a"), 
              by = "v_a_id") %>%
    filter(vaccine %in% d_v_dt$vaccine) %>%
    select(country, v_a_id, vaccine, year, age, fvps)

  # Full factorial country-year-age grid (we join results to this)
  full_dt = expand_grid(
    country = all_countries(),
    year    = o$data_years, 
    age     = o$data_ages) %>%
    as.data.table()
  
  # ---- xxxx ----
  
  # TODO: Allow each d_v_a to be 'targeted' or 'non-targeted'
  
  total_list = list()
  
  # Iterate through these diseases
  for (disease in diseases) {
    
    message(" - ", disease)
    
    vaccine_info = d_v_dt %>% 
      filter(disease == !!disease) %>% 
      select(vaccine, schedule)
    
    # ---- Step 1: collate coverage over sources
    
    has_booster = any(vaccine_info$schedule == "booster")
    
    if (has_booster == TRUE) {
      
      browser()
      
      coverage_dt %>% 
        inner_join(y  = vaccine_info, 
                   by = "vaccine") %>%
        select(country, year, age, schedule, fvps) %>%
        pivot_wider(names_from  = schedule, 
                    values_from = fvps)
      
      
    }
    
    browser()
    
    
    
    coverage_dt
    
    # Initiate list for effective results
    effective_list = list()
    
    # Iterate through these vaccines
    for (vaccine in vaccine_info$vaccine) {
      
      # Effiacy profile for this vaccine
      efficacy = efficacy_dt[[vaccine]]
      
      # Whether vaccine is primary series or booster dose
      schedule = ifelse(
        test = grepl("_BX$", vaccine), 
        yes  = "booster", 
        no   = "primary")
      
      # Apply effective coverage function to each coverage entry
      effective_list[[vaccine]] = coverage_dt %>%
        filter(vaccine == !!vaccine) %>%
        dtapply(efficacy_fn, efficacy) %>%
        rbindlist() %>%
        # First summarise for effective FVPs...
        group_by(country, year, age) %>%
        summarise(effective_fvps = sum(value)) %>%   # Assumes targeted SIA
        # summarise(effective_fvps = 1 - prod(1 - value)) %>%  # Assumes non-targeted SIA
        ungroup() %>%
        # Join with full factorial grid...
        full_join(y  = full_dt, 
                  by = names(full_dt)) %>%
        replace_na(list(effective_fvps = 0)) %>%
        # Append vaccine details...
        mutate(schedule = schedule) %>%
        arrange(country, year, age) %>%
        as.data.table()
    }
    
    
    browser()
    
    effective_dt = rbindlist(effective_list) %>%
      pivot_wider(names_from  = schedule, 
                  values_from = effective_fvps) # %>%
      # mutate(primary, booster)
    
    # First summarise for cumulative total coverage...
    group_by(country, year, age) %>%
      # summarise(effective_coverage = 1 - prod(1 - value)) %>%  # Assumes non-targeted vaccination
      summarise(effective_coverage = min(sum(value), 1)) %>%   # Assumes targeted vaccination
      ungroup()
    
    
    
    
    
    total_list[[id]] = effective_coverage_dt %>%
      mutate(d_v_a_id = id, .after = 1)
    
    browser()
    
    relative_risk_fn = get("gbd_rr")

    # Load either deaths_averted or deaths_disease for this strata
    rr_dt = table("gbd_estimates") %>%
      filter(d_v_a_id == id) %>%
      select(-d_v_a_id) %>%
      # Append all-cause deaths...
      full_join(y  = table("wpp_death"),
                by = c("country", "year", "age")) %>%
      rename(deaths_allcause = death) %>%
      # Append total coverage...
      inner_join(y  = effective_coverage_dt,
                 by = c("country", "year", "age")) %>%
      arrange(country, year, age) # %>%
      # Calculate relative risk...
      # relative_risk_fn()
    
    # Use this to estimate deaths without a vaccine...
    mutate(deaths_without = deaths_disease / (1 - efficacy * effective_coverage), 
           deaths_without = ifelse(effective_coverage == 0, NA, deaths_without))



    # ---- Use deaths averted to calculate DALYs averted ----

    browser()

    run_dalys()
  }
  
  # Save total coverage datatable to file
  total_dt = rbindlist(total_list)
  save_rds(total_dt, "non_modelled", "effective_coverage")
  
  # ---- Impact visualisation plots ----
  
  plot_effective_coverage()
}

# ---------------------------------------------------------
# Calculate relative risk for GBD disease
# ---------------------------------------------------------
gbd_rr = function(strata_dt) {
  
  # Append vaccine efficacy
  strata_dt %<>%
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_id") %>%
    # Append vaccine efficacy...
    left_join(y  = table("vaccine_efficacy"), 
              by = c("disease", "vaccine")) %>%
    select(all_of(names(strata_dt)), efficacy) %>%
    # Use this to estimate deaths without a vaccine...
    mutate(deaths_without = deaths_disease / (1 - efficacy * effective_coverage), 
           deaths_without = ifelse(effective_coverage == 0, NA, deaths_without))
  
  # Shorthand variable for readability
  #
  # Note that w - d equivalent to a (ie deaths averted from vimc_rr function)
  o = strata_dt$deaths_allcause  # Deaths [o]bserved
  d = strata_dt$deaths_disease   # Deaths attributable to this [d]isease
  w = strata_dt$deaths_without   # Deaths estimated [w]ithout a vaccine
  e = strata_dt$efficacy         # Vaccine [e]fficacy
  
  # Calculate relative risk
  rr_dt = strata_dt %>%
    mutate(rr = (o - d + (1 - e) * w) / (o - d + w)) %>%
    # Remove redundant variables...
    select(-efficacy, -deaths_without)
  
  return(rr_dt)
}

# ---------------------------------------------------------
# Calculate relative risk for VIMC disease
# ---------------------------------------------------------
vimc_rr = function(strata_dt) {
  
  browser()
  
  # Shorthand variable for readability
  o = strata_dt$deaths_allcause  # Deaths [o]bserved
  a = strata_dt$deaths_averted   # Deaths [a]verted from vaccination
  c = strata_dt$effective_coverage   # Lifetime vaccine [c]overage
  
  # Calculate relative risk
  rr_dt = strata_dt %>%
    mutate(rr = (o - (a * (1 - c) / c)) / (o + a)) %>%
    # Tidy up trivial values...
    mutate(rr = ifelse(effective_coverage == 0, NA, rr), 
           rr = ifelse(is.infinite(rr), NA, rr))
  
  return(rr_dt)
}

# ---------------------------------------------------------
# Extract averted deaths from relative risk equation
# ---------------------------------------------------------
calculate_averted_deaths = function(deaths_obs, coverage, rr) {
  
  # rr = (o - (a * (1 - c) / c)) / (o + a)
  #
  # => a = (o * (1 - rr) * c) / (1 + (rr - 1) * c)
  #
  # where:
  #   o = deaths observed
  #   a = deaths averted from vaccine
  #   c = coverage
  
  # Estimate deaths averted given relative risk
  #
  # NOTE: Equation derived by solving vimc_rr for a
  averted_deaths <- deaths_obs * (coverage * (1 - rr) / 
                                    (1 - coverage * (1 - rr)))
  
  return(averted_deaths)
}

# ---------------------------------------------------------
# xxxxxxxxxxx
# ---------------------------------------------------------
efficacy_fn = function(data, efficacy) {
  
  # Indicies for years and ages
  year_idx = match(data$year, o$data_years) : length(o$data_years)
  age_idx  = match(data$age,  o$data_ages)  : length(o$data_ages)
  
  # Index upto only the smallest of these two vectors
  vec_idx = 1 : min(length(year_idx), length(age_idx))
  
  # Effective FVPs: initial efficacy and immunity decay
  effective_fvps = data$fvps * efficacy
  
  # These form the only non-trivial entries
  effective_dt = data.table(
    country = data$country, 
    year    = o$data_years[year_idx[vec_idx]],
    age     = o$data_ages[age_idx[vec_idx]], 
    value   = effective_fvps[vec_idx])
  
  return(effective_dt)
}

