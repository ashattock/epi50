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
  
  d_v_a_dt = table("disease") %>%
    filter(source == "gbd") %>%
    left_join(y  = table("d_v_a"), 
              by = "disease") %>%
    # TEMP...
    filter(activity == "routine") %>%
    select(d_v_a_id, disease, vaccine, activity)
  
  total_list = list()
  
  for (id in d_v_a_dt$d_v_a_id) {
    
    # Details of this d_v_a
    d_v_a = d_v_a_dt[d_v_a_id == id, ]
    
    total_coverage_dt = table("coverage") %>%
      left_join(y  = table("v_a"), 
                by = "v_a_id") %>%
      # Coverage is per vaccine, not per strata...
      filter(vaccine  == d_v_a$vaccine, 
             activity == d_v_a$activity) %>%
      # Calculate coverage by cohort...
      total_coverage(d_v_a)  # See coverage.R
    
    total_list[[id]] = total_coverage_dt %>%
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
      inner_join(y  = total_coverage_dt,
                 by = c("country", "year", "age")) %>%
      arrange(country, year, age) %>%
      # Calculate relative risk...
      relative_risk_fn()



    # ---- Use deaths averted to calculate DALYs averted ----

    browser()

    run_dalys()
  }
  
  # Save total coverage datatable to file
  total_dt = rbindlist(total_list)
  save_rds(total_dt, "non_modelled", "total_coverage")
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
    mutate(deaths_without = deaths_disease / (1 - efficacy * total_coverage), 
           deaths_without = ifelse(total_coverage == 0, NA, deaths_without))
  
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
  c = strata_dt$total_coverage   # Lifetime vaccine [c]overage
  
  # Calculate relative risk
  rr_dt = strata_dt %>%
    mutate(rr = (o - (a * (1 - c) / c)) / (o + a)) %>%
    # Tidy up trivial values...
    mutate(rr = ifelse(total_coverage == 0, NA, rr), 
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
# Perform sanity checks on relative risk calculations
# ---------------------------------------------------------
check_rr = function(rr_dt) {
  
  # Throw error if no valid relative risk results
  if (nrow(rr_dt[rr > 0 & rr < 1]) == 0)
    stop("Error calculating relative risk")
  
  browser()
  
  # Look for missing coverage where deaths averted are non-zero
  if (nrow(rr_dt[coverage == 0 & deaths_averted > 0]) > 0) {
    
    # Proportion of cases
    prop <- round(nrow(rr_dt[coverage == 0 & deaths_averted > 0]) /
                    nrow(rr_dt[deaths_averted > 0]) * 100, 2)
    
    # Throw warning
    warning("Missing coverage in ", prop,
            "% of country-age-years with deaths averted")
  }
  
  # Check for non-sensical numbers
  if (any(range(rr_dt[!is.na(rr)]$rr) < 0 | range(rr_dt[!is.na(rr)]$rr) > 1)) {
    
    # Proportion of cases
    prop <- round(nrow(rr_dt[coverage > 0 & (rr < 0 | rr > 1)]) /
                    nrow(rr_dt) * 100, 2)
    
    # Throw warning
    warning("Over 1 or less than 0 mortality reduction in ", prop,
            "% of country-age-years")
  }
}

