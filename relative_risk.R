###########################################################
# RELATIVE RISK
#
# Calculate relative risk given coverage details and deaths 
# averted. Then impute relative risk for missing countries
# using GBD covariates.
#
###########################################################

# Variables:
#  - deaths_obs: GBD-estimated all cause number of deaths
#  - deaths_disease: GBD-estimated deaths attributable to a disease
#  - deaths_averted: model-estimated deaths averted
#  - predict_averted: predicted deaths averted using stats model (for countries not in VIMC ??)

# ---------------------------------------------------------
# Impute relative risk for a each d-v-a (aka strata)
# ---------------------------------------------------------
run_relative_risk = function() {
  
  # Only continue if specified by do_step
  if (!is.element(2, o$do_step)) return()
  
  message("* Calculating relative risk")
  
  # Initiate list for relative risk results
  rr_list = list()
  
  # All strata to run
  all_strata = table("d_v_a") %>%
    left_join(y = table("disease"), 
              by = "disease") %>%
    select(d_v_a_id, disease, vaccine, activity, source)
  
  
  
  
  all_strata %<>%
    filter(disease == "Measles")
  
  
  
  
  # Repeat this process for each strata
  for (id in all_strata$d_v_a_id) {
    
    # Details of this strata
    strata = all_strata[d_v_a_id == id, ]
    
    # Display progress message to user
    msg = unlist(strata[, !"d_v_a_id"])
    message(paste0(" - ", paste(msg, collapse = ", ")))
    
    # Prepare for relative risk calculation
    dt = calculate_rr(strata)
    
    browser()
    
    # Append country-specific covariates (eg social-economic index) from GBD
    dt = merge_rr_covariates(dt)
    
    # Place within above function(s)...
    
    # Throw error if no valid relative risk results
    if (nrow(dt[rr > 0 & rr < 1]) == 0)
      stop("Error calculating relative risk")
    
    # ---- Fit a model to predict relative risk ----
    
    # @Austin: Why binomial when predictor is continuous?
    
    age_knots = table("age_knots") %>%
      filter(disease == strata$disease) %>%
      pull(age_knots) %>%
      eval_str()
    
    # Fit a binomial GLM for relative risk using GBD covariates
    fit = glm(formula = rr ~ haqi + sdi + year + mx + 
                splines::bs(age, knots = age_knots),
              data    = dt[rr > 0 & rr < 1],
              family  = "binomial")
    
    # Predict relative risk
    dt[, pred := predict(fit, dt)]
    
    # Transform these predictions (?? Not sure what's happening here ??)
    # For floating point precision
    dt[, pred_rr := ifelse(pred > 0, 1 / (1 + exp(-pred)), exp(pred) / (exp(pred) + 1))]
    
    # Remove covariates
    dt[, c("pred", "haqi", "sdi", "mx") := NULL]
    dt[, d_v_a_id := strata$d_v_a_id]
    
    browser()  
    
    # This function is the inverse of vimc_rr() ??
    # Is this only for non-VIMC countries?
    
    dt[, predict_averted := get_averted_deaths(deaths_obs, coverage, pred_rr)]
    
    # Store relative risk in list
    rr_list[[id]] = dt
    
    # ... and also save result to file
    # save_file(dt, o$pth$relative_risk, paste1("strata", id))
  }
  
  # ---- Combine all results ----
  
  message(" - Concatenating relative risk results")
  
  # Bind all stratas into single datatable
  rr_dt = rbindlist(rr_list, fill = TRUE)
  
  # Save relative risk calculations and predictions to file
  save_file(rr_dt, o$pth$relative_risk, "relative_risk")
  
  # ---- Plot diagnostics ----
  
  # Check flag
  if (o$plot_diagnostics) {
    
    message(" - Plotting diagnostics")
    
    browser()
    
    plot_strata_fit(rr_dt)
  }
}

# ---------------------------------------------------------
# Calculate relative risk for a given strata
# ---------------------------------------------------------
calculate_rr = function(strata) {
  
  if (strata$source == "vimc") use_table = "vimc_impact" # "vimc_impact" or "vimc_yov"
  if (strata$source == "gbd")  use_table = "gbd_estimates"
  
  # Total vaccine coverage by cohort
  total_coverage_dt = table("coverage") %>%
    left_join(y  = table("v_a"), 
              by = "v_a_id") %>%
    # Coverage is per vaccine, not per strata...
    filter(vaccine  == strata$vaccine, 
           activity == strata$activity, 
           year %in% o$analysis_years) %>%  # TODO: Move this 'year' filter to prepare.R
    # Calculate coverage by cohort...
    total_coverage()  # See coverage.R
  
  relative_risk_fn = get(paste1(strata$source, "rr"))
  
  browser()
  
  # Load either deaths_averted or deaths_disease for this strata
  strata_dt = table(use_table) %>%
    filter(d_v_a_id == strata$d_v_a_id) %>%
    # Append all-cause deaths...
    inner_join(y  = table("wpp_death"), 
               by = c("country", "year", "age")) %>%
    rename(deaths_allcause = death) %>%
    arrange(country, year, age) %>%
    # Append total coverage...
    inner_join(y  = total_coverage_dt, 
               by = c("country", "d_v_a_id", "year", "age")) %>%
    arrange(country, year, age) %>%
    # Calculate relative risk...
    relative_risk_fn()
  
  # Perform sanity checks on relative risk calculations
  # check_rr(strata_dt)
  
  return(strata_dt)
}

# ---------------------------------------------------------
# Calculate relative risk for VIMC disease
# ---------------------------------------------------------
vimc_rr = function(strata_dt) {
  
  # Relative risk tends to 1 as coverage tends to 1, could be negative
  #
  # Equation: rr = (o - (a * (1 - c) / c)) / (o + a)

  # Shorthand variable for readability
  o = strata_dt$deaths_allcause  # Deaths [o]bserved
  a = strata_dt$deaths_averted   # Deaths [a]verted from vaccination
  c = strata_dt$total_coverage   # Lifetime vaccine [c]overage
  
  # Calculate relative risk
  relative_risk_dt = strata_dt %>%
    mutate(rr = (o - (a * (1 - c) / c)) / (o + a)) %>%
    # Tidy up trivial values...
    mutate(rr = ifelse(total_coverage == 0, NA, rr), 
           rr = ifelse(is.infinite(rr), NA, rr))

  # browser()
  
  return(relative_risk_dt)
}

# ---------------------------------------------------------
# Calculate relative risk for GBD disease
# ---------------------------------------------------------
gbd_rr = function(x) {
  
  # Relative risk tends to 1 as coverage tends to 1, could be negative
  #
  # Equation: rr = (o - d + w * (1 - e)) / (o - d + w)
  #
  # where:
  #   o = deaths observed
  #   d = deaths attributable to this disease
  #   w = deaths estimated without a vaccine
  #   e = vaccine efficacy
  # 
  # Note that w - d equivalent to a (ie deaths averted from vimc_rr function)
  
  browser()
  
  strata = table("d_v_a") %>%
    filter(d_v_a_id == unique(x$d_v_a_id))
  
  # Vaccine efficacy
  strata_efficacy = efficacy %>%
    filter(disease == strata$disease, 
           vaccine == strata$vaccine) %>%
    pull(mean)
  
  # Number of deaths expected without vaccination ??
  out_dt <- copy(x)
  out_dt[coverage > 0, deaths_no := deaths_disease /
           (1 - strata_efficacy * coverage)]
  
  # Calculate relative risk
  out_dt[coverage > 0, rr :=
           (deaths_obs - deaths_disease + (1 - strata_efficacy) * deaths_no) /
           (deaths_obs - deaths_disease + deaths_no)]
  
  # Remove now-redundant parameter
  out_dt[, deaths_no := NULL]
  
  return(out_dt[])
}

# ---------------------------------------------------------
# Perform sanity checks on relative risk calculations
# ---------------------------------------------------------
check_rr = function(dt) {
  
  # Look for missing coverage where deaths averted are non-zero
  if (nrow(dt[coverage == 0 & deaths_averted > 0]) > 0) {
    
    # Proportion of cases
    prop <- round(nrow(dt[coverage == 0 & deaths_averted > 0]) /
                    nrow(dt[deaths_averted > 0]) * 100, 2)
    
    # Throw warning
    warning("Missing coverage in ", prop,
            "% of country-age-years with deaths averted")
  }
  
  # Check for non-sensical numbers
  if (any(range(dt[!is.na(rr)]$rr) < 0 | range(dt[!is.na(rr)]$rr) > 1)) {
    
    # Proportion of cases
    prop <- round(nrow(dt[coverage > 0 & (rr < 0 | rr > 1)]) /
                    nrow(dt) * 100, 2)
    
    # Throw warning
    warning("Over 1 or less than 0 mortality reduction in ", prop,
            "% of country-age-years")
  }
}

# ---------------------------------------------------------
# Append country-specific covariates (eg social-economic index) from GBD
# ---------------------------------------------------------
merge_rr_covariates = function(dt) {
  
  # NOTES:
  #  1) WWP stands for 'World Population Prospects'
  #     - nx := number of people
  #     - mx := mortality rate
  #     - fx := fertility rate
  #     - mig := migration rate
  #  2) Within the GBD (Global Burden of Disease) data
  #  a) HAQi stands for 'Healthcare Access and Quality index'
  #     See: https://www.healthdata.org/research-article/healthcare-access-and-quality-index-based-mortality-causes-amenable-personal-health
  #  b) SDI stands for 'Socio-demographic Index'
  #     See: https://www.healthdata.org/taxonomy/glossary/socio-demographic-index-sdi
  
  # All countries, years, and ages
  dt = expand_grid(
    country = unique(country_table$country),
    age     = o$data_ages,
    year    = 2000 : 2095) %>% # ?? Why 2095?
    as.data.table() %>%
    merge(dt, by = c("country", "age", "year"), all.x = T)
  
  # Load WPP data
  wpp_pop = table("wpp_pop")
  
  browser() # wpp_pop already summarised over gender
  
  # Add mortality
  wpp_input = wpp_pop
  mx_dt <- wpp_input[, .(mx = mean(mx)), by = .(country, year, age)]
  dt <- merge(dt, mx_dt, by = c("country", "age", "year"), all.x = T)
  
  browser()
  
  # Add GBD covariates (SDI and HAQi)
  dt <- merge(dt, gbd_covariates, by = c("country", "year"), all.x = T)
  
  return(dt)
}

# ---------------------------------------------------------
# Extract averted deaths from relative risk equation
# ---------------------------------------------------------
get_averted_deaths = function(deaths_obs, coverage, rr) {
  
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

