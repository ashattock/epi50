###########################################################
# RELATIVE RISK
#
# Calculate relative risk given coverage details and deaths 
# averted. Then impute relative risk for missing countries
# using GBD covariates.
#
# NOTES:
# 1) HAQi stands for 'Healthcare Access and Quality index'
#    See: https://www.healthdata.org/research-article/healthcare-access-and-quality-index-based-mortality-causes-amenable-personal-health
# 2) SDI stands for 'Socio-demographic Index'
#    See: https://www.healthdata.org/taxonomy/glossary/socio-demographic-index-sdi
#
###########################################################

# Variables:
#  - deaths_allcause: GBD-estimated all cause number of deaths
#  - deaths_disease: GBD-estimated deaths attributable to a disease
#  - deaths_averted: VIMC model-estimated deaths averted
#  - predict_averted: predicted deaths averted using statistical model

# ---------------------------------------------------------
# Impute relative risk for a each d-v-a (aka strata)
# ---------------------------------------------------------
run_relative_risk = function() {
  
  # Only continue if specified by do_step
  if (!is.element(2, o$do_step)) return()
  
  message("* Calculating relative risk")
  
  # All strata to run
  all_strata = table("d_v_a") %>%
    left_join(y = table("disease"), 
              by = "disease") %>%
    select(d_v_a_id, disease, vaccine, activity, source)
  
  
  
  
  
  all_strata %<>%
    filter(disease == "Measles") # Dip Measles
  
  
  
  
  # Initiate list for predicted relative risk results
  rr_list = list()
  
  # Repeat this process for each strata
  for (id in all_strata$d_v_a_id) {
    
    # Details of this strata
    strata = all_strata[d_v_a_id == id, ]
    
    # Display progress message to user
    msg = unlist(strata[, !"d_v_a_id"])
    message(paste0(" - ", paste(msg, collapse = ", ")))
    
    # ---- Total coverage ----
    
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
    
    # ---- Calculate relative risk ----
    
    # Load table and rr function relevant for this strata source
    estimates_table  = paste1(strata$source, "estimates") %>% table()
    relative_risk_fn = paste1(strata$source, "rr") %>% get()
    
    # Load either deaths_averted or deaths_disease for this strata
    rr_dt = estimates_table %>%
      filter(d_v_a_id == strata$d_v_a_id) %>%
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
    
    # ---- Append covariates ----
    
    # GBD-derived covariates by country and age (see above notes)
    covariate_dt = table("gbd_covariates")
    
    # TODO: Is death rate actually helping??
    
    # All-cause death rate by country, year, and age
    # death_rate_dt = table("wpp_pop") %>%
    #   inner_join(y  = table("wpp_death"),
    #              by = c("country", "year", "age")) %>%
    #   mutate(death_rate = death / pop) %>%
    #   select(country, year, age, death_rate)
    
    # Data used to fit statistical model
    rr_data_dt = rr_dt %>%
      # Append covariates...
      inner_join(y  = covariate_dt,
                 by = c("country", "year")) %>%
      # Extract what we'll use for rr predictions...
      filter(rr > 0 & rr < 1) %>%
      select(year, age, sdi, haqi, rr)
    
    # All points to predict
    rr_pred_dt = rr_dt %>%
      # Append covariates...
      inner_join(y  = covariate_dt,
                 by = c("country", "year")) %>%
      # Extract only points we want to predict for...
      filter(total_coverage > 0) %>%
      select(country, year, age, sdi, haqi, rr)
    
    # plot_dt = rr_data_dt %>%
    #   select(sdi, haqi, rr) %>%
    #   mutate(rr_type = ifelse(is.na(rr), "pred", "data"), 
    #          rr      = ifelse(is.na(rr), 0, rr)) %>%
    #   pivot_longer(cols = -c(rr, rr_type), 
    #                names_to = "var") %>%
    #   as.data.table()
    # 
    # g = ggplot(plot_dt) +
    #   aes(x = value, y = rr) + 
    #   geom_bin_2d() +
    #   facet_grid(rr_type~var) 
    
    # ---- Normalise ----
    
    transform_rr = function(mat, t0, t1)
      dt = as.data.table(t((mat - t0) / (t1 - t0)))
    
    data_mat = t(as.matrix(rr_data_dt))
    pred_mat = t(as.matrix(rr_pred_dt[, !"country"]))
    
    t0 = rowMins(data_mat)
    t1 = rowMaxs(data_mat)
    
    norm_data_dt = transform_rr(data_mat, t0, t1)
    norm_pred_dt = transform_rr(pred_mat, t0, t1) %>%
      cbind(country = rr_pred_dt$country)
    
    # ---- Fit a model to predict relative risk ----
    
    # # Parse spline age knots for this pathogen
    # age_knots = table("age_knots") %>%
    #   filter(disease == strata$disease) %>%
    #   pull(age_knots) %>%
    #   eval_str()
    # 
    # # Fit a spline spanning all ages
    # age_spline_dt = bs(rr_fit_dt$age, knots = age_knots) %>%
    #   as_named_dt(paste1("age", 1 : ncol(.)))
    
    # Fit a binomial GLM for relative risk using all covariates
    rr_fit = glm(formula = rr ~ year + age + sdi + haqi,
                 data    = norm_data_dt) #,
    # family  = "binomial")
    
    rr_predict_dt = norm_pred_dt %>%
      mutate(rr_predict = predict(rr_fit, .), 
             rr_predict2 = ifelse(rr_predict > 0, 
                                  1 / (1 + exp(-rr_predict)), 
                                  exp(rr_predict) / (exp(rr_predict) + 1)))
    
    plot_dt = rr_predict_dt %>%
      # filter(rr > 0 & rr < 1) %>%
      pivot_longer(cols = c(rr_predict, rr_predict2),
                   names_to  = "predict_type",
                   values_to = "prediction")

    g = ggplot(plot_dt, aes(x = rr, y = prediction)) +
      geom_point() + # aes(colour = year)) +
      facet_wrap(~predict_type, scales = "free_y") +
      xlim(0, 1) +
      ylim(0, 1)
    
    browser()
    
    # ---- xxxxxx ----
    
    browser()
    
    # predict_dt[, predict_averted := calculate_averted_deaths(deaths_obs, coverage, pred_rr)]
    
    rr_list[[id]] = predict_dt
    
    browser()
  }
  
  # ---- Combine all results ----
  
  message(" - Concatenating relative risk results")
  
  browser() # Why do we still need to fill?
  
  # Bind all stratas into single datatable
  rr_dt = rbindlist(rr_list, fill = TRUE)
  
  # Perform sanity checks on relative risk calculations
  check_rr(rr_dt)
  
  # Save relative risk calculations and predictions to file
  save_rds(rr_dt, "relative_risk", "relative_risk")
  
  # ---- Plot diagnostics ----
  
  # Check flag
  if (o$plot_diagnostics) {
    
    message(" - Plotting diagnostics")
    
    browser()
    
    plot_strata_fit(rr_dt)
  }
}

# ---------------------------------------------------------
# Calculate relative risk for VIMC disease
# ---------------------------------------------------------
vimc_rr = function(strata_dt) {
  
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
# Calculate relative risk for GBD disease
# ---------------------------------------------------------
gbd_rr = function(strata_dt) {
  
  # Append vaccine efficacy
  strata_dt %<>%
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_id") %>%
    # Append vaccine efficacy...
    left_join(y  = table("gbd_efficacy"), 
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

