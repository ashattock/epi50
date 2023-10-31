###########################################################
# IMPUTE
#
# Impute deaths averted for VIMC pathogens for countries
# not modelled by VIMC. Uses several GBD covariates, amongst
# Other predictors.
#
# NOTES:
# 1) HAQi stands for 'Healthcare Access and Quality index'
#    See: https://www.healthdata.org/research-article/healthcare-access-and-quality-index-based-mortality-causes-amenable-personal-health
# 2) SDI stands for 'Socio-demographic Index'
#    See: https://www.healthdata.org/taxonomy/glossary/socio-demographic-index-sdi
#
###########################################################

# ---------------------------------------------------------
# Parent function for imputing missing VIMC countries
# ---------------------------------------------------------
run_impute = function() {
  
  # Only continue if specified by do_step
  if (!is.element(2, o$do_step)) return()
  
  message("* Running country imputation")
  
  # Load target to fit to: impact per FVP
  target_dt = get_target()
  
  # Calculate relative-risk for all d-v-a combinations 
  impute_dt = table("d_v_a") %>%
    pluck("d_v_a_id") %>%
    lapply(do_impute, target = target_dt) %>%
    rbindlist()
  
  # predict_dt[, predict_averted := calculate_averted_deaths(deaths_obs, coverage, pred_rr)]
  
  # Save relative risk calculations and predictions to file
  save_rds(impute_dt, "impute", "impute_result")
  
  # Check flag
  if (o$plot_diagnostics) {
    
    message(" - Plotting diagnostics")
    
    plot_target()
    
    plot_covariates()
    
    plot_impute_fit()
    
    browser()
    
    # plot_impute_countries()
  }
}

# ---------------------------------------------------------
# Perform imputation
# ---------------------------------------------------------
do_impute = function(d_v_a_id, target) {
  
  # Details of this d_v_a
  d_v_a = table("d_v_a") %>%
    filter(d_v_a_id == !!d_v_a_id) %>%
    left_join(y  = table("disease"), 
              by = "disease") %>%
    select(disease, vaccine, activity, source)
  
  # Display progress message to user
  message(paste0(" - ", paste(d_v_a, collapse = ", ")))
  
  # ---- Append covariates ----
  
  # All-cause death rate by country, year, and age
  infant_mortality_dt = table("wpp_pop") %>%
    filter(age == 0) %>%  # TODO: Instead filter by table("wiise_vaccine") age group
    inner_join(y  = table("wpp_death"),
               by = c("country", "year", "age")) %>%
    # Calculate infant mortality rate...
    mutate(imr = death / pop) %>%
    select(country, year, imr)
  
  # Append covariates to target
  target_dt = target %>%
    filter(d_v_a_id == !!d_v_a_id) %>%
    inner_join(y  = table("gbd_covariates"),
               by = c("country", "year")) %>%
    inner_join(y  = infant_mortality_dt,
               by = c("country", "year")) %>%
    # Calculate n years of estimates...
    mutate(n_years = 1) %>%
    group_by(country) %>%
    mutate(n_years = cumsum(n_years)) %>%
    ungroup() %>%
    as.data.table()
  
  # Data used to fit statistical model
  data_dt = target_dt %>%
    filter(!is.na(target)) %>%
    select(target, n_years, coverage, sdi, haqi, imr) %>%
    # Remove target outliers for better normalisation...
    mutate(lower = mean(target) - 3 * sd(target), 
           upper = mean(target) + 3 * sd(target), 
           outlier = target < lower | target > upper) %>%
    filter(outlier == FALSE) %>%
    select(-outlier, -lower, -upper)
  
  # Values to predict for (including data used for fitting)
  pred_dt = target_dt %>%
    select(all_of(names(data_dt)))
  
  # ---- Check for trivial case ----
  
  # Return out if no data available
  if (nrow(data_dt) == 0) {
    
    message(" !! Insufficient data for imputation !!")
    
    # We'll store a blank datatable in this case
    fit = list(data = data_dt)
    
    # Save to file
    save_rds(fit, "impute", "impute", d_v_a_id)
    
    return()
  }
  
  # ---- Normalise predictors and response ----
  
  transform_fn = function(x, a, b)
    y = t((x - a) / (b - a)) %>% as.data.table()
  
  retransform_fn = function(y, a, b)
    x = y * (b["target"] - a["target"]) + a["target"]
  
  data_mat = t(as.matrix(data_dt))
  pred_mat = t(as.matrix(pred_dt))
  
  a = rowMins(data_mat)
  b = rowMaxs(data_mat)
  
  norm_data_dt = transform_fn(data_mat, a, b)
  norm_pred_dt = transform_fn(pred_mat, a, b)
  
  # ---- Fit a model to predict impact per FVP ----
  
  # Fit a binomial GLM for relative risk using all covariates
  fit_model = glm(
    formula = target ~ n_years + coverage + sdi + haqi + imr, 
    data    = norm_data_dt)
  
  # Use fitted model to predict 
  result_dt = target_dt %>%
    select(country, d_v_a_id, year) %>%
    # Predict impact per FVP...
    cbind(norm_pred_dt) %>%
    mutate(predict = predict(fit_model, .), 
           predict = pmax(predict, 0)) %>%  # Do not predict negative
    # Remove predictors...
    select(country, d_v_a_id, year, target, predict) %>%
    # Back-transform target and prediction...
    mutate(target  = retransform_fn(target,  a, b), 
           predict = retransform_fn(predict, a, b))
  
  fit = list(
    model  = fit_model, 
    data   = norm_data_dt, 
    result = result_dt)
  
  # Save to file
  save_rds(fit, "impute", "impute", d_v_a_id)
  
  return(result_dt)
}

# ---------------------------------------------------------
# Load/calculate target variable: impact per FVP
# ---------------------------------------------------------
get_target = function() {
  
  # Wrangle VIMC impact estimates
  impact_dt = table("vimc_estimates") %>%
    # Sum impact over age...
    group_by(country, d_v_a_id, year) %>%
    summarise(impact_abs = sum(deaths_averted)) %>%
    ungroup() %>%
    filter(impact_abs > 0) %>%
    # Cumulative sum impact...
    arrange(country, d_v_a_id, year) %>%
    group_by(country, d_v_a_id) %>%
    mutate(impact_cum = cumsum(impact_abs)) %>%
    ungroup() %>%
    # Append d_v_a details...
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_id") %>%
    as.data.table()
  
  # Extract coverage
  coverage_dt = table("coverage") %>%
    # Summarise over age...
    group_by(country, v_a_id, year) %>%
    summarise(coverage = mean(coverage), 
              fvps_abs = sum(fvps)) %>%
    ungroup() %>%
    # Cumulative sum FVPs...
    arrange(country, v_a_id, year) %>%
    group_by(country, v_a_id) %>%
    mutate(fvps_cum = cumsum(fvps_abs)) %>%
    ungroup() %>%
    # Append v_a details...
    left_join(y  = table("v_a"), 
              by = "v_a_id") %>%
    filter(vaccine %in% unique(impact_dt$vaccine)) %>%
    left_join(y  = table("d_v_a"), 
              by = c("vaccine", "activity")) %>%
    select(country, d_v_a_id, year, coverage, 
           fvps_abs, fvps_cum) %>%
    as.data.table()
  
  # Combine into single datatable
  target_dt = coverage_dt %>%
    left_join(y  = impact_dt, 
              by = c("country", "d_v_a_id", "year")) %>%
    select(-disease, -vaccine, -activity) %>%
    # Impact per FVP...
    mutate(target = impact_cum / fvps_cum) # %>%
    # select(country, d_v_a_id, year, target, coverage)
  
  # Save this datatable to file for plotting purposes
  save_rds(target_dt, "impute", "target")
  
  return(target_dt)
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
  
  browser() # TODO: Needs to be redone
  
  # Estimate deaths averted given relative risk
  #
  # NOTE: Equation derived by solving vimc_rr for a
  averted_deaths <- deaths_obs * (coverage * (1 - rr) / 
                                    (1 - coverage * (1 - rr)))
  
  return(averted_deaths)
}

