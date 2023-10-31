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
  
  target_dt = get_target()
  
  # plot_target(target_dt)
  
  # Calculate relative-risk for all d-v-a combinations 
  impute_dt = table("d_v_a") %>%
    pluck("d_v_a_id") %>%
    lapply(do_impute, target = target_dt) %>%
    rbindlist()
  
  browser()
  
  # predict_dt[, predict_averted := calculate_averted_deaths(deaths_obs, coverage, pred_rr)]
  
  # Save relative risk calculations and predictions to file
  save_rds(impute_dt, "impute", "impute")
  
  browser()
  
  # Check flag
  if (o$plot_diagnostics) {
    
    message(" - Plotting diagnostics")
    
    browser()
    
    # plot_covariates()
    plot_impute_fit()
    
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
  
  plot_dt = data_dt %>%
    pivot_longer(cols = -target, 
                 names_to = "covariate") %>%
    arrange(covariate, target) %>%
    as.data.table()
  
  g1 = ggplot(plot_dt) +
    aes(x = target, y = value, colour = covariate) +
    geom_point(alpha = 0.2) +
    facet_wrap(~covariate, scales = "free_y")
  
  # ---- Normalise ----
  
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
  
  if (d_v_a_id == 4)
    browser()
  
  # Fit a binomial GLM for relative risk using all covariates
  fit_model = glm(
    formula = target ~ n_years + coverage + sdi + haqi + imr, 
    data    = norm_data_dt)
  
  # Use fitted model to predict 
  result_dt = target_dt %>%
    select(country, d_v_a_id, year) %>%
    # Predict impact per FVP...
    cbind(norm_pred_dt) %>%
    mutate(predict = predict(fit_model, .)) %>%
    # Remove predictors...
    select(country, d_v_a_id, year, target, predict) %>%
    # Back-transform target and prediction...
    mutate(target  = retransform_fn(target,  a, b), 
           predict = retransform_fn(predict, a, b))
    
  plot_dt = result_dt %>%
    filter(!is.na(target))
  
  g2 = ggplot(plot_dt, aes(x = target, y = predict)) +
    geom_point(alpha = 0.5) +
    geom_abline(colour = "red")
  
  fit = list(
    model  = fit_model, 
    data   = data_dt, 
    result = result_dt)
  
  # Save relative risk calculations and predictions to file
  save_rds(fit, "impute", "impute", d_v_a_id)
  
  return(result_dt)
}

# ---------------------------------------------------------
# Load/calculate target variable: impact per FVP
# ---------------------------------------------------------
get_target = function() {
  
  # First load population size of each country over time
  # pop_dt = table("wpp_pop") %>%
  #   group_by(country, year) %>%
  #   summarise(pop = sum(pop)) %>%
  #   ungroup() %>%
  #   as.data.table()
  
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
    # Cumulative relative to 100k people...
    # left_join(y  = pop_dt, 
    #           by = c("country", "year")) %>%
    # mutate(impact_rel = o$per_person * impact_cum / pop) %>%
    # select(country, d_v_a_id, year, impact = impact_rel) %>%
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
    # Cumulative relative to 100k people...
    # left_join(y  = pop_dt, 
    #           by = c("country", "year")) %>%
    # mutate(fvps_rel = o$per_person * fvps_cum / pop) %>%
    # select(country, v_a_id, year, coverage, fvps = fvps_rel) %>%
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

