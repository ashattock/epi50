###########################################################
# IMPUTE
#
# Impute deaths averted for VIMC pathogens for countries
# not modelled by VIMC.
#
###########################################################

# ---------------------------------------------------------
# Parent function for imputing missing VIMC countries
# ---------------------------------------------------------
run_impute = function() {
  
  # Only continue if specified by do_step
  if (!is.element(4, o$do_step)) return()
  # if (!is.element(7, o$do_step)) return()
  
  message("* Running country imputation")
  
  # TODO: Extend to allow fitting to temporally extrapolated data 
  #       => useful messages on determinants of vaccine impact over time 
  
  # TODO: Repeat process for DALYs
  
  # Create set of models to evaluate
  list[models, covars] = get_models()
  
  # Load response variable (impact per FVP)
  target_dt = get_impute_data()
  
  # Return out if no training data identified
  if (nrow(target_dt) == 0)
    return()
  
  # ---- Perform regression ----
  
  # TEMP: Use basic or full imputation method
  #
  # OPTIONS:
  #  basic_impute - IA2030 method using GBD covariates
  #  perform_impute1 - Helen's time series regression method
  #  perform_impute2 - Same as Helen's method, refactored code
  method = "perform_impute2"
  
  # TEMP: Ignoring problematic cases for now
  ignore = c(7, 11, 14)
  
  # Call country imputation function
  predict_dt = table("d_v_a") %>%
    filter(source == "vimc") %>%
    # TEMP: Ignoring problematic cases for now...
    filter(!d_v_a_id %in% ignore) %>%  
    # Apply geographical imputation model...
    pull(d_v_a_id) %>%
    lapply(FUN = get(method), 
           models = models, 
           covars = covars, 
           target = target_dt) %>%
    rbindlist() %>%
    select(d_v_a_id, country, year, impact_impute)
  
  # Apply imputations where needed
  impute_dt = target_dt %>%
    select(d_v_a_id, country, year, fvps_cum, impact_cum) %>%
    # Append predictions...
    left_join(y  = predict_dt, 
              by = c("d_v_a_id", "country", "year")) %>%
    # Apply predictions where imputation is needed...
    mutate(impact = ifelse(
      test = is.na(impact_cum),
      yes  = impact_impute,
      no   = impact_cum)) %>%
    select(d_v_a_id, country, year, 
           fvps = fvps_cum, impact) %>%
    # Assume any missing values are zero impact...
    replace_na(list(impact = 0)) %>%
    # TEMP: Bound impact above by 1...
    mutate(impact = pmin(impact, 1))
    
  # Save imputed results to file
  save_rds(impute_dt, "impute", "impute_result")
  
  # ---- Plot results ----
  
  # Plot predictor-response relationships
  # plot_covariates()
  
  # Plot model choice
  # plot_model_choice()
  
  # Plot predictive performance for each country
  # plot_impute_perform()
  
  # Plot fitted imputation model for each country
  # plot_impute_fit()
  
  # Plot imputation quality of fit for all countries
  # plot_impute_quality()
  
  # Plot tornado plots of predictors by d_v_a
  # plot_tornado_d_v_a()
  
  # Plot tornado plots of predictors by region
  # plot_tornado_region()
}

# ---------------------------------------------------------
# Create set of regression models to evaluate
# ---------------------------------------------------------
get_models = function() {
  
  # List of available covariates
  covars = list(
    c0 = "log(coverage)",
    c1 = "log(coverage_minus_1)", 
    c2 = "log(coverage_minus_2)", 
    c3 = "log(coverage_minus_3)", 
    c4 = "log(coverage_minus_4)", 
    p  = "pop_0to14", 
    g  = "gini", 
    h  = "HDI", 
    a  = "attended_births")
  
  # Define models (using shorthand covariate references)
  models = list(
    m1  = "c0", 
    m2  = "c0 + c1", 
    m3  = "c0 + c1 + c2", 
    m4  = "c0 + c1 + c2 + c3", 
    m5  = "c0 + c1 + c2 + c3 + c4", 
    m6  = "c0 + c1 + c2", 
    m7  = "c0 + c1 + c2 + c3 + c4 + p + g", 
    m8  = "c0 + c1 + c2 + c3 + c4 + p + g + a", 
    m9  = "c0 + c1 + c2 + c3 + c4 + p + g + a", 
    m10 = "c0 + c1 + c2 + c3 + h + p + g", 
    m11 = "c0 + c1 + c2 + h + p + g", 
    m12 = "c0 + c1 + h + p + g", 
    m13 = "c0 + h + p + g")
  
  return(list(models, covars))
}

# ---------------------------------------------------------
# Load/calculate target (impact per FVP) for modelled pathogens
# ---------------------------------------------------------
get_impute_data = function() {
  
  # Population size of each country over time
  pop_dt = table("wpp_pop") %>%
    lazy_dt() %>%
    group_by(country, year) %>%
    summarise(pop = sum(pop)) %>%
    ungroup() %>%
    as.data.table()
  
  # Wrangle VIMC impact estimates
  impact_dt = table("vimc_estimates") %>%
    lazy_dt() %>%
    # Sum impact over age...
    group_by(d_v_a_id, country, year) %>%
    summarise(impact_abs = sum(deaths_averted)) %>%
    ungroup() %>%
    mutate(impact_abs = pmax(impact_abs, 0)) %>%
    # Scale results to per capita...
    left_join(y  = pop_dt, 
              by = c("country", "year")) %>%
    mutate(impact_rel = impact_abs / pop) %>%
    select(-pop) %>%
    # Cumulative sum impact...
    arrange(d_v_a_id, country, year) %>%
    group_by(d_v_a_id, country) %>%
    mutate(impact_cum = cumsum(impact_rel)) %>%
    ungroup() %>%
    as.data.table()
  
  # Extract FVPs
  fvps_dt = table("coverage") %>%
    lazy_dt() %>%
    # Only impute pathogens and years for which we've VIMC estimates...
    filter(d_v_a_id %in% unique(impact_dt$d_v_a_id), 
           year     %in% unique(impact_dt$year)) %>%
    # Summarise over age...
    group_by(d_v_a_id, country, year) %>%
    summarise(fvps_abs = sum(fvps)) %>%
    ungroup() %>%
    # Scale results to per capita...
    left_join(y  = pop_dt, 
              by = c("country", "year")) %>%
    mutate(fvps_rel = fvps_abs / pop) %>%
    select(-pop) %>%
    # Cumulative sum FVPs...
    arrange(d_v_a_id, country, year) %>%
    group_by(d_v_a_id, country) %>%
    mutate(fvps_cum = cumsum(fvps_rel)) %>%
    ungroup() %>%
    as.data.table()
  
  # Combine into single datatable
  target_dt = fvps_dt %>%
    left_join(y  = impact_dt, 
              by = c("d_v_a_id", "country", "year")) %>%
    # Impact per FVP...
    mutate(target = impact_cum / fvps_cum)
  
  # Save this datatable to file for plotting purposes
  save_rds(target_dt, "impute", "target")
  
  # Throw a warning if no target data identified
  if (nrow(target_dt) == 0)
    warning("No imputation training data identified")
  
  return(target_dt)
}

# ---------------------------------------------------------
# Perform imputation
# ---------------------------------------------------------
perform_impute2 = function(d_v_a_id, models, covars, target) {
  
  # Stepwise regression
  # TODO: Update to lasso regularisation for optimal predictor selection
  
  # Extract name of this d-v-a
  d_v_a_name = table("d_v_a") %>%
    filter(d_v_a_id == !!d_v_a_id) %>%
    pull(d_v_a_name)
  
  # Display progress message to user
  message(" - ", d_v_a_name)
  
  # ---- Append covariates ----
  
  # Full list of covariates used by specified models
  retain_covars = which_covars(models, covars)
  
  # Summarise vaccination coverage by country and year
  coverage_dt = table("coverage") %>%
    lazy_dt() %>%
    filter(d_v_a_id == !!d_v_a_id) %>%
    # Recalculate for whole population...
    group_by(country, year) %>%
    summarise(fvps   = sum(fvps),
              cohort = sum(cohort)) %>%
    ungroup() %>%
    mutate(coverage = fvps / cohort) %>%
    as.data.table()
  
  # TODO: Are GBD covariates still used/needed?
  
  # Create time-series tibble with all covariates
  target_ts = target %>%
    lazy_dt() %>%
    # Data for this d-v-a...
    filter(d_v_a_id == !!d_v_a_id) %>%
    select(-d_v_a_id) %>%
    # Append vaccination coverage...
    full_join(y = coverage_dt,  
              by = c("country", "year")) %>%
    # Append covariates: GBD...
    full_join(y  = table("gbd_covariates"),
              by = c("country", "year")) %>%
    # Append covariates: UNICEF...
    # full_join(y  = table("unicef_covariates"),
    #           by = c("country", "year"),
    #           relationship = "many-to-many") %>%
    # Append covariates: Gapminder...
    full_join(y  = table("gapminder_covariates"),  # TODO: multiple entries for COD(Congo, Kinshasa)
              by = c("country", "year"), 
              relationship = "many-to-many") %>%
    # Append period...
    # HCJ - have I interpreted this correctly? period is a potential predictor?
    left_join(y  = get_period(), 
              by = "year") %>%
    # Summarise to single row for each country per year...
    arrange(country, year) %>%
    group_by(country, year) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    # Create dummy variables for historic coverage (NA prior to our data)
    mutate(coverage_minus_1 = lag(coverage, 1),
           coverage_minus_2 = lag(coverage, 2),
           coverage_minus_3 = lag(coverage, 3),
           coverage_minus_4 = lag(coverage, 4),
           coverage_minus_5 = lag(coverage, 5),
           coverage_minus_6 = lag(coverage, 6),
           coverage_minus_7 = lag(coverage, 7),
           coverage_minus_8 = lag(coverage, 8)) %>%
    # Create dummy variables for historic health_spending (NA prior to our data)
    mutate(health_spending_minus_1 = lag(health_spending, 1),
           health_spending_minus_2 = lag(health_spending, 2),
           health_spending_minus_3 = lag(health_spending, 3),
           health_spending_minus_4 = lag(health_spending, 4),
           health_spending_minus_5 = lag(health_spending, 5),
           health_spending_minus_6 = lag(health_spending, 6),
           health_spending_minus_7 = lag(health_spending, 7),
           health_spending_minus_8 = lag(health_spending, 8)) %>%
    # Only retain covariates defined in models...
    select(country, year, target, 
           all_of(retain_covars)) %>%
    # Execute dtplyr...
    as.data.table() %>%
    # Convert to time-series tibble...
    as_tsibble(index = year, 
               key   = country) 
  
  # ---- Evaluate all user-defined models ----
  
  message("  > Evaluating models")
  
  # Subset training data (which we have impact estimates for)
  data_ts = target_ts %>%
    # Remove zeros to allow for log transformation...
    filter(target > 0) %>%
    # Remove country if insufficient data points for fitting...
    group_by(country) %>%
    filter(n() >= o$min_data_requirement) %>% 
    ungroup()
  
  # Evaluate all models in parallel
  if (o$parallel$impute)
    model_list = mclapply(
      X   = names(models),
      FUN = evaluate_model, 
      models = models, 
      covars = covars, 
      data   = data_ts,
      mc.cores = o$n_cores,
      mc.preschedule = FALSE)
  
  # Evaluate all models consecutively
  if (!o$parallel$impute)
    model_list = lapply(
      X   = names(models),
      FUN = evaluate_model, 
      models = models, 
      covars = covars, 
      data   = data_ts)
  
  # ---- Model selection ----
  
  message("  > Model selection")

  # For each country, select the model with the best AICc
  model_choice = model_list %>%
    lapply(report) %>%
    rbindlist() %>%
    # Remove null models...
    filter(!is.infinite(AICc)) %>% 
    # Retain only best fit model (if equal, keep the first)...
    group_by(country) %>%
    slice_min(AICc, with_ties = FALSE) %>%
    unique() %>%
    # Reappend best model...
    left_join(y  = rbindlist(model_list), 
              by = c("country", "model_id")) %>%
    # Reduce down to keep only model and AICc... 
    select(country, model_id, tslm, AICc) %>%
    mutate(d_v_a_id = d_v_a_id, 
           .before = 1) %>%
    # Convert to mable class...
    as_mable(key   = country, 
             model = tslm) %>%
    suppressWarnings()
  
  # Extract parameters of best fitting model for each country
  model_fit = tidy(model_choice) %>%
    select(d_v_a_id, country, model_id, term, 
           estimate, std.error, p.value) %>%
    as.data.table()
  
  # Evaluate models on the training data
  #
  # @HCJ - why is this needed if model 13 does the predictions?
  # So we can see how well the 'best' models fits the target?
  # model_eval = augment(model_choice) %>%
  #   select(country, year, target, 
  #          prediction = .fitted) %>%
  #   as.data.table()
  
  # ---- Predictions ----
  
  message("  > Model predictions")
  
  # Prediction are performed at regional level
  region_dt = table("country") %>%
    select(country, region)
  
  # TODO: Determine results by region AND period
  
  # TEMP: Use model 13 (most commonly selected) to predict unseen countries, for now
  # TODO: Generalise: allow model selection for imputed country by e.g. region. For now, model 13 works well in general
  # TODO: Choose most appropriate method for selecting coefficients e.g. nearest neighbours
  model_13_region = model_list[[13]] %>%
    tidy() %>%
    # Append region...
    left_join(y  = region_dt, 
              by = "country") %>%
    # Median coefficient by region to avoid outliers...
    group_by(region, term) %>%
    summarise(estimate = median(estimate, na.rm = TRUE)) %>%
    ungroup() %>%
    # Spread to wide format...
    pivot_wider(
      names_from  = term,
      names_glue  = "{term}_coefficient",
      values_from = estimate) %>%
    as.data.table()
  
  # evaluate_predictions()
  
  quote = function(x, q = '"') 
    paste0(q, x, q)
  
  predict_covars = models[13] %>%
    interpret_covars(covars) %>%
    str_remove_all(" ") %>%
    str_split("\\+") %>%
    pluck(1)
  
  predict_str = predict_covars %>%
    paste1("coefficient") %>%
    quote("`") %>%
    paste(predict_covars, sep = " * ") %>%
    paste(collapse = " + ")
  
  predict_fn   = paste0("prediction = exp(", predict_str, ")")
  predict_eval = paste0("mutate(predictors, ", predict_fn, ")")
  
  # Select countries for imputation
  predictors = target_ts %>% 
    inner_join(y  = region_dt, 
               by = "country") %>% 
    left_join(y  = model_13_region, 
              by = "region")
  
  # Impute target values using coefficients from WHO region
  predict_dt = eval_str(predict_eval) %>%
    select(country, year, prediction) %>% 
    as.data.table()
  
  # ---- Format output ----
  
  # Apply predictions to impute missing impact estimates
  result_dt = target %>%
    filter(d_v_a_id == !!d_v_a_id) %>%
    # Append predictions...
    left_join(y  = predict_dt, 
              by = c("country", "year")) %>%
    select(d_v_a_id, country, year, fvps_cum, 
           impact_cum, target, prediction) %>%
    # Multiply through to obtain cumulative impact over time...
    mutate(impact_impute = fvps_cum * prediction, 
           .after = impact_cum)
  
  # Also format predictors for use in plotting
  data_dt = data_ts %>%
    mutate(d_v_a_id = d_v_a_id, 
           .before = 1) %>%
    as.data.table()
  
  # Store the data used, fitted model, and result
  fit = list(
    model  = model_choice,  # NOTE: Only for non-imputed
    report = model_fit,
    result = result_dt, 
    data   = data_dt)
  
  # Save to file
  save_rds(fit, "impute", "impute", d_v_a_id)
  
  return(result_dt)
}

# ---------------------------------------------------------
# Evaluate given user-specified model
# ---------------------------------------------------------
evaluate_model = function(id, models, covars, data) {
  
  # Interpret covariate references
  model_str = interpret_covars(models[[id]], covars)
  
  # Construct full model string to be evaluated
  model_fn   = paste0("TSLM(log(target) ~ ", model_str, ")")
  model_eval = paste0("model(data, tslm = ", model_fn, ")")
  
  # Evaluate model and append model reference
  model_mab = eval_str(model_eval) %>%
    mutate(model_id = id, 
           .before  = tslm) %>%
    suppressWarnings()
  
  return(model_mab)
}

# ---------------------------------------------------------
# xxxxxxxxxxx
# ---------------------------------------------------------
evaluate_predictions = function() {
  
  browser() 
}

# ---------------------------------------------------------
# Evaluate given user-specified model
# ---------------------------------------------------------
interpret_covars = function(model, covars) {
  
  # Interpret shorthand references
  for (covar in names(covars))
    model = str_replace_all(
      string  = model, 
      pattern = paste0("\\b", covar, "\\b"), 
      replacement = covars[[covar]])
  
  return(model)
}

# ---------------------------------------------------------
# Extract covariates used by this model(s)
# ---------------------------------------------------------
which_covars = function(models, covars) {
  
  # Shorthand covariates used in specified models
  covars_used = unlist(models) %>%
    paste(collapse = " + ") %>%
    str_split_1(pattern = " \\+ ") %>%
    unique()
  
  # Associated names of covariate columns 
  covar_names = covars[covars_used] %>%
    unlist(covars) %>%
    str_remove("^.*\\(+") %>%
    str_remove("\\)+$")
  
  return(covar_names)
}

# ---------------------------------------------------------
# Easily convert between year and period
# ---------------------------------------------------------
get_period = function() {
  
  # Indices of period change
  year_idx = seq(
    from = o$period_length, 
    to   = length(o$years), 
    by   = o$period_length) + 1
  
  # Format into full year-period datatable
  period_dt = tibble(year = o$years[year_idx]) %>%
    mutate(period = 1 : n()) %>%
    full_join(y  = tibble(year = o$year), 
              by = "year") %>%
    arrange(year) %>%
    fill(period, .direction = "updown") %>%
    as.data.table()
  
  return(period_dt)
}

