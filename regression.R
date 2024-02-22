###########################################################
# REGRESSION
#
# Fit a user-defined series of models to cumulative impact 
# per FVP. This process has two use cases:
#  1) Impute impact for countries not modelled by VIMC
#  2) To infer drivers of impact (on complete set of estimates)
#
###########################################################

# ---------------------------------------------------------
# Parent function for regression modelling
# ---------------------------------------------------------
run_regression = function(case) {
  
  # Only continue if specified by do_step
  if (case == "impute" & !is.element(4, o$do_step)) return()
  if (case == "infer"  & !is.element(7, o$do_step)) return()
  
  message("* Running regression process: ", case)
  
  # TODO: Repeat process for DALYs
  
  # Create set of models to evaluate
  list[models, covars] = define_models(case)
  
  # Load response variable (impact per FVP)
  target = get_regression_data(case)
  
  # Return out if no training data identified
  if (nrow(target) == 0)
    return()
  
  # ---- Perform regression ----
  
  # TEMP: Use basic or full imputation method
  #
  # OPTIONS:
  #  basic_regression - IA2030 method using GBD covariates
  #  perform_regression1 - Helen's time series regression method
  #  perform_regression2 - Same as Helen's method, refactored code
  method = "perform_regression2"
  
  # Which sources of public health impact are to be modelled
  if (case == "impute") use_sources = qc(vimc)
  if (case == "infer")  use_sources = qc(vimc, static, extern)
  
  # TEMP: Ignoring problematic cases for now
  ignore = c(7, 11, 14,   # Non-routine, small numbers
             16, 18, 22,  # DTP boosters causing problems
             25)          # Polio results not ready yet
             
  # Call country imputation function
  predict_dt = table("d_v_a") %>%
    filter(source %in% use_sources) %>%
    # TEMP: Ignoring problematic cases for now...
    filter(!d_v_a_id %in% ignore) %>%  
    # Apply geographical imputation model...
    pull(d_v_a_id) %>%
    lapply(FUN = get(method), 
           models = models, 
           covars = covars, 
           target = target, 
           case   = case) %>%
    rbindlist() %>%
    select(d_v_a_id, country, year, impact_impute)
  
  # Apply imputations where needed
  # 
  # NOTE: A trivial process in the infer case
  impute_dt = target %>%
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
    mutate(impact = pmin(impact, 1))  # HCJ: This can be removed when fitting looks good
  
  # Save imputed results to file
  save_rds(impute_dt, case, case, "result")
  
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
# Define set of regression models to evaluate
# ---------------------------------------------------------
define_models = function(case) {
  
  # List of available covariates
  covars = list(
    cov0  = "log(coverage)",
    cov1  = "log(coverage_minus_1)", 
    cov2  = "log(coverage_minus_2)", 
    cov3  = "log(coverage_minus_3)", 
    cov4  = "log(coverage_minus_4)", 
    pop14 = "pop_0to14", 
    gini  = "gini", 
    hdi   = "HDI", 
    ab    = "attended_births", 
    hs0   = "health_spending", 
    hs1   = "health_spending_minus_1", 
    hs2   = "health_spending_minus_2")
  
  # Define models (using shorthand covariate references)
  models = list(
    
    # Models for imputing missing countries
    impute = list(
      m1  = "cov0", 
      m2  = "cov0 + cov1", 
      m3  = "cov0 + cov1 + cov2", 
      m4  = "cov0 + cov1 + cov2 + cov3", 
      m5  = "cov0 + cov1 + cov2 + cov3 + cov4", 
      # m6  = "cov0 + cov1 + cov2",
      m7  = "cov0 + cov1 + cov2 + cov3 + cov4 + pop14 + gini", 
      m8  = "cov0 + cov1 + cov2 + cov3 + cov4 + pop14 + gini + ab", 
      # m9  = "cov0 + cov1 + cov2 + cov3 + cov4 + pop14 + gini + ab",
      m10 = "cov0 + cov1 + cov2 + cov3 + hdi + pop14 + gini", 
      m11 = "cov0 + cov1 + cov2 + hdi + pop14 + gini", 
      m12 = "cov0 + cov1 + hdi + pop14 + gini", 
      m13 = "cov0 + hdi + pop14 + gini"), 
    
    # Models for inferring key drivers of impact 
    infer = list(
      x1  = "cov0", 
      x8  = "cov0 + cov1 + cov2 + cov3 + cov4 + pop14 + gini + ab", 
      x13 = "cov0 + hdi + pop14 + gini"))
  
  return(list(models[[case]], covars))
}

# ---------------------------------------------------------
# Load/calculate target variable (impact per FVP)
# ---------------------------------------------------------
get_regression_data = function(case) {
  
  # Population size of each country over time
  pop_dt = table("wpp_pop") %>%
    lazy_dt() %>%
    group_by(country, year) %>%
    summarise(pop = sum(pop)) %>%
    ungroup() %>%
    as.data.table()
  
  # Impact estimates in imputation case: VIMC pathogens, VIMC countries
  if (case == "impute") {
    outcomes_dt = table("vimc_estimates") %>%
      rename(impact = deaths_averted)
  }
  
  # Impact estimates in imputation case: all modelled results
  if (case == "infer")
    outcomes_dt = read_rds("history", "deaths_averted")
  
  # TODO: Probably best to remove previously imputed estimates in the infer case
  
  # Convert estimates to cumulative form
  impact_dt = outcomes_dt %>%
    lazy_dt() %>%
    # Sum impact over age...
    group_by(d_v_a_id, country, year) %>%
    summarise(impact_abs = sum(impact)) %>%
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
    # Subset pathogens and years...
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
  save_rds(target_dt, case, case, "target")
  
  # Throw a warning if no target data identified
  if (nrow(target_dt) == 0)
    warning("No training data identified")
  
  return(target_dt)
}

# ---------------------------------------------------------
# Perform regression
# ---------------------------------------------------------
perform_regression2 = function(d_v_a_id, models, covars, target, case) {
  
  # Stepwise regression
  # TODO: Update to lasso regularisation for optimal predictor selection
  
  # Extract name of this d-v-a
  d_v_a_name = table("d_v_a") %>%
    filter(d_v_a_id == !!d_v_a_id) %>%
    pull(d_v_a_name)
  
  # Display progress message to user
  message(" - ", d_v_a_name)
  
  # Append all required covariates - see separate function
  target_ts = append_covariates(d_v_a_id, models, covars, target)
  
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
  
  # ---- Predictions ----
  
  message("  > Model predictions")
  
  # Impute case: predict missing data using regional best models 
  if (case == "impute") {
    
    # TODO: Generalise: allow model selection for imputed country by e.g. region. For now, model 13 works well in general
    # TODO: Choose most appropriate method for selecting coefficients e.g. nearest neighbours
    
    # TEMP: Use model 13 (a commonly selected model) for predictions
    use_model = "m13"
    
    # Evaluate this model - see seperate function
    predict_dt = evaluate_predictions(
      id         = use_model, 
      model_list = model_list, 
      target     = target_ts, 
      case       = case)
  }
  
  # Infer case: just a case of evaluating on the training data
  if (case == "infer") {
    
    # Evaluate models on the training data
    predict_dt = augment(model_choice) %>%
      select(country, year, prediction = .fitted) %>%
      as.data.table()
  }
  
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
  save_rds(fit, case, case, d_v_a_id)
  
  return(result_dt)
}

# ---------------------------------------------------------
# Append all required covariates
# ---------------------------------------------------------
append_covariates = function(d_v_a_id, models, covars, target) {
  
  # TODO: Are GBD covariates still used/needed?
  
  # ---- Identify covariates from model specification ----
  
  # Shorthand covariates used in specified models
  covars_used = unlist(models) %>%
    paste(collapse = " + ") %>%
    str_split_1(pattern = " \\+ ") %>%
    unique()
  
  # Associated names of covariate columns 
  covars_retain = covars[covars_used] %>%
    unlist(covars) %>%
    str_remove("^.*\\(+") %>%
    str_remove("\\)+$")
  
  # ---- Define covariates to be lagged ----
  
  # Details of covariates we wish to lag
  lag_dt = covars_retain %>%
    str_split("_minus_", simplify = TRUE) %>%
    as_named_dt(c("covar", "idx")) %>%
    filter(idx > 0)
  
  # Extract all to-be-lagged covariates
  covars_lag  = unique(lag_dt$covar)
  n_lag_years = max(as.numeric(lag_dt$idx))
  
  # Small function to apply lag to given covariate
  covar_lag_fn = function(dt) {
    
    # Iterate through covariates and years to lag
    for (i in covars_lag) {
      for (j in seq_len(n_lag_years)) {
        
        # Incrementally offset by one year
        dt[[paste1(i, "minus", j)]] = lag(dt[[i]], j)
      }
    }
    
    return(dt)
  }
  
  # ---- Format coverage (a key predictor) ----
  
  # Summarise vaccination coverage by country and year
  coverage_dt = table("coverage") %>%
    lazy_dt() %>%
    filter(d_v_a_id == !!d_v_a_id) %>%
    # Summarise over age groups...
    group_by(country, year) %>%
    summarise(fvps   = sum(fvps),
              cohort = sum(cohort)) %>%
    ungroup() %>%
    # Recalculate for whole population...
    mutate(coverage = fvps / cohort) %>%
    select(country, year, coverage) %>%
    as.data.table()
  
  # ---- Append all other covariates ----
  
  # Create time-series tibble with all covariates
  target_ts = target %>%
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
    # TODO: This shouldn't be necessary
    arrange(country, year) %>%
    group_by(country, year) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    # Lag any necessary covariates...
    group_by(country) %>%
    covar_lag_fn() %>%
    ungroup() %>%
    # Only retain covariates defined in models...
    select(country, year, target, 
           all_of(covars_retain)) %>%
    # Convert to time-series tibble...
    as_tsibble(index = year, 
               key   = country) 
  
  return(target_ts)
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
# Evalulate chosen model for all settings
# ---------------------------------------------------------
evaluate_predictions = function(id, model_list, target, case) {
  
  # TODO: Set up predictions also by period
  
  # Full set of models available - we'll subset for this modal ID
  list[models, covars] = define_models(case)
  
  # Index of this model in model list
  idx = which(names(models) == id)
  
  # ---- Determine predictors by region and period ----
  
  # Prediction are performed at regional level
  region_dt = table("country") %>%
    select(country, region)
  
  # Take the median coefficient across each region
  coefficient_dt = model_list %>%
    pluck(idx) %>%
    tidy() %>%
    lazy_dt() %>%
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
  
  # ---- Construct predictor function call ----
  
  # Small function to wrap a string in quotes
  quote = function(x, q = '"') 
    paste0(q, x, q)
  
  # Column names of predictors
  predict_covars = models[[id]] %>%
    interpret_covars(covars) %>%
    str_remove_all(" ") %>%
    str_split("\\+") %>%
    pluck(1)
  
  # Construct linear product of predictors and coefficients
  predict_str = predict_covars %>%
    paste1("coefficient") %>%
    quote("`") %>%
    paste(predict_covars, sep = " * ") %>%
    paste(collapse = " + ")
  
  # Construct complete function call to be evaluated
  predict_fn   = paste0("prediction = exp(", predict_str, ")")
  predict_eval = paste0("predictors %>% mutate(", predict_fn, ")")
  
  # Append coefficients to predictors
  predictors = target %>% 
    inner_join(y  = region_dt, 
               by = "country") %>% 
    left_join(y  = coefficient_dt, 
              by = "region") %>% 
    as.data.table()
  
  # Evaluate function call to predict all target values
  predict_dt = eval_str(predict_eval) %>%
    select(country, year, prediction)
  
  return(predict_dt)
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

