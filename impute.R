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
  if (!is.element(3, o$do_step)) return()
  
  message("* Running country imputation")
  
  # TODO: Repeat process for DALYs
  
  # Load target to fit to (impact per FVP)
  target_dt = get_impute_data()
  
  # ---- Perform imputation ----
  
  browser() # Check d_v_a table
  
  # Impute missing countries for all d-v-a combinations
  impute_dt = table("d_v_a") %>%
    # Filter for VIMC pathogens only...
    left_join(y  = table("disease"),
              by = "disease") %>%
    filter(source == "vimc") %>%
    # Apply geographical imputation model...
    pull(d_v_a_id) %>%
    lapply(perform_impute, target = target_dt) %>%
    rbindlist() %>%
    # Merge VIMC estimates with those just imputed...
    mutate(impact = ifelse(
      test = is.na(target),
      yes  = impact_impute,
      no   = impact_cum)) %>%
    select(country, d_v_a_id, year,
           fvps = fvps_cum, impact)
  
  # Save imputed results to file
  save_rds(impute_dt, "impute", "impute_result")
  
  # ---- Plot results ----
  
  # NOTE: All plotting functionality lives in plotting.R
  
  # Plot predictor-response relationships
  plot_covariates()
  
  # Plot imputation quality of fit
  plot_impute_quality()
  
  # Plot train-predict countries
  plot_impute_countries()
}

# ---------------------------------------------------------
# Perform imputation
# ---------------------------------------------------------
perform_impute = function(d_v_a_id, target) {
  
  # Details of this d_v_a
  d_v_a_name = data.table(d_v_a_id = d_v_a_id) %>%
    append_d_v_a_name() %>%
    pull(d_v_a_name)
  
  # Display progress message to user
  message(" - ", d_v_a_name)
  
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
    # Append GBD indices and infant mortality...
    left_join(y  = table("gbd_covariates"),
              by = c("country", "year")) %>%
    left_join(y  = infant_mortality_dt,
              by = c("country", "year")) %>%
    # Calculate n years of estimates...
    mutate(n_years = 1) %>%
    group_by(country) %>%
    mutate(n_years = cumsum(n_years)) %>%
    ungroup() %>%
    as.data.table()
  
  # TODO: We can either include or exclude zero here - arguably with zero is better...
  
  # Data used to fit statistical model
  data_dt = target_dt %>%
    filter(!is.na(target)) %>%
    # filter(target > 0) %>%
    select(target, n_years, sdi, haqi, imr) %>%
    # Remove target outliers for better normalisation...
    mutate(lower = mean(target) - 3 * sd(target), 
           upper = mean(target) + 3 * sd(target), 
           outlier = target < lower | target > upper) %>%
    filter(outlier == FALSE) %>%
    select(-outlier, -lower, -upper)
  
  # Sanity check that we have no NAs here
  if (any(is.na(data_dt)))
    stop("NA values identified in predictors")
  
  # Values to predict for (including data used for fitting)
  pred_dt = target_dt %>%
    select(all_of(names(data_dt)))
  
  # ---- Check for trivial case ----
  
  # Return out if no data available
  if (nrow(data_dt[target > 0]) < 10) {
    
    message(" !! Insufficient data for imputation !!")
    
    # Store trivial outcomes
    fit = list(data = data_dt, result = NULL)
    
    # Save to file
    save_rds(fit, "impute", "impute", d_v_a_id)
    
    return()
  }
  
  # ---- Normalise predictors and response ----
  
  # Function to normalise ready for fitting
  transform_fn = function(x, a, b)
    y = t((x - a) / (b - a)) %>% as.data.table()
  
  # Function to back transform to original scale
  retransform_fn = function(y, a, b)
    x = y * (b["target"] - a["target"]) + a["target"]
  
  # Matrices of points to fit with and points to predict for
  data_mat = t(as.matrix(data_dt))
  pred_mat = t(as.matrix(pred_dt))
  
  # Min and max in data used for fitting
  a = rowMins(data_mat)
  b = rowMaxs(data_mat)
  
  # Use these min ana max values to normalise
  norm_data_dt = transform_fn(data_mat, a, b)
  norm_pred_dt = transform_fn(pred_mat, a, b)
  
  # ---- Fit a model to predict impact per FVP ----
  
  # Fit a GLM for impact per FVP using all covariates
  fit_model = glm(
    formula = target ~ n_years + sdi + haqi + imr, 
    data    = norm_data_dt)
  
  # Use fitted model to predict 
  result_dt = target_dt %>%
    select(country, d_v_a_id, year, fvps_cum, impact_cum) %>%
    # Predict impact per FVP...
    cbind(norm_pred_dt) %>%
    mutate(predict = predict(fit_model, .), 
           predict = pmax(predict, 0)) %>%
    # Remove predictors...
    select(country, d_v_a_id, year, fvps_cum, impact_cum, 
           target, predict) %>%
    # Back-transform target and prediction...
    mutate(target  = retransform_fn(target,  a, b), 
           predict = retransform_fn(predict, a, b)) %>%
    # Multiply through to obtain cumulative impact over time...
    mutate(impact_impute = fvps_cum * predict, 
           .after = impact_cum)
  
  # Sanity check that all predicted values are legitimate
  if (any(is.na(result_dt$predict)))
    stop("NA values identified in predicted impact")
  
  # Store the fitted model, the data used, and the result
  fit = list(
    model   = fit_model, 
    data    = norm_data_dt, 
    result  = result_dt)
  
  # Save to file
  save_rds(fit, "impute", "impute", d_v_a_id)
  
  return(result_dt)
}

# ---------------------------------------------------------
# Load/calculate target (impact per FVP) for modelled pathogens
# ---------------------------------------------------------
get_impute_data = function() {
  
  # Population size of each country over time
  pop_dt = table("wpp_pop") %>%
    group_by(country, year) %>%
    summarise(pop = sum(pop)) %>%
    ungroup() %>%
    as.data.table()
  
  # Wrangle VIMC impact estimates
  impact_dt = table("vimc_estimates") %>%
    # Sum impact over age...
    group_by(country, d_v_a_id, year) %>%
    summarise(impact_abs = sum(deaths_averted)) %>%
    ungroup() %>%
    mutate(impact_abs = pmax(impact_abs, 0)) %>%
    # Scale results to per capita...
    left_join(y  = pop_dt, 
              by = c("country", "year")) %>%
    mutate(impact_rel = impact_abs / pop) %>%
    select(-pop) %>%
    # Cumulative sum impact...
    arrange(country, d_v_a_id, year) %>%
    group_by(country, d_v_a_id) %>%
    mutate(impact_cum = cumsum(impact_rel)) %>%
    ungroup() %>%
    as.data.table()
  
  browser() # Check d_v_a table
  
  # Extract FVPs
  fvps_dt = table("coverage") %>%
    # Append d_v_a details...
    left_join(y  = table("v_a"),
              by = "v_a_id") %>%
    left_join(y  = table("d_v_a"),
              by = c("vaccine", "activity")) %>%
    # Only impute pathogens and years for which we've VIMC estimates...
    filter(d_v_a_id %in% unique(impact_dt$d_v_a_id), 
           year     %in% unique(impact_dt$year)) %>%
    # Summarise over age...
    group_by(country, d_v_a_id, year) %>%
    summarise(fvps_abs = sum(fvps)) %>%
    ungroup() %>%
    # Scale results to per capita...
    left_join(y  = pop_dt, 
              by = c("country", "year")) %>%
    mutate(fvps_rel = fvps_abs / pop) %>%
    select(-pop) %>%
    # Cumulative sum FVPs...
    arrange(country, d_v_a_id, year) %>%
    group_by(country, d_v_a_id) %>%
    mutate(fvps_cum = cumsum(fvps_rel)) %>%
    ungroup() %>%
    as.data.table()
  
  # Combine into single datatable
  target_dt = fvps_dt %>%
    left_join(y  = impact_dt, 
              by = c("country", "d_v_a_id", "year")) %>%
    # Impact per FVP...
    mutate(target = impact_cum / fvps_cum)
  
  # Save this datatable to file for plotting purposes
  save_rds(target_dt, "impute", "target")
  
  return(target_dt)
}

