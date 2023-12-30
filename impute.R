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
  # Summarise vaccination coverage by country, by year
  coverage_dt = table("coverage") %>%
                as.data.frame() %>%
                  select(-coverage) %>% # This is coverage by age group. We recalculate for whole population
                  group_by(country, v_a_id, year) %>%
                  mutate(total_fvps = sum(fvps),
                         total_pop = sum(cohort)) %>%
                  summarise(total_fvps = first(total_fvps),
                            total_pop = first(total_pop)) %>%
                  mutate(coverage = total_fvps / total_pop) %>%
                  ungroup() %>%
                  select(-c(total_fvps, total_pop))
  
 # browser()
    # Append covariates to target
  target_dt = target %>%
    filter(d_v_a_id == !!d_v_a_id) %>%
    # Append vaccination coverage, GBD covariates and Gapminder covariates
        full_join(y = coverage_dt,  
              by = c("country", "year")) %>%
        full_join(y  = table("gbd_covariates"),
              by = c("country", "year")) %>%
        inner_join(y  = table("gapminder_covariates"),    #TODO: multiple entries for COD(Congo, Kinshasa)
              by = c("country", "year"), relationship = "many-to-many") %>%
    arrange(country, year) %>%
    
    # Summarise to single row for each d_v_a_id per country per year (accounting for new v_a_id functionality)
    group_by(country, d_v_a_id, year) %>%
    filter(row_number() == 1) %>%
    select(-v_a_id) %>% 
    ungroup () %>%
    
    # Create dummy variables for historic coverage (NA prior to our data)
    mutate(coverage_minus_1 = lag(coverage, 1),
           coverage_minus_2 = lag(coverage, 2),
           coverage_minus_3 = lag(coverage, 3),
           coverage_minus_4 = lag(coverage, 4),
           coverage_minus_5 = lag(coverage, 5),
           coverage_minus_6 = lag(coverage, 6),
           coverage_minus_7 = lag(coverage, 7),
           coverage_minus_8 = lag(coverage, 8),
           coverage_minus_9 = lag(coverage, 9)) %>%
  
    # Create dummy variables for historic health_spending (NA prior to our data)
    mutate(health_spending_minus_1 = lag(health_spending, 1),
           health_spending_minus_2 = lag(health_spending, 2),
           health_spending_minus_3 = lag(health_spending, 3),
           health_spending_minus_4 = lag(health_spending, 4),
           health_spending_minus_5 = lag(health_spending, 5),
           health_spending_minus_6 = lag(health_spending, 6),
           health_spending_minus_7 = lag(health_spending, 7),
           health_spending_minus_8 = lag(health_spending, 8),
           health_spending_minus_9 = lag(health_spending, 9)) %>%
    
      as.data.table()
  
#browser()  
   # Convert to tsibble format for time series regression by country
   data_dt = target_dt %>%
    filter(!is.na(target)) %>% 
     # Replace zero impact with minimal impact for log transformation
     mutate(target = ifelse(target==0, 1e-20, target)) %>%
     as_tsibble(index = year, key = c(country,d_v_a_id))
   
  # TODO Automate model comparison for each d_v_a
  # Ad-hoc model comparison (AICc) was conducted here....
   impact_model = data_dt %>%
     filter(country %in% c("ECU", "FJI", "GEO")) %>%   # TODO: generalise to avoid failure to fit models with paucity of data
     model(tslm = TSLM(log(target) ~ log(coverage) +
                                     log(coverage_minus_1) +
                                     log(coverage_minus_2) +
                                     HDI +
                                     pop_0to14 +
                                     attended_births +
                                     gini)) 
   
   #report(impact_model) # Quick check summary
  #browser() 
   # Arrange in tidy format for easy access of estimates, p-values etc.
   model_fit = impact_model %>% tidy() 
   
   # Link model output back to WHO region
   regions = as.data.frame(data_dt) %>% select(country, region_short) %>% unique()
   model_fit = inner_join(x=regions, y=model_fit, by="country")
   
   # Plot estimates of regression predictors by region
   plot_model_fit_df = model_fit %>%
                        #filter(p.value <= 0.05) %>%
                        select(country, region_short, term, estimate) %>%
                        spread(term, estimate)# %>%
                       # rename(check = 'log(coverage)')
#browser()   
   
  # Plot correlation between two predictor variables
  ggplot(plot_model_fit_df, aes(gini, HDI, colour = region_short)) +
                     geom_point()
   
  # Plot correlation between multiple predictor variables   
  plot_model_fit_df %>% select(gini, HDI, `log(coverage)`, pop_0to14) %>% ggpairs(upper = list(continuous = wrap("cor", method = "spearman")))
  
  
  # Plot model fit for a single country 
  plot_df = augment(impact_model) %>%
     filter(country == "ECU")
  
  ggplot(data = plot_df, aes(x = target, y = .fitted)) +
     geom_point() +
     labs(
       y = "Fitted (predicted values)",
       x = "Data (actual values)",
       title = paste("Vaccine impact of", d_v_a_name, "in", plot_df$country)
     ) +
     geom_abline(intercept = 0, slope = 1)

  
  # Plot model fit for a single country 
  plot_df = augment(impact_model) %>%
    filter(country == "GEO")
  
  ggplot(data = plot_df, aes(x = year)) +
     geom_line(aes(y = target, colour = "Data")) +
     geom_line(aes(y = .fitted, colour = "Fitted")) +
     labs(y = NULL,
          title = paste("Vaccine impact of", d_v_a_name, "in", plot_df$country)
     ) +
     scale_colour_manual(values=c(Data="black",Fitted="#D55E00")) +
     guides(colour = guide_legend(title = NULL))
   
   browser()
  # Manually explore associations between predictor variables for different geographical regions and time points
  #  explore_dt =  data_dt %>% as.data.table() %>% # Transform to data table to remove country as categorical variable
  #                  filter(#year > 2000 & year <= 2020 &
                    #region_short == "AFR" &
  #                  target > 2e-20
  #                  ) %>%
  #                select(-country) %>%
  #                select(target, gini, health_spending, coverage, coverage_minus_1,coverage_minus_2,coverage_minus_3,sdi) 

   # explore_dt %>%   ggpairs(upper = list(continuous = wrap("cor", method = "spearman"))) # Use Spearman rank correlation to account for outliers
  
   # Reset data used to fit statistical model
   data_dt = target_dt %>%
     filter(!is.na(target)) %>%
     select(target, coverage, sdi)
            
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
    formula = target ~ coverage + sdi, 
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
 # if (any(is.na(result_dt$predict)))
#    stop("NA values identified in predicted impact")
  
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

