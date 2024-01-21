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
  if (!is.element(4, o$do_step)) return()
  
  message("* Running country imputation")
  
  # TODO: Repeat process for DALYs
  
  # Load target to fit to (impact per FVP)
  target_dt = get_impute_data()
  
  # ---- Perform imputation ----
  
  # Impute missing countries for all d-v-a combinations
  impute_dt = table("d_v_a") %>%
    
  # TODO: For debugging only!!  
    filter(d_v_a_id %in% c(1)) %>%
    
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
  
  message(" - ", table("d_v_a")[d_v_a_id, d_v_a_name])
  
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

    # Append covariates to target
  target_dt = target %>%
    filter(d_v_a_id == !!d_v_a_id) %>%
    # Append vaccination coverage, GBD covariates and Gapminder covariates
        full_join(y = coverage_dt,  
              by = c("country", "year")) %>%
        full_join(y  = table("gbd_covariates"),
              by = c("country", "year")) %>%
        full_join(y  = table("gapminder_covariates"),    #TODO: multiple entries for COD(Congo, Kinshasa)
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
           coverage_minus_3 = lag(coverage, 3)) %>%
  
    # Create dummy variables for historic health_spending (NA prior to our data)
    mutate(health_spending_minus_1 = lag(health_spending, 1),
           health_spending_minus_2 = lag(health_spending, 2),
           health_spending_minus_3 = lag(health_spending, 3),
           health_spending_minus_4 = lag(health_spending, 4),
           health_spending_minus_5 = lag(health_spending, 5),
           health_spending_minus_6 = lag(health_spending, 6),
           health_spending_minus_7 = lag(health_spending, 7),
           health_spending_minus_8 = lag(health_spending, 8)) %>%
    
    as_tsibble(index = year, key = c(country,d_v_a_id)) 
  
   # Convert to tsibble format for time series regression by country
   data_dt = target_dt %>%
    filter(!is.na(target)) %>% 
     
     # Remove zeros to allow for log transformation
     filter(target != 0) %>%
     
     # TODO: Activate this functionality (will need to select appropriate model downstream for comparing time periods)
     # Please extend to allow fitting to the temporally extrapolated data -> useful messages on determinants of vaccine impact over time 
     
     # Set time periods
     mutate(period = case_when(year < 1984 ~ 1,
                               year >= 1984 & year < 1994 ~ 2,
                               year >= 1994 & year < 2004 ~ 3,
                               year >= 2004 & year < 2014 ~ 4,
                               year >= 2014 & year < 2024 ~ 5)) %>%
     
     # Prepare to fit a model to each country / d_v_a_id combination
     as_tsibble(index = year, key = c(country,d_v_a_id)) %>%
     
     # TODO: Here is the option to fit models for each decade, for each country, for each d_v_a. There may be issues with data completeness for some combinations
     #as_tsibble(index = year, key = c(country,d_v_a_id, period)) %>%
     
     group_by(country) %>%
     filter(n() > 4) %>% # remove if fewer than 4 non-zero values for a given country (insufficient for fitting)
     ungroup()
   

   
  # Stepwise regression
  # TODO LATER: Update to lasso regularisation for optimal predictor selection
  # TODO LATER: Call sets of predictors from table and not write out this long way (see broom)

  # Model 1: log(coverage)
  model_1 =  data_dt %>%
              model(tslm = TSLM(log(target) ~ log(coverage))) %>%
              mutate(model_number = 1) 
   
  model_report_1 = report(model_1) %>% # Glance of fit stats inc R-squared, p-value, AIC, AICc etc.
                  mutate(model_number = 1) %>%
    select(-c(r_squared, adj_r_squared, sigma2, statistic, AIC, BIC, CV, deviance, df.residual, rank))
  
  model_fit_1 = model_1 %>% 
                  tidy() %>% # Arrange in tidy format for easy access of estimates, p-values etc.
                  mutate(model_number = 1) 
    
  # Model 2: log(coverage), log(coverage_minus_1)
  model_2 = data_dt %>%
    model(tslm = TSLM(log(target) ~ log(coverage) +
                        log(coverage_minus_1) 
                         )) %>%
    mutate(model_number = 2)
  
  model_report_2 = report(model_2) %>% # Glance of fit stats inc R-squared, p-value, AIC, AICc etc.
    mutate(model_number = 2) %>%
    select(-c(r_squared, adj_r_squared, sigma2, statistic, AIC, BIC, CV, deviance, df.residual, rank))
  
  model_fit_2 = model_2 %>% 
                 tidy() %>% # Arrange in tidy format for easy access of estimates, p-values etc.
                 mutate(model_number = 2) 
  
  # Model 3: log(coverage), log(coverage_minus_1), log(coverage_minus_2)
  model_3 = data_dt %>%
    model(tslm = TSLM(log(target) ~ log(coverage) +
                        log(coverage_minus_1) +
                        log(coverage_minus_2) 
    )) %>%
    mutate(model_number = 3)
  
  model_report_3 = report(model_3) %>% # Glance of fit stats inc R-squared, p-value, AIC, AICc etc.
    mutate(model_number = 3) %>%
    select(-c(r_squared, adj_r_squared, sigma2, statistic, AIC, BIC, CV, deviance, df.residual, rank))
  
  model_fit_3 = model_3 %>% 
                tidy() %>% # Arrange in tidy format for easy access of estimates, p-values etc.
                mutate(model_number = 3) 
  
  # Model 4: log(coverage), log(coverage_minus_1), log(coverage_minus_2), log(coverage_minus_3)
  model_4 = data_dt %>%
    model(tslm = TSLM(log(target) ~ log(coverage) +
                        log(coverage_minus_1) +
                        log(coverage_minus_2) +
                        log(coverage_minus_3) 
                        )) %>%
    mutate(model_number = 4)
  
  model_report_4 = report(model_4) %>% # Glance of fit stats inc R-squared, p-value, AIC, AICc etc.
    mutate(model_number = 4) %>%
    select(-c(r_squared, adj_r_squared, sigma2, statistic, AIC, BIC, CV, deviance, df.residual, rank))
  
  
  model_fit_4 = model_4 %>% 
                 tidy() %>% # Arrange in tidy format for easy access of estimates, p-values etc.
                 mutate(model_number = 4) 
  
  # Model 5: log(coverage), log(coverage_minus_1), log(coverage_minus_2), log(coverage_minus_3), log(coverage_minus_4)
  model_5 = data_dt %>%
    model(tslm = TSLM(log(target) ~ log(coverage) +
                        log(coverage_minus_1) +
                        log(coverage_minus_2) +
                        log(coverage_minus_3) +
                        log(coverage_minus_4) 
    )) %>%
    mutate(model_number = 5)
  
  model_report_5 = report(model_5) %>% # Glance of fit stats inc R-squared, p-value, AIC, AICc etc.
    mutate(model_number = 5) %>%
    select(-c(r_squared, adj_r_squared, sigma2, statistic, AIC, BIC, CV, deviance, df.residual, rank))
  
  model_fit_5 = model_5 %>% 
                tidy() %>% # Arrange in tidy format for easy access of estimates, p-values etc.
                mutate(model_number = 5) 
  
  # Model 6: log(coverage), log(coverage_minus_1), log(coverage_minus_2), log(coverage_minus_3), log(coverage_minus_4), HDI
  model_6 = data_dt %>%
    model(tslm = TSLM(log(target) ~ log(coverage) +
                        log(coverage_minus_1) +
                        log(coverage_minus_2) +
                        log(coverage_minus_3) +
                        log(coverage_minus_4) +
                        HDI 
    )) %>%
    mutate(model_number = 6)
  
  model_report_6 = report(model_6) %>% # Glance of fit stats inc R-squared, p-value, AIC, AICc etc.
    mutate(model_number = 6) %>%
    select(-c(r_squared, adj_r_squared, sigma2, statistic, AIC, BIC, CV, deviance, df.residual, rank))
  
  model_fit_6 = model_6 %>% 
                tidy() %>% # Arrange in tidy format for easy access of estimates, p-values etc.
                mutate(model_number = 6) 
  
  # Model 7: log(coverage), log(coverage_minus_1), log(coverage_minus_2), log(coverage_minus_3), log(coverage_minus_4), HDI, pop_0to14
  model_7 = data_dt %>%
    model(tslm = TSLM(log(target) ~ log(coverage) +
                        log(coverage_minus_1) +
                        log(coverage_minus_2) +
                        log(coverage_minus_3) +
                        log(coverage_minus_4) +
                        HDI +
                        pop_0to14 
    )) %>%
    mutate(model_number = 7)
  
  model_report_7 = report(model_7) %>% # Glance of fit stats inc R-squared, p-value, AIC, AICc etc.
    mutate(model_number = 7) %>%
    select(-c(r_squared, adj_r_squared, sigma2, statistic, AIC, BIC, CV, deviance, df.residual, rank))
  
  model_fit_7 = model_7 %>% 
                  tidy() %>% # Arrange in tidy format for easy access of estimates, p-values etc.
                  mutate(model_number = 7) 
  
  # Model 8: log(coverage), log(coverage_minus_1), log(coverage_minus_2), log(coverage_minus_3), log(coverage_minus_4), HDI, pop_0to14, gini
  model_8 = data_dt %>%
    model(tslm = TSLM(log(target) ~ log(coverage) +
                        log(coverage_minus_1) +
                        log(coverage_minus_2) +
                        log(coverage_minus_3) +
                        log(coverage_minus_4) +
                        HDI +
                        pop_0to14 +
                        gini
    )) %>%
    mutate(model_number = 8)
  
  model_report_8 = report(model_8) %>% # Glance of fit stats inc R-squared, p-value, AIC, AICc etc.
    mutate(model_number = 8) %>%
    select(-c(r_squared, adj_r_squared, sigma2, statistic, AIC, BIC, CV, deviance, df.residual, rank))
  
  model_fit_8 = model_8 %>% 
                 tidy() %>% # Arrange in tidy format for easy access of estimates, p-values etc.
                 mutate(model_number = 8) 
  
  # Model 9: log(coverage), log(coverage_minus_1), log(coverage_minus_2), log(coverage_minus_3), log(coverage_minus_4), HDI, pop_0to14, attended_births, gini
  model_9 = data_dt %>%
    model(tslm = TSLM(log(target) ~ log(coverage) +
                        log(coverage_minus_1) +
                        log(coverage_minus_2) +
                        log(coverage_minus_3) +
                        log(coverage_minus_4) +
                        HDI +
                        pop_0to14 +
                        attended_births +
                        gini
    )) %>%
    mutate(model_number = 9)
  
  model_report_9 = report(model_9) %>% # Glance of fit stats inc R-squared, p-value, AIC, AICc etc.
    mutate(model_number = 9) %>%
    select(-c(r_squared, adj_r_squared, sigma2, statistic, AIC, BIC, CV, deviance, df.residual, rank))
  
  model_fit_9 = model_9 %>% 
                tidy() %>% # Arrange in tidy format for easy access of estimates, p-values etc.
                mutate(model_number = 9) 
  
  # Model 10: log(coverage), log(coverage_minus_1), log(coverage_minus_2), log(coverage_minus_3), HDI, pop_0to14,  gini
  model_10 = data_dt %>%
    model(tslm = TSLM(log(target) ~ log(coverage) +
                        log(coverage_minus_1) +
                        log(coverage_minus_2) +
                        log(coverage_minus_3) +
                        HDI +
                        pop_0to14 +
                      #  attended_births +
                        gini
    )) %>%
    mutate(model_number = 10)
  
  model_report_10 = report(model_10) %>% # Glance of fit stats inc R-squared, p-value, AIC, AICc etc.
    mutate(model_number = 10) %>%
    select(-c(r_squared, adj_r_squared, sigma2, statistic, AIC, BIC, CV, deviance, df.residual, rank))
  
  model_fit_10 = model_10 %>%
                    tidy() %>%# Arrange in tidy format for easy access of estimates, p-values etc.
                    mutate(model_number =10)
  
  # Model 11: log(coverage), log(coverage_minus_1), log(coverage_minus_2), HDI, pop_0to14,  gini
  model_11 = data_dt %>%
    model(tslm = TSLM(log(target) ~ log(coverage) +
                        log(coverage_minus_1) +
                        log(coverage_minus_2) +
                        HDI +
                        pop_0to14 +
                        #attended_births +
                        gini
    )) %>%
    mutate(model_number = 11)
  
  model_report_11 = report(model_11) %>% # Glance of fit stats inc R-squared, p-value, AIC, AICc etc.
    mutate(model_number = 11) %>%
    select(-c(r_squared, adj_r_squared, sigma2, statistic, AIC, BIC, CV, deviance, df.residual, rank))
  
  model_fit_11 = model_11 %>% 
                   tidy() %>%# Arrange in tidy format for easy access of estimates, p-values etc.
                   mutate(model_number = 11)
  
  # Model 12: log(coverage), log(coverage_minus_1), HDI, pop_0to14, gini
  model_12 = data_dt %>%
    model(tslm = TSLM(log(target) ~ log(coverage) +
                        log(coverage_minus_1) +
                        HDI +
                        pop_0to14 +
                        gini
    )) %>%
    mutate(model_number = 12)
  
  model_report_12 = report(model_12) %>% # Glance of fit stats inc R-squared, p-value, AIC, AICc etc.
    mutate(model_number = 12) %>%
    select(-c(r_squared, adj_r_squared, sigma2, statistic, AIC, BIC, CV, deviance, df.residual, rank))
  
  model_fit_12 = model_12 %>% 
    tidy() %>%# Arrange in tidy format for easy access of estimates, p-values etc.
    mutate(model_number = 12)
  
  # Model 13: log(coverage), HDI, pop_0to14, gini
  model_13 = data_dt %>%
    model(tslm = TSLM(log(target) ~ log(coverage) +
                        HDI +
                        pop_0to14 +
                        gini
    )) %>%
    mutate(model_number = 13)
  
  model_report_13 = report(model_13) %>% # Glance of fit stats inc R-squared, p-value, AIC, AICc etc.
    mutate(model_number = 13) %>%
    select(-c(r_squared, adj_r_squared, sigma2, statistic, AIC, BIC, CV, deviance, df.residual, rank))
  
  model_fit_13 = model_13 %>%
                  tidy() %>%# Arrange in tidy format for easy access of estimates, p-values etc.
                  mutate(model_number = 13)
 
  # For each country, select the model with the best AICc (lowest number)
  model_choice = rbind(model_report_1,
                       model_report_2,
                       model_report_3,
                       model_report_4,
                       model_report_5,
                       model_report_6,
                       model_report_7,
                       model_report_8,
                       model_report_9,
                       model_report_10,
                       model_report_11,
                       model_report_12,
                       model_report_13) %>%
                 arrange(country) %>%
                 group_by(country) %>%
                 filter(!is.infinite(AICc)) %>% # remove null models
                 slice_min(AICc, with_ties = FALSE) %>% # if two models are equally the best, keep the first
                 select(-c(.model, df, log_lik, p_value))
   
   # Extract parameters of best fitting model for each country (according to AICc)
   model_fit = rbind(model_fit_1,
                     model_fit_2,
                     model_fit_3,
                     model_fit_4,
                     model_fit_5,
                     model_fit_6,
                     model_fit_7,
                     model_fit_8,
                     model_fit_9,
                     model_fit_10,
                     model_fit_11,
                     model_fit_12,
                     model_fit_13)
   
   model_fit = left_join(model_choice, model_fit, by=c("country", "model_number"))
 

   # Extract best fitting model for each country (according to AICc)
   all_models = rbind(model_1,
                     model_2,
                     model_3,
                     model_4,
                     model_5,
                     model_6,
                     model_7,
                     model_8,
                     model_9,
                     model_10,
                     model_11,
                     model_12,
                     model_13)
   
   best_model = left_join(model_choice, all_models, by=c("country", "model_number", "d_v_a_id")) %>%
                 as_mable(key = c(country, d_v_a_id), model = tslm)
   
   
   # Link model output back to WHO region
   regions = as.data.frame(target_dt) %>% select(country, region_short) %>% unique()
   model_fit = inner_join(x=regions, y=model_fit, by="country")
   model_choice = inner_join(x=regions, y=model_choice, by='country')
   
 
   ## TODO: Plot to move to plotting.R
   # Explore model selection by region
   ggplot(data = model_choice, aes(x = model_number)) +
     geom_histogram(binwidth = 1) +
        facet_wrap(~region_short)
   
   #TODO: Plot to move to plotting.R, fix for d_v_a_name, prettify
   # Plot data vs. fitted for all countries (model fit)
   plot_df = augment(best_model) 
   
   ggplot(data = plot_df, aes(x = target, y = .fitted)) +
     geom_point() +
    labs(
       y = "Fitted (predicted values)",
      x = "Data (actual values)",
       title = paste("Vaccine impact of", d_v_a_name)
     ) +
    geom_abline(intercept = 0, slope = 1) +
     facet_wrap(~country, ncol = 21)
   
   #TODO: Plot to move to plotting.R, fix for d_v_a_name, prettify
     # Plot model fit for a single country 
     plot_df = augment(best_model) %>%
       filter(country == "HND")
     
     ggplot(data = plot_df, aes(x = year)) +
       geom_point(aes(y = target, colour = "Data")) +
       geom_line(aes(y = .fitted, colour = "Fitted")) +
       labs(y = NULL,
           title = paste("Vaccine impact of", d_v_a_name, "in", plot_df$country)
     ) +
    scale_colour_manual(values=c(Data="black",Fitted="#D55E00")) +
     guides(colour = guide_legend(title = NULL))
   
     #TODO: Plot to move to plotting.R, fix for d_v_a_name, prettify
     # Plot model fit for all countries
   plot_df = augment(best_model)
   
   ggplot(data = plot_df, aes(x = year)) +
     geom_point(aes(y = target, colour = "Data")) +
    geom_line(aes(y = .fitted, colour = "Fitted")) +
     labs(y = NULL,
         title = paste("Vaccine impact of", d_v_a_name)
     ) +
    scale_colour_manual(values=c(Data="black",Fitted="#D55E00")) +
     guides(colour = guide_legend(title = NULL))  +
    facet_wrap(~country, ncol = 21)
   
 
  # Select countries for imputation
  impute_dt = target_dt %>% filter(! country %in% data_dt$country)
  
  # TODO LATER: Generalise. Allow model selection for imputed country by e.g. region. For now, model 13 works well in general
  # Model 13 is the most commonly selected, take median coefficient by WHO region (helps to avoid outliers)
  model_fit_13 = inner_join(x=regions, y=model_fit_13, by="country")
  
  model_13_region = model_fit_13 %>%
                      group_by(region_short, term) %>%
                      summarise(estimate = median(estimate, na.rm = TRUE), n = n()) %>%
                      pivot_wider(names_from = term,
                                  names_glue = "{term}_coefficient",
                                  values_from = estimate)
  
  # TODO LATER: Choose most appropriate method for selecting coefficients e.g. nearest neighbours
  # Impute target values using coefficients from WHO region
  impute_dt = left_join(impute_dt, model_13_region, by = c("region_short")) %>%
                mutate(estimate = exp((`log(coverage)_coefficient` * log(coverage)) +
                                       (HDI_coefficient * HDI) +
                                       (pop_0to14_coefficient * pop_0to14)+
                                       (gini_coefficient * gini))) %>% 
                select(-c(names(model_13_region))) %>%
                inner_join(y = regions, by = "country") 
  
  # Plot model fit for all countries
  # Store fitted values for VIMC countries
  best_model_output = augment(best_model) %>%
                        mutate(estimate = .fitted) %>%
                        select(-c(model_number, AICc, .model, .resid, .innov, .fitted))
                       
  
  impute_output = impute_dt %>% select(country, d_v_a_id, year, estimate) %>%
                                 mutate(target = NA) %>%
                                 filter(!is.na(d_v_a_id))
  
  estimate_dt = bind_rows(best_model_output, impute_output)
  
  #TODO: Plot to move to plotting.R, fix for d_v_a_name, prettify
  plot_df = estimate_dt
  
  ggplot(data = plot_df, aes(x = year)) +
    geom_point(aes(y = target, colour = "Data")) +
    geom_line(aes(y = .fitted, colour = "Fitted")) +
    labs(y = NULL,
         title = paste("Vaccine impact of", d_v_a_name)
    ) +
    scale_colour_manual(values=c(Data="black",Fitted="#D55E00")) +
    guides(colour = guide_legend(title = NULL))  +
    facet_wrap(~country, ncol = 21)
  
  # Recombine estimated impact with predictor data
  recombine_dt = full_join(data_dt, best_model_output, by = c("country", "d_v_a_id", "year")) %>% 
                  rename(target = target.x)
               
  # Combine original and imputed values 
  result_dt = bind_rows(recombine_dt, impute_dt) %>%
      # Remove predictors...
      select(country, d_v_a_id, year, fvps_cum, impact_cum, target, estimate) %>%
      # Multiply through to obtain cumulative impact over time...
  mutate(impact_impute = fvps_cum * estimate, .after = impact_cum)

  
  covariates_dt = data_dt %>% 
                    as.data.table() %>%
                    select(c(d_v_a_id, target, coverage, HDI, gini, pop_0to14))

  
  # Store the fitted model, the data used, and the result
  fit = list(
    model   = best_model, # N.B> Only for non-imputed
    report  = model_fit,
    choice  = model_choice,
    data    = covariates_dt, 
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

