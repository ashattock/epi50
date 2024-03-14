###########################################################
# HISTORY
#
# Bring everything together: evaluate non-linear impact 
# functions to determine vaccine impact. Back project using
# ratio of first few years of impact estimates.
#
###########################################################

# ---------------------------------------------------------
# Use impact functions to calculate historical impact
# ---------------------------------------------------------
run_history = function(metric) {
  
  # Only continue if specified by do_step
  if (!is.element(6, o$do_step)) return()
  
  message("* Calculating impact of historical coverage: ", metric)
  
  # ---- Extract FVP data to be evaluated ----
  
  message(" > Preparing historical coverage")
  
  # Population size of each country over time
  pop_dt = table("wpp_pop") %>%
    lazy_dt() %>%
    group_by(country, year) %>%
    summarise(pop = sum(pop)) %>%
    ungroup() %>%
    as.data.table()
  
  # Extract FVPs over time
  #
  # NOTE: We do not need cumulative FVPs yet, these are only for evaluating
  #       impact fns, and we first need to subset what is to be evaluated.
  fvps_dt = table("coverage") %>%
    lazy_dt() %>%
    # Summarise over age...
    group_by(d_v_a_id, country, year) %>%
    summarise(fvps_abs = sum(fvps)) %>%
    ungroup() %>%
    # Scale results to per capita...
    left_join(y  = pop_dt, 
              by = c("country", "year")) %>%
    select(d_v_a_id, country, year, pop, fvps_abs) %>%
    mutate(fvps_rel = fvps_abs / pop) %>%
    as.data.table()
  
  # From which years have impact functions been fit from
  start_fit_dt = read_rds("impact", "impact", metric, "data") %>%
    lazy_dt() %>%
    select(d_v_a_id, country, year) %>%
    group_by(d_v_a_id, country) %>%
    slice_min(year, with_ties = FALSE) %>%
    ungroup() %>%
    rename(start_year = year) %>%
    as.data.table()
  
  # FVPs to evaulate using impact functions 
  eval_dt = 
    # Begin for full factorial of points...
    expand_grid(
      d_v_a_id  = table("d_v_a")$d_v_a_id, 
      country   = all_countries(), 
      eval_year = o$years) %>%
    # Remove years prior to fit start year...
    inner_join(y  = start_fit_dt, 
               by = c("d_v_a_id", "country")) %>%
    filter(eval_year >= start_year) %>%
    select(d_v_a_id, country, year = eval_year) %>%
    # Append FVPs...
    inner_join(y  = fvps_dt, 
               by = c("d_v_a_id", "country", "year")) %>%
    # Cumulative sum FVPs...
    arrange(d_v_a_id, country, year) %>%
    group_by(d_v_a_id, country) %>%
    mutate(fvps_cum = cumsum(fvps_rel)) %>%
    ungroup() %>%
    as.data.table()
  
  # ---- Evaluate impact functions -----
  
  message(" > Evaluating impact functions")
  
  # Evaluate impact of relevant FVPs
  result_fit_dt = eval_dt %>%
    # Rename solely for use in evaluation fn...
    rename(fvps = fvps_cum) %>%
    evaluate_impact_function(metric = metric) %>%
    # Convert back to meaningful names...
    rename(fvps_cum   = fvps, 
           impact_cum = impact) %>%
    # Reverse cumsum to derive annual relative impact...
    lazy_dt() %>%
    group_by(d_v_a_id, country) %>%
    mutate(impact_rel = rev_cumsum(impact_cum)) %>%
    ungroup() %>%
    # Rescale back to population scale...
    mutate(impact_abs = impact_rel * pop) %>%
    as.data.table()
  
  # ---- Back project using initial impact ratios -----
  
  message(" > Back projecting")
  
  # Initial impact per FVPs - used to back project
  #
  # NOTE: Idea behind init_impact_years is to smooth out any 
  #       initially extreme or jumpy impact ratios
  initial_ratio_dt = result_fit_dt %>%
    lazy_dt() %>%
    group_by(d_v_a_id, country) %>%
    # Take the first init_impact_years years...
    slice_min(order_by  = year, 
              n         = o$init_impact_years, 
              with_ties = FALSE) %>%
    # Take the mean over these first years...
    summarise(impact_mean = mean(impact_abs), 
              fvps_mean   = mean(fvps_abs)) %>%
    ungroup() %>%
    # Calculate the mean initial ratio...
    mutate(initial_ratio = impact_mean / fvps_mean) %>%
    select(d_v_a_id, country, initial_ratio) %>%
    as.data.table()
  
  # Save initial ratio to file for diagnostic plotting
  save_rds(initial_ratio_dt, "history", "initial_ratio", metric) 
  
  # Back project by applying initial ratio 
  back_project_dt = result_fit_dt %>%
    lazy_dt() %>%
    select(d_v_a_id, country, year, impact_abs) %>%
    # Join with full FVPs data...
    full_join(y  = fvps_dt, 
              by = c("d_v_a_id", "country", "year")) %>%
    select(d_v_a_id, country, year, 
           fvps   = fvps_abs, 
           impact = impact_abs) %>%
    arrange(country, d_v_a_id, year) %>%
    # Append impact ratio...
    left_join(y  = initial_ratio_dt, 
              by = c("d_v_a_id", "country")) %>%
    replace_na(list(initial_ratio = 0)) %>%
    # Ratio of past impact assumed consistent with initial years...
    mutate(impact = ifelse(
      test = is.na(impact), 
      yes  = fvps * initial_ratio, 
      no   = impact)) %>%
    select(-initial_ratio) %>%
    as.data.table()
  
  # ---- Append external models ----
  
  message(" > Appending external models")
  
  # Load up results from external models
  extern_dt = table("extern_estimates") %>%
    lazy_dt() %>%
    filter(d_v_a_id %in% table("d_v_a")$d_v_a_id) %>%
    rename(impact = !!paste1(metric, "averted")) %>%
    # Summarise results over age...
    group_by(d_v_a_id, country, year) %>%
    summarise(impact = sum(impact)) %>%
    ungroup() %>%
    # TODO: Update placeholder with actual values
    mutate(fvps = 1, .before = impact) %>%
    as.data.table()
  
  # Concatenate results from external models
  result_dt = back_project_dt %>%
    rbind(extern_dt)
  
  # Save results to file
  save_rds(result_dt, "history", "burden_averted", metric) 
  
  # ---- Convert deaths to YLL ----
  
  # Only approporiate/necessary when computing deaths averted
  if (metric == "deaths") {
    
    message(" > Calculating years of life lost")
    
    # Apply age structure and calculate YLL from deaths
    yll_dt = result_dt %>%
      # Apply age structure...
      expand_grid(impact_age_multiplier()) %>%
      mutate(deaths = impact * scaler) %>%
      select(d_v_a_id, country, year, age, deaths) %>%
      # Append life expectency...
      mutate(life_exp = 80) %>%  # TODO: Use proper WPP data
      # Calculate years of life lost...
      mutate(yll = deaths * pmax(0, life_exp - age)) %>%
      group_by(d_v_a_id, country, year) %>%
      summarise(impact = sum(yll)) %>%
      ungroup() %>%
      # Reappend FVPs for consistent formatting...
      left_join(y  = fvps_dt, 
                by = c("d_v_a_id", "country", "year")) %>%
      rename(fvps = fvps_abs) %>%
      select(all_names(result_dt)) %>%
      as.data.table()
      
    # Save results to file
    save_rds(yll_dt, "history", "burden_averted_yll") 
  }
  
  # ---- Plot outcomes ----
  
  # Plot inital impact ratios used to back project
  plot_impact_fvps(metric, scope = "initial")
}

# ---------------------------------------------------------
# Evaluate impact function given FVPs
# ---------------------------------------------------------
evaluate_impact_function = function(data_dt = NULL, metric = NULL) {
  
  # ---- Load stuff ----
  
  # TEMP: For now take mean, but later sample from posterior
  coef_dt = read_rds("impact", "posteriors", metric) %>% 
    lazy_dt() %>%
    group_by(d_v_a_id, country, fn, param) %>% 
    summarise(value = mean(value)) %>% 
    ungroup() %>% 
    arrange(d_v_a_id, country, param) %>% 
    as.data.table()
  
  # ---- Interpret trivial argument ----
  
  # Countries to evaluate
  if (is.null(data_dt)) {
    
    # Default points at which to evaluate
    x_eval = seq(0, o$eval_x_scale, length.out = 101)
    
    # Construct evaluation datatable
    data_dt = coef_dt %>%
      select(d_v_a_id, country) %>%
      unique() %>%
      expand_grid(fvps = x_eval) %>%
      as.data.table()
  }
  
  # ---- Evaluate best fit model and coefficients ----
  
  # Function to valuate best coefficients
  eval_fn = function(i, ids, fns, data, coef) {
    
    # Index d-v-a country ID
    id = ids[i, ]
    
    # message(paste(id, collapse = ", "))
    
    # Exract FVPs for this ID
    data %<>% inner_join(id, by = c("d_v_a_id", "country"))
    
    # Fitted function and parameters
    coef %<>% inner_join(id, by = c("d_v_a_id", "country"))
    
    # Load fitted function
    fn = fns[[unique(coef$fn)]]
    
    # Call function with fitted coefficients
    impact = fn(x = data$fvps, p = coef$value)
    
    # Output in datatable form
    result = cbind(data, impact)
    
    return(result)
  }
  
  # Load set of functions that may be evaluated
  fns = fn_set()
  
  # All country - dva combos to evaluate
  ids = data_dt %>%
    select(d_v_a_id, country) %>%
    unique() %>%
    inner_join(y  = coef_dt,
               by = c("d_v_a_id", "country")) %>%
    select(d_v_a_id, country) %>%
    unique()
  
  # Apply evaluations in parallel
  if (o$parallel$history)
    result_list = mclapply(
      X    = seq_row(ids), 
      FUN  = eval_fn, 
      ids  = ids,
      fns  = fns,
      data = data_dt, 
      coef = coef_dt,
      mc.cores = o$n_cores,
      mc.preschedule = FALSE)
  
  # Apply evaluations consecutively
  if (!o$parallel$history)
    result_list = lapply(
      X   = seq_row(ids), 
      FUN = eval_fn,
      ids  = ids,
      fns  = fns,
      data = data_dt, 
      coef = coef_dt)
  
  # Squash results into single datatable
  result_dt = rbindlist(result_list) %>%
    # Transform impact to real scale...
    mutate(impact = impact / o$impact_scaler)
  
  return(result_dt)
}

# ---------------------------------------------------------
# Calulate child mortality rates in vaccine and no vaccine scenarios
# ---------------------------------------------------------
mortality_rates = function(age_bound = 0, grouping = "none") {
  
  # NOTE: Options for 'grouping' argument: "none", "region", or "income"
  
  # Construct grouping datatable to be joined to results
  grouping_dt = table("country") %>%
    left_join(y  = table("income_status"), 
              by = "country") %>%
    mutate(none = "none") %>%
    select(country, year, group = !!grouping)
  
  # ---- Demography ----
  
  # Population as per WPP - needed to convert to rates
  pop_dt = table("wpp_pop") %>%
    filter(age <= age_bound) %>%
    left_join(y  = grouping_dt, 
              by = c("country", "year")) %>%
    group_by(group, year) %>%
    summarise(pop = sum(pop)) %>%
    ungroup() %>%
    arrange(group, year) %>%
    as.data.table()
  
  # Child deaths as recorded by WPP
  deaths_dt = table("wpp_deaths") %>%
    filter(age <= age_bound) %>%
    left_join(y  = grouping_dt, 
              by = c("country", "year")) %>%
    group_by(group, year) %>%
    summarise(deaths = sum(deaths)) %>%
    ungroup() %>%
    arrange(group, year) %>%
    as.data.table()
  
  # ---- Age-structured lives and births saved ----
  
  # Vaccine impact disaggregated by age
  age_effect = impact_age_multiplier() # %>%
    # filter(age <= age_bound)
  
  # Estimated child deaths averted by vaccination
  impact_dt = read_rds("history", "burden_averted_deaths") %>%
    expand_grid(age_effect) %>%
    # Summarise over d-v-a...
    lazy_dt() %>%
    group_by(country, year, age) %>%
    summarise(lives_saved = sum(impact * scaler)) %>%
    ungroup() %>%
    # Tidy up...
    arrange(country, year) %>%
    as.data.table()
  
  # fertility_dt  = table("wpp_fertility_total") %>%
  #   select(country, year, fertility_total) %>%
  #   left_join(y  = table("wpp_fertility_age"), 
  #             by = c("country", "year")) %>%
  #   # Total fertility per person (assuming 50% women)...
  #   mutate(fertility = fertility_total * fertility_age / 2) %>%
  #   select(-fertility_total, -fertility_age)
  # 
  # saved_dt = expand_grid(
  #   country = all_countries(), 
  #   year    = o$years,
  #   age     = o$ages,
  #   count   = o$ages) %>%
  #   # xxx...
  #   filter(age < count) %>%
  #   mutate(year_alt = 1 + max(o$ages) - count + year, 
  #          age_alt  = 1 + max(o$ages) - count + age) %>%
  #   select(-count) %>%
  #   filter(year_alt <= max(o$years)) %>%
  #   arrange(country, year, age, age_alt) %>%
  #   # xxx...
  #   left_join(y  = impact_dt, 
  #             by = c("country", "year", "age")) %>%
  #   replace_na(list(lives_saved = 0)) %>%
  #   select(country, year = year_alt, age = age_alt, lives_saved) %>%
  #   # Further, these would not contribute to fertility...
  #   left_join(y  = fertility_dt, 
  #             by = c("country", "year", "age")) %>%
  #   replace_na(list(fertility = 0)) %>%
  #   mutate(births_saved = lives_saved * fertility) %>%
  #   select(-fertility) %>%
  #   # xxx...
  #   mutate(lives_saved = ifelse(age <= age_bound, lives_saved, 0)) %>%
  #   # Apply and summarise over grouping...
  #   left_join(y  = grouping_dt,
  #             by = c("country", "year")) %>%
  #   lazy_dt() %>%
  #   group_by(group, year) %>%
  #   summarise(lives_saved  = sum(lives_saved), 
  #             births_saved = sum(births_saved)) %>%
  #   ungroup() %>%
  #   as.data.table()
  
  # ---- Mortality rates ----
  
  # All child deaths by year
  averted_dt = impact_dt %>%
    left_join(y  = grouping_dt,
              by = c("country", "year")) %>%
    lazy_dt() %>%
    group_by(group, year) %>%
    summarise(deaths_averted = sum(lives_saved)) %>%
    ungroup() %>%
    as.data.table()
  
  # Calculate mortality rates in each scenario
  mortality_dt = pop_dt %>%
    # Join child death estimates...
    left_join(y  = deaths_dt, 
              by = c("group", "year")) %>%
    left_join(y  = averted_dt, 
              by = c("group", "year")) %>%
    mutate(deaths_alt1 = deaths + deaths_averted) %>%
    # Calculate child mortality rates...
    group_by(group) %>%
    mutate(rate      = deaths      / pop, 
           rate_alt1 = deaths_alt1 / pop,
           rate_alt1 = pmin(rate_alt1, rate_alt1[1]),
           rate_alt2 = rate_alt1[1]) %>%
    ungroup() %>%
    mutate(deaths_alt2 = pop * rate_alt2, 
           .after = deaths_alt1) %>%
    select(-pop, -deaths_averted) %>%
    as.data.table()
  
  return(mortality_dt)
}

# ---------------------------------------------------------
# Scaler to estimate vaccination impact disaggregated by age
# ---------------------------------------------------------
impact_age_multiplier = function() {
  
  # TODO: Extract impact distribution by age for each d-v-a
  #       using original impact estimates
  
  # TEMP: Assume impact is highest for infants and decays exponentially
  age_effect = data.table(age = o$ages) %>%
    mutate(scaler = exp(-sqrt(2 * age)), 
           scaler = scaler / sum(scaler)) 
  
  return(age_effect)
}

