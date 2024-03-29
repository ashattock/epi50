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
  
  # Eausy reference external d-v-a IDs
  extern_id = table("d_v_a")[source == "extern", d_v_a_id]
  
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
    arrange(d_v_a_id, country, year) %>%
    filter(!d_v_a_id %in% extern_id) %>%
    # Append impact ratio...
    left_join(y  = initial_ratio_dt, 
              by = c("d_v_a_id", "country")) %>%
    replace_na(list(initial_ratio = 0)) %>%
    # Ratio of past impact assumed consistent with initial years...
    mutate(impact = ifelse(
      test = is.na(impact), 
      yes  = fvps * initial_ratio, 
      no   = impact)) %>%
    select(d_v_a_id, country, year, impact) %>%
    as.data.table()
  
  # ---- Append external models ----
  
  message(" > Appending external models")
  
  # Load up results from external models
  extern_dt = table("extern_estimates") %>%
    lazy_dt() %>%
    rename(impact = !!paste1(metric, "averted")) %>%
    # Summarise results over age...
    group_by(d_v_a_id, country, year) %>%
    summarise(impact = sum(impact)) %>%
    ungroup() %>%
    as.data.table()
  
  # Concatenate results from external models
  result_dt = back_project_dt %>%
    rbind(extern_dt) %>%
    arrange(d_v_a_id, country, year)
  
  # Save results to file
  save_rds(result_dt, "history", "burden_averted", metric) 
  
  # ---- Supporting results ----
  
  # Only approporiate/necessary when computing deaths averted
  if (metric == "deaths") {
    
    # ---- Age structure of impact ----
    
    message(" > Summarising age structure")
    
    # Load all results and summarise age at impact
    age_impact_dt = table("vimc_estimates") %>%
      # Append all results...
      bind_rows(table("static_estimates")) %>%
      bind_rows(table("extern_estimates")) %>%
      select(d_v_a_id, country, year, age, impact = deaths_averted) %>%
      # Pregnancy vaccines are a special case...
      left_join(y  = table("d_v_a"), 
                by = "d_v_a_id") %>%
      mutate(impact = ifelse(
        test = grepl("_px$", vaccine) & age > 0, 
        yes  = 0, 
        no   = impact)) %>%
      # Attribute all impact to infants...
      mutate(impact = ifelse(
        test = grepl("_px$", vaccine) & age == 0, 
        yes  = 1, 
        no   = impact)) %>%
      # Summarise deaths over space and time...
      lazy_dt() %>%
      group_by(d_v_a_id, year, age) %>%
      summarise(value = sum(abs(impact))) %>%
      ungroup() %>%
      # Normalise absolute numbers...
      group_by(d_v_a_id, year) %>%
      mutate(scaler = value / sum(value)) %>%
      ungroup() %>%
      replace_na(list(scaler = 0)) %>%
      # Mean over time...
      group_by(d_v_a_id, age) %>%
      summarise(scaler = mean(scaler)) %>%
      ungroup() %>%
      # Normalise means...
      group_by(d_v_a_id) %>%
      mutate(scaler = scaler / sum(scaler)) %>%
      ungroup() %>%
      as.data.table()
    
    # Expand for each age
    age_effect_dt = table("d_v_a") %>%
      select(d_v_a_id) %>%
      expand_grid(age = o$ages) %>%
      left_join(y  = age_impact_dt, 
                by = c("d_v_a_id", "age")) %>%
      replace_na(list(scaler = 0)) %>%
      as.data.table()
    
    # g = age_effect_dt %>%
    #   group_by(d_v_a_id) %>%
    #   mutate(scaler = cumsum(scaler)) %>%
    #   ungroup() %>%
    #   format_d_v_a_name() %>%
    #   ggplot() +
    #   aes(x = age,
    #       y = scaler,
    #       colour = d_v_a_name) +
    #   geom_line() +
    #   facet_wrap(~d_v_a_name) +
    #   xlim(0, 50)
    
    # Save to tables cache
    save_table(age_effect_dt, "impact_age_multiplier")
    
    # ---- Convert deaths to YLL ----
    
    message(" > Calculating years of life lost")
    
    # Apply age structure and calculate YLL from deaths
    yll_dt = result_dt %>%
      lazy_dt() %>%
      # Append life expectancy...
      left_join(y  = table("wpp_life_exp"), 
                by = c("country", "year")) %>%
      select(-age) %>%
      # Apply age structure...
      left_join(y  = table("impact_age_multiplier"), 
                by = "d_v_a_id", 
                relationship = "many-to-many") %>%
      mutate(deaths = impact * scaler) %>%
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
      arrange(d_v_a_id, country, year) %>%
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
  
  # ---- Demography ----
  
  # Population as per WPP - needed to convert to rates
  pop_dt = table("wpp_pop") %>%
    filter(age <= age_bound) %>%
    group_by(country, year) %>%
    summarise(pop = sum(pop)) %>%
    ungroup() %>%
    as.data.table()
  
  # Child deaths as recorded by WPP
  deaths_dt = table("wpp_deaths") %>%
    filter(age <= age_bound) %>%
    group_by(country, year) %>%
    summarise(deaths = sum(deaths)) %>%
    ungroup() %>%
    as.data.table()
  
  # ---- Age-structured deaths averted ----

  # Vaccine impact disaggregated by age
  age_effect = table("impact_age_multiplier") %>%
    filter(age <= age_bound + 1) %>%
    group_by(d_v_a_id) %>%
    summarise(scaler = sum(scaler)) %>%
    ungroup() %>%
    as.data.table()
  
  # g = ggplot(age_effect %>%
  #              format_d_v_a_name()) +
  #   aes(fill = d_v_a_name,
  #       x    = d_v_a_name,
  #       y    = scaler) +
  #   geom_col() +
  #   theme(axis.text.x = element_text(
  #     angle = 50, hjust = 1))
  
  # Estimated child deaths averted by vaccination
  averted_dt = read_rds("history", "burden_averted_deaths") %>%
    left_join(y  = age_effect, 
              by = "d_v_a_id") %>%
    mutate(impact_age = impact * scaler) %>%
    # Summarise over d-v-a...
    group_by(country, year) %>%
    summarise(averted = sum(impact_age)) %>%
    ungroup() %>%
    as.data.table()
  
  # ---- Mortality rates ----
  
  # Construct grouping datatable to be joined to results
  grouping_dt = table("country") %>%
    left_join(y  = table("income_status"), 
              by = "country") %>%
    mutate(none = "none") %>%
    select(group = !!grouping, country, year)
  
  # Calculate mortality rates in each scenario
  mortality_dt = grouping_dt %>%
    left_join(y  = pop_dt, 
              by = c("country", "year")) %>%
    # Join child death estimates...
    left_join(y  = deaths_dt, 
              by = c("country", "year")) %>%
    left_join(y  = averted_dt, 
              by = c("country", "year")) %>%
    # Summarise over group...
    group_by(group, year) %>%
    summarise(pop     = sum(pop), 
              deaths  = sum(deaths), 
              averted = sum(averted)) %>%
    # Calculate child mortality rates...
    mutate(deaths_alt1 = deaths + averted, 
           rate        = deaths      / pop, 
           rate_alt1   = deaths_alt1 / pop) %>%
    group_by(group) %>%
    # Alternative scenario: no improvement in mortality rate...
    mutate(rate_alt2 = rate_alt1[1]) %>%
    ungroup() %>%
    mutate(deaths_alt2 = pop * rate_alt2, 
           .after = deaths_alt1) %>%
    select(-pop, -averted) %>%
    as.data.table()
  
  # rate_alt1 = pmin(rate_alt1, rate_alt1[1])
  
  # ---- Format output ----
  
  # Use more descriptive scenario names
  scenarios = qc(vaccine, no_vaccine, no_other)
  col_names = c("group", "year", scenarios)
  
  # Function to format output datatable
  mortality_format_fn = function(dt, metric) {

    # Selection of metrics and melt to long format
    formated_dt = dt %>%
      select(group, year, starts_with(metric)) %>%
      as_named_dt(col_names) %>%
      pivot_longer(cols = -c(group, year), 
                   names_to = "scenario") %>%
      select(scenario, group, year, value) %>%
      arrange(scenario, group, year) %>%
      as.data.table()
    
    return(formated_dt)
  }
  
  # Store results in list
  mortality = list(
    rate  = mortality_format_fn(mortality_dt, "rate"),
    value = mortality_format_fn(mortality_dt, "deaths"))
  
  return(mortality)
}

