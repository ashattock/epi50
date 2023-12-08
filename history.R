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
run_history = function() {
  
  # Only continue if specified by do_step
  if (!is.element(5, o$do_step)) return()
  
  message("* Calculating impact of historical coverage")
  
  # ---- Extract FVP data to be evaluated ----
  
  message(" - Preparing historical coverage")
  
  # Population size of each country over time
  pop_dt = table("wpp_pop") %>%
    group_by(country, year) %>%
    summarise(pop = sum(pop)) %>%
    ungroup() %>%
    as.data.table()
  
  # Extract FVPs over time
  #
  # NOTE: We do not need cumulative FVPs yet, these are only for evaluating
  #       impact fns, and we first need to subset what is to be evalated.
  fvps_dt = table("coverage") %>%
    # Append d_v_a details...
    left_join(y  = table("v_a"),
              by = "v_a_id") %>%
    left_join(y  = table("d_v_a"),
              by = c("vaccine", "activity")) %>%
    # Summarise over age...
    group_by(country, d_v_a_id, year) %>%
    summarise(fvps_abs = sum(fvps)) %>%
    ungroup() %>%
    # Scale results to per capita...
    left_join(y  = pop_dt, 
              by = c("country", "year")) %>%
    select(country, d_v_a_id, year, pop, fvps_abs) %>%
    mutate(fvps_rel = fvps_abs / pop) %>%
    as.data.table()
  
  # From which years have impact functions been fit from
  start_fit_dt = read_rds("impact", "data") %>%
    select(country, d_v_a_id, year) %>%
    group_by(country, d_v_a_id) %>%
    slice_min(year, with_ties = FALSE) %>%
    ungroup() %>%
    rename(start_year = year) %>%
    as.data.table()
  
  # FVPs to evaulate using impact functions 
  eval_dt = 
    # Begin for full factorial of points...
    expand_grid(
      country   = all_countries(), 
      d_v_a_id  = table("d_v_a")$d_v_a_id, 
      eval_year = o$years) %>%
    # Remove years prior to fit start year...
    inner_join(y  = start_fit_dt, 
               by = c("country", "d_v_a_id")) %>%
    filter(eval_year >= start_year) %>%
    select(country, d_v_a_id, year = eval_year) %>%
    # Append FVPs...
    inner_join(y  = fvps_dt, 
               by = c("country", "d_v_a_id", "year")) %>%
    # Cumulative sum FVPs...
    arrange(country, d_v_a_id, year) %>%
    group_by(country, d_v_a_id) %>%
    mutate(fvps_cum = cumsum(fvps_rel)) %>%
    ungroup() %>%
    as.data.table()
  
  # ---- Evaluate impact and back-project -----
  
  message(" - Evaluating impact functions")
  
  # Evaluate impact of relevant FVPs
  result_fit_dt = eval_dt %>%
    # Rename solely for use in evaluation fn...
    rename(fvps = fvps_cum) %>%
    evaluate_impact_function() %>%
    # Convert back to meaningful names...
    rename(fvps_cum   = fvps, 
           impact_cum = impact) %>%
    # Reverse cumsum to derive annual relative impact...
    group_by(country, d_v_a_id) %>%
    mutate(impact_rel = rev_cumsum(impact_cum)) %>%
    ungroup() %>%
    # Rescale back to population scale...
    mutate(impact_abs = impact_rel * pop) %>%
    as.data.table()
  
  # Initial impact per FVPs - used to back project
  #
  # NOTE: Idea behind init_impact_years is to smooth out any 
  #       initially extreme or jumpy impact ratios
  init_ratio_dt = result_fit_dt %>%
    group_by(country, d_v_a_id) %>%
    # Take the first init_impact_years years...
    slice_min(order_by  = year, 
              n         = o$init_impact_years, 
              with_ties = FALSE) %>%
    # Take the mean over these first years...
    summarise(impact_mean = mean(impact_abs), 
              fvps_mean   = mean(fvps_abs)) %>%
    ungroup() %>%
    # Calculate the mean initial ratio...
    mutate(init_ratio = impact_mean / fvps_mean) %>%
    select(country, d_v_a_id, init_ratio) %>%
    as.data.table()
  
  # init_ratio_dt %>%
  #   slice_max(init_ratio, n = 10)
  # 
  # plot_dt = init_ratio_dt %>%
  #   append_d_v_a_name() %>%
  #   filter(init_ratio > 0)
  # 
  # g = ggplot(plot_dt) + 
  #   aes(x = init_ratio, 
  #       y = after_stat(scaled), 
  #       colour = d_v_a_name, 
  #       fill   = d_v_a_name) + 
  #   geom_density(alpha = 0.2) + 
  #   facet_wrap(~d_v_a_name) +
  #   scale_x_continuous(trans  = "log10", 
  #                      labels = scientific)
  
  # Apply the 
  result_dt = result_fit_dt %>%
    select(country, d_v_a_id, year, impact_abs) %>%
    # Join with full FVPs data...
    full_join(y  = fvps_dt, 
              by = c("country", "d_v_a_id", "year")) %>%
    select(country, d_v_a_id, year, 
           fvps   = fvps_abs, 
           impact = impact_abs) %>%
    arrange(country, d_v_a_id, year) %>%
    # Append impact ratio...
    left_join(y  = init_ratio_dt, 
              by = c("country", "d_v_a_id")) %>%
    replace_na(list(init_ratio = 0)) %>%
    # Ratio of past impact assumed consistent with initial years...
    mutate(impact = ifelse(
      test = is.na(impact), 
      yes  = fvps * init_ratio, 
      no   = impact)) %>%
    select(-init_ratio)
  
  # Save results to file
  save_rds(result_dt, "results", "results") 
  
  # ---- Plot results ----
  
  # Plot historical impact over time
  plot_history()
}

# ---------------------------------------------------------
# # Evaluate impact function given FVPs
# ---------------------------------------------------------
evaluate_impact_function = function(eval_dt = NULL) {
  
  # ---- Load stuff ----
  
  # Load results: best fit functions and associtaed coefficients
  best_dt = read_rds("impact", "best_model")
  coef_dt = read_rds("impact", "coef")
  
  # Best coefficients
  best_coef = coef_dt %>%
    inner_join(y  = best_dt, 
               by = c("country", "d_v_a_id", "fn")) %>%
    mutate(par = as.list(setNames(value, coef))) %>% 
    group_by(country, d_v_a_id, fn) %>% 
    summarise(par = list(par)) %>% 
    ungroup() %>% 
    arrange(country, d_v_a_id) %>% 
    as.data.table()
  
  # ---- Interpret trivial argument ----
  
  # Countrues to evaluate
  if (is.null(eval_dt)) {
    
    # Default points at which to evaluate
    x_eval = seq(0, o$eval_x_scale, length.out = 101)
    
    # Construct evaluation datatable
    eval_dt = best_coef %>%
      select(country, d_v_a_id) %>%
      expand_grid(fvps = x_eval) %>%
      as.data.table()
  }
  
  # ---- Evaluate best fit model and coefficients ----
  
  # Function to valuate best coefficients
  eval_fn = function(i) {
    
    # message(paste(c_d_v_a[i, ], collapse = ", "))
    
    # Exract FVPs for this c_d_v_a
    data = c_d_v_a[i, ] %>%
      inner_join(y  = eval_dt, 
                 by = c("country", "d_v_a_id"))
    
    # Fitted function and parameters
    pars = c_d_v_a[i, ] %>%
      inner_join(y  = best_coef, 
                 by = c("country", "d_v_a_id"))
    
    # Load fitted function
    fn = fn_set()[[pars$fn]]
    
    # Call function with fitted coefficients
    impact = do.call(fn, c(list(x = data$fvps), pars$par[[1]]))
    
    # Output result in datatable form
    result = cbind(data, impact)
    
    return(result)
  }
  
  # All country - dva combos to evaluate
  c_d_v_a = eval_dt %>%
    inner_join(y  = best_coef,
               by = c("country", "d_v_a_id")) %>%
    select(country, d_v_a_id) %>%
    unique()
  
  # Apply the evaluation funtion
  result_dt = seq_row(c_d_v_a) %>%
    lapply(eval_fn) %>%
    rbindlist() %>%
    # Transform impact to real scale...
    mutate(impact = impact / o$impact_scaler)
  
  return(result_dt)
}

