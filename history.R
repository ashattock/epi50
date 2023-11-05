###########################################################
# HISTORY
#
# xxxxxxxxxxxxxxxxx
#
###########################################################

# ---------------------------------------------------------
# Parent function for xxxxxxxxxxx
# ---------------------------------------------------------
run_history = function() {
  
  # Only continue if specified by do_step
  if (!is.element(5, o$do_step)) return()
  
  message("* Calculating impact of historical coverage")
  
  # STEPS
  # - Use linear function for cases with insufficient data
  # - Assume linear impact for any coverage PRIOR to VIMC estimates (<2000)
  # - Assume linear coverage up to start of WIISE data (<1980)
  # - Apply impact functions to historical coverage data
  # - Job done
  
  # ---- Apply impact functions to historical data ----
  
  # Population size of each country over time
  pop_dt = table("wpp_pop") %>%
    group_by(country, year) %>%
    summarise(pop = sum(pop)) %>%
    ungroup() %>%
    as.data.table()
  
  eval_dt = table("coverage") %>%
    group_by(country, v_a_id, year) %>%
    summarise(fvps_temporal = sum(fvps)) %>%
    mutate(fvps_cum = cumsum(fvps_temporal)) %>%
    ungroup() %>%
    left_join(y  = pop_dt, 
              by = c("country", "year")) %>%
    mutate(fvps = fvps_cum / pop) %>%
    left_join(y  = table("v_a"), 
              by = "v_a_id") %>%
    left_join(y  = table("d_v_a"), 
              by = c("vaccine", "activity"), 
              relationship = "many-to-many") %>%
    select(country, d_v_a_id, year, pop,
           fvps_temporal, fvps_cum, fvps) %>%
    as.data.table()
  
  result_dt = evaluate_impact_function(eval_dt) %>%
    # Population weight for population-level impact...
    mutate(impact_cum = impact * pop) %>%
    # Revere cumsum to derive annual impact...
    group_by(country, d_v_a_id) %>%
    mutate(impact_temporal = rev_cumsum(impact_cum)) %>%
    ungroup() %>%
    # Tidy up...
    select(country, d_v_a_id, year, 
           fvps   = fvps_temporal, 
           impact = impact_temporal) %>%
    as.data.table()
  
  # Save results to file
  save_rds(result_dt, "results", "results") 
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
  best_coef = best_dt %>%
    inner_join(y  = coef_dt, 
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
  
  c_d_v_a = eval_dt %>%
    inner_join(y  = best_coef,
               by = c("country", "d_v_a_id")) %>%
    select(country, d_v_a_id) %>%
    unique()
  
  result_dt = nrow(c_d_v_a) %>%
    seq_len() %>%
    lapply(eval_fn) %>%
    rbindlist() %>%
    # Transform impact to real scale...
    mutate(impact = impact / o$impact_scaler)
    
  return(result_dt)
}

