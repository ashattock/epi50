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
  # - Load impact function data
  # - Use linear function for cases with insufficient data
  # - Assume linear impact for any coverage PRIOR to VIMC estimates (<2000)
  # - Assume linear coverage up to start of WIISE data (<1980)
  # - Apply impact functions to historical coverage data
  # - Job done
  
  # ---- Load impact function data ----

  # coef_dt = read_rds("impact", "coef") 

  # best_dt = read_rds("impact", "best_model")
  
  browser()
  
  
  
  
  # ---- Apply impact functions to historical data ----
  
  browser()
  
  history_dt
  
  evaluate_impact_function(x = history_dt$fvps)
  

}

# ---------------------------------------------------------
# # Evaluate impact function given FVPs
# ---------------------------------------------------------
evaluate_impact_function = function(country = NULL, d_v_a_id = NULL, x = NULL) {
  
  # NOTE: x = FVPs, y = impact (ie deaths averted)
  
  # ---- Load stuff ----
  
  # Load results: best fit functions and associtaed coefficients
  best_dt = read_rds("impact", "best_model")
  coef_dt = read_rds("impact", "coef")
  
  # ---- Interpret trivial arguments ----
  
  # Countrues to evaluate
  if (is.null(country))
    country = unique(best_dt$country)
  
  # D-V-A combinations to evaluate
  if (is.null(d_v_a_id))
    d_v_a_id = unique(best_dt$d_v_a_id)
  
  # Points at which to evaluate
  if (is.null(x)) {
    x_pts = seq(0, o$per_person * o$eval_x_scale, length.out = 101)
    x = data.table(expand_grid(country, d_v_a_id, fvps = x_pts))
  }
  
  # ---- Evaluate best fit model and coefficients ----
  
  # Function to valuate best coefficients
  eval_fn = function(i) {
    
    # Datatable row of interest
    a = best_coef[i, ]
    
    # Extract x from input datatable
    x = x[country == a$country & d_v_a_id == a$d_v_a_id, fvps]
    
    # Call function of interest with these coefficients
    y = do.call(fn_set()[[a$fn]], c(list(x = x), a$par[[1]]))
    
    # Output in datatable form
    out = cbind(a, x = x, y = y)
  }
  
  # Best coefficients
  best_coef = coef_dt %>%
    filter(country %in% !!country, 
           d_v_a_id   %in% !!d_v_a_id) %>%
    inner_join(y  = best_dt, 
               by = c("country", "d_v_a_id", "fn")) %>%
    mutate(par = as.list(setNames(value, coef))) %>% 
    group_by(country, d_v_a_id, fn) %>% 
    summarise(par = list(par)) %>% 
    ungroup() %>% 
    as.data.table()
  
  # Evaluate x to get best fit y for each row in best_coef
  best_fit = nrow(best_coef) %>%
    seq_len() %>%
    lapply(eval_fn) %>%
    rbindlist() %>%
    # Transform impact to real scale...
    mutate(y = y / o$impact_scaler) %>%
    select(country, d_v_a_id, fn,
           fvps   = x, 
           impact = y)
  
  return(best_fit)
}

