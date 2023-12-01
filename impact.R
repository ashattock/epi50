###########################################################
# IMPACT
#
# An alternative to impact factor for representing impact as
# a function of vaccine coverage.
#
# TODO: Think about age age here... can we justify summarising?
# I think we can, but d-v-a must consistently have targeted
# the same age groups over time.
#
###########################################################

# ---------------------------------------------------------
# Parent function for calculating non-linear impact
# ---------------------------------------------------------
run_impact = function() {
  
  # Only continue if specified by do_step
  if (!is.element(4, o$do_step)) return()
  
  message("* Fitting impact functions")
  
  # ---- FVPs and impact estimates ----
  
  message(" - Preparing FVP-impact data")
  
  # Prepare impact-FVP data to fit to
  data_dt = get_impact_data()
  
  # Exploratory plots of data used to fit impact functions
  plot_impact_data()
  
  # ---- Determine best fitting model ----
  
  message(" - Evaluating impact functions")
  
  # Country-disease-vaccine-activity combinations
  c_d_v_a = data_dt %>%
    select(country, d_v_a_id) %>%
    unique()
  
  # Number of rows in this table
  n_row = nrow(c_d_v_a)
  
  # Functions we'll attempt to fit with
  fns = fn_set()
  
  # Preallocate lists to store results
  coef = aic = r2 = vector('list', n_row)
  
  # Initiate progress bar
  pb = start_progress_bar(n_row)
  
  # Iterate through as instances
  for (i in seq_len(n_row)) {
    x = c_d_v_a[i, ]
    
    # Attempt to fit all fns and determine most suitable
    result = get_best_model(fns, data_dt, x$country, x$d_v_a_id)
    
    # Store results
    coef[[i]] = result$coef
    aic[[i]]  = result$aic
    r2[[i]]   = result$r2
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  
  # Close progress bar
  close(pb)
  
  # ---- Store results ----
  
  message(" - Selecting best functions")
  
  # Collapse all fn coefficients into single datatable
  coef_dt = rbindlist(coef)
  
  # Save to file
  save_rds(coef_dt, "impact", "coef") 
  
  # Collapse suitability metrics into single datatable
  aic_dt = rbindlist(aic, fill = TRUE)
  r2_dt  = rbindlist(r2,  fill = TRUE)
  
  # Save to file
  save_rds(aic_dt, "impact", "aic")
  save_rds(r2_dt,  "impact", "r2")
  
  if (o$model_metric == "aicc") {
    
    # Extract best fitting function based on key metrics
    best_dt =
      # Bind AICc and R-squared results
      rbind(aic_dt %>% pivot_longer(cols = names(fns)),
            r2_dt  %>% pivot_longer(cols = names(fns))) %>%
      pivot_wider(names_from  = c(name, metric),
                  values_from = value) %>%
      # If we hit thresholds, trivialise result...
      # mutate(lin0_aicc = ifelse(lin0_r2 > o$r2_threshold0, -Inf, lin0_aicc),
      #        log3_aicc = ifelse(log3_r2 < o$r2_threshold1, Inf, log3_aicc)) %>%
      select(country, d_v_a_id, lin0 = lin0_aicc, log3 = log3_aicc) %>%
      # Select function with lowest AICc...
      pivot_longer(cols = names(fns),
                   names_to = "fn") %>%
      group_by(country, d_v_a_id) %>%
      slice_min(value, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      # Tidy up...
      select(country, d_v_a_id, fn) %>%
      as.data.table()
  }
  
  if (o$model_metric == "r2") {
    
    # Extract best fitting function based on key metrics
    best_dt = 
      # Bind AICc and R-squared results
      rbind(aic_dt %>% pivot_longer(cols = names(fns)), 
            r2_dt  %>% pivot_longer(cols = names(fns))) %>%
      pivot_wider(names_from  = c(name, metric), 
                  values_from = value) %>%
      select(country, d_v_a_id, lin0 = lin0_r2, log3 = log3_r2) %>%
      # Select function with highest R-squared...
      pivot_longer(cols = names(fns), 
                   names_to = "fn") %>%
      group_by(country, d_v_a_id) %>%
      slice_max(value, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      # Tidy up...
      select(country, d_v_a_id, fn) %>%
      as.data.table()
  }
  
  # Save to file
  save_rds(best_dt, "impact", "best_model")
  
  # ---- Plot results ----
  
  # NOTE: All plotting functionality lives in plotting.R
  
  # Plot function selection statistics
  plot_model_selection()
  
  # Plot impact function evaluation
  plot_model_fits()
}

# ---------------------------------------------------------
# xxxxxxxx
# ---------------------------------------------------------
get_impact_data = function() {
  
  # Population size of each country over time
  pop_dt = table("wpp_pop") %>%
    group_by(country, year) %>%
    summarise(pop = sum(pop)) %>%
    ungroup() %>%
    as.data.table()
  
  # Load impact estimates from VIMC (inlcuding imputed)
  #
  # NOTE: Result of imputation is already in cumulative form
  vimc_dt = read_rds("impute", "impute_result") %>%
    # Scale results to per capita...
    left_join(y  = pop_dt, 
              by = c("country", "year")) %>%
    mutate(fvps   = fvps / pop, 
           impact = impact / pop) %>%
    select(-pop)
  
  # Load non-modelled impact estimates
  gbd_dt = read_rds("non_modelled", "deaths_averted_vaccine") %>%
    # Scale results to per capita...
    left_join(y  = pop_dt, 
              by = c("country", "year")) %>%
    mutate(fvps   = fvps / pop, 
           impact = impact / pop) %>%
    select(-pop) %>%
    # Convert to cumulative FVP and impact...
    group_by(country, d_v_a_id) %>%
    mutate(fvps_cum   = cumsum(fvps), 
           impact_cum = cumsum(impact)) %>%
    ungroup() %>%
    # Tidy up...
    select(country, d_v_a_id, year,
           fvps   = fvps_cum,
           impact = impact_cum) %>%
    as.data.table()
  
  # Impact estimates per capita from all sources
  data_dt = rbind(vimc_dt, gbd_dt) %>%
    # Impact per FVP...
    mutate(impact_fvp = impact / fvps) %>%
    # Ingore cases with a single data point...
    add_count(country, d_v_a_id) %>%
    filter(n > 1) %>%
    select(-n)
  
  # Save to file
  save_rds(data_dt, "impact", "data") 
  
  return(data_dt)
}

# ---------------------------------------------------------
# Set of functions to fit - we'll determine the 'best'
# ---------------------------------------------------------
fn_set = function(dict = FALSE) {
  
  # Set of statistical models / functions we want to test
  out = list(
    lin0 = function(x, a)       y = a*x,
    # lin  = function(x, a, b)    y = a*x + b,
    # quad = function(x, a, b, c) y = a*x^2 + b*x + c,
    log3 = function(x, a, b, c) y = logistic(x, a, b, upper = c))
  
  # Alternative functionality - return dictionary
  if (dict == TRUE)
    out = c(
      lin0 = "Gradiant", 
      # lin  = "Linear", 
      # quad = "Quadratic", 
      log3 = "Logistic")
  
  return(out)
}

# ---------------------------------------------------------
# Parent function to determine best fitting function
# ---------------------------------------------------------
get_best_model = function(fns, data_dt, country, d_v_a_id) {
  
  # Reduce data down to what we're interested in
  c_d_v_a_dt = data_dt %>%
    filter(country  == !!country, 
           d_v_a_id == !!d_v_a_id) %>%
    select(x = fvps,
           y = impact) %>%
    # Multiply impact for more consistent x-y scales...
    mutate(y = y * o$impact_scaler)
  
  # Do not fit if insufficient data
  if (nrow(c_d_v_a_dt) <= 3)
    return()
  
  # Declare x and y values for which we want to determine a relationship
  x = c_d_v_a_dt$x
  y = c_d_v_a_dt$y
  
  # Use optim algorithm to get good starting point for MLE
  start = credible_start(fns, x, y)
  
  # Run MLE from this starting point
  fit = run_mle(fns, start, x, y)
  
  # Determine AICc value for model suitability
  result = model_quality(fns, fit, country, d_v_a_id, x, y)
  
  return(result)
}

# ---------------------------------------------------------
# Determine credible starting points for MLE - it needs it
# ---------------------------------------------------------
credible_start = function(fns, x, y, plot = FALSE) {
  
  # Initialise starting point for sigma in likelihood function
  s0 = 1  # This is essentially a placeholder until run_mle 
  
  # Let's start with any old points
  start = list(
    lin0 = list(s = s0, a = 1),
    # lin  = list(s = s0, a = 1, b = 1),
    log3 = list(s = s0, a = 1, b = 1, c = 1))
  
  # Define an objective function to minimise - sum of squares
  obj_fn = function(fn, ...) {
    
    # Squared difference
    diff_sq = (y - fns[[fn]](x, ...)) ^ 2
    
    # The sum of the squared difference
    obj_val = list(y = sum(diff_sq))
    
    return(obj_val)
  }
  
  # Define model-specific calls to objective function
  asd_fn = list(
    lin0 = function(p, args) obj_fn("lin0", a = p[1]),
    # lin  = function(p, args) obj_fn("lin",  a = p[1], b = p[2]),
    log3 = function(p, args) obj_fn("log3", a = p[1], b = p[2], c = p[3]))
  
  # Initiate list for storing plotting datatables
  plot_list = list()
  
  # Points to evaluate when plotting
  x_eval = seq(0, max(x), length.out = 101)
  
  # Iterate through stats models
  for (fn in names(fns)) {
    
    # Number of parameters for this model
    n_pars = length(start[[fn]]) - 1
    
    # Run ASD optimisation algorithm
    optim = asd(
      fn   = asd_fn[[fn]],
      x0   = rep(1, n_pars),
      lb   = rep(0, n_pars), 
      ub   = rep(Inf, n_pars),
      max_iters = 10000)
    
    # Overwrite starting point with optimal parameters
    start[[fn]] = c(s0, optim$x) %>%
      setNames(names(start[[fn]])) %>%
      as.list()
    
    # Store diagnostic plotting data if desired
    if (plot == TRUE) {
      
      # Construct plotting datatable for this model
      plot_args = c(list(x_eval), as.list(optim$x))
      plot_list[[fn]] = data.table(
        x = x_eval,
        y = do.call(fns[[fn]], as.list(plot_args)), 
        fn = fn)
    }
  }
  
  # Create diagnostic plot if desired
  if (plot == TRUE) {
    
    # Plot the quality of fit for each model starting point
    g_asd = ggplot(rbindlist(plot_list)) +
      aes(x = x, y = y) +
      geom_line(aes(colour = fn)) + 
      geom_point(data = data.table(x = x, y = y),
                 colour = "black")
  }
  
  return(start)
}

# ---------------------------------------------------------
# Fit MLE for each fn using prevriously determined start point
# ---------------------------------------------------------
run_mle = function(fns, start, x, y, plot = FALSE) {
  
  # Log likelihood function to maximise
  likelihood = function(s, y_pred)
    ll = -sum(dnorm(x = y, mean = y_pred, sd = s, log = TRUE))
  
  # Annoyingly we need to name likelihood inputs, hence wrapper functions
  model = list(
    lin0 = function(s, a)       likelihood(s, fns$lin0(x, a)),
    # lin  = function(s, a, b)    likelihood(s, fns$lin(x, a, b)),
    log3 = function(s, a, b, c) likelihood(s, fns$log3(x, a, b, c)))
  
  # Initiate list for storing plotting datatables
  fit = plot_list = list()
  
  # Points to evaluate when plotting
  x_eval = seq(0, max(x), length.out = 101)
  
  # Iterate through stats models
  for (fn in names(fns)) {
    
    # Vector of lower bounds
    n_pars = length(start[[fn]]) - 1  # Zero for function coefficients
    l_bnds = c(1e-1, rep(0, n_pars))  # Small but non-zero for sigma
    
    # Attempt to fit MLE model using prevriously determined start point
    fit_result = tryCatch(
      expr  = mle(minuslogl = model[[fn]], 
                  start = start[[fn]], 
                  lower = l_bnds,
                  nobs  = length(y)),
      
      # Return trivial if MLE failed
      error   = function(e) NULL,
      warning = function(w) NULL)
    
    # Store results if we've had success
    if (!is.null(fit_result)) {
      
      # Store fit and extract coefficients for plotting
      fit[[fn]] = fit_result
      coef      = fit_result@coef[-1]
      
      # Store diagnostic plotting data if desired
      if (plot == TRUE) {
        
        # Construct plotting datatable for this model
        plot_args = c(list(x_eval), as.list(coef))
        plot_list[[fn]] = data.table(
          x = x_eval,
          y = do.call(fns[[fn]], as.list(plot_args)), 
          fn = fn)
      }
    }
  }
  
  # Create diagnostic plot if desired
  if (plot == TRUE) {
    
    # Plot the quality of fit for each model
    g_mle = ggplot(rbindlist(plot_list)) +
      aes(x = x, y = y) +
      geom_line(aes(colour = fn)) + 
      geom_point(data = data.table(x = x, y = y),
                 colour = "black")
  }
  
  return(fit)
}

# ---------------------------------------------------------
# Determine model quality - primarily this is via AICc
# ---------------------------------------------------------
model_quality = function(fns, fit, country, d_v_a_id, x, y) {
  
  # Return out if no fits succesful 
  if (length(fit) == 0)
    return()
  
  # ---- Extract coefficients ----
  
  # Coefficients for each successful model
  coef = unlist(lapply(fit, function(a) a@coef[-1]))
  coef = data.table(var   = names(coef), 
                    value = coef) %>%
    separate(var, c("fn", "coef")) %>%
    mutate(country = country, 
           d_v_a_id     = d_v_a_id, 
           .before = 1)
  
  # ---- AICc ----
  
  # Calculate AIC - adjusted for sample size
  aic = sapply(fit, AICc) %>%
    as.list() %>%
    as.data.table() %>%
    mutate(country = country, 
           d_v_a_id   = d_v_a_id, 
           metric  = "aicc", 
           .before = 1)
  
  # ---- R squared ----
  
  # Function to compute adjusted R-squared
  r2_fn = function(a) {
    
    # Extract optimal coefficients in list format
    args = as.list(fit[[a]]@coef[-1])
    
    # Evaluate at all points in x
    y_eval = do.call(fns[[a]], c(list(x), args))
    
    # Sample size and number of parameters
    n = fit[[a]]@nobs
    k = length(fit[[a]]@details$par)  - 1
    
    # Calculate R-squared (coefficient of determination)
    r2 = suppressWarnings(cor(y, y_eval) ^ 2)
    
    # Adjusted R-squared - considers number of parameters
    r2_adjusted = 1 - ((1 - r2) * (n - 1) / (n - k - 1))
    
    return(r2_adjusted)
  }
  
  # Apply function to succeful models
  r2 = names(fit) %>%
    lapply(function(a) r2_fn(a)) %>%
    setNames(names(fit)) %>%
    as.data.table() %>%
    mutate(country = country, 
           d_v_a_id   = d_v_a_id, 
           metric  = "r2", 
           .before = 1)
  
  return(list(coef = coef, aic = aic, r2 = r2))
}

# ---------------------------------------------------------
# Use AICc rather than AIC to reduce overfitting
# ---------------------------------------------------------
AICc = function(x) {
  
  # See en.wikipedia.org/wiki/Akaike_information_criterion
  
  # Sample size
  n = x@nobs
  
  # Number of parameters
  k = length(x@details$par) - 1
  
  # The usual AIC term
  aic_term = stats::AIC(x)
  
  # An additional penalty term for small sample size
  pen_term = (2*k^2 + 2*k) / (n - k - 1)
  
  # Sum these terms
  aicc = aic_term # + pen_term
  
  return(aicc)
}

