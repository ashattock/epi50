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

  message("* Calculating impact per FVP")
  
  # ---- FVPs and impact estimates ----
  
  # Population size of each country over time
  pop_dt = table("wpp_pop") %>%
    group_by(country, year) %>%
    summarise(pop = sum(pop)) %>%
    ungroup() %>%
    as.data.table()
  
  # Impact estimates (VIMC & imputed) per capita
  impact_dt = read_rds("impute", "impute_result") %>%
    left_join(y  = pop_dt, 
              by = c("country", "year")) %>%
    mutate(fvps   = o$per_person * fvps / pop, 
           impact = o$per_person * impact / pop) %>%
    select(-pop) %>%
    # Impact per FVP...
    mutate(impact_fvp = impact / fvps) %>%
    # Ingore cases with a single data point...
    add_count(country, d_v_a_id) %>%
    filter(n > 1) %>%
    select(-n)
  
  # ---- Exploratory plots ----
  
  browser()
  
  # Impact per FVP over time
  g1 = (ggplot(impact_dt) +
          aes(x = year, y = impact_fvp, colour = country) +
          geom_line(show.legend = FALSE) +
          facet_wrap(~d_v_a_id, scales = "free_y")) %>%
    prettify1(save = c("Year", "impact", "FVP"))
  
  # Cumulative FVPs vs cumulative deaths averted
  g2 = (ggplot(impact_dt) + 
           aes(x = fvps, y = impact, colour = country) +
           geom_line(show.legend = FALSE) +
           facet_wrap(~d_v_a_id, scales = "free")) %>%
    prettify1(save = c("FVP", "impact"))
  
  # ---- Determine best fitting model ----
  
  # Country-disease-vaccine-activity combinations
  c_d_v_a = vimc_dt %>%
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
    result = get_best_model(fns, vimc_dt, x$country, x$d_v_a_id)
    
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
get_best_model = function(fns, vimc_dt, country, d_v_a_id) {
  
  # Reduce data down to what we're interested in
  data_dt = prep_data(vimc_dt, country, d_v_a_id)
  
  # Do not fit if insufficient data
  if (nrow(data_dt) <= 3)
    return()
  
  # Declare x and y values for which we want to determine a relationship
  x = data_dt$x
  y = data_dt$y
  
  # Use optim algorithm to get good starting point for MLE
  start = prep_start(fns, x, y)
  
  # Run MLE from this starting point
  fit = run_mle(fns, start, x, y)
  
  # Determine AICc value for model suitability
  result = model_quality(fns, fit, country, d_v_a_id, x, y)
  
  return(result)
}

# ---------------------------------------------------------
# Prepare data for fitting - for this C, D, V, and A
# ---------------------------------------------------------
prep_data = function(vimc_dt, country, d_v_a_id) {
  
  # Reduce to data of interest
  data_dt = vimc_dt %>%
    filter(country  == !!country, 
           d_v_a_id == !!d_v_a_id) %>%
    select(x = fvps_rel,
           y = impact_rel) %>%
    # Multiply impact for more consistent x-y scales...
    mutate(y = y * o$impact_scaler)
  
  return(data_dt)
}

# ---------------------------------------------------------
# Determine credible starting points for MLE - it needs it
# ---------------------------------------------------------
prep_start = function(fns, x, y) {
  
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
      args = NULL,
      lb   = rep(0, n_pars), 
      ub   = rep(Inf, n_pars),
      max_iters  = 10000, 
      plot_iters = NULL, 
      verbose    = FALSE)
    
    # Construct plotting datatable for this model
    plot_args = c(list(x_eval), as.list(optim$x))
    plot_list[[fn]] = data.table(
      x = x_eval,
      y = do.call(fns[[fn]], as.list(plot_args)), 
      fn = fn)
    
    # Overwrite starting point with optimal parameters
    start[[fn]] = c(s0, optim$x) %>%
      setNames(names(start[[fn]])) %>%
      as.list()
  }
  
  # Squash all models into one datatable
  plot_dt = rbindlist(plot_list)
  
  # Plot the quality of fit for each model starting point
  g_asd = ggplot(plot_dt, aes(x = x, y = y)) +
    geom_line(aes(colour = fn)) + 
    geom_point(data = data.table(x = x, y = y),
               colour = "black")
  
  return(start)
}

# ---------------------------------------------------------
# Fit MLE for each fn using prevriously determined start point
# ---------------------------------------------------------
run_mle = function(fns, start, x, y) {
  
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
      
      # Construct plotting datatable for this model
      plot_args = c(list(x_eval), as.list(coef))
      plot_list[[fn]] = data.table(
        x = x_eval,
        y = do.call(fns[[fn]], as.list(plot_args)), 
        fn = fn)
    }
  }
  
  # Squash all models into one datatable
  plot_dt = rbindlist(plot_list)
  
  # Plot the quality of fit for each model
  g_mle = ggplot(plot_dt, aes(x = x, y = y)) +
    geom_line(aes(colour = fn)) + 
    geom_point(data = data.table(x = x, y = y),
               colour = "black")
  
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

# ---------------------------------------------------------
# Evaluate best fit model
# ---------------------------------------------------------
evaluate_best_model = function(country = NULL, d_v_a_id = NULL, x = NULL) {
  
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
    x = data.table(expand_grid(country, d_v_a_id, fvps_rel = x_pts))
  }
  
  # ---- Evaluate best fit model and coefficients ----
  
  # NOTE: x = FVPs, y = impact (ie deaths averted)
  
  # Function to valuate best coefficients
  eval_fn = function(i) {
    
    # Datatable row of interest
    a = best_coef[i, ]
    
    # Extract x from input datatable
    x = x[country == a$country & d_v_a_id == a$d_v_a_id, fvps_rel]
    
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
           fvps_rel   = x, 
           impact_rel = y)
  
  return(best_fit)
}

# ---------------------------------------------------------
# Plot occurancces of each 'best' model
# ---------------------------------------------------------
plot_model_counts = function(focus = "log3") {
  
  # Load stuff: best fit functions and associtaed coefficients
  best_dt = read_rds("impact", "best_model")
  
  # ---- Plot function count ----
  
  # Simple plotting function with a few features
  plot_count = function(var, fig = "count", ord = "n") {
    
    # Determine order - with 'focus' function first
    fn_ord  = c(focus, setdiff(unique(best_dt$fn), focus))
    fn_dict = fn_set(dict = TRUE)
    
    # Number of times each model is optimal
    count_dt = best_dt %>% 
      rename(var = !!var) %>% 
      # Number and proportion of each fn...
      count(var, fn) %>%
      group_by(var) %>%
      mutate(total = sum(n)) %>%
      ungroup() %>%
      mutate(p = n / total) %>%
      # Set appropriate plotting order...
      rename(val = !!ord) %>%
      select(var, fn, val) %>%
      pivot_wider(names_from  = fn, 
                  values_from = val, 
                  values_fill = 0) %>%
      arrange_at(fn_ord) %>%
      # Final formatting...
      pivot_longer(cols = -var, 
                   names_to  = "fn", 
                   values_to = "val") %>%
      mutate(fn  = recode(fn, !!!fn_dict), 
             fn  = fct_inorder(fn),
             var = fct_inorder(var)) %>%
      as.data.table()
    
    # Check figure type flag
    if (fig == "count") {
      
      # Number of occurances
      g = (ggplot(count_dt[val > 0]) + 
             aes(x = var, y = val, fill = fn) + 
             geom_col() + 
             coord_flip()) %>%
        prettify2(save = c("Count", focus, var, ord))
    }
    
    # Check figure type flag
    if (fig == "density") {
      
      # Density of occurances
      g = ggplot(count_dt[fn == fn_dict[focus]]) + 
        aes(x = val) +
        geom_bar()
    }
    
    return(g)
  }
  
  # ---- A variety of plots ----
  
  # Plot by disease-vaccine-activity
  g1 = plot_count("d_v_a_id", ord = "n")
  g2 = plot_count("d_v_a_id", ord = "p")
  
  # Plot by country
  g3 = plot_count("country", fig = "count")
  g4 = plot_count("country", fig = "density")
  
  # Save the last figure
  save_fig(g4, "Density", focus)
}

# ---------------------------------------------------------
# Plot occurancces of each 'best' model
# ---------------------------------------------------------
plot_model_fits = function(focus, zoom = TRUE) {
  
  # Evaluate best fit model
  best_fit = evaluate_best_model()
  
  # Also load the data - we'll plot fits against this
  vimc_dt = read_rds("impact", "vimc_dt")
  
  # ---- Plot fitted FVPs vs impact ----
  
  # Plot of all succesful fits - not just the best
  g1 = (ggplot(best_fit[fn == focus]) + 
          aes(x = fvps_rel, 
              y = impact_rel, 
              colour = country, 
              group  = country) + 
          geom_line(data   = best_fit[fn != focus],
                    colour = "grey40", 
                    alpha  = 0.2, 
                    show.legend = FALSE) + 
          geom_line(size = 1, show.legend = FALSE) + 
          facet_wrap(~d_v_a_id, scales = "free_y")) %>% 
    prettify1(save = c("Best model", focus))
  
  # ---- Plot fitted against data ----
  
  # ALl combinations of C-D-V-A for this focus
  c_d_v_a = best_fit %>%
    filter(fn == focus) %>%
    select(country, d_v_a_id) %>%
    unique()
  
  # Data associated with all focus model fits
  focus_data = vimc_dt %>%
    inner_join(y  = c_d_v_a,
               by = c("country", "d_v_a_id")) %>%
    select(country, d_v_a_id, fvps_rel, impact_rel)
  
  # Determine max data
  max_data = focus_data %>%
    group_by(d_v_a_id) %>%
    slice_max(fvps_rel, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(d_v_a_id, x_max = fvps_rel) %>%
    mutate(x_max = ceiling(x_max * 10) / 10) %>%
    as.data.table()
  
  # Trivialise panel zoom if flag turned off
  if (zoom == FALSE)
    max_data$x_max = max(max_data$x_max, o$eval_x_scale)
  
  # Now shift focus soley to one function
  focus_fit = best_fit %>%
    filter(fn == focus) %>%
    select(-fn) %>%
    left_join(y  = max_data,
              by = "d_v_a_id") %>%
    filter(fvps_rel - 1e-4 <= x_max) %>%
    select(-x_max)
  
  # Plot characteristics depends on zoom
  scales = ifelse(zoom, "free", "free_y")
  name   = ifelse(zoom, "zoom", "no zoom")
  
  # Plot best fits of focus function against the data
  g2 = (ggplot(focus_fit) +
          aes(x = fvps_rel, y = impact_rel, colour = country) +
          geom_point(data = focus_data,
                     size = 0.75,
                     alpha = 0.5,
                     show.legend = FALSE) +
          geom_line(show.legend = FALSE) +
          facet_wrap(~d_v_a_id, scales = scales)) %>%
    prettify1(save = c("Fit data", focus, name))
}

# ---------------------------------------------------------
# Apply colour scheme and tidy up axes - impact plots
# ---------------------------------------------------------
prettify1 = function(g, save = NULL) {
  
  # Colour map to sample from
  map = "pals::kovesi.rainbow"
  
  # Axes label dictionary
  lab_dict = c(
    coverage    = "Coverage of target population",
    year        = "Year",
    fvps        = "Fully vaccinated persons (FVPs) in one year",
    fvps_100k   = "Fully vaccinated persons (FVPs) per 100k people",
    fvps_cum    = "Cumulative fully vaccinated persons (FVPs)",
    fvps_rel    = "Cumulative fully vaccinated persons (FVPs) per population-person",
    impact      = "Deaths averted",
    impact_100k = "Deaths averted per 100k population",
    impact_cum  = "Cumulative deaths averted", 
    impact_rel  = "Cumulative deaths averted per population-person", 
    impact_fvp  = "Deaths averted per fully vaccinated person (FVP)")
  
  # Extract info from the plot
  g_info = ggplot_build(g)
  
  # Number of colours to generates - one per country
  all_country  = table("country")$country
  plot_country = unique(g_info$plot$data$country)
  
  # Construct colours from map
  all_cols  = colour_scheme(map, n = length(all_country)) 
  plot_cols = all_cols[all_country %in% plot_country]
  
  # Apply the colours
  g = g + scale_colour_manual(values = plot_cols)
  
  # Prettyify axes
  g = g + 
    scale_x_continuous(
      name   = lab_dict[g_info$plot$labels$x], 
      expand = expansion(mult = c(0, 0.05)),  
      labels = comma) +
    scale_y_continuous(
      name   = lab_dict[g_info$plot$labels$y], 
      expand = expansion(mult = c(0, 0.05)),
      labels = comma)
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(strip.text    = element_text(size = 18), #10),
          axis.title    = element_text(size = 22), #15),
          axis.text     = element_text(size = 12, angle = 30, hjust = 1),
          # axis.text     = element_text(size = 7, angle = 30, hjust = 1),
          axis.line     = element_blank(),
          panel.border  = element_rect(linewidth = 1, colour = "black", fill = NA),
          panel.spacing = unit(0.5, "lines"),
          strip.background = element_blank()) 
  
  # Save plots to file
  if (!is.null(save))
    save_fig(g, save)
  
  return(g)
}

# ---------------------------------------------------------
# Apply colour scheme and tidy up axes - count plots
# ---------------------------------------------------------
prettify2 = function(g, save = NULL) {
  
  # Construct manual colour scheme
  cols = c("grey60", "dodgerblue1")
  
  # Apply the colours
  g = g + scale_fill_manual(name = "Best model", values = cols)
  
  # Prettyify axes
  g = g + scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                             breaks = pretty_breaks(), 
                             name   = "Count")
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title.y  = element_blank(),
          axis.title.x  = element_text(size = 18),
          axis.text     = element_text(size = 10),
          axis.line     = element_blank(),
          panel.border  = element_rect(linewidth = 1, colour = "black", fill = NA),
          panel.spacing = unit(0.5, "lines"),
          legend.title  = element_text(size = 14),
          legend.text   = element_text(size = 12),
          legend.key    = element_blank(),
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines")) 
  
  # Save plots to file
  if (!is.null(save))
    save_fig(g, save)
  
  return(g)
}

