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
  if (!is.element(5, o$do_step)) return()
  
  message("* Fitting impact functions")
  
  # ---- FVPs and impact estimates ----

  message(" - Preparing FVP-impact data")

  # Prepare impact-FVP data to fit to
  data_dt = get_impact_data()

  # Exploratory plots of data used to fit impact functions
  plot_impact_data()
  
  # Plot all-time impact per FVPs
  plot_impact_fvps(scope = "all_time")

  # ---- Model fitting ----

  message(" - Evaluating impact functions")

  # Detect cores to use (only 1 if not running in parallel)
  n_cores = ifelse(o$parallel, detectCores(), 1)

  message("  > Using ", n_cores, " cores")

  # Country-disease-vaccine-activity combinations
  run_dt = data_dt %>%
    select(country, d_v_a_id) %>%
    unique() %>%
    mutate(d_v_a_str = str_pad(d_v_a_id, 2, pad = "0"),
           run_id    = paste1(country, d_v_a_str)) %>%
    select(-d_v_a_str)

  # Initiate progress bar
  pb = progress_init(nrow(run_dt))

  # Run get_best_model in parallel
  if (o$parallel)
    mclapply(X    = run_dt$run_id,
             FUN  = get_best_model,
             run  = run_dt,
             data = data_dt,
             pb   = pb,
             mc.cores = n_cores,
             mc.preschedule = FALSE)

  # Run get_best_model consecutively
  if (!o$parallel)
    lapply(X    = run_dt$run_id,
           FUN  = get_best_model,
           run  = run_dt,
           data = data_dt,
           pb   = pb)

  # Close progress bar
  progress_close(pb)

  # ---- Model selection ----

  # Compile all results
  compile_results(run_dt)

  # Select best function for each country-d_v_a combination
  model_selection()

  # ---- Plot results ----

  # Plot function selection statistics
  plot_model_selection()
  
  # Plot impact vs coverage by vaccine, income, and decade 
  # plot_impact_coverage()
  
  # Plot impact function evaluation
  plot_model_fits()
}

# ---------------------------------------------------------
# Prepare impact-FVP data to fit to
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
  vimc_dt = read_rds("impute", "impute_result")
  
  # Load static model impact estimates
  gbd_dt = read_rds("static", "deaths_averted_vaccine") %>%
    # Scale results to per capita...
    left_join(y  = pop_dt, 
              by = c("country", "year")) %>%
    mutate(fvps   = fvps / pop, 
           impact = impact / pop) %>%
    select(-pop) %>%
    # Convert to cumulative FVP and impact...
    group_by(country, d_v_a_id) %>%
    mutate(fvps   = cumsum(fvps), 
           impact = cumsum(impact)) %>%
    ungroup() %>%
    as.data.table()
  
  # Impact estimates per capita from all sources
  data_dt = rbind(vimc_dt, gbd_dt) %>%
    mutate(impact_fvp = impact / fvps)
  
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
    linear1   = function(x, a)       y = a * x,
    logistic2 = function(x, a, b)    y = a / (1 + exp(-b * (x - a/2))), 
    logistic3 = function(x, a, b, c) y = logistic(x, a, b, upper = c))
  
  # Alternative functionality - return dictionary
  if (dict == TRUE)
    out = c(
      linear1   = "Gradiant (1 parameter)", 
      logistic2 = "Logistic (2 parameters)",
      logistic3 = "Sigmodal (3 parameters)")
  
  return(out)
}

# ---------------------------------------------------------
# Parent function to determine best fitting function
# ---------------------------------------------------------
get_best_model = function(id, run, data, pb) {
  
  # Details of this run
  this_run = run[run_id == id]
  
  # Reduce data down to what we're interested in
  fit_data_dt = data %>%
    filter(country  == this_run$country,
           d_v_a_id == this_run$d_v_a_id) %>%
    select(x = fvps,
           y = impact) %>%
    # Multiply impact for more consistent x-y scales...
    mutate(y = y * o$impact_scaler)
  
  # Append the origin (zero vaccine, zero impact)
  fit_data_dt = data.table(x = 1e-12, y = 0) %>%
    rbind(fit_data_dt)
  
  # Do not fit if insufficient data
  if (nrow(fit_data_dt) > 1) {
    
    # Declare x and y values for which we want to determine a relationship
    x = fit_data_dt$x
    y = fit_data_dt$y
    
    # Functions we'll attempt to fit with
    fns = fn_set()
    
    # Use optim algorithm to get good starting point for MLE
    start = credible_start(fns, x, y)
    
    # Run MLE from this starting point
    fit = run_mle(fns, start, x, y)
    
    # # Determine AICc value for model suitability
    model_quality(fns, fit, x, y, id)
  }
  
  # Update progress bar
  progress_idx = which(run$run_id == id)
  progress_update(pb, progress_idx)
}

# ---------------------------------------------------------
# Determine credible starting points for MLE - it needs it
# ---------------------------------------------------------
credible_start = function(fns, x, y, plot = FALSE) {
  
  # Initialise starting point for sigma in likelihood function
  s0 = 1  # This is essentially a placeholder until run_mle 
  
  # Let's start with any old points
  start = list(
    linear1   = list(s = s0, a = 1),
    logistic2 = list(s = s0, a = 1, b = 1), 
    logistic3 = list(s = s0, a = 1, b = 1, c = 1))
  
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
    linear1   = function(p, args) obj_fn("linear1",   a = p[1]),
    logistic2 = function(p, args) obj_fn("logistic2", a = p[1], b = p[2]),
    logistic3 = function(p, args) obj_fn("logistic3", a = p[1], b = p[2], c = p[3]))
  
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
      lb   = rep(1e-10, n_pars),
      ub   = rep(Inf,   n_pars),
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
# Fit MLE for each fn using previously determined start point
# ---------------------------------------------------------
run_mle = function(fns, start, x, y, plot = FALSE) {
  
  # Log likelihood function to maximise
  likelihood = function(s, y_pred)
    ll = -sum(dnorm(x = y, mean = y_pred, sd = s, log = TRUE))
  
  # Annoyingly we need to name likelihood inputs, hence wrapper functions
  model = list(
    linear1   = function(s, a)       likelihood(s, fns$linear1(x, a)),
    logistic2 = function(s, a, b)    likelihood(s, fns$logistic2(x, a, b)), 
    logistic3 = function(s, a, b, c) likelihood(s, fns$logistic3(x, a, b, c)))
  
  # Initiate list for storing plotting datatables
  fit = plot_list = list()
  
  # Points to evaluate when plotting
  x_eval = seq(0, max(x), length.out = 101)
  
  # Iterate through stats models
  for (fn in names(fns)) {
    
    # Vector of lower bounds
    n_pars = length(start[[fn]]) - 1  # Zero for function coefficients
    l_bnds = c(1e-1, rep(0, n_pars))  # Small but non-zero for sigma
    
    # Catch errors when running in parallel
    if (!o$parallel) {
      
      # Attempt to fit MLE model using prevriously determined start point
      fit_result = tryCatch(
        expr = mle(
          minuslogl = model[[fn]],
          start     = start[[fn]],
          lower     = l_bnds,
          nobs      = length(y)),
        
        # Return trivial if MLE failed
        error   = function(e) return(NULL),
        warning = function(w) return(NULL))
    }
    
    # Do not catch errors when running in parallel
    if (o$parallel) {
      
      # Attempt to fit MLE model using prevriously determined start point
      fit_result = mle(
        minuslogl = model[[fn]], 
        start     = start[[fn]], 
        lower     = l_bnds,
        nobs      = length(y))
    }
    
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
model_quality = function(fns, fit, x, y, run_id) {
  
  # Return out if no fits succesful 
  if (length(fit) == 0)
    return()
  
  # Coefficients for each successful model
  coef = unlist(lapply(fit, function(a) a@coef[-1]))
  coef = data.table(var   = names(coef), 
                    value = coef) %>%
    separate(var, c("fn", "coef")) %>%
    mutate(run_id = run_id, 
           .before = 1)
  
  # Save to file
  save_name = paste1(run_id, "coef")
  save_rds(coef, "impact_runs", save_name) 
  
  # Calculate AIC - adjusted for sample size
  aicc = sapply(fit, AICc) %>%
    as.list() %>%
    as.data.table() %>%
    mutate(run_id = run_id) %>%
    pivot_longer(cols = -run_id, 
                 names_to  = "fn", 
                 values_to = "aicc") %>%
    as.data.table()
  
  # Save to file
  save_name = paste1(run_id, "aicc")
  save_rds(aicc, "impact_runs", save_name) 
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
# Compile all results into full datatables
# ---------------------------------------------------------
compile_results = function(run_dt) {
  
  message(" - Compiling all results")
  
  # Repeat for each type of output
  for (type in c("coef", "aicc")) {
    
    # All file names for this type of result
    names = paste1(run_dt$run_id, type)
    files = paste0(o$pth$impact_runs, names, ".rds")
    
    # Files that exist 
    #
    # NOTE: Files with insufficient data are not saved
    result_files = files[file.exists(files)]
    
    # Load all results and squash into single datatable
    results_dt = result_files %>%
      lapply(readRDS) %>%
      rbindlist(fill = TRUE)
    
    # Append run details
    compiled_dt = run_dt %>%
      left_join(y  = results_dt, 
                by = "run_id") %>%
      select(-run_id)
    
    # Save to file
    save_rds(compiled_dt, "impact", type)
  }
}

# ---------------------------------------------------------
# Select best function considering complexity
# ---------------------------------------------------------
model_selection = function() {
  
  message(" - Selecting best functions")
  
  # Extract best fitting function based on AICc
  best_dt = read_rds("impact", "aicc") %>%
    # Select function with lowest AICc...
    group_by(country, d_v_a_id) %>%
    slice_min(aicc, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    # Tidy up...
    select(country, d_v_a_id, fn) %>%
    as.data.table()
  
  # Save to file
  save_rds(best_dt, "impact", "best_model")
}

