###########################################################
# PLOTTING
#
# All plotting functionality in one place.
#
###########################################################

# ---------------------------------------------------------
# Plot total number of FVP over time for each d_v_a
# ---------------------------------------------------------
plot_total_fvps = function() {
  
  message("  > Plotting total number of FVP")
  
  # Flag for whether to plot FVPs cumulatively over time
  cumulative = TRUE
  
  # Number of FVPs by source of data
  source_dt = table("coverage") %>%
    # Summarise over countries and age...
    group_by(v_a_id, source, year) %>%
    summarise(fvps = sum(fvps)) %>%
    ungroup() %>%
    # Cumulative FVPs...
    group_by(v_a_id, source) %>%
    mutate(fvps_cum = cumsum(fvps)) %>%
    ungroup() %>%
    # Report for each d_v_a...
    left_join(y  = table("v_a"), 
              by = "v_a_id") %>%
    full_join(y  = table("d_v_a"), 
              by = c("vaccine", "activity"), 
              relationship = "many-to-many") %>%
    # Tidy up...
    select(d_v_a_id, source, year, fvps, fvps_cum) %>%
    arrange(d_v_a_id, source, year) %>%
    append_d_v_a_name() %>%
    as.data.table()
  
  # Total FVPs (sum of all sources)
  total_dt = source_dt %>%
    # Summarise over source...
    group_by(d_v_a_id, year) %>%
    summarise(fvps = sum(fvps)) %>%
    ungroup() %>%
    # Cumulative FVPs...
    group_by(d_v_a_id) %>%
    mutate(fvps_cum = cumsum(fvps)) %>%
    ungroup() %>%
    # Tidy up...
    append_d_v_a_name() %>%
    as.data.table()
  
  # Metric to use for y axis
  y = ifelse(cumulative, "fvps_cum", "fvps")
  
  # Plot FVPs over time for each d_v_a
  g = ggplot(source_dt) + 
    aes(x = year, y = !!sym(y)) + 
    geom_line(aes(colour = source)) + 
    geom_line(data = total_dt, 
              linetype = "dashed",
              colour   = "black") + 
    facet_wrap(~d_v_a_name)
  
  # Save to file
  save_fig(g, "FVPs over time", dir = "data_visualisation")
}

# ---------------------------------------------------------
# Plot age targets as defined by WIISE and VIMC coverage data
# ---------------------------------------------------------
plot_coverage_age_density = function() {
  
  message("  > Plotting coverage data density by age")
  
  # Construct plotting datatable
  plot_dt = table("coverage") %>%
    mutate(trans_age = pmax(age, 1), .after = age) %>%
    append_v_a_name()
  
  # Plot age density of coverage data by source
  g = ggplot(plot_dt) +
    aes(x = trans_age, 
        y = after_stat(scaled), 
        colour = source,
        fill   = source) +
    geom_density(alpha = 0.2) +
    facet_wrap(~v_a_name)
  
  # Apply meaningful scale
  g = g + scale_x_continuous(
    trans  = "log2", 
    limits = c(1, 2^7))
  
  # Save to file
  save_fig(g, "Coverage density by age", dir = "data_visualisation")
}

# ---------------------------------------------------------
# Plot impact-FVP relationships prior to imputation
# ---------------------------------------------------------
plot_target = function() {
  
  message("  > Plotting impact-FVP relationships")
  
  # Stat to plot for both types of figures
  #
  # OPTIONS: "abs" or "cum"
  stat_1 = "cum"  # Timing figure
  stat_2 = "cum"  # Ratio figure
  
  # Number of sets to divide countries into
  n_sets = 8
  
  # ---- Plot set up ----
  
  # Divide countries into 'sets'
  set_dt = table("country") %>%
    select(country) %>%
    mutate(set = rep(1 : n_sets, nrow(.))[1 : nrow(.)], 
           set = factor(set, levels = 1 : n_sets))
  
  # Targets for VIMC diseases split by d_v_a
  plot_list = read_rds("impute", "target") %>%
    filter(!is.na(target)) %>%
    inner_join(y  = set_dt, 
               by = "country") %>%
    append_d_v_a_name() %>%
    split(f = .$d_v_a_name)
  
  # Save figures for each d_v_a
  for (d_v_a in names(plot_list)) {
    
    # ---- Plot 1: time lag between coverage and impact ----
    
    # Construct plotting datatable
    plot_dt = plot_list %>%
      pluck(d_v_a) %>%
      # Select metrics of choice...
      select(country, set, year, 
             fvps   = !!paste1("fvps", stat_1),
             impact = !!paste1("impact", stat_1)) %>%
      # Remove countries that have only one data point...
      add_count(country) %>%
      filter(n > 1) %>%
      select(-n) %>%
      # Normalise for better visuals...
      pivot_longer(cols = c(fvps, impact),
                   names_to = "variable") %>%
      group_by(country, variable) %>% 
      mutate(norm = value / max(value)) %>%
      ungroup() %>%
      filter(!is.na(norm)) %>%
      as.data.table()
    
    # Plot (any) time lag between deployment and impact
    g1 = ggplot(plot_dt) +
      aes(x = year, y = norm, colour = country) +
      geom_line(show.legend = FALSE) +
      facet_grid(variable ~ set, scales = "free_y")
    
    # ---- Plot 2: time-varying difference between FVP and impact ----
    
    # Construct plotting datatable
    plot_dt = plot_list %>%
      pluck(d_v_a) %>%
      # Select metrics of choice...
      select(country, set, year, 
             fvps   = !!paste1("fvps", stat_2),
             impact = !!paste1("impact", stat_2)) %>%
      # Calculate FVP - impact rate...
      filter(impact > 0) %>%
      mutate(impact_fvp = fvps / impact) %>%
      # Remove countries that have only one data point...
      add_count(country) %>%
      filter(n > 1) %>%
      select(-n)
    
    # Plot (any) time-varying difference between FVP and impact
    g2 = ggplot(plot_dt) +
      aes(x = year, y = impact_fvp, colour = country) +
      geom_line(show.legend = FALSE) +
      facet_wrap(~set, scales = "free_y", nrow = 2)
    
    # Save in nested directory
    dir = c("imputation", "target")
    
    # Save figure to file
    save_fig(g1, "VIMC impact-FVP timing", d_v_a, dir = dir)
    save_fig(g2, "VIMC impact-FVP ratio",  d_v_a, dir = dir)
  }
}

# ---------------------------------------------------------
# Plot (any) correlation between covariates and imputation target
# ---------------------------------------------------------
plot_covariates = function() {
  
  message("  > Plotting covariate-target relationships")
  
  # ---- Load data used for fitting ----
  
  # Function to load imputation data
  load_data_fn = function(id) {
    
    # Load imputation result for this d_v_a
    impute = read_rds("impute", "impute", id)
    
    # If result is trivial, return trivial data
    if (is.null(impute$result)) {
      data = NULL 
      
    } else {  # Otherwise extract data
      data = impute$data %>%
        mutate(d_v_a_id = id, .before = 1)
    }
    
    return(data)
  }
  
  # Load imputation data for all d-v-a
  data_dt = table("d_v_a") %>%
    pluck("d_v_a_id") %>%
    lapply(load_data_fn) %>%
    rbindlist()
  
  # ---- Produce plot ----
  
  # Construct tidy plotting datatable
  plot_dt = data_dt %>%
    pivot_longer(cols = -c(d_v_a_id, target), 
                 names_to = "covariate") %>%
    arrange(d_v_a_id, covariate, target) %>%
    append_d_v_a_name() %>%
    as.data.table()
  
  # Plot covariates vs imputation target
  g = ggplot(plot_dt) +
    aes(x = target, y = value, colour = covariate) +
    geom_point(alpha = 0.2, shape = 16, show.legend = FALSE) +
    facet_grid(covariate~d_v_a_name, scales = "free")
  
  # Save figure to file
  save_fig(g, "Covariate relationships", dir = "imputation")
}

# ---------------------------------------------------------
# Plot truth vs predicted for imputation training data
# ---------------------------------------------------------
plot_impute_fit = function() {
  
  message("  > Plotting imputation quality of fit")
  
  # ---- Load results from fitting ----
  
  # Function to load imputation results
  load_results_fn = function(id)
    result = read_rds("impute", "impute", id)$result
  
  # Load imputation results for all d-v-a
  results_dt = table("d_v_a") %>%
    pluck("d_v_a_id") %>%
    lapply(load_results_fn) %>%
    rbindlist()
  
  # ---- Construct plotting datatables ----
  
  # Prepare datatable for plotting
  plot_dt = results_dt %>%
    filter(!is.na(target)) %>%
    select(-country) %>%
    append_d_v_a_name() %>%
    # Remove target outliers for better normalisation...
    group_by(d_v_a_name) %>%
    mutate(lower = mean(target) - 3 * sd(target), 
           upper = mean(target) + 3 * sd(target), 
           outlier = target < lower | target > upper) %>%
    ungroup() %>%
    filter(outlier == FALSE) %>%
    select(-outlier, -lower, -upper) %>%
    as.data.table()
  
  # Maximum value in each facet (target or predict)
  blank_dt = plot_dt %>%
    mutate(max_value = pmax(target, predict)) %>%
    group_by(d_v_a_name) %>%
    summarise(max_value = max(max_value)) %>%
    ungroup() %>%
    expand_grid(type = c("target", "predict")) %>%
    pivot_wider(names_from  = type, 
                values_from = max_value) %>%
    as.data.table()
  
  # ---- Produce plot ----
  
  # Single plot with multiple facets
  g = ggplot(plot_dt) +
    aes(x = target, y = predict, color = d_v_a_name) +
    geom_point(alpha = 0.5, shape = 16, show.legend = FALSE) +
    geom_blank(data = blank_dt) +    # For square axes
    geom_abline(colour = "black") +  # To see quality of predict vs target
    facet_wrap(~d_v_a_name, scales = "free")
  
  # Save figure to file
  save_fig(g, "Imputation fit", dir = "imputation")
}

# ---------------------------------------------------------
# Plot country-aggregated imputation errors
# ---------------------------------------------------------
plot_impute_countries = function() {
  
  message("  > Plotting country-aggregated imputation errors")
  
  # Plot on log10 scale (second figure only)
  log_scale = TRUE
  
  # ---- Load results from fitting ----
  
  # Function to load imputation results
  load_results_fn = function(id)
    result = read_rds("impute", "impute", id)$result
  
  # Load imputation results for all d-v-a
  results_dt = table("d_v_a") %>%
    pluck("d_v_a_id") %>%
    lapply(load_results_fn) %>%
    rbindlist() %>%
    append_d_v_a_name()
  
  # ---- Plot 1: annual error by country ----
  
  # Truth vs predicted over time for training data
  annual_dt = results_dt %>%
    select(country, d_v_a_name, year, 
           vimc   = impact_cum, 
           impute = impact_impute) %>%
    filter(!is.na(vimc)) %>%
    mutate(lower = pmin(vimc, impute), 
           upper = pmax(vimc, impute))
  
  # Plot annual errors by country
  g = ggplot(annual_dt, aes(x = year)) +
    geom_ribbon(aes(ymin = lower, 
                    ymax = upper, 
                    fill = country),
                alpha = 0.3) +
    geom_line(aes(y = vimc, colour = country), 
              linewidth = 0.5) +
    geom_line(aes(y = impute, colour = country), 
              linewidth = 0.5, 
              linetype  = "dashed") +
    facet_wrap(~d_v_a_name, scales = "free_y")
  
  # Remove legend
  g = g + theme(legend.position = "none")
  
  # Save figure to file
  save_fig(g, "Imputation error annual", dir = "imputation")
  
  # ---- Plot 2: total error by country ----
  
  # Where imputed countries lie in terms of magnitude
  total_dt = results_dt %>%
    # Take cumulative values for each country...
    group_by(country, d_v_a_name) %>%
    summarise(truth   = max(impact_cum), 
              predict = max(impact_impute)) %>%
    ungroup() %>%
    # VIMC as truth-predict scatter, imputed along diagonal...
    mutate(source = ifelse(is.na(truth), "impute", "vimc"), 
           truth  = ifelse(is.na(truth), predict, truth)) %>%
    # Set lower bound to stay finite if transforming scale...
    mutate(truth   = pmax(truth,   1), 
           predict = pmax(predict, 1)) %>%
    arrange(d_v_a_name, desc(source), country) %>%
    as.data.table()
  
  # Maximum value in each facet (target or predict)
  blank_dt = total_dt %>%
    mutate(max_value = pmax(truth, predict)) %>%
    group_by(d_v_a_name) %>%
    summarise(max_value = max(max_value)) %>%
    ungroup() %>%
    expand_grid(type = c("truth", "predict")) %>%
    pivot_wider(names_from  = type, 
                values_from = max_value) %>%
    as.data.table()
  
  # Plot truth vs predicted along with imputed countries
  g = ggplot(total_dt, aes(x = truth, y = predict)) +
    geom_abline(colour = "black") +  # To see quality of truth vs predict
    geom_blank(data = blank_dt) +    # For square axes
    geom_point(aes(colour = source)) + 
    facet_wrap(~d_v_a_name, scales = "free")
  
  # Check scaling flag
  if (log_scale == TRUE) {
    
    # Transform both axes
    g = g + 
      scale_x_continuous(trans = "log10") + 
      scale_y_continuous(trans = "log10")
  }
  
  # Construct file name to save figure
  save_name = "Imputation error total"
  if (log_scale == TRUE)
    save_name = c(save_name, "log scale")
  
  # Save figure to file
  save_fig(g, save_name, dir = "imputation")
}

# ---------------------------------------------------------
# Exploratory plots of data used to fit impact functions
# ---------------------------------------------------------
plot_impact_data = function() {
  
  message("  > Plotting impact function fitting data")
  
  # Load data used for impact function fitting
  data_dt = read_rds("impact", "data") %>%
    append_d_v_a_name()
  
  # Impact per FVP over time
  g1 = ggplot(data_dt) +
    aes(x = year, y = impact_fvp, colour = country) +
    geom_line(show.legend = FALSE) +
    facet_wrap(~d_v_a_name, scales = "free_y")
  # prettify1(save = c("Year", "impact", "FVP"))
  
  # Cumulative FVPs vs cumulative deaths averted
  g2 = ggplot(data_dt) + 
    aes(x = fvps, y = impact, colour = country) +
    geom_line(show.legend = FALSE) +
    facet_wrap(~d_v_a_name, scales = "free")
  # prettify1(save = c("FVP", "impact"))
  
  # Figure sub-directory to save to
  dir = "impact_functions"
  
  # Save figures to file
  save_fig(g1, "Data - impact ratio", dir = dir)
  save_fig(g2, "Data - cumulative FVP vs impact", dir = dir)
}

# ---------------------------------------------------------
# Plot function selection statistics
# ---------------------------------------------------------
plot_model_selection = function() {
  
  message("  > Plotting impact function counts")
  
  # Load stuff: best fit functions and associtaed coefficients
  best_dt = read_rds("impact", "best_model") %>%
    append_d_v_a_name()
  
  # ---- Plot function count ----
  
  # Simple plotting function with a few features
  plot_selection = function(var, type = "count", stat = "n") {
    
    # Determine order - with 'focus' function first
    fn_dict = fn_set(dict = TRUE)
    
    # Number of times each model is optimal
    selection_dt = best_dt %>% 
      rename(var = !!var) %>% 
      # Number and proportion of each fn...
      count(var, fn) %>%
      group_by(var) %>%
      mutate(total = sum(n)) %>%
      ungroup() %>%
      mutate(p = n / total) %>%
      # Set appropriate plotting order...
      rename(val = !!stat) %>%
      select(var, fn, val) %>%
      pivot_wider(names_from  = fn, 
                  values_from = val, 
                  values_fill = 0) %>%
      arrange_at(rev(names(fn_dict))) %>%
      # Final formatting...
      pivot_longer(cols = -var, 
                   names_to  = "fn", 
                   values_to = "val") %>%
      mutate(fn  = recode(fn, !!!fn_dict), 
             fn  = fct_inorder(fn),
             var = fct_inorder(var)) %>%
      as.data.table()
    
    # Check figure type flag
    if (type == "count") {
      
      # Number of occurances
      g = ggplot(selection_dt[val > 0]) + 
        aes(x = var, y = val, fill = fn) + 
        geom_col() + 
        coord_flip()
      # prettify2(save = c("Count", focus, var, stat))
    }
    
    # Check figure type flag
    if (type == "density") {
      
      # Density of occurances
      g = ggplot(selection_dt) + 
        aes(x = val, fill = fn) +
        geom_bar()
    }
    
    return(g)
  }
  
  # ---- A variety of plots ----
  
  # Figure sub-directory to save to 
  save_dir = "impact_functions"
  
  # Plot by disease-vaccine-activity
  g1 = plot_selection("d_v_a_name", stat = "n")
  g2 = plot_selection("d_v_a_name", stat = "p")
  
  save_fig(g1, "Selection", "pathogen", "number",     dir = save_dir)
  save_fig(g2, "Selection", "pathogen", "proportion", dir = save_dir)
  
  # Plot by country
  g3 = plot_selection("country", type = "count")
  g4 = plot_selection("country", type = "density")
  
  # Save the last figure
  save_fig(g3, "Selection", "country", dir = save_dir)
  save_fig(g4, "Density",   "country", dir = save_dir)
}

# ---------------------------------------------------------
# Plot impact function evaluation
# ---------------------------------------------------------
plot_model_fits = function() {
  
  message("  > Plotting impact function fits")
  
  # Load data used for impact function fitting
  data_dt = read_rds("impact", "data") %>%
    append_d_v_a_name()
  
  # Evaluate only as far as we have data
  max_data = data_dt %>%
    group_by(country, d_v_a_name) %>%
    slice_max(fvps, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(country, d_v_a_name, x_max = fvps) %>%
    # Increment up one so we plot slightly past the data
    mutate(x_max = x_max + o$eval_x_scale / 100) %>%
    as.data.table()
  
  # Evaluate selected impact function
  best_fit = evaluate_impact_function() %>%
    append_d_v_a_name()
  
  # Apply max_data so we only plot up to data (or just past)
  plot_dt = best_fit %>%
    left_join(y  = max_data,
              by = c("country", "d_v_a_name")) %>%
    filter(fvps < x_max) %>%
    append_d_v_a_name() %>%
    select(country, d_v_a_name, fvps, impact)
  
  # Plot function evaluation against the data
  g = ggplot(plot_dt) +
    aes(x = fvps, y = impact, colour = country) +
    geom_point(data = data_dt,
               size = 0.75,
               alpha = 0.5,
               show.legend = FALSE) +
    geom_line(show.legend = FALSE) +
    facet_wrap(~d_v_a_name, scales = "free")
  # prettify1(save = c("Fit data", focus, name))
  
  save_fig(g, "Impact function evaluation", dir = "impact_functions")
}

# ---------------------------------------------------------
# Plot annual totals - diagnostic figure to check alignment of means
# ---------------------------------------------------------
plot_annual_total = function() {
  
  message("  > Plotting annual totals")
  
  browser() # Needs updating for EPI50 pipeline...
  
  # Load modelled total deaths per year
  scenario_total = try_load(o$pth$impact_factors, "scenario_total")
  
  # Load uncertainty draws - we want to see the means of these the same as above
  draws_dt = try_load(o$pth$uncertainty, "draws")
  
  # Uncertainty per year (all diseases and countries) across draws
  annual_dt = draws_dt %>%
    # Melt to long format...
    pivot_longer(cols = starts_with("draw"), 
                 names_to = "draw") %>%
    # Total deaths averted per year (per draw)...
    group_by(year, draw) %>%
    summarise(value = sum(value)) %>%
    ungroup() %>%
    # Mean and bounds across draws...
    group_by(year) %>%
    summarise(mean =  mean(value), 
              lower = quantile(value, 0.05), 
              upper = quantile(value, 0.95)) %>%
    ungroup() %>%
    as.data.table()
  
  # Plot mean and bounds from draws, and overlay modelled mean
  g = ggplot(annual_dt, aes(x = year)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                colour = "red", fill = "red", alpha = 0.5) +
    geom_line(aes(y = mean), colour = "red", linewidth = 2) +
    geom_point(data    = scenario_total,
               mapping = aes(y = deaths_averted),
               colour  = "black")
  
  # Prettify plot
  g %<>% ggpretty(
    x_lab = "Year", 
    y_lab = "Deaths averted")
  
  # Save figure to file
  save_name = "Uncertainty bounds - Annual total"
  save_fig(g, save_name, dir = "uncertainty")
}

# ---------------------------------------------------------
# Plot uncertainty draws
# ---------------------------------------------------------
plot_uncertainty_draws = function() {
  
  message("  > Plotting uncertainty draws")
  
  # Flag for transforming y axis to log10 scale
  y_transform = FALSE
  
  # ---- Load and format plot datatables ----
  
  # Load draws from file
  draws_dt = try_load(o$pth$uncertainty, "draws")
  
  # Append source and collapse to d_v_at
  details_dt = draws_dt %>%
    left_join(disease_table, by = "disease") %>%
    unite("d_v_at", disease, vaccine, activity_type)
  
  # Mean deaths averted (not from draws)
  mean_dt = details_dt %>%
    group_by(d_v_at, source) %>%
    summarise(deaths_averted = sum(deaths_averted)) %>%
    ungroup() %>%
    as.data.table()
  
  # Sampled deaths averted
  samples_dt = details_dt %>%
    pivot_longer(cols = starts_with("draw"), 
                 names_to = "draw") %>%
    group_by(d_v_at, source, draw) %>%
    summarise(deaths_averted = sum(value)) %>%
    ungroup() %>%
    as.data.table()
  
  # ---- Produce plot ----
  
  # Plot draws and then mean on top
  g = ggplot(samples_dt, aes(x = d_v_at, y = deaths_averted)) + 
    geom_violin(aes(colour = source, fill = source), 
                alpha = 0.2) + 
    geom_point(data = mean_dt, colour = "black", size = 2)
  
  # Prettify plot
  g %<>% ggpretty(
    x_lab = "Disease - vaccine - activity", 
    y_lab = "Deaths averted",
    x_rotate = TRUE, 
    y_pretty = FALSE)
  
  # Transform y axis to log10 scale if desired
  scale_fn = ifelse(y_transform, "scale_y_log10", "scale_y_continuous")
  g = g + get(scale_fn)(labels = label_comma(), 
                        expand = expansion(mult = c(0, 0.05)))
  
  # Save figure to file
  save_name = "Uncertainty draws - All diseases"
  save_fig(g, save_name, dir = "uncertainty")
}

# ---------------------------------------------------------
# Plot parameters of fitted beta distribution to vaccine efficacy
# ---------------------------------------------------------
plot_gbd_uncertainty_dist = function() {
  
  message("  > Plotting GBD uncertainty distribution")
  
  # Points over which to evaluate beta distribution
  eval_pts = seq(0, 1, length.out = 100)
  
  # ---- Load and format plot datatables ----
  
  browser() # gbd_efficacy is now efficacy
  
  # Load actual vaccine efficacy for GBD diseases and collapse disease-vaccine
  efficacy_dt = gbd_efficacy %>%
    left_join(y  = disease_table, 
              by = "disease") %>%
    unite("d_v", disease_name, vaccine) %>%
    select(d_v, mean, lower, upper)
  
  # Load fitted parameters and collapse disease-vaccine
  beta_pars = try_load(o$pth$uncertainty, "gbd_beta_pars") %>%
    left_join(y  = disease_table, 
              by = "disease") %>%
    unite("d_v", disease_name, vaccine) %>%
    select(d_v, p1, p2)
  
  # Mean and 90% CI of fitted beta distribution
  beta_summary_dt = beta_pars %>%
    mutate(par_mean  = p1 / (p1 + p2),       # Mean of a beta distribution
           par_lower = qbeta(0.05, p1, p2),  # Lower bound
           par_upper = qbeta(0.95, p1, p2))  # Upper bound
  
  # Function to evaluate beta distribution
  beta_fn = function(p, x)
    dbeta(x = eval_pts, p[1], p[2])
  
  # Evaluate beta for fitted parameters for each disease-vaccine
  beta_dist_dt = beta_pars %>%
    select(p1, p2) %>%
    apply(1, beta_fn, x = eval_pts) %>%
    # Convert to tidy datatable...
    as_named_dt(beta_pars$d_v) %>%
    mutate(x = eval_pts) %>%
    pivot_longer(cols = - x, 
                 names_to = "d_v") %>%
    # Normalise probability distributions...
    group_by(d_v) %>%
    mutate(value = value / max(value)) %>%
    ungroup() %>%
    arrange(d_v, x) %>%
    as.data.table()
  
  # ---- Produce plot ----
  
  # Plot actual vaccine efficacy for each disease-vaccine
  g = ggplot(efficacy_dt) + 
    geom_errorbar(aes(y = 1.1, xmin = lower, xmax = upper), 
                  width = 0.1, size = 1.5, colour = "black") + 
    geom_point(aes(y = 1.1, x = mean), 
               size = 3, colour = "darkred") +
    facet_wrap(~d_v, nrow = 1)
  
  # Plot fitted beta distribution
  g = g + geom_line(data = beta_dist_dt, aes(x = x, y = value), 
                    size = 2, colour = "darkblue") 
  
  # Function for adding vlines to plot
  add_vline = function(g, var, col)
    g = g + geom_vline(data     = beta_summary_dt, 
                       mapping  = aes_string(xintercept = var), 
                       colour   = col, 
                       linetype = "dashed")
  
  # Append mean and 90% CI from fitted beta distribution
  g %<>% add_vline("par_mean", "darkred") %>%
    add_vline("par_lower", "black") %>%
    add_vline("par_upper", "black")
  
  # Prettify plot
  g %<>% ggpretty(
    x_lab = "Efficacy", 
    y_lab = "Probability density")
  
  # Final prettifying touches
  g = g + theme(axis.text.y = element_blank())
  
  # Save figure to file
  save_name = "Uncertainty distributions - GBD diseases"
  save_fig(g, save_name, dir = "uncertainty")
}

# ---------------------------------------------------------
# Plot optimisation perfomance of dist fit to vaccine efficacy
# ---------------------------------------------------------
plot_gbd_uncertainty_fit = function() {
  
  message("  > Plotting GBD uncertainty fit")
  
  # Grid points for diagnostic plot
  n_grid = 100  # n_grid^2 total evaluations per disease
  
  # Padding around heatmap
  padding = 0.05
  
  # Colours for good and bad objective values
  colours = list(
    good = "blue", 
    bad  = "white", 
    best = "darkred")  # Best fit parameter value(s)
  
  # ---- Evaluate objective function for grid of points ----
  
  # Points to evaluate across each parameter
  eval_pts = seq(o$par_lower, o$par_upper, length.out = n_grid)
  
  # Create grid of all points to evaluate
  grid_dt = expand_grid(
    p1  = eval_pts,
    p2  = eval_pts,
    obj = NA) %>%
    as.data.table()
  
  # Initate list of grids
  grid_list = list()
  
  browser() # gbd_efficacy is now efficacy
  
  # Loop through diseases
  for (i in 1 : nrow(gbd_efficacy)) {
    
    # Vaccine efficacy details
    v = gbd_efficacy[i, .(mean, lower, upper)]
    
    # To assess performance, evaluate every point in a grid
    #
    # NOTE: See uncertainty.R for function gbd_obj_fn()
    for (j in 1 : nrow(grid_dt))
      grid_dt$obj[[j]] = gbd_obj_fn(unlist(grid_dt[j, .(p1, p2)]), v)
    
    grid_list[[i]] = grid_dt %>%
      mutate(disease = gbd_efficacy[i, disease], 
             vaccine = gbd_efficacy[i, vaccine])
  }
  
  # Bind datatables and collapse disease-vaccine
  plot_dt = rbindlist(grid_list) %>%
    left_join(y  = disease_table, 
              by = "disease") %>%
    unite("d_v", disease_name, vaccine) %>%
    select(d_v, p1, p2, obj)
  
  # Also load best fit parameters (see uncertainty.R)
  fitted_pars = try_load(o$pth$uncertainty, "gbd_beta_pars") %>%
    left_join(y  = disease_table, 
              by = "disease") %>%
    unite("d_v", disease_name, vaccine) %>%
    mutate(p1 = log(p1), p2 = log(p2)) %>%
    select(d_v, p1, p2)
  
  # ---- Produce plot ----
  
  # Heat map of objective function for each disease-vaccine
  g = ggplot(plot_dt, aes(x = p1, y = p2)) +
    geom_tile(aes(fill = obj), colour = NA) +
    facet_wrap(~d_v) + 
    scale_fill_gradient(low  = colours$good, 
                        high = colours$bad)
  
  # Plot pre-determined optimal value on top
  g = g + geom_point(data   = fitted_pars, 
                     colour = colours$best, 
                     size   = 3)
  
  # Prettify plot
  g %<>% ggpretty(
    title = "Optimisation performance",
    x_lab = "Shape parameter 1", 
    y_lab = "Shape parameter 2", 
    x_pretty = FALSE, 
    y_pretty = FALSE)
  
  # Remove padding and text from heatmap
  g = g + scale_x_continuous(expand = c(padding, padding)) + 
    scale_y_continuous(expand = c(padding, padding)) + 
    theme(panel.border = element_blank(), 
          axis.text.x  = element_blank(), 
          axis.text.y  = element_blank(),
          axis.ticks   = element_blank())
  
  # Save figure to file
  save_name = "Uncertainty fit - GBD diseases"
  save_fig(g, save_name, dir = "uncertainty")
}

# ---------------------------------------------------------
# Convert d_v_a_id into human-readable sting
# ---------------------------------------------------------
append_d_v_a_name = function(id_dt) {
  
  # Append d_v_a description
  name_dt = id_dt %>%
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_id") %>%
    mutate(d_v_a_name = paste0(disease, " (", 
                               vaccine, "): ", 
                               activity), 
           .after = d_v_a_id) %>%
    select(-disease, -vaccine, -activity)
  
  return(name_dt)
}

# ---------------------------------------------------------
# Convert v_a_id into human-readable sting
# ---------------------------------------------------------
append_v_a_name = function(id_dt) {
  
  # Append v_a description
  name_dt = id_dt %>%
    left_join(y  = table("v_a"), 
              by = "v_a_id") %>%
    mutate(v_a_name = paste0(vaccine, ": ", 
                             activity), 
           .after = v_a_id) %>%
    select(-vaccine, -activity)
  
  return(name_dt)
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

# ---------------------------------------------------------
# Save a ggplot figure to file with default settings
# ---------------------------------------------------------
save_fig = function(g, ..., dir = NULL) {
  
  # Collapse inputs into vector of strings
  fig_name_parts = unlist(list(...))
  
  # Construct file name to concatenate with file path
  save_name = paste(fig_name_parts, collapse = " - ")
  
  # Construct path to save file to
  save_path = o$pth$figures
  if (!is.null(dir)) {
    ext_path  = paste(dir, collapse = file_sep())
    save_path = paste0(save_path, ext_path, file_sep())
  }
  
  # Create directory if it exists
  if (!dir.exists(save_path))
    dir.create(save_path)
  
  # Repeat the saving process for each image format in figure_format
  for (fig_format in o$figure_format) {
    full_path = paste0(save_path, save_name, ".", fig_format)
    
    # Save figure (size specified in options.R)
    ggsave(full_path, 
           plot   = g, 
           device = fig_format, 
           dpi    = o$save_resolution, 
           width  = o$save_width, 
           height = o$save_height, 
           units  = o$save_units)
  }
}

