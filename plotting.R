###########################################################
# PLOTTING
#
# All plotting functionality in one place.
#
###########################################################

#' Plot coverage for a specific country
#' @param dt A data.table with coverage data
#' @return A plot of the coverage data
#' @method plot coverage
#' @export
plot.coverage <- function(x, ...) {
  gg <- ggplot2::ggplot(
    x,
    ggplot2::aes(
      x = year,
      y = value,
      color = vaccine_short,
      linetype = as.factor(gender)
    )
  ) +
    ggplot2::geom_line()
  
  return(gg)
}

## Scatter against covariates
scatter_rr <- function(dt, x_var) {
  for (a in 0:9) {
    gg <- ggplot(dt[age == a], aes(x = get(x_var), y = rr)) +
      geom_point(size = 0.1, alpha = 0.2) +
      facet_wrap(~vaccine_short, scales = "free_y") +
      theme_bw() +
      ggtitle(paste(x_var, "vs mortality reduction by vaccine: Age", a)) +
      xlab(x_var)
    print(gg)
  }
}

#' Make a map showing presence or absence of an indicator
#' 
#' @param countries A character vector of iso3 codes for countries
#' @param title A string with the title of the plot
#' @returns A ggplot object with a world map
#' @examples 
#' map_countries(country_table$country, "All countries")
#' @export
map_locations <- function(countries, title) {
  ggplot2::theme_set(ggplot2::theme_bw())
  world <- rnaturalearth::ne_countries(
    scale = "medium",
    continent = c(
      "north america", "africa", "south america",
      "europe", "asia", "oceania"
    ),
    returnclass = "sf"
  )
  
  world$present <- ifelse(world$iso_a3 %in% countries, 1, NA)
  gg <- ggplot2::ggplot(data = world) +
    ggplot2::geom_sf(ggplot2::aes(fill = as.factor(present))) +
    ggplot2::ggtitle(
      title,
      subtitle = paste0("(", length(countries), " countries)")
    ) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_fill_discrete(na.value = "gray95")
  return(gg)
}

plot_age_year <- function(dt, log_transform = F, value_name = "") {
  dt <- unique(dt[, .(year, age, value)])
  
  gg <- ggplot(dt, aes(x = year, y = age, fill = value)) + 
    geom_tile() +
    xlab("Year") + ylab("Age") + labs(fill = value_name) +
    viridis::scale_fill_viridis(
      option = "viridis", 
      direction = -1, 
      trans = ifelse(log_transform, "log10", "identity")) + 
    theme_minimal() +
    coord_fixed() +
    theme(text = element_text(size = 20))
  print(gg)
}

# ---------------------------------------------------------
# xxxxxxx
# ---------------------------------------------------------
plot_strata_fit <- function(rr_dt) {
  
  # Prepare datatable for plotting
  plot_dt = rr_dt %>%
    rename(truth   = deaths_averted, 
           predict = averted) %>%
    filter(truth   > 0, 
           predict > 0) %>%
    left_join(d_v_at_table, by = "d_v_at_id") %>%
    mutate(strata = paste0(vaccine, ": ", activity_type))
  
  # Maximum value in each strata (truth or predict)
  blank_dt = plot_dt %>%
    mutate(max_value = ceiling(pmax(truth, predict))) %>%
    group_by(strata) %>%
    summarise(max_value = max(max_value)) %>%
    ungroup() %>%
    expand_grid(type = c("truth", "predict")) %>%
    pivot_wider(names_from  = type, 
                values_from = max_value) %>%
    as.data.table()
  
  # Single plot with multiple facets
  #
  # NOTE: Removed the log10 scaling
  g = ggplot(plot_dt, aes(x = truth, y = predict)) +
    geom_point(aes(color = age), size = 0.5, alpha = 0.5) +
    geom_blank(data = blank_dt) +  # For square axes
    geom_abline(slope = 1) +  # To see quality of predict = truth
    facet_wrap(~strata, scales = "free")
  
  # Use a nice colour scheme
  g = g + viridis::scale_color_viridis(option = "viridis")
  
  # Prettify
  g = g + theme_bw() +
    xlab("Truth") + 
    ylab("Predicted")
  
  # Save figure to file
  save_fig(g, "Strata fit", dir = "diagnostics")
}

# ---------------------------------------------------------
# Plot uncertainty draws
# ---------------------------------------------------------
plot_draws = function(fig_name) {
  
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
  save_fig(g, fig_name, dir = "diagnostics")
}

# ---------------------------------------------------------
# Plot annual totals - diagnostic figure to check alignment of means
# ---------------------------------------------------------
plot_annual_total = function(fig_name) {
  
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
  save_fig(g, fig_name, dir = "diagnostics")
}

# ---------------------------------------------------------
# Plot parameters of fitted beta distribution to vaccine efficacy
# ---------------------------------------------------------
plot_gbd_uncertainty_dist = function(fig_name) {
  
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
  save_fig(g, fig_name, dir = "diagnostics")
}

# ---------------------------------------------------------
# Plot optimisation perfomance of dist fit to vaccine efficacy
# ---------------------------------------------------------
plot_gbd_uncertainty_fit = function(fig_name) {
  
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
  save_fig(g, fig_name, dir = "diagnostics")
}

# ---------------------------------------------------------
# Prettify ggplot figure
# ---------------------------------------------------------
ggpretty = function(g, cols = NULL, colour = NULL, fill = NULL, title = NULL, x_rotate = FALSE,
                    x_lab = NULL, y_lab = NULL, x_pretty = TRUE, y_pretty = TRUE) {
  
  # Set line colours if specified
  if (!is.null(colour))
    g = g + scale_colour_manual(values = cols[[colour]], 
                                name   = first_cap(colour))
  
  # Set patch colours if specified
  if (!is.null(fill))
    g = g + scale_fill_manual(values = cols[[fill]], 
                              name   = first_cap(fill))
  
  # Prettify x axis
  if (x_pretty && !x_rotate)
    g = g + scale_x_continuous(breaks = scales::pretty_breaks())
  
  # Prettify y axis
  if (y_pretty)
    g = g + scale_y_continuous(breaks = scales::pretty_breaks(), 
                               expand = expansion(mult = c(0, 0.05)))
  
  # Set labels
  g = g + ggtitle(title) +
    xlab(x_lab) + 
    ylab(y_lab)
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(plot.title    = element_text(size = o$font_size[1], hjust = 0.5),
          axis.title    = element_text(size = o$font_size[2]),
          axis.text     = element_text(size = o$font_size[3]),
          axis.text.x   = element_text(angle = ifelse(x_rotate, 50, 0), 
                                       hjust = ifelse(x_rotate, 1, 0.5)), 
          axis.line     = element_blank(),
          panel.border  = element_rect(linewidth = 1, colour = "black", fill = NA),
          panel.spacing = unit(1, "lines"),
          strip.text    = element_text(size = o$font_size[4]),
          strip.background = element_blank(),
          legend.title  = element_text(size = o$font_size[5]),
          legend.text   = element_text(size = o$font_size[6]))
  
  return(g)
}

# ---------------------------------------------------------
# Save a ggplot figure to file with default settings
# ---------------------------------------------------------
save_fig = function(g, ..., dir = "figures") {
  
  # Collapse inputs into vector of strings
  fig_name_parts = unlist(list(...))
  
  # Construct file name to concatenate with file path
  save_name = paste(fig_name_parts, collapse = " - ")
  
  # Repeat the saving process for each image format in figure_format
  for (fig_format in o$figure_format) {
    save_pth  = paste0(o$pth[[dir]], save_name, ".", fig_format)
    
    # Save figure (size specified in options.R)
    ggsave(save_pth, 
           plot   = g, 
           device = fig_format, 
           dpi    = o$save_resolution, 
           width  = o$save_width, 
           height = o$save_height, 
           units  = o$save_units)
  }
}

