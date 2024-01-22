###########################################################
# PLOTTING
#
# All plotting functionality in one place.
#
###########################################################

# ---------------------------------------------------------
# Plot methodology figure to be used in paper
# ---------------------------------------------------------
plot_scope = function() {
  
  message("  > Plotting country-disease scope")
  
  # Manually set tidy y axis limit
  y_max = 625  # In millions
  
  # Linear time interpolation for less pixilated figures
  smoothness = 10  # Higher value for smoother plot
  
  # Dictionary for full impact source descriptions
  impact_dict = c(
    extern = "Dynamic modelling (external to VIMC)",
    vimc   = "Dynamic modelling (contributing to VIMC)", 
    static = "Static modelling", 
    impute = "Geographic imputation model", 
    extrap = "Temporal extrapolation model")
  
  # Associated colours
  impact_colours = c("#EB7D5B", "#FED23F", "#B5D33D", "#6CA2EA", "#442288")
  
  # ---- Number of FVPs by pathogen ----
  
  browser()
  
  # Number of FVPs over time
  fvps_dt = table("coverage") %>%
    filter(coverage > 0) %>%
    # Append disease details...
    left_join(y  = table("v_a"), 
              by = "v_a_id") %>%
    left_join(y  = table("d_v_a"), 
              by = c("vaccine", "activity")) %>%
    # Concept here is FVP, so remove birth dose & boosters...
    filter(!grepl("_BX$", vaccine), 
           !grepl("_BD$", vaccine)) %>% 
    # Summarise over age...
    group_by(disease, country, year) %>%
    summarise(fvps = sum(fvps)) %>%
    ungroup() %>%
    as.data.table()
  
  # TEMP: A placeholder for polio (copy of measles FVP)
  polio_dt = fvps_dt %>%
    filter(disease == "Measles") %>%
    mutate(disease = "Polio", 
           fvps = fvps * 1.02)
  
  # TEMP: Append these polio placeholder values
  fvps_dt %<>% rbind(polio_dt)
  
  # ---- Source of impact estimates ----
  
  # Static model approach
  static_dt = table("gbd_estimates") %>%
    select(disease, country, year) %>%
    unique() %>%
    mutate(class = "static")
  
  # VIMC approach
  vimc_dt = table("vimc_estimates") %>%
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_id") %>%
    select(disease, country, year) %>%
    unique() %>%
    mutate(class = "vimc")
  
  # VIMC country imputation method
  impute_dt = table("disease") %>%
    filter(source == "vimc") %>%
    select(disease) %>%
    expand_grid(
      country = all_countries(), 
      year    = unique(vimc_dt$year)) %>%
    left_join(y  = vimc_dt, 
              by = c("disease", "country", "year")) %>%
    filter(is.na(class)) %>%
    mutate(class = "impute") %>%
    as.data.table()
  
  # Other modelled pathogens
  #
  # TEMP: Read in polio impact table when ready
  extern_dt = fvps_dt %>%
    filter(disease == "Polio") %>%
    select(disease, country, year) %>%
    mutate(class = "extern")
  
  # ---- Construct plotting datatable ----
  
  # Combine all impact sources
  all_dt = rbind(static_dt, vimc_dt, impute_dt, extern_dt) %>%
    # Append FVPs...
    right_join(y  = fvps_dt, 
               by = c("disease", "country", "year")) %>%
    # Anything not yet specified is time-exrapolated...
    mutate(class = ifelse(is.na(class), "extrap", class)) %>%
    group_by(disease, class, year) %>%
    summarise(fvps = sum(fvps)) %>%
    ungroup() %>%
    as.data.table()
  
  # Year range of analysis
  year1 = min(o$years)
  year2 = max(o$years)
  
  # Smoothen over non-trivial years
  plot_dt = expand_grid(
    disease = unique(fvps_dt$disease), 
    class   = names(impact_dict), 
    year    = seq(year1, year2, by = 1 / smoothness)) %>%
    # Append results and source of impact...
    left_join(y  = all_dt, 
              by = c("disease", "class", "year")) %>%
    # Interpolate annual values for smoother plot...
    mutate(year_int = floor(year)) %>%
    group_by(disease, class, year_int) %>%
    mutate(n = sum(fvps, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(n > 0) %>%
    select(-year_int, -n) %>%
    # Interpolate annual values for smoother plot...
    group_by(disease, class) %>%
    mutate(fvps = na_interpolation(fvps)) %>%
    ungroup() %>%
    # Convert to units of millions...
    mutate(fvps = fvps / 1e6) %>%
    # filter(fvps > 0) %>%
    # Use impact source descriptions...
    mutate(class = recode(class, !!!impact_dict), 
           class = factor(class, rev(impact_dict))) %>%
    # Use full disease names...
    left_join(y  = table("disease"), 
              by = "disease") %>%
    replace_na(list(disease_name = "Poliomyelitis")) %>%
    # Total number of FVP over time...
    group_by(disease) %>%
    mutate(total = sum(fvps)) %>%
    ungroup() %>%
    # And use this to set disease order...
    arrange(desc(total), disease, class, year) %>%
    select(disease = disease_name, class, year, fvps, total) %>%
    mutate(disease = fct_inorder(disease)) %>%
    as.data.table()
  
  # ---- Construct label datatable ----
  
  # Label description string
  label_str = paste0("Total (", year1, "-", year2, "): ")
  
  # Construct labels: total FVPs over analysis timeframe
  label_dt = plot_dt %>%
    select(disease, total) %>%
    unique() %>%
    mutate(total = round(total / (1e3 * smoothness), 2), 
           label = paste0(label_str, total, " billion")) %>%
    # Set coordinates...
    mutate(year = year1 + 0.01 * (year2 - year1), 
           fvps = y_max * 0.9)
  
  # ---- Produce plot ----
  
  # Plot FVP over time per pathogen and impact source
  g = ggplot(plot_dt) +
    aes(x = year, 
        y = fvps) +
    geom_bar(
      mapping  = aes(fill = class), 
      stat     = "identity", 
      position = "stack",  
      width    = 1 / smoothness) +
    # Add total labels...
    geom_text(
      data    = label_dt, 
      mapping = aes(label = label), 
      size    = 3.5,
      hjust   = 0, 
      vjust   = 1) + 
    # Some intricate faceting...
    facet_rep_wrap(
      facets = ~disease, 
      ncol   = 1, 
      labeller = label_wrap_gen(width = 20), 
      strip.position = "right", 
      repeat.tick.labels = FALSE) + 
    # Set colours and legend title...
    scale_fill_manual(
      values = impact_colours, 
      name   = "Source of impact estimates") +
    guides(fill = guide_legend(reverse = TRUE)) +
    # Prettify x axis...
    scale_x_continuous(
      limits = c(year1 - 1 / smoothness, 
                 year2 + 1 / smoothness), 
      expand = expansion(mult = c(0, 0)), 
      breaks = seq(year1, year2, by = 5)) +  
    # Prettify y axis...
    scale_y_continuous(
      name   = "Number of people receiving final primary dose (in millions)", 
      limits = c(0, y_max), 
      labels = comma,
      expand = expansion(mult = c(0, 0)))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(strip.text    = element_text(size = 14),
          strip.text.y  = element_text(angle = 0, hjust = 0),
          strip.background = element_blank(), 
          axis.title.x  = element_blank(),
          axis.title.y  = element_text(
            size = 20, margin = margin(l = 10, r = 20)),
          axis.text     = element_text(size = 8),
          axis.ticks    = element_blank(), 
          axis.line     = element_line(linewidth = 0.25),
          panel.spacing = unit(-0.6, "lines"), 
          panel.grid.major.y = element_line(linewidth = 0.25),
          legend.title  = element_text(size = 14),
          legend.text   = element_text(size = 11),
          legend.key    = element_blank(),
          legend.position = "right", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save to file
  save_fig(g, "Analysis scope", dir = "data_visualisation")
  save_fig(g, "Figure 1", dir = "manuscript")
}

# ---------------------------------------------------------
# Plot total number of FVP over time for each d_v_a
# ---------------------------------------------------------
plot_total_fvps = function() {
  
  message("  > Plotting total number of FVP")
  
  # ---- Plot 1: by d_v_a ----
  
  # Flag for whether to plot FVPs cumulatively over time
  cumulative = TRUE
  
  # Number of FVPs by source of data
  source_dt = table("coverage_source") %>%
    # Summarise over countries and age...
    group_by(v_a_id, source, year) %>%
    summarise(fvps = sum(fvps) / 1e9) %>%
    ungroup() %>%
    # Cumulative FVPs...
    group_by(v_a_id, source) %>%
    mutate(fvps_cum = cumsum(fvps)) %>%
    ungroup() %>%
    # Report for each d_v_a...
    left_join(y  = table("v_a"), 
              by = "v_a_id") %>%
    full_join(y  = table("d_v_a"), 
              by = c("vaccine", "activity")) %>%
    select(d_v_a_name, source, year, fvps, fvps_cum) %>%
    # Tidy up...
    format_d_v_a_name() %>%
    arrange(d_v_a_name, source, year) %>%
    as.data.table()
  
  # Total FVPs (sum of all sources)
  #
  # NOTE: Not necessarily equal to sum of all sources as
  #       SIA are assumed to be only partially targeted
  total_dt = table("coverage") %>%
    # Summarise over countries and age...
    group_by(v_a_id, year) %>%
    summarise(fvps = sum(fvps) / 1e9) %>%
    ungroup() %>%
    # Cumulative FVPs...
    group_by(v_a_id) %>%
    mutate(fvps_cum = cumsum(fvps)) %>%
    ungroup() %>%
    # Report for each d_v_a...
    left_join(y  = table("v_a"), 
              by = "v_a_id") %>%
    full_join(y  = table("d_v_a"), 
              by = c("vaccine", "activity")) %>%
    select(d_v_a_id, year, fvps, fvps_cum) %>%
    # Tidy up...
    format_d_v_a_name() %>%
    arrange(d_v_a_id, year) %>%
    as.data.table()
  
  # Metric to use for y axis
  y = ifelse(cumulative, "fvps_cum", "fvps")
  
  # Plot FVPs over time for each d_v_a
  g = ggplot(source_dt) + 
    aes(x = year, y = !!sym(y)) + 
    geom_line(
      mapping   = aes(colour = source), 
      linewidth = 1.5) + 
    geom_line(
      data      = total_dt, 
      linetype  = "dashed",
      colour    = "black", 
      linewidth = 1.5) + 
    # Facet with strip text wrapping...
    facet_wrap(
      facets   = vars(d_v_a_name), 
      labeller = label_wrap_gen(width = 24), 
      scales   = "free_y") + 
    # Prettify x axis...
    scale_x_continuous(
      limits = c(min(o$years), max(o$years)), 
      expand = expansion(mult = c(0, 0)), 
      breaks = seq(min(o$years), max(o$years), by = 10)) +  
    # Prettify y axis...
    scale_y_continuous(
      name   = "Total receiving full schedule (in billions)", 
      labels = comma,
      expand = expansion(mult = c(0, NA)))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title.x  = element_blank(),
          axis.title.y  = element_text(
            size = 20, margin = margin(l = 10, r = 20)),
          axis.text     = element_text(size = 10),
          axis.text.x   = element_text(hjust = 1, angle = 50), 
          axis.line     = element_blank(),
          strip.text    = element_text(size = 12),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"),
          panel.grid.major.y = element_line(linewidth = 0.5),
          legend.title  = element_blank(),
          legend.text   = element_text(size = 14),
          legend.key    = element_blank(),
          legend.position = "right", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save to file
  save_dir = "data_visualisation"
  save_fig(g, "FVPs by source", dir = save_dir)
  
  # ---- Plot 2: by disease & vaccine ----
  
  # Produce plot for both disease and vaccine
  for (d_v in c("disease", "vaccine")) {
    
    # Sum FVPs for each disease or vaccine
    d_v_dt = source_dt %>%
      left_join(y  = table("d_v_a"), 
                by = "d_v_a_name") %>%
      select(d_v = !!d_v, activity, year, fvps) %>%
      group_by(d_v, activity, year) %>%
      summarise(fvps = sum(fvps)) %>%
      ungroup() %>%
      group_by(d_v, activity) %>%
      mutate(fvps_cum = cumsum(fvps)) %>%
      ungroup() %>%
      as.data.table()
    
    # Extend back/forth for attractive plotting
    plot_dt = 
      expand_grid(
        d_v      = unique(d_v_dt$d_v), 
        activity = unique(d_v_dt$activity), 
        year     = o$years) %>%
      left_join(y  = d_v_dt, 
                by = c("d_v", "activity", "year")) %>%
      group_by(d_v, activity) %>%
      fill(fvps_cum, .direction = "down")  %>%
      replace_na(list(fvps_cum = 0)) %>%
      ungroup() %>%
      as.data.table()
    
    # Plot FVPs over time for each disease or vaccine
    g = ggplot(plot_dt) + 
      aes(x = year, 
          y = !!sym(y), 
          fill = activity) + 
      geom_area() + 
      facet_wrap(~d_v) + 
      # Prettify x axis...
      scale_x_continuous(
        limits = c(min(o$years), max(o$years)), 
        expand = expansion(mult = c(0, 0)), 
        breaks = seq(min(o$years), max(o$years), by = 10)) +  
      # Prettify y axis...
      scale_y_continuous(
        name   = "Total receiving full primary or booster schedule (in millions)", 
        labels = comma,
        expand = expansion(mult = c(0, NA)))
    
    # Prettify theme
    g = g + theme_classic() + 
      theme(axis.title.x  = element_blank(),
            axis.title.y  = element_text(
              size = 20, margin = margin(l = 10, r = 20)),
            axis.text     = element_text(size = 9),
            axis.text.x   = element_text(hjust = 1, angle = 50), 
            axis.line     = element_blank(),
            strip.text    = element_text(size = 12),
            strip.background = element_blank(), 
            panel.border  = element_rect(
              linewidth = 0.5, fill = NA),
            panel.spacing = unit(1, "lines"),
            panel.grid.major.y = element_line(linewidth = 0.5),
            legend.title  = element_blank(),
            legend.text   = element_text(size = 12),
            legend.key    = element_blank(),
            legend.position = "right", 
            legend.key.height = unit(2, "lines"),
            legend.key.width  = unit(2, "lines"))
    
    # Save to file
    save_fig(g, paste0("FVPs by ", d_v), dir = save_dir)
  }
}

# ---------------------------------------------------------
# Plot smoothed FVP for static model pathogens
# ---------------------------------------------------------
plot_smooth_fvps = function() {
  
  message("  > Plotting smoothed FVPs (static model pathogens)")
  
  # ---- Plot 1: data vs smoothing ----
  
  # Smoothing data - group by country and age
  smooth_dt = table("smoothed_fvps") %>%
    mutate(country_age = paste1(country, age), 
           fvps_smooth = fvps_smooth / 1e6, 
           fvps        = fvps / 1e6) %>%
    append_v_a_name() %>%
    select(v_a_name, country, country_age, 
           year, fvps, fvps_smooth)
  
  # Plot the data with associated smoothing
  g1 = ggplot(smooth_dt) +
    aes(x = year, 
        colour = country, 
        group  = country_age) +
    # Plot data...
    geom_point(
      mapping = aes(y = fvps),
      show.legend = FALSE) +
    # Plot smoothing...
    geom_line(
      mapping = aes(y = fvps_smooth),
      show.legend = FALSE) +
    # Simple faceting...
    facet_wrap(~v_a_name) + 
    # Prettify x axis...
    scale_x_continuous(
      expand = expansion(mult = c(0, 0)), 
      breaks = pretty_breaks()) +  
    # Prettify y axis...
    scale_y_continuous(
      name   = "Fully vaccinated people (in millions)", 
      labels = comma,
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05)))
  
  # Prettify theme
  g1 = g1 + theme_classic() + 
    theme(axis.title.x  = element_blank(),
          axis.title.y  = element_text(
            size = 20, margin = margin(l = 10, r = 20)),
          axis.text     = element_text(size = 10),
          # axis.text.x   = element_text(hjust = 1, angle = 50), 
          axis.line     = element_blank(),
          strip.text    = element_text(size = 14),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"))
  
  # ---- Plot 2: smoothing error ----
  
  # Total error for each vaccine
  diagnostic_dt = smooth_dt %>%
    # Summarise for each vaccine...
    group_by(v_a_name) %>%
    summarise(fvps = sum(fvps),
              fvps_smooth = sum(fvps_smooth)) %>%
    ungroup() %>%
    # Take the absolute difference...
    mutate(diff = fvps - fvps_smooth,
           abs  = abs(diff)) %>%
    # Calculate total error...
    mutate(err = round(100 * abs / fvps, 2), 
           err = paste0("Error: ", err, "%")) %>%
    # Append error to vaccine description...
    mutate(v_a_name = paste0(v_a_name, "\n", err)) %>%
    # Set plotting order by abs diff...
    arrange(abs) %>%
    mutate(v_a_name = fct_inorder(v_a_name)) %>%
    as.data.table()
  
  # Plot total smoothing errors
  g2 = ggplot(diagnostic_dt) +
    aes(x = v_a_name, 
        y = diff, 
        fill = abs) +
    geom_col(show.legend = FALSE) +
    coord_flip() + 
    # Prettify x axis (noting coord_flip)...
    scale_y_continuous(
      name   = "Total difference in FVPs after smoothing (in millions)", 
      labels = comma)
  
  # Prettify theme
  g2 = g2 + theme_classic() + 
    theme(axis.text     = element_text(size = 12),
          axis.title.x  = element_text(
            size = 20, margin = margin(b = 10, t = 20)),
          axis.title.y  = element_blank(),
          axis.line     = element_blank(),
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA))
  
  # ---- Save diagnostic plots ----
  
  # Directory to save to
  save_dir = "data_visualisation"
  
  # Save figures to file
  save_fig(g1, "Data smoothing outcomes",   dir = save_dir)
  save_fig(g2, "Data smoothing difference", dir = save_dir)
}

# ---------------------------------------------------------
# Plot coverage density by disease
# ---------------------------------------------------------
plot_coverage = function() {
  
  # All coverage values by disease, vaccine, and age class
  all_coverage_dt = table("coverage") %>%
    # Classify by age group...
    mutate(age_group = ifelse(age == 0, "infant", "other")) %>%
    # Classify by disease...
    left_join(y  = table("v_a"), 
              by = "v_a_id") %>%
    left_join(y  = table("d_v_a"), 
              by = c("vaccine", "activity")) %>%
    # Tidy up...
    select(disease, vaccine, age_group, coverage) %>%
    arrange(disease, vaccine, age_group)
  
  # Details for file destination
  save_name = "Coverage value density"
  save_dir  = "data_visualisation"
  
  # Produce plot for both disease and vaccine
  for (d_v in c("disease", "vaccine")) {
    
    # Select appropriate column: disease or vaccine
    plot_dt = all_coverage_dt %>%
      select(d_v = !!d_v, age_group, coverage)
    
    # Plot coverage value density
    g = ggplot(plot_dt) + 
      aes(x = coverage, 
          y = after_stat(scaled), 
          colour = age_group,
          fill   = age_group) + 
      geom_density(alpha = 0.2) + 
      facet_wrap(~d_v) + 
      # Prettify x axis...
      scale_x_continuous(
        name   = "Coverage",
        labels = percent,
        limits = c(0, 1), 
        expand = expansion(mult = c(0, 0)), 
        breaks = pretty_breaks()) +  
      # Prettify y axis...
      scale_y_continuous(
        name   = "Density", 
        limits = c(0, 1), 
        expand = expansion(mult = c(0, 0.1)), 
        breaks = pretty_breaks())
    
    # Prettify theme
    g = g + theme_classic() + 
      theme(axis.text     = element_text(size = 10),
            axis.title.x  = element_text(
              size = 20, margin = margin(b = 10, t = 20)),
            axis.title.y  = element_text(
              size = 20, margin = margin(l = 10, r = 20)),
            axis.line     = element_blank(),
            strip.text    = element_text(size = 14),
            strip.background = element_blank(), 
            panel.border  = element_rect(
              linewidth = 0.5, fill = NA),
            panel.spacing = unit(1, "lines"),
            legend.title  = element_blank(),
            legend.text   = element_text(size = 14),
            legend.key    = element_blank(),
            legend.position = "right", 
            legend.key.height = unit(2, "lines"),
            legend.key.width  = unit(2, "lines"))
    
    # Save to file
    save_fig(g, save_name, d_v, dir = save_dir)
  }
}

# ---------------------------------------------------------
# Plot age targets as defined by WIISE and VIMC coverage data
# ---------------------------------------------------------
plot_coverage_age_density = function() {
  
  message("  > Plotting coverage data density by age")
  
  # Construct plotting datatable
  plot_dt = table("coverage_source") %>%
    mutate(trans_age = pmax(age, 1), .after = age) %>%
    append_v_a_name()
  
  # Plot age density of coverage data by source
  g = ggplot(plot_dt) +
    aes(x = trans_age, 
        y = after_stat(scaled), 
        colour = source,
        fill   = source) +
    geom_density(alpha = 0.2) +
    # Facet with strip text wrapping...
    facet_wrap(
      facets   = vars(v_a_name), 
      labeller = label_wrap_gen(width = 24)) + 
    # Prettify x axis...
    scale_x_continuous(
      name   = "Age (log2 scale)",
      trans  = "log2", 
      limits = c(1, 2^7), 
      expand = c(0, 0), 
      breaks = 2^(0:7)) +  
    # Prettify y axis...
    scale_y_continuous(
      name   = "Density", 
      limits = c(0, 1), 
      expand = expansion(mult = c(0, 0.1)), 
      breaks = pretty_breaks())
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.text     = element_text(size = 10),
          axis.title.x  = element_text(
            size = 20, margin = margin(b = 10, t = 20)),
          axis.title.y  = element_text(
            size = 20, margin = margin(l = 10, r = 20)),
          axis.line     = element_blank(),
          strip.text    = element_text(size = 12),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"),
          legend.title  = element_blank(),
          legend.text   = element_text(size = 12),
          legend.key    = element_blank(),
          legend.position = "right", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save to file
  save_fig(g, "Coverage density by age", dir = "data_visualisation")
}

# ---------------------------------------------------------
# Plot Global Burden of Disease death estimates by age
# ---------------------------------------------------------
plot_gbd_estimates = function() {
  
  message("  > Plotting GBD death estimates by age")
  
  # Define age groups and associated upper bounds
  age_groups = c(
    "Infants"   = 1,
    "1-4 years" = 5,
    "5-14"      = 15,
    "15-49"     = 50,
    "50-69"     = 70,
    "70+ years" = max(o$ages))
  
  # Map each age bin to respective age group
  age_group_dt = data.table(age = o$ages) %>%
    mutate(group_idx = match(age, age_groups),
           group = names(age_groups[group_idx])) %>%
    fill(group, .direction = "up") %>%
    select(age, age_group = group)
  
  # Load GBD estimates and categorise into age groups
  plot_dt = table("gbd_estimates") %>%
    append_d_v_t_name() %>%
    left_join(y  = age_group_dt,
              by = "age") %>%
    # Summarise for broad age groups...
    group_by(disease, year, age_group) %>%
    summarise(deaths = sum(deaths_disease)) %>%
    ungroup() %>%
    # Set factors for meaningful plotting order...
    mutate(age_group = factor(age_group, names(age_groups))) %>%
    as.data.table()
  
  # Plot deaths over time by age group
  g = ggplot(plot_dt) +
    aes(x = year, 
        y = deaths, 
        fill = age_group) +
    geom_bar(stat = "identity") +
    facet_wrap(~disease, scales = "free_y") + 
    # Set colour scheme...
    scale_fill_manual(
      name   = "Age group",
      values = colour_scheme(
        map = "brewer::paired", 
        n   = length(age_groups))) +
    # Prettify y axis...
    scale_y_continuous(
      name   = "GBD-estimated number of deaths", 
      labels = comma,
      expand = expansion(mult = c(0, 0.05)), 
      breaks = pretty_breaks())
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title.x  = element_blank(),
          axis.title.y  = element_text(
            size = 20, margin = margin(l = 10, r = 20)),
          axis.text     = element_text(size = 10),
          axis.line     = element_blank(),
          strip.text    = element_text(size = 14),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"),
          legend.title  = element_blank(),
          legend.text   = element_text(size = 14),
          legend.key    = element_blank(),
          legend.position = "right", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save to file
  save_fig(g, "GBD deaths by age group", dir = "data_visualisation")
}

# ---------------------------------------------------------
# Plot static model pathogen vaccine efficacy with immunity decay
# ---------------------------------------------------------
plot_vaccine_efficacy = function() {
  
  message("  > Plotting vaccine efficacy profiles")
  
  # Function to group similar vaccines but split by dose
  schedule_fn = function(dt) {
    
    # Append descriptive columns
    shedule_dt = dt %>%
      # Primary schedule or booster dose...
      mutate(schedule = ifelse(
        !str_detect(vaccine, "_BX$"), "primary", "booster")) %>%
      mutate(schedule = factor(schedule, c("primary", "booster"))) %>%
      # Append disease-vaccine name
      mutate(vaccine = str_remove(vaccine, "[0-9]|(_BX)"), 
             d_v = paste0(disease, " (", vaccine, ")"), 
             d_v = fct_inorder(d_v))
    
    return(shedule_dt)
  }
  
  # Load vaccine efficacy profiles
  plot_dt = table("vaccine_efficacy_profiles") %>%
    schedule_fn() %>%
    select(d_v, schedule, time, profile)
  
  # Load data used to calculate these profiles
  data_dt = table("vaccine_efficacy") %>%
    left_join(y  = table("d_v"), 
              by = "vaccine") %>%
    schedule_fn() %>%
    select(d_v, schedule, 
           init     = efficacy, 
           time     = decay_x, 
           halflife = decay_y) %>%
    pivot_longer(cols = c(init, halflife), 
                 values_to = "profile") %>%
    mutate(time = ifelse(name == "init", 0, time)) %>%
    filter(!is.na(time), 
           !is.na(profile)) %>%
    select(all_of(names(plot_dt))) %>%
    as.data.table()
  
  # Plot vaccine efficacy with waning immunity (if any)
  g = ggplot(plot_dt) + 
    aes(x = time, y = profile, colour = d_v) + 
    geom_line(linewidth = 2) + 
    geom_point(data = data_dt, 
               size = 5) + 
    facet_grid(~schedule) +
    # Prettify x axis...
    scale_x_continuous(
      name   = "Years after completion of full schedule",
      expand = expansion(mult = c(0.05, 0.05)), 
      breaks = pretty_breaks()) +
    # Prettify y axis...
    scale_y_continuous(
      name   = "Vaccine efficacy (death reduction)", 
      labels = percent,
      limits = c(0, 1), 
      expand = expansion(mult = c(0, 0.05)), 
      breaks = pretty_breaks()) + 
    # Prettify legend (needed for y spacing to take effect)...
    guides(fill = guide_legend(byrow = TRUE))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.text     = element_text(size = 11),
          axis.title.x  = element_text(
            size = 20, margin = margin(b = 10, t = 20)),
          axis.title.y  = element_text(
            size = 20, margin = margin(l = 10, r = 20)),
          axis.line     = element_blank(),
          strip.text    = element_text(size = 16),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"),
          panel.grid.major.y = element_line(linewidth = 0.25),
          legend.title  = element_blank(),
          legend.text   = element_text(size = 11),
          legend.key    = element_blank(),
          legend.position = "right", 
          legend.spacing.y  = unit(1, "lines"),
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save figure to file
  save_fig(g, "Vaccine efficacy profiles", dir = "static_models")
}

# ---------------------------------------------------------
# Plot effective coverage with waning immunity for static model pathogens
# ---------------------------------------------------------
plot_effective_coverage = function() {
  
  message("  > Plotting effective coverage by year and age")
  
  # Plot only up to a certain age
  age_max = 50
  
  # Manually define appropriate number of colours 
  colours = colour_scheme("pals::brewer.blues", n = 8)
  
  # Repeat for disease and vaccine type
  for (by in c("disease", "type")) {
    
    # Load previously calculated total coverage file
    effective_dt = read_rds("static", "effective_coverage", by)
    
    # Population weight over all countries
    plot_dt = effective_dt %>%
      append_d_v_t_name() %>%
      select(country, by = !!by, year, age, coverage) %>%
      filter(age <= age_max) %>%
      # Append population size...
      left_join(y  = table("wpp_pop"),
                by = c("country", "year", "age")) %>%
      mutate(n = pop * coverage) %>%
      # Population weighted coverage...
      group_by(by, year, age) %>%
      summarise(effective_coverage = sum(n / sum(pop))) %>%
      ungroup() %>%
      as.data.table()
    
    # Plot each pathogen by year ana age
    g = ggplot(plot_dt) + 
      aes(x = year, y = age, fill = effective_coverage) + 
      geom_tile() +
      facet_wrap(~by) + 
      # Set continuous colour bar...
      scale_fill_gradientn(
        colours = colours, 
        limits  = c(0, 1), 
        breaks  = pretty_breaks(), 
        label   = percent,
        guide   = guide_colourbar(
          title = "Effective coverage")) + 
      # Prettify x axis...
      scale_x_continuous(
        expand = c(0, 0), 
        breaks = seq(min(o$years), max(o$years), by = 5)) +
      # Prettify y axis...
      scale_y_continuous(
        name   = "Vaccine efficacy (death reduction)", 
        expand = c(0, 0), 
        breaks = pretty_breaks())
    
    # Prettify theme
    g = g + theme_classic() + 
      theme(axis.title.x  = element_blank(),
            axis.title.y  = element_text(
              size = 20, margin = margin(l = 10, r = 20)),
            axis.text     = element_text(size = 10),
            axis.text.x   = element_text(hjust = 1, angle = 50), 
            axis.line     = element_blank(),
            strip.text    = element_text(size = 14),
            strip.background = element_blank(), 
            panel.border  = element_rect(
              linewidth = 0.5, fill = NA),
            panel.spacing = unit(1, "lines"),
            legend.position    = "top", 
            legend.title  = element_text(
              size = 14, margin = margin(r = 20)),
            legend.text   = element_text(size = 10),
            legend.key.height = unit(1.5, "lines"),
            legend.key.width  = unit(6,   "lines"))
    
    # Save figure to file
    save_name = "Effective coverage by year and age"
    save_fig(g, save_name, by, dir = "static_models")
  }
}

# ---------------------------------------------------------
# Plot deaths and DALYs averted for static model pathogens
# ---------------------------------------------------------
plot_static = function() {
  
  message("  > Plotting static model impact results")
  
  # Deaths disease/averted dictionary
  metric_dict = c(
    deaths_disease = "Estimated disease-specific deaths (GBD 2019)", 
    deaths_averted = "Estimated deaths averted deaths from static model")
  
  # Associated colours
  metric_colours = c("darkred", "navyblue")
  
  # ---- Plot by disease ----
  
  # Load previously calculated total coverage file
  averted_dt = read_rds("static", "deaths_averted_disease")
  
  # Summarise results over country and age
  disease_dt = averted_dt %>%
    pivot_longer(cols = starts_with("deaths"), 
                 names_to = "metric") %>%
    # Summarise over countries...
    group_by(disease, year, metric) %>%
    summarise(value = sum(value) / 1e6) %>%
    ungroup() %>%
    arrange(metric, disease, year) %>%
    # Recode deaths disease/averted...
    append_d_v_t_name() %>%
    mutate(metric = recode(metric, !!!metric_dict), 
           metric = factor(metric, metric_dict)) %>%
    as.data.table()
  
  # Plot deaths and deaths averted by disease
  g = ggplot(disease_dt) + 
    aes(x = year, 
        y = value, 
        colour = metric) + 
    geom_line(linewidth = 2) + 
    facet_wrap(~disease) + 
    # Set colours and prettify legend...
    scale_colour_manual(values = metric_colours) + 
    guides(color = guide_legend(
      byrow = TRUE, ncol = 1)) +
    # Prettify x axis...
    scale_x_continuous(
      # limits = c(min(o$years), max(o$years)), 
      expand = expansion(mult = c(0, 0)), 
      breaks = pretty_breaks()) +  
    # Prettify y axis...
    scale_y_continuous(
      name   = "Number of people (in millions)", 
      labels = comma,
      expand = expansion(mult = c(0, 0.05)))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title.x  = element_blank(),
          axis.title.y  = element_text(
            size = 20, margin = margin(l = 10, r = 20)),
          axis.text     = element_text(size = 10),
          axis.text.x   = element_text(hjust = 1, angle = 50), 
          axis.line     = element_blank(),
          strip.text    = element_text(size = 14),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"),
          panel.grid.major.y = element_line(linewidth = 0.5),
          legend.title  = element_blank(),
          legend.text   = element_text(size = 14),
          legend.key    = element_blank(),
          legend.position = "bottom", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(3, "lines"))
  
  # Save figure to file
  save_fig(g, "Deaths averted by disease", dir = "static_models")
  
  # ---- Plot by vaccine ----
  
  # Load previously calculated total coverage file
  averted_dt = read_rds("static", "deaths_averted_vaccine")
  
  # Summarise results over country
  vaccine_dt = averted_dt %>%
    # Convert from d_v_a to v_a...
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_id") %>%
    left_join(y  = table("vaccine"), 
              by = "vaccine") %>%
    # Summarise over countries...
    group_by(vaccine_name, year) %>%
    summarise(deaths_averted = sum(impact)) %>%
    ungroup() %>%
    arrange(vaccine_name, year) %>%
    as.data.table()
  
  # Plot deaths and deaths averted by disease
  g = ggplot(vaccine_dt) + 
    aes(x = year, 
        y = deaths_averted, 
        colour = vaccine_name) + 
    geom_line(linewidth   = 2, 
              show.legend = FALSE) + 
    facet_wrap(~vaccine_name,
               scales = "free_y") + 
    # Prettify x axis...
    scale_x_continuous(
      # limits = c(min(o$years), max(o$years)), 
      expand = expansion(mult = c(0, 0)), 
      breaks = pretty_breaks()) +  
    # Prettify y axis...
    scale_y_continuous(
      name   = "Deaths averted by vaccine", 
      labels = comma,
      limits = c(0, NA), 
      expand = expansion(mult = c(0, 0.05)))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title.x  = element_blank(),
          axis.title.y  = element_text(
            size = 18, margin = margin(l = 10, r = 20)),
          axis.text     = element_text(size = 9),
          axis.text.x   = element_text(hjust = 1, angle = 50), 
          axis.line     = element_blank(),
          strip.text    = element_text(size = 12),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"),
          panel.grid.major.y = element_line(linewidth = 0.5))
  
  # Save figure to file
  save_fig(g, "Deaths averted by vaccine", dir = "static_models")
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
    left_join(y  = table("disease"), 
              by = "disease") %>%
    filter(source == "vimc") %>%
    pull(d_v_a_id) %>%
    lapply(load_data_fn) %>%
    rbindlist()
  
  # ---- Produce plot ----
  
  # Construct tidy plotting datatable
  plot_dt = data_dt %>%
    pivot_longer(cols = -c(d_v_a_id, target), 
                 names_to = "covariate") %>%
    arrange(d_v_a_id, covariate, target) %>%
    format_d_v_a_name() %>%
    select(d_v_a_name, target, covariate, value) %>%
    as.data.table()
  
  # Plot covariates vs imputation target
  g = ggplot(plot_dt) +
    aes(x = target, 
        y = value, 
        colour = covariate) +
    geom_point(
      alpha = 0.2, 
      shape = 16, 
      show.legend = FALSE) +
    # Facet in 2 dims with label wrapping...
    facet_grid(
      rows     = vars(covariate), 
      cols     = vars(d_v_a_name),
      labeller = label_wrap_gen(width = 20), 
      scales   = "free") + 
    # Set colour scheme...
    scale_colour_manual(
      values = colour_scheme(
        map = "pals::kovesi.rainbow", 
        n   = n_unique(plot_dt$covariate))) +
    # Prettify x axis...
    scale_x_continuous(
      name   = "Normalised predictor", 
      limits = c(0, 1), 
      expand = c(0, 0), 
      breaks = pretty_breaks()) +  
    # Prettify y axis...
    scale_y_continuous(
      name   = "Normalised response (impact per FVP)", 
      limits = c(0, 1), 
      expand = c(0, 0), 
      breaks = pretty_breaks())
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.text     = element_text(size = 6),
          axis.text.x   = element_text(hjust = 1, angle = 50),
          axis.title.x  = element_text(
            size = 16, margin = margin(l = 4, r = 8)),
          axis.title.y  = element_text(
            size = 16, margin = margin(b = 4, t = 8)),
          axis.line     = element_blank(),
          strip.text    = element_text(size = 8),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(0.4, "lines"))
  
  # Save figure to file
  save_fig(g, "Covariate relationships", dir = "imputation")
}

# ---------------------------------------------------------
# Plot truth vs predicted for imputation training data
# ---------------------------------------------------------
plot_impute_quality = function() {
  
  message("  > Plotting imputation quality of fit")
  
  # ---- Load results from fitting ----
  
  # Function to load imputation results
  load_results_fn = function(id)
    result = read_rds("impute", "impute", id)$result
  
  # Load imputation results for all d-v-a
  results_dt = table("d_v_a") %>%
    left_join(y  = table("disease"), 
              by = "disease") %>%
    filter(source == "vimc") %>%
    pull(d_v_a_id) %>%
    lapply(load_results_fn) %>%
    rbindlist()
  
  # ---- Construct plotting datatables ----
  
  # Prepare datatable for plotting
  plot_dt = results_dt %>%
    filter(!is.na(target)) %>%
    select(-country) %>%
    format_d_v_a_name() %>%
    # Remove target outliers for better normalisation...
    group_by(d_v_a_name) %>%
    mutate(lower = mean(target) - 3 * sd(target), 
           upper = mean(target) + 3 * sd(target), 
           outlier = target < lower | target > upper) %>%
    ungroup() %>%
    filter(outlier == FALSE) %>%
    select(-outlier, -lower, -upper, -d_v_a_id) %>%
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
    aes(x = target, 
        y = predict, 
        color = d_v_a_name) +
    # Plot truth vs predicted...
    geom_point(alpha = 0.35, 
               shape = 16, 
               show.legend = FALSE) +
    # For square axes...
    geom_blank(data = blank_dt) +
    # To see quality of predict vs target...
    geom_abline(colour = "black") + 
    # Simple faceting with wrap labelling...
    facet_wrap(
      facets   = vars(d_v_a_name), 
      labeller = label_wrap_gen(width = 30), 
      scales   = "free") + 
    # Prettify x axis...
    scale_x_continuous(
      name   = "Imputation target", 
      labels = scientific,
      limits = c(0, NA), 
      expand = c(0, 0), 
      breaks = pretty_breaks()) +  
    # Prettify y axis...
    scale_y_continuous(
      name   = "Imputation prediction", 
      labels = scientific,
      limits = c(0, NA),
      expand = c(0, 0), 
      breaks = pretty_breaks())
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.text     = element_text(size = 8),
          axis.text.x   = element_text(hjust = 1, angle = 50),
          axis.title.x  = element_text(
            size = 18, margin = margin(l = 10, r = 20)),
          axis.title.y  = element_text(
            size = 18, margin = margin(b = 10, t = 20)),
          axis.line     = element_blank(),
          strip.text    = element_text(size = 12),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(0.5, "lines"))
  
  # Save figure to file
  save_fig(g, "Imputation fit", dir = "imputation")
}

# ---------------------------------------------------------
# Plot country-aggregated imputation errors
# ---------------------------------------------------------
plot_impute_countries = function() {
  
  message("  > Plotting country-aggregated imputation errors")
  
  # ---- Load results from fitting ----
  
  # Function to load imputation results
  load_results_fn = function(id)
    result = read_rds("impute", "impute", id)$result
  
  # Load imputation results for all d-v-a
  results_dt = table("d_v_a") %>%
    left_join(y  = table("disease"), 
              by = "disease") %>%
    filter(source == "vimc") %>%
    pull(d_v_a_id) %>%
    lapply(load_results_fn) %>%
    rbindlist() %>%
    format_d_v_a_name()
  
  # ---- Plot 1: annual error by country ----
  
  # Truth vs predicted over time for training data
  annual_dt = results_dt %>%
    select(d_v_a_name, country, year, 
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
    # Country truth...
    geom_line(aes(y = vimc, colour = country), 
              linewidth = 0.5) +
    # Country predicted...
    geom_line(aes(y = impute, colour = country), 
              linewidth = 0.5, 
              linetype  = "dashed") +
    facet_wrap(~d_v_a_name, scales = "free_y") + 
    # Remove legend...
    theme(legend.position = "none")
  
  # Save figure to file
  save_fig(g, "Imputation error annual", dir = "imputation")
  
  # ---- Plot 2: total error by country ----
  
  # Where imputed countries lie in terms of magnitude
  total_dt = results_dt %>%
    # Take cumulative values for each country...
    group_by(d_v_a_name, country) %>%
    summarise(truth   = max(impact_cum), 
              predict = max(impact_impute)) %>%
    ungroup() %>%
    # VIMC as truth-predict scatter, imputed along diagonal...
    mutate(source = ifelse(is.na(truth), "impute", "vimc"), 
           truth  = ifelse(is.na(truth), predict, truth)) %>%
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
  
  # Save figure to file
  save_fig(g, "Imputation error total", dir = "imputation")
}

# ---------------------------------------------------------
# Exploratory plots of data used to fit impact functions
# ---------------------------------------------------------
plot_impact_data = function() {
  
  message("  > Plotting impact function fitting data")
  
  # Load data used for impact function fitting
  data_dt = read_rds("impact", "data") %>%
    format_d_v_a_name()
  
  # ---- Plot 1: impact per FVP over time ----
  
  # Impact per FVP over time
  g1 = ggplot(data_dt) +
    aes(x = year, 
        y = impact_fvp, 
        colour = country) +
    geom_line(show.legend = FALSE) +
    # Faceting with wrap labelling...
    facet_wrap(
      facets   = vars(d_v_a_name), 
      labeller = label_wrap_gen(width = 30), 
      scales   = "free_y") + 
    # Set colour scheme...
    scale_colour_manual(
      values = colour_scheme(
        map = "viridis::viridis", 
        n   = n_unique(all_countries()))) +
    # Prettify x axis...
    scale_x_continuous(
      expand = c(0, 0), 
      breaks = pretty_breaks()) +  
    # Prettify y axis...
    scale_y_continuous(
      name   = "Impact per fully vaccinated", 
      labels = scientific,
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05)), 
      breaks = pretty_breaks())
  
  # Prettify theme
  g1 = g1 + theme_classic() + 
    theme(axis.text     = element_text(size = 8),
          axis.text.x   = element_text(hjust = 1, angle = 50),
          axis.title.x  = element_text(
            size = 16, margin = margin(l = 10, r = 20)),
          axis.title.y  = element_text(
            size = 16, margin = margin(b = 10, t = 20)),
          axis.line     = element_blank(),
          strip.text    = element_text(size = 10),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(0.5, "lines"))
  
  # Save figure to file
  save_fig(g1, "Data - impact ratio", dir = "impact_functions")
  
  # ---- Plot 2: cumulative impact vs cumulative FVP ----
  
  # Cumulative FVPs vs cumulative deaths averted
  g2 = ggplot(data_dt) + 
    aes(x = fvps, 
        y = impact, 
        colour = country) +
    geom_line(show.legend = FALSE) +
    # Faceting with wrap labelling...
    facet_wrap(
      facets   = vars(d_v_a_name), 
      labeller = label_wrap_gen(width = 30), 
      scales   = "free") + 
    # Set colour scheme...
    scale_colour_manual(
      values = colour_scheme(
        map = "viridis::viridis", 
        n   = n_unique(all_countries()))) +
    # Prettify x axis...
    scale_x_continuous(
      name   = "Cumulative FVPs per capita (including new birth cohorts)", 
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05)), 
      breaks = pretty_breaks()) +  
    # Prettify y axis...
    scale_y_continuous(
      name   = "Cumulative impact per capita (including new birth cohorts)", 
      labels = scientific,
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05)), 
      breaks = pretty_breaks())
  
  # Prettify theme
  g2 = g2 + theme_classic() + 
    theme(axis.text     = element_text(size = 8),
          axis.title.x  = element_text(
            size = 16, margin = margin(l = 10, r = 20)),
          axis.title.y  = element_text(
            size = 16, margin = margin(b = 10, t = 20)),
          axis.line     = element_blank(),
          strip.text    = element_text(size = 10),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(0.5, "lines"))
  
  # Save figure to file
  save_fig(g2, "Data - cumulative FVP vs impact", 
           dir = "impact_functions")
}

# ---------------------------------------------------------
# Plot function selection statistics
# ---------------------------------------------------------
plot_model_selection = function() {
  
  message("  > Plotting impact model selection")
  
  # Load stuff: best fit functions and associtaed coefficients
  best_dt = read_rds("impact", "best_model") %>%
    format_d_v_a_name()
  
  # ---- Plot function count ----
  
  # Simple plotting function with a few features
  plot_selection = function(var, type = "count", stat = "number") {
    
    # Determine order - with 'focus' function first
    fn_dict = fn_set(dict = TRUE)
    
    # Number of times each model is optimal
    selection_dt = best_dt %>% 
      rename(var = !!var) %>% 
      # Number and proportion of each fn...
      count(var, fn, name = "number") %>%
      group_by(var) %>%
      mutate(total = sum(number)) %>%
      ungroup() %>%
      mutate(proportion = number / total) %>%
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
      g = ggplot(selection_dt) + 
        aes(x    = var, 
            y    = val, 
            fill = fn) + 
        geom_col() + 
        coord_flip() + 
        # Prettify x axis (noting coord_flip)...
        scale_y_continuous(
          name   = paste(first_cap(stat), "of countries"), 
          expand = expansion(mult = c(0, 0.05)), 
          breaks = pretty_breaks())
      
      # Prettify theme
      g = g + theme_classic() + 
        theme(axis.text     = element_text(size = 12),
              axis.title.x  = element_text(
                size = 20, margin = margin(b = 10, t = 20)),
              axis.title.y  = element_blank(),
              axis.line     = element_blank(),
              panel.border  = element_rect(
                linewidth = 0.5, fill = NA),
              legend.title  = element_blank(),
              legend.text   = element_text(size = 12),
              legend.key    = element_blank(),
              legend.position = "right", 
              legend.key.height = unit(2, "lines"),
              legend.key.width  = unit(2, "lines"))
    }
    
    # Check figure type flag
    if (type == "density") {
      
      # Density of occurances
      g = ggplot(selection_dt) + 
        aes(x    = val, 
            fill = fn) +
        geom_bar() + 
        # Prettify x axis...
        scale_x_continuous(
          name   = "Number of model selections for a country", 
          breaks = pretty_breaks()) +  
        # Prettify y axis...
        scale_y_continuous(
          name   = paste(first_cap(stat), "of countries"), 
          expand = expansion(mult = c(0, 0.05)), 
          breaks = pretty_breaks())
      
      # Prettify theme
      g = g + theme_classic() + 
        theme(axis.text     = element_text(size = 12),
              axis.title.x  = element_text(
                size = 20, margin = margin(b = 10, t = 20)),
              axis.title.y  = element_text(
                size = 20, margin = margin(l = 10, r = 20)),
              axis.line     = element_blank(),
              panel.border  = element_rect(
                linewidth = 0.5, fill = NA),
              legend.title  = element_blank(),
              legend.text   = element_text(size = 12),
              legend.key    = element_blank(),
              legend.position = "right", 
              legend.key.height = unit(2, "lines"),
              legend.key.width  = unit(2, "lines"))
    }
    
    return(g)
  }
  
  # ---- A variety of plots ----
  
  # Figure sub-directory to save to 
  save_dir = "impact_functions"
  
  # Plot by disease-vaccine-activity
  g1 = plot_selection("d_v_a_name", stat = "number")
  g2 = plot_selection("d_v_a_name", stat = "proportion")
  
  # Filename stem
  save_name = "Selection by pathogen"
  
  # Save both figures
  save_fig(g1, save_name, "number",     dir = save_dir)
  save_fig(g2, save_name, "proportion", dir = save_dir)
  
  # Plot by country
  g3 = plot_selection("country", type = "density")
  
  # Save the last figure
  save_name = "Selection density by country"
  save_fig(g3, save_name, dir = save_dir)
}

# ---------------------------------------------------------
# Plot impact function evaluation
# ---------------------------------------------------------
plot_model_fits = function() {
  
  message("  > Plotting impact function fits")
  
  # Load data used for impact function fitting
  data_dt = read_rds("impact", "data") %>%
    format_d_v_a_name()
  
  # Evaluate only as far as we have data
  max_data = data_dt %>%
    group_by(d_v_a_name, country) %>%
    slice_max(fvps, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(d_v_a_name, country, x_max = fvps) %>%
    # Increment up one so we plot slightly past the data
    mutate(x_max = x_max + o$eval_x_scale / 100) %>%
    as.data.table()
  
  # Evaluate selected impact function
  best_fit = evaluate_impact_function() %>%
    format_d_v_a_name()
  
  # Apply max_data so we only plot up to data (or just past)
  plot_dt = best_fit %>%
    left_join(y  = max_data,
              by = c("d_v_a_name", "country")) %>%
    filter(fvps < x_max) %>%
    format_d_v_a_name() %>%
    select(d_v_a_name, country, fvps, impact)
  
  # Plot function evaluation against the data
  g = ggplot(plot_dt) +
    aes(x = fvps, 
        y = impact, 
        colour = country) +
    # Plot data, then fit on top...
    geom_point(
      data = data_dt,
      size = 0.75,
      alpha = 0.5,
      show.legend = FALSE) +
    geom_line(show.legend = FALSE) +
    # Faceting with wrap labelling...
    facet_wrap(
      facets   = vars(d_v_a_name), 
      labeller = label_wrap_gen(width = 30), 
      scales   = "free") + 
    # Set colour scheme...
    scale_colour_manual(
      values = colour_scheme(
        map = "viridis::viridis", 
        n   = n_unique(all_countries()))) +
    # Prettify x axis...
    scale_x_continuous(
      name   = "Cumulative FVPs per capita (including new birth cohorts)", 
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05)), 
      breaks = pretty_breaks()) +  
    # Prettify y axis...
    scale_y_continuous(
      name   = "Cumulative impact per capita (including new birth cohorts)", 
      labels = scientific,
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05)), 
      breaks = pretty_breaks())
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.text     = element_text(size = 8),
          axis.text.x   = element_text(hjust = 1, angle = 50),
          axis.title.x  = element_text(
            size = 16, margin = margin(b = 10, t = 20)),
          axis.title.y  = element_text(
            size = 16, margin = margin(l = 10, r = 20)),
          axis.line     = element_blank(),
          strip.text    = element_text(size = 10),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(0.5, "lines"))
  
  save_fig(g, "Impact function evaluation", dir = "impact_functions")
}

# ---------------------------------------------------------
# Plot impact ratios - either all or initial only
# ---------------------------------------------------------
plot_impact_fvps = function(scope) {
  
  width = 0.3
  
  # Function for averaging (mean or median)
  avg_fn = get("mean")
  
  # ---- Load data based on scope ----
  
  # All time plot
  if (scope == "all_time") {
    
    message("  > Plotting all-time impact per FVP")
    
    # Load initial ratio data
    impact_dt = read_rds("impact", "data")
    
    # Set a descriptive y-axis title
    y_lab = "Impact per fully vaccinated person (log10 scale)"
  }
  
  # Initial year plot
  if (scope == "initial") {
    
    message("  > Plotting initial impact per FVP")
    
    # Load initial ratio data
    impact_dt = read_rds("impact", "initial_ratio") %>%
      rename(impact_fvp = initial_ratio)
    
    # Set a descriptive y-axis title
    y_lab = "Initial impact per FVP used for back projection (log10 scale)"
  }
  
  # ---- Classify by income status ----
  
  # Load income status dictionary
  income_dict = table("income_dict")
  
  # Load income status of each country
  income_dt = table("income_status") %>%
    filter(year == max(year)) %>%
    left_join(y  = income_dict, 
              by = "income") %>%
    select(country, income = income_name)
  
  # Classify by income group
  data_dt = impact_dt %>%
    filter(impact_fvp > 0) %>%
    # Append d_v_a description...
    format_d_v_a_name() %>%
    # Append income status description...
    left_join(y  = income_dt, 
              by = "country") %>%
    mutate(income = factor(
      x      = income, 
      levels = rev(income_dict$income_name))) %>%
    select(d_v_a = d_v_a_name, income, impact_fvp) %>%
    unique()
  
  # ---- Plotting coordinates ----
  
  # Set x values for each d_v_a
  x_major = data_dt %>%
    group_by(d_v_a) %>%
    summarise(order = avg_fn(impact_fvp)) %>%
    ungroup() %>%
    arrange(desc(order)) %>%
    mutate(x_major = 1 : n()) %>%
    select(-order) %>%
    as.data.table()
  
  # Offset each income group
  x_minor = data_dt %>%
    select(income) %>%
    unique() %>%
    arrange(income) %>%
    mutate(x_minor = seq(
      from = -width / 2, 
      to   = width / 2, 
      length.out = n()))
  
  # Width of offset of points within income group
  x_spray = width / (nrow(income_dict) * 2)
  
  # Bring all x and y values together
  plot_dt = data_dt %>%
    left_join(y  = x_major, 
              by = "d_v_a") %>%
    left_join(y  = x_minor, 
              by = "income") %>%
    mutate(x = x_major + x_minor) %>%
    select(d_v_a, income, x, y = impact_fvp) %>%
    arrange(x)
  
  # Extract bounds (after transformation)
  lb = floor(min(log10(data_dt$impact_fvp)))
  ub = ceiling(max(log10(data_dt$impact_fvp)))
  
  # ---- Vaccine and income status averages ----
  
  # Average by vaccine and income group
  income_average_dt = plot_dt %>%
    group_by(income, x) %>%
    summarise(y = avg_fn(y)) %>%
    ungroup() %>%
    arrange(x) %>%
    as.data.table()
  
  # Average by vaccine (over all income groups)
  vaccine_average_dt = data_dt %>%
    group_by(d_v_a) %>%
    summarise(y = avg_fn(impact_fvp)) %>%
    ungroup() %>%
    left_join(y  = x_major,
              by = "d_v_a") %>%
    select(d_v_a, x = x_major, y) %>%
    arrange(x) %>%
    as.data.table()
  
  # ---- Produce primary plot ----
  
  # Plot impact per FVP
  g = ggplot(plot_dt) +
    aes(x = x, y = y) +
    # Plot all points by vaccine and income...
    geom_point(
      mapping  = aes(colour = income),
      alpha    = 0.1,
      shape    = 16,
      stroke   = 0,
      position = position_jitter(
        width = x_spray,
        seed  = 1)) + 
    # Average of each vaccine...
    geom_point(
      data   = vaccine_average_dt,
      colour = "black",
      shape = 95, # Horizontal lines
      size  = 20) +
    # Average of each income group...
    geom_point(
      data    = income_average_dt, 
      mapping = aes(fill = income), 
      color   = "black", 
      shape   = 23, 
      size    = 3) + 
    # Prettify x axis...
    scale_x_continuous(
      breaks = x_major$x_major, 
      labels = x_major$d_v_a) +
    # Prettify y axis...
    scale_y_continuous(
      name   = y_lab,
      trans  = "log10",
      labels = scientific, 
      limits = c(10 ^ lb, 10 ^ ub),
      expand = c(0, 0),
      breaks = 10 ^ rev(ub : lb))
  
  # ---- Prettify and save ----
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.text     = element_text(size = 10),
          axis.text.x   = element_text(hjust = 1, angle = 50),
          axis.title.x  = element_blank(),
          axis.title.y  = element_text(
            size = 16, margin = margin(l = 10, r = 20)),
          axis.line     = element_blank(),
          panel.grid.major.y = element_line(linewidth = 0.25),
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"))
  
  # Filename to save to depending on scope
  save_stem  = "Density of impact per FVP"
  save_scope = str_replace(scope, "_", " ")
  
  # Save these figures to file
  save_fig(g, save_stem, save_scope, dir = "impact_functions")
  
  # Save all time plot as key manuscript figure
  if (scope == "all_time")
    save_fig(g, "Figure 3", dir = "manuscript")
}

# ---------------------------------------------------------
# Plot impact vs coverage by vaccine, income, and decade 
# ---------------------------------------------------------
plot_impact_coverage = function() {
  
  message("  > Plotting impact against coverage")
  
  browser()
  
  # Function for averaging (mean or median)
  avg_fn = get("mean")
  
  # ---- Classify by income status ----
  
  # Load initial ratio data
  impact_dt = read_rds("impact", "data")
  
  # Load income status dictionary
  income_dict = table("income_dict")
  
  # Load income status of each country
  income_dt = table("income_status") %>%
    filter(year == max(year)) %>%
    left_join(y  = income_dict, 
              by = "income") %>%
    select(country, income = income_name)
  
  browser()
  
  # Classify by income group
  data_dt = impact_dt %>%
    filter(impact_fvp > 0) %>%
    # Append income status description...
    left_join(y  = income_dt, 
              by = "country") %>%
    mutate(income = factor(
      x      = income, 
      levels = rev(income_dict$income_name))) %>%
    select(d_v_a = d_v_a_name, income, impact_fvp) %>%
    unique()
  
  # ---- Plotting coordinates ----
  
  # Set x values for each d_v_a
  x_major = data_dt %>%
    group_by(d_v_a) %>%
    summarise(order = avg_fn(impact_fvp)) %>%
    ungroup() %>%
    arrange(desc(order)) %>%
    mutate(x_major = 1 : n()) %>%
    select(-order) %>%
    as.data.table()
  
  # Offset each income group
  x_minor = data_dt %>%
    select(income) %>%
    unique() %>%
    arrange(income) %>%
    mutate(x_minor = seq(
      from = -width / 2, 
      to   = width / 2, 
      length.out = n()))
  
  # Width of offset of points within income group
  x_spray = width / (nrow(income_dict) * 2)
  
  # Bring all x and y values together
  plot_dt = data_dt %>%
    left_join(y  = x_major, 
              by = "d_v_a") %>%
    left_join(y  = x_minor, 
              by = "income") %>%
    mutate(x = x_major + x_minor) %>%
    select(d_v_a, income, x, y = impact_fvp) %>%
    arrange(x)
  
  # Extract bounds (after transformation)
  lb = floor(min(log10(data_dt$impact_fvp)))
  ub = ceiling(max(log10(data_dt$impact_fvp)))
  
  # ---- Vaccine and income status averages ----
  
  income_average_dt = plot_dt %>%
    group_by(income, x) %>%
    summarise(y = avg_fn(y)) %>%
    ungroup() %>%
    arrange(x) %>%
    as.data.table()
  
  vaccine_average_dt = data_dt %>%
    group_by(d_v_a) %>%
    summarise(y = avg_fn(impact_fvp)) %>%
    ungroup() %>%
    left_join(y  = x_major,
              by = "d_v_a") %>%
    select(d_v_a, x = x_major, y) %>%
    arrange(x) %>%
    as.data.table()
  
  # ---- Produce primary plot ----
  
  # Plot impact per FVP
  g = ggplot(plot_dt) +
    aes(x = x, y = y) +
    # Plot all points by vaccine and income...
    geom_point(
      mapping  = aes(colour = income),
      alpha    = 0.1,
      shape    = 16,
      stroke   = 0,
      position = position_jitter(
        width = x_spray,
        seed  = 1)) + 
    # Average of each vaccine...
    geom_point(
      data   = vaccine_average_dt,
      colour = "black",
      shape = 95, # Horizontal lines
      size  = 20) +
    # Average of each income group...
    geom_point(
      data    = income_average_dt, 
      mapping = aes(fill = income), 
      color   = "black", 
      shape   = 23, 
      size    = 3) + 
    # Prettify x axis...
    scale_x_continuous(
      breaks = x_major$x_major, 
      labels = x_major$d_v_a) +
    # Prettify y axis...
    scale_y_continuous(
      name   = "Impact per fully vaccinated person (log10 scale)",
      trans  = "log10",
      labels = scientific, 
      limits = c(10 ^ lb, 10 ^ ub),
      expand = c(0, 0),
      breaks = 10 ^ rev(ub : lb))
  
  # ---- Prettify and save ----
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.text     = element_text(size = 10),
          axis.text.x   = element_text(hjust = 1, angle = 50),
          axis.title.x  = element_blank(),
          axis.title.y  = element_text(
            size = 16, margin = margin(l = 10, r = 20)),
          axis.line     = element_blank(),
          panel.grid.major.y = element_line(linewidth = 0.25),
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"))
  
  # Filename to save to depending on scope
  save_stem  = "Impact - coverage by income and decade"
  save_scope = str_replace(scope, "_", " ")
  
  # Save these figures to file
  save_fig(g, save_stem, save_scope, dir = "impact_functions")
  
  # Save all time plot as key manuscript figure
  # if (scope == "all_time")
  #   save_fig(g, "Figure 3", dir = "manuscript")
}

# ---------------------------------------------------------
# Main results plot - historical impact over time
# ---------------------------------------------------------
plot_historical_impact = function() {
  
  message("  > Plotting historical impact")
  
  # Dictionary for temporal and cumulative subplots
  impact_dict = c(
    impact_cum = "Cumulative deaths averted (in millions)", 
    impact     = "Deaths averted per year (in millions)")
  
  # ---- Construct plotting data ----
  
  # Prepare final results
  results_dt = read_rds("results", "deaths_averted") %>%
    # Append full disease names...
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_id") %>%
    left_join(y  = table("disease"), 
              by = "disease") %>%
    # Cumulative results for each disease...
    group_by(disease_name, year) %>%
    summarise(impact = sum(impact) / 1e6) %>%
    mutate(impact_cum = cumsum(impact)) %>%
    ungroup() %>%
    rename(disease = disease_name) %>%
    # Tidy format for single plot...
    pivot_longer(cols = c(impact, impact_cum), 
                 names_to = "metric") %>%
    mutate(metric = recode(metric, !!!impact_dict), 
           metric = factor(metric, impact_dict)) %>%
    arrange(metric, disease, year) %>%
    as.data.table()
  
  # Construct labels: total FVPs over analysis timeframe
  label_dt = results_dt %>%
    filter(metric == impact_dict[["impact_cum"]]) %>%
    group_by(disease) %>%
    summarise(total = round(max(value), 1)) %>%
    ungroup() %>%
    mutate(total = paste0("Total: ", total, " million"), 
           label = paste0(disease, "\n", total)) %>%
    select(disease, disease_label = label) %>%
    as.data.table()
  
  # Append total labels to plotting data
  plot_dt = results_dt %>%
    left_join(y  = label_dt, 
              by = "disease") %>%
    select(disease_label, year, metric, value)
  
  # ---- Produce plot ----
  
  # Stacked yearly bar plot
  g = ggplot(plot_dt) +
    aes(x = year, 
        y = value, 
        fill = disease_label) + 
    geom_col() +
    # Facet by temporal-cumulative metric...
    facet_wrap(~metric, scales = "free_y") +
    # Set colours...
    scale_fill_manual(
      values = get_palette("disease")) + 
    # Prettify y axis...
    scale_y_continuous(
      labels = comma, 
      breaks = pretty_breaks(),
      expand = expansion(mult = c(0, 0.05))) +
    # Prettiy x axis...
    scale_x_continuous(
      limits = c(min(o$years) - 1, 
                 max(o$years) + 1), 
      expand = expansion(mult = c(0, 0)), 
      breaks = seq(min(o$years), max(o$years), by = 5)) +
    # Prettify legend (needed for y spacing to take effect)...
    guides(fill = guide_legend(byrow = TRUE))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title    = element_blank(),
          axis.text     = element_text(size = 11),
          axis.text.x   = element_text(hjust = 1, angle = 50), 
          axis.line     = element_blank(),
          strip.text    = element_text(size = 14),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"),
          panel.grid.major.y = element_line(linewidth = 0.25),
          legend.title  = element_blank(),
          legend.text   = element_text(size = 11),
          legend.key    = element_blank(),
          legend.position = "right", 
          legend.spacing.y  = unit(1, "lines"),
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save these figures to file
  save_fig(g, "Historical impact", dir = "historical_impact")
  save_fig(g, "Figure 2", dir = "manuscript")
}

# ---------------------------------------------------------
# Plot change in child mortality rates over time
# ---------------------------------------------------------
plot_child_mortality = function() {
  
  message("  > Plotting child mortality rates")
  
  # ---- Figure properties ----
  
  metric_dict = list(
    avert = "Deaths of children under 5 averted", 
    rate  = "Mortality rate for children under 5", 
    exp   = "Global life expectancy", 
    cov   = "Global vaccine coverage")
  
  # Dictionary for each vaccine case
  case_dict = list( 
    vaccine    = "Vaccination as obserevd", 
    no_vaccine = "No historical vaccination")
  
  # ---- Construct plotting datatables ----
  
  mortality_dt = mortality_rates(grouping = "none")
  
  # Cumulative deaths averted
  averted_dt = mortality_dt %>%
    mutate(total_averted = cumsum(averted)) %>%
    mutate(metric = "avert", 
           case   = "vaccine") %>%
    select(metric, case, year, value = total_averted)
  
  # Mortality rate over time
  rate_dt = mortality_dt %>%
    select(year, vaccine, no_vaccine) %>%
    pivot_longer(cols = -year, 
                 names_to = "case") %>%
    mutate(metric = "rate") %>%
    select(metric, case, year, value) %>%
    as.data.table()
  
  # Life expectancy over time
  life_exp_dt = NULL
  
  # Vaccine coverage over time - global average
  coverage_dt = table("coverage_global") %>%
    filter(coverage > 0, 
           !grepl("_BD$", vaccine)) %>%
    # Combine DTP3 estimates...
    mutate(vaccine = ifelse(
      test = vaccine %in% qc(Dip3, Tet3, aPer3), 
      yes  = "DTP3", 
      no   = vaccine)) %>%
    group_by(vaccine, year) %>%
    summarise(coverage = mean(coverage)) %>%
    ungroup() %>%
    # Append descriptive vaccine names...
    left_join(y  = table("vaccine"), 
              by = "vaccine") %>%
    replace_na(list(vaccine_name = "DTP third dose")) %>%
    mutate(metric = "cov", 
           case   = "vaccine") %>%
    # Tidy up...
    select(metric, case, vaccine_name, 
           year, value = coverage) %>%
    as.data.table()
  
  # Combine into single plotting datatable
  lines_dt = averted_dt %>%
    bind_rows(rate_dt) %>%
    bind_rows(life_exp_dt) %>%
    bind_rows(coverage_dt) %>%
    replace_na(list(vaccine_name = "-")) %>%
    mutate(metric = recode(metric, !!!metric_dict),
           metric = factor(metric, metric_dict), 
           case   = recode(case,   !!!case_dict), 
           case   = factor(case,   case_dict))
  
  area_dt = averted_dt %>%
    rbind(rate_dt) %>%
    rbind(life_exp_dt) %>%
    pivot_wider(names_from = case) %>%
    replace_na(list(no_vaccine = 0)) %>%
    mutate(y_min = pmin(vaccine, no_vaccine), 
           y_max = pmax(vaccine, no_vaccine)) %>%
    mutate(metric = recode(metric, !!!metric_dict),
           metric = factor(metric, metric_dict)) %>%
    select(metric, year, y_min, y_max) %>%
    as.data.table()
  
  # ---- Set y axis properties for each facet ----
  
  fct_fn = function(x) factor(metric_dict[[x]], metric_dict)
  
  y1 = scale_y_continuous(
    labels = comma, 
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05)), 
    breaks = pretty_breaks())
  
  y2 = scale_y_continuous(
    labels = percent, 
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05)), 
    breaks = pretty_breaks())

  y3 = scale_y_continuous(
    labels = percent, 
    limits = c(0, 1),
    expand = c(0, 0), 
    breaks = pretty_breaks())

  y_list = list(
    metric == fct_fn("avert") ~ y1,
    metric == fct_fn("rate")  ~ y2, 
    metric == fct_fn("exp")   ~ y1, 
    metric == fct_fn("cov")   ~ y3)
  
  n_cols  = n_unique(coverage_dt$vaccine_name)
  colours = colour_scheme("viridis::viridis", n = n_cols)
  
  # ---- Produce plot ----
  
  # Plot all metrics over time
  g = ggplot(area_dt) +
    aes(x = year) +
    # Plot area...
    geom_ribbon(
      mapping = aes(
        ymin = y_min, 
        ymax = y_max), 
      fill  = "grey", 
      alpha = 0.7) +
    # Plot lines...
    geom_line(
      data    = lines_dt,
      mapping = aes(
        y = value,
        linetype = case, 
        colour   = vaccine_name), 
      linewidth = 1.1) +
    # Facet by metric...
    facet_wrap(
      facets = vars(metric), 
      scales = "free_y", 
      ncol   = 1) +
    # Set colour scheme...
    scale_color_manual(
      values = c("black", colours)) + 
    # Prettify y axis...
    facetted_pos_scales(y = y_list) + 
    # Prettiy x axis...
    scale_x_continuous(
      limits = c(min(o$years), max(o$years)), 
      expand = expansion(mult = c(0, 0)), 
      breaks = seq(min(o$years), max(o$years), by = 5)) +
    # Prettify legend...
    guides(colour = "none")
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title    = element_blank(),
          axis.text     = element_text(size = 10),
          axis.text.x   = element_text(hjust = 1, angle = 50), 
          axis.line     = element_blank(),
          strip.text    = element_text(size = 16),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"),
          legend.title  = element_blank(),
          legend.text   = element_text(size = 14),
          legend.position = "right", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save figure to file
  save_fig(g, "Child mortality rates", dir = "historical_impact")
  
  # Also save as main manuscript figure
  save_fig(g, "Figure 4a", dir = "manuscript")
}

# ---------------------------------------------------------
# Plot regional differences in child mortality changes
# ---------------------------------------------------------
plot_mortality_change = function() {
  
  message("  > Plotting regional changes in child mortality")
  
  # Description of metric and year scope
  metric_str = "under 5 mortality rate"
  years_str  = paste0("(", paste(range(o$years), collapse = "-"), ")")
  
  # Dictionary for metric type
  type_dict = list(
    diff    = paste("Absolute decrease in", metric_str, years_str),
    rel     = paste("Relative decrease in", metric_str, years_str), 
    contrib = paste("Contribution of vaccination to decrease in", 
                    metric_str, years_str))
  
  # Metric type colour scheme (+1 for global average)
  colours = c("#EBAA2D", "#DF721F", "#71C2A9", "#808080")
  
  # Mortality rates - by region and global average
  region_dt = mortality_rates(grouping = "region")
  world_dt  = mortality_rates(grouping = "none") %>%
    mutate(group = "World")
  
  # Contribution of vaccination to decrease in child mortality
  contrib_dt = rbind(region_dt, world_dt) %>%
    select(group, year, 
           end0 = no_vaccine, 
           end  = vaccine) %>%
    # Start (year 1) and end (year n) values...
    group_by(group) %>%
    mutate(start = end0[year == min(year)],
           diff0 = start - end0,
           diff  = start - end) %>%
    ungroup() %>%
    # We're interested in the value come the final year...
    filter(year == max(year)) %>%
    select(-year) %>%
    # Relative decrease and contribution of vaccination...
    mutate(rel     = diff / start, 
           contrib = 1 - diff0 / diff) %>%
    select(-end0, -diff0) %>%
    as.data.table()
  
  # Order regions by absolute decrease in mortality rates
  order = contrib_dt %>%
    filter(group != "World") %>%
    arrange(desc(diff)) %>%
    pull(group) %>%
    c("World")
  
  # Labels bars for clarity
  label_dt = contrib_dt %>%
    # Convert to percentage...
    mutate(across(
      .cols = c(diff, start, end), 
      .fns  = ~ . * 100)) %>%
    # Construct string explaining start and end values...
    mutate(label = sprintf(
      "%.1f%%\n(%.1f%% to %.1f%%)", 
      diff, start, end)) %>%
    # Only required for absolute bars...
    mutate(type = "diff") %>%
    select(group, type, label)
  
  # Construct plotting datatable
  plot_dt = contrib_dt %>%
    # Retain only what we want to plot...
    select(group, diff, rel, contrib) %>%
    pivot_longer(cols = -group, 
                 names_to = "type") %>%
    # Labels: simple and expansive versions...
    mutate(str = sprintf("%.0f%%", value * 100)) %>%
    left_join(y  = label_dt, 
              by = c("group", "type")) %>%
    mutate(label = ifelse(is.na(label), str, label)) %>%
    # Position according to metric type...
    group_by(type) %>%
    mutate(nudge   = max(value) * 0.05, 
           lab_pos = value + nudge) %>%
    ungroup() %>%
    select(-str, -nudge) %>%
    # Define colour groups (world is special case)...
    mutate(colour = match(type, names(type_dict)), 
           colour = ifelse(
             test = group == "World", 
             yes  = max(colour) + 1, 
             no   = colour), 
           colour = as.character(colour)) %>%
    # Set factors for desired ordering...
    mutate(type  = recode(type, !!!type_dict), 
           type  = factor(type, type_dict), 
           group = factor(group, order)) %>%
    as.data.table()
  
  # Plot regional results as bars
  g = ggplot(plot_dt) +
    aes(x = group, 
        y = value,  
        fill = colour) +
    # Plot bars...
    geom_bar(
      stat  = "identity", 
      width = 0.6) +
    # Add descriptive text...
    geom_text(
      mapping = aes(
        label = label, 
        y     = lab_pos), 
      vjust   = 0, 
      size    = 3.5) + 
    # Facet by type...
    facet_wrap(
      facets = vars(type), 
      ncol   = 1, 
      scales = "free") + 
    # Set colours...
    scale_fill_manual(values = colours) +
    # Prettify y axis...
    scale_y_continuous(
      labels = percent, 
      limits = c(0, NA), 
      expand = expansion(mult = c(0, 0.2)), 
      breaks = pretty_breaks())
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title    = element_blank(),
          axis.text.x   = element_text(size = 14),
          axis.text.y   = element_text(size = 10),
          axis.line     = element_blank(),
          strip.text    = element_text(size = 20),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(2, "lines"),
          panel.grid.major.y = element_line(linewidth = 0.25),
          legend.position = "none")
  
  # Save figure to file
  save_fig(g, "Child mortality change by region", 
           dir = "historical_impact")
  
  # Also save as main manuscript figure
  save_fig(g, "Figure 4b", dir = "manuscript")
}

# ---------------------------------------------------------
# xxxxxxxxxx
# ---------------------------------------------------------
plot_prob_death_age = function() {
  
  message("  > Plotting probability of death by age")
  
  age_bins = 4
  
  # Dictionary for each vaccine case
  #
  # NOTE: Using reverse order here for prettier plot
  case_dict = list( 
    no_vaccine = "No historical vaccination",
    vaccine    = "Vaccination as obserevd")
  
  # Estimated child deaths averted by vaccination
  averted_dt = read_rds("results", "deaths_averted") %>%
    filter(year == max(o$years)) %>%
    # Append region...
    left_join(y  = table("country"), 
              by = "country") %>%
    # xxx ...
    group_by(region) %>%
    summarise(total_averted = sum(impact)) %>%
    ungroup() %>%
    # xxx ...
    expand_grid(impact_age_multiplier()) %>%
    mutate(age     = age + 1, 
           averted = total_averted * scaler) %>%
    select(region, age, averted) %>%
    as.data.table()
  
  # Probability of death by age
  plot_dt = table("wpp_pop") %>%
    filter(year == max(o$years)) %>%
    left_join(y  = table("wpp_deaths"), 
              by = c("country", "year", "age")) %>%
    # xxx ...
    mutate(age = age + 1) %>%
    filter(age %in% 2 ^ (0 : age_bins)) %>%
    # Append region...
    left_join(y  = table("country"), 
              by = "country") %>%
    # Average prob of death by region and year...
    group_by(region, age) %>%
    summarise(pop    = sum(pop), 
              deaths = sum(deaths)) %>%
    ungroup() %>%
    # Append deaths averted...
    left_join(y  = averted_dt, 
              by = c("region", "age")) %>%
    replace_na(list(averted = 0)) %>%
    mutate(vaccine    = deaths / pop, 
           no_vaccine = (deaths + averted) / pop) %>%
    # Melt to long format ready for plotting...
    select(region, age, vaccine, no_vaccine) %>%
    pivot_longer(cols = c(vaccine, no_vaccine), 
                 names_to = "case") %>%
    # Set regional order...
    arrange(desc(value)) %>%
    mutate(region = fct_inorder(region)) %>%
    # Set vaccine case order...
    mutate(case = recode(case, !!!case_dict), 
           case = factor(case, case_dict)) %>%
    arrange(region, case, age) %>%
    as.data.table()
  
  g = ggplot(plot_dt) +
    aes(x = age, 
        y = value, 
        fill = case) +
    geom_col(position = "dodge") +
    # Facet by region...
    facet_wrap(
      facets = vars(region), 
      scales = "free_y") +
    # Set colour scheme...
    scale_fill_manual(
      values = o$palette_who) + 
    # Prettify x axis...
    scale_x_continuous(
      name   = "Age (log2 scale)",
      trans  = "log2", 
      breaks = 2 ^ (0 : age_bins)) +
    # Prettify y axis...
    scale_y_continuous(
      name   = paste("Probability of death before", 
                     "n^th birthday in", max(o$years)),
      labels = percent, 
      limits = c(0, NA), 
      expand = expansion(mult = c(0, 0.05)), 
      breaks = pretty_breaks())
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title   = element_text(size = 20),
          axis.title.x = element_text(
            margin = margin(t = 20, b = 10)),
          axis.title.y = element_text(
            margin = margin(l = 10, r = 20)),
          axis.text    = element_text(size = 10),
          axis.ticks   = element_blank(), 
          axis.line    = element_line(linewidth = 0.25),
          strip.text   = element_text(size = 18),
          strip.background = element_blank(), 
          panel.spacing = unit(1, "lines"), 
          panel.grid.major.y = element_line(linewidth = 0.25),
          legend.title = element_blank(),
          legend.text  = element_text(size = 14),
          legend.key   = element_blank(),
          legend.position = "bottom", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save figure to file
  save_fig(g, "Probability of death by age", dir = "historical_impact")
  
  # Also save as main manuscript figure
  save_fig(g, "Figure 7", dir = "manuscript")
}

# ---------------------------------------------------------
# xxxxxxxxxx
# ---------------------------------------------------------
plot_survival_increase = function() {
  
  message("  > Plotting increase in childhood survival")
  
  age_bound = 50
  
  title = "Historical vaccination compared to hypothetical no vaccination"
  
  # Estimated child deaths averted by vaccination
  averted_dt = read_rds("results", "deaths_averted") %>%
    # ...
    filter(year == max(o$years)) %>%
    mutate(year = as.character(year)) %>%
    # Append region...
    left_join(y  = table("country"), 
              by = "country") %>%
    # ...
    group_by(region) %>%
    summarise(total_averted = sum(impact)) %>%
    ungroup() %>%
    # ...
    expand_grid(impact_age_multiplier()) %>%
    mutate(averted = total_averted * scaler) %>%
    select(region, age, averted) %>%
    as.data.table()
  
  survival_dt = table("wpp_pop") %>%
    left_join(y  = table("wpp_deaths"), 
              by = c("country", "year", "age")) %>%
    # xxx ...
    filter(year == max(o$years), 
           age  <= age_bound) %>%
    # Append region...
    left_join(y  = table("country"), 
              by = "country") %>%
    # Average prob of death by region and year...
    group_by(region, age) %>%
    summarise(pop    = sum(pop), 
              deaths = sum(deaths)) %>%
    ungroup() %>%
    # xxx ...
    left_join(y  = averted_dt, 
              by = c("region", "age")) %>%
    # xxx...
    group_by(region) %>%
    mutate(c_pop   = cumsum(pop), 
           c_death = cumsum(deaths), 
           c_avert = cumsum(averted)) %>%
    ungroup() %>%
    # xxx ...
    mutate(vaccine    = c_death / c_pop, 
           no_vaccine = (c_death + c_avert) / c_pop, 
           relative   = 1 - vaccine / no_vaccine) %>%
    select(region, age, pop, value = relative) %>%
    as.data.table()
  
  world_dt = survival_dt %>%
    group_by(age) %>%
    mutate(weight = pop / sum(pop)) %>%
    summarise(value = sum(value * weight)) %>%
    ungroup() %>%
    mutate(region = "World", 
           global = 1) %>%
    bind_rows(survival_dt) %>%
    select(-pop) %>%
    replace_na(list(global = 0)) %>%
    arrange(global, region) %>%
    mutate(region = fct_inorder(region)) %>%
    as.data.table()
  
  fvps_dt = table("coverage") %>%
    left_join(y  = table("country"), 
              by = "country") %>%
    group_by(age) %>%
    summarise(fvps = sum(fvps)) %>%
    ungroup() %>%
    mutate(value = fvps / sum(fvps)) %>%
    filter(age <= age_bound) %>%
    as.data.table()
  
  # ---- Produce plot ----
  
  colours = c(get_palette("region"), "#000000")
  
  g = ggplot(fvps_dt) +
    aes(x = age, 
        y = value) +
    geom_col(alpha = 0.4) +
    geom_line(
      data    = world_dt, 
      mapping = aes(
        colour    = region,
        linetype  = fct_rev(factor(global)),
        linewidth = global)) +
    # Set colour scheme
    scale_colour_manual(values = colours) +
    scale_linewidth(range = c(1.2, 3)) +
    # Add a figure title
    ggtitle(title) + 
    # Prettiy x axis...
    scale_x_continuous(
      name   = "Age in years", 
      expand = expansion(add = 1),
      breaks = seq(0, age_bound, by = 5)) +
    # Prettify y axis...
    scale_y_continuous(
      name   = "Relative increase in survival probability", 
      labels = percent, 
      expand = expansion(mult = c(0, 0.05)), 
      breaks = pretty_breaks(), 
      sec.axis = sec_axis(
        trans = ~ .,
        name  = "Age at vaccination (all vaccines)",
        labels = percent, 
        breaks = pretty_breaks())) +
    # Prettify legend...
    guides(linetype  = "none", 
           linewidth = "none", 
           colour = guide_legend(
             nrow = 1, 
             override.aes = list(linewidth = 2)))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title   = element_text(size = 20),
          axis.title.x = element_text(
            margin = margin(t = 20, b = 10)),
          axis.title.y.left = element_text(
            margin = margin(l = 10, r = 20)),
          axis.title.y.right = element_text(
            margin = margin(l = 20, r = 10), color = "grey"),
          axis.text = element_text(size = 10),
          axis.text.y.right = element_text(color = "grey"),
          axis.line  = element_blank(),
          strip.text = element_text(size = 16),
          strip.background = element_blank(), 
          plot.title    = element_text(
            margin = margin(t = 10, b = 20), 
            size   = 24,
            hjust  = 0.5), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"),
          legend.title  = element_blank(),
          legend.text   = element_text(size = 14),
          legend.position = "bottom", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save figure to file
  save_fig(g, "Increase in survival", dir = "historical_impact")
  
  # Also save as main manuscript figure
  save_fig(g, "Figure 6", dir = "manuscript")
}

# ---------------------------------------------------------
# Plot measles deaths in context of all cause deaths
# ---------------------------------------------------------
plot_measles_in_context = function() {
  
  message("  > Plotting measles in context")
  
  # Upper age bound for estimates
  age_bound = 5
  
  # Metrics to plot
  metric_dict = list(
    cum  = "Cumulative number of deaths in children under %i",
    abs  = "Annual number of deaths in children under %i",
    norm = "Normalised cause of death for children under %i",
    cov  = "Number of vaccine doses delivered (%i-%i)")
  
  # Cause of death
  cause_dict = list(
    all_cause = "All cause deaths",
    measles   = "Measles as cause of death", 
    other     = "Other cause of death")
  
  # Define vaccine scenarios
  case_dict = list( 
    vaccine    = "Vaccination as obserevd", 
    no_vaccine = "No historical vaccination")
  
  # Name of d-v-a to plot
  measles_id = "xxxMeasles"
  
  # ID of coverage values to plot
  coverage_ids = c(4, 5) # 6
  
  # Colour schemes
  colours = list(
    area  = c("#B0B0B0", "#0189C5"), 
    line  = c("#000000", "#0189C5"), 
    doses = c("#C70A0A", "#901616", "#400808"))
  
  # ---- Deaths by cause ----
  
  # Death data - WPP and external measles models
  deaths = list(
    wpp    = table("wpp_deaths"), 
    extern = table("extern_deaths"))
  
  # Measles deaths with and without vaccine
  measles_deaths = deaths$extern %>%
    format_d_v_a_name() %>%
    # Pathogen and age group of interest...
    filter(d_v_a_name == measles_id, 
           age        <= age_bound) %>%
    pivot_longer(cols = c(vaccine, no_vaccine), 
                 names_to = "case") %>%
    # Summarise over age bins of interest...
    group_by(case, year) %>%
    summarise(measles = sum(value)) %>%
    ungroup() %>%
    as.data.table()
  
  # All cause deaths from WPP
  all_deaths = deaths$wpp %>%
    # Age group of interest...
    filter(age <= age_bound) %>%
    group_by(year) %>%
    summarise(deaths = sum(deaths)) %>%
    ungroup() %>%
    # Repeat for each vaccine scenario...
    expand_grid(case = names(case_dict)) %>%
    select(case, year, deaths) %>%
    # Append number of deaths averted from vaccination...
    left_join(y  = measles_deaths, 
              by = c("case", "year")) %>%
    group_by(year) %>%
    mutate(averted = measles - measles[case == "vaccine"]) %>%
    ungroup() %>%
    # Increment all cause deaths in scenario of no vaccination...
    mutate(all_cause = deaths + averted) %>%
    select(case, year, all_cause) %>%
    arrange(case, year) %>%
    as.data.table()
  
  # Combine causes of death for measles in context
  context_dt = all_deaths %>%
    left_join(y  = measles_deaths, 
              by = c("case", "year")) %>%
    mutate(other = all_cause - measles) %>% 
    # Melt cause of death to long format...
    pivot_longer(cols = names(cause_dict), 
                 names_to  = "cause", 
                 values_to = "abs") %>%
    # Calculate normalised cause of death...
    group_by(case, year) %>%
    mutate(norm = abs / abs[cause == "all_cause"]) %>%
    ungroup() %>%
    # Calculate cumulative values...
    group_by(case, cause) %>%
    mutate(cum = cumsum(abs)) %>%
    ungroup() %>%
    # Convert to tidy format...
    pivot_longer(cols = any_of(names(metric_dict)), 
                 names_to = "metric") %>%
    # Set appropriate order...
    select(case, metric, cause, year, value) %>%
    arrange(case, metric, cause, year) %>%
    as.data.table()
  
  # ---- Construct plotting datatables ----
  
  # First update metric dictionary with values
  metric_dict = list(
    cum  = sprintf(metric_dict$cum,  age_bound),
    abs  = sprintf(metric_dict$abs,  age_bound),
    norm = sprintf(metric_dict$norm, age_bound),
    cov  = sprintf(metric_dict$cov,  min(o$years), max(o$years)))

  # Lines for measles and all cause
  line_dt = context_dt %>%
    filter(cause != "other") %>%
    # Don't plot normalised all cause as always one...
    filter(!(metric == "norm" & 
               cause == "all_cause")) %>%
    # Recode variables and set ordering...
    mutate(case   = recode(case,   !!!case_dict), 
           case   = factor(case,   case_dict), 
           metric = recode(metric, !!!metric_dict), 
           metric = factor(metric, metric_dict), 
           cause  = recode(cause,  !!!cause_dict),
           cause  = factor(cause,  cause_dict))
  
  # Area for measles and other cause in vaccine case only
  area_dt = context_dt %>%
    filter(cause != "all_cause", 
           case  == "vaccine") %>%
    select(-case) %>%
    # Recode variables and set ordering...
    mutate(metric = recode(metric, !!!metric_dict), 
           metric = factor(metric, metric_dict), 
           cause  = recode(cause,  !!!cause_dict),
           cause  = factor(cause,  cause_dict))
  
  # Gloabl measles coverage over time
  coverage_dt = table("coverage_global") %>%
    select(-vaccine) %>%
    filter(v_a_id %in% coverage_ids) %>%
    # Append metric details...
    append_v_a_name() %>%
    mutate(metric = metric_dict$cov, 
           metric = factor(metric, metric_dict)) %>%
    # Tidy up...
    select(metric, v_a_name, year, value = coverage) %>%
    as.data.table()
  
  # ---- Produce plot ----
  
  # Plot measles in context of all cause deaths
  g = ggplot(area_dt) +
    aes(x = year, 
        y = value) +
    # Plot cause of death...
    geom_area(
      mapping = aes(fill = fct_rev(cause)),
      alpha   = 0.5,
      show.legend = FALSE) +
    # Plot all cause deaths...
    geom_line(
      data    = line_dt,
      mapping = aes(
        colour = cause, 
        linetype = case), 
      linewidth = 1.1) +
    # Set primary colour scheme...
    scale_colour_manual(
      name   = "Cause of death", 
      values = colours$line) + 
    scale_fill_manual(
      values = colours$area) +
    # Prettify primary legends...
    guides(colour = guide_legend(
      order   = 2, 
      reverse = TRUE)) + 
    # Plot number of doses...
    ggnewscale::new_scale_color() + 
    geom_line(
      data    = coverage_dt,
      mapping = aes(colour = v_a_name)) +
    # Set doses colour scheme...
    scale_colour_manual(
      name   = "Vaccine", 
      values = colours$doses) + 
    # Facet by metric...
    facet_wrap(
      facets = vars(metric),
      scales = "free_y",
      ncol   = 1) +
    # Prettiy x axis...
    scale_x_continuous(
      limits = c(min(o$years), max(o$years)), 
      expand = expansion(mult = c(0, 0)), 
      breaks = seq(min(o$years), max(o$years), by = 5)) +
    # Prettify y axis...
    expand_limits(y = c(0, 1)) + 
    scale_y_continuous(
      labels = comma, 
      expand = expansion(mult = c(0, 0)), 
      breaks = pretty_breaks()) +
    # Prettify primary legends...
    guides(
      colour    = guide_legend(
        order   = 3), 
      linetype  = guide_legend(
        order   = 1, 
        title   = "Vaccine scenario"))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title    = element_blank(),
          axis.text     = element_text(size = 10),
          axis.text.x   = element_text(hjust = 1, angle = 50), 
          axis.line     = element_blank(),
          strip.text    = element_text(size = 16),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"),
          legend.title  = element_text(size = 14),
          legend.text   = element_text(size = 12),
          legend.margin = margin(t = 40, b = 40, l = 10),
          legend.position = "right", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save figure to file
  save_fig(g, "Measles in context", dir = "historical_impact")
  
  # Also save as main manuscript figure
  save_fig(g, "Figure 5", dir = "manuscript")
}

# ---------------------------------------------------------
# Helen's exploratory figures
# ---------------------------------------------------------
plot_exploratory = function() {
  
  browser() # To be fully integrated...
  
  # Plot estimates of regression predictors by region
  plot_model_fit_df = model_fit %>%
    select(country, region_short, term, estimate) %>%
    pivot_wider(names_from = term,
                values_from = estimate) %>%
    group_by(region_short) %>%
    as.data.table()
  
  # Plot Spearman rank correlation between coefficients of regression
  plot_model_fit_df %>% 
    select(`log(coverage)`,
           `log(coverage_minus_1)`,
           `log(coverage_minus_2)`,
           `log(coverage_minus_3)`,
           `log(coverage_minus_4)`,
           gini,
           HDI,
           `(Intercept)`,
           attended_births) %>%
    ggpairs()
  
  ggpairs(plot_model_fit_df, columns = c(3,5,6), ggplot2::aes(colour=region_short))
  
  # Plot density of coefficients of predictors by region
  ggpairs(plot_model_fit_df,
          columns = c(3,6:10),
          aes(colour = region_short,
              alpha = 0.5))
  
  # Plot data vs fitted for a single country
  plot_df = augment(model_1) %>%
    filter(country == "THA")
  
  ggplot(data = plot_df, aes(x = target, y = .fitted)) +
    geom_point() +
    labs(
      y = "Fitted (predicted values)",
      x = "Data (actual values)",
      title = paste("Vaccine impact of", d_v_a_name, "in", plot_df$country)
    ) +
    geom_abline(intercept = 0, slope = 1)
  
  # Plot data vs. fitted for all countries
  plot_df = augment(best_model)
  
  ggplot(data = plot_df, aes(x = target, y = .fitted)) +
    geom_point() +
    labs(
      y = "Fitted (predicted values)",
      x = "Data (actual values)",
      title = paste("Vaccine impact of", d_v_a_name)
    ) +
    geom_abline(intercept = 0, slope = 1) +
    facet_wrap(~country, ncol = 21)
  
  # Plot model fit for a single country
  plot_df = augment(best_model) %>%
    filter(country == "AGO")
  
  ggplot(data = plot_df, aes(x = year)) +
    geom_point(aes(y = target, colour = "Data")) +
    geom_line(aes(y = .fitted, colour = "Fitted")) +
    labs(y = NULL,
         title = paste("Vaccine impact of", d_v_a_name, "in", plot_df$country)
    ) +
    scale_colour_manual(values=c(Data="black",Fitted="#D55E00")) +
    guides(colour = guide_legend(title = NULL))
  
  # Plot model fit for all countries
  plot_df = augment(best_model)
  
  ggplot(data = plot_df, aes(x = year)) +
    geom_point(aes(y = target, colour = "Data")) +
    geom_line(aes(y = .fitted, colour = "Fitted")) +
    labs(y = NULL,
         title = paste("Vaccine impact of", d_v_a_name)
    ) +
    scale_colour_manual(values=c(Data="black",Fitted="#D55E00")) +
    guides(colour = guide_legend(title = NULL))  +
    facet_wrap(~country, ncol = 21)
  
  # Manually explore associations between predictor variables for different geographical regions and time points
  explore_dt =  data_dt %>% as.data.table() %>% # Transform to data table to remove country as categorical variable
    filter(#year > 2000 & year <= 2020 &
      region_short == "AFR" &
        target > 2e-20) %>%
    select(-country) %>%
    select(target, gini, health_spending, coverage, 
           coverage_minus_1, coverage_minus_2, coverage_minus_3, sdi)
  
  explore_dt %>% 
    ggpairs(upper = list(continuous = wrap("cor", method = "spearman"))) # Use Spearman rank correlation to account for outliers
  
  # Explore model selection by region
  ggplot(data = model_choice, aes(x = model_number)) +
    geom_histogram() +
    facet_wrap(~region_short)
  
  # Explore prob. density of coefficients of predictors
  plot_df = model_fit %>% filter(term == "pop_0to14")
  
  ggplot(data = plot_df, aes(x=estimate) ) +
    geom_density() +
    facet_wrap(~region_short, ncol=1)
}

# ---------------------------------------------------------
# Helen's tornado plot of predictor coefficients
# ---------------------------------------------------------
plot_tornado = function() {
  
  # TODO: Split by decade, facet by region OR d_v_a_id OR predictor
  
  # Explore density of coefficients of predictors
  plot_dt = impute_1$model_fit %>%
    select(-c(country, d_v_a_id.x, d_v_a_id.y, model_number, .model, AICc)) %>%
    filter(#!term == "HDI" &
      !region_short == "NA" &
        p.value <= 0.05) %>%
    mutate(model = region_short)
  
  g = dwplot(plot_dt) + 
    geom_vline(xintercept = 0, linetype = 2) +
    ggtitle(paste0("Predicting ", d_v_a_name, " vaccine impact"))
  
  g = g + theme_classic()
}

# ---------------------------------------------------------
# Convert v_a_id into human-readable sting
# ---------------------------------------------------------
append_v_a_name = function(id_dt) {
  
  # Vaccine descriptive names
  vaccine_dt = table("vaccine") %>%
    select(vaccine, .v = vaccine_name)
  
  # Append v_a description
  name_dt = id_dt %>%
    left_join(y  = table("v_a"), 
              by = "v_a_id") %>%
    rename(.a = activity) %>%
    # Append descriptive names...
    left_join(y  = vaccine_dt, 
              by = "vaccine") %>%
    select(-vaccine) %>%
    # Construct full names...
    mutate(v_a_name = paste0(.v, ": ", .a), 
           .after = v_a_id) %>%
    # Set plotting order...
    mutate(v_a_name = fct_inorder(v_a_name)) %>%
    select(-.v, -.a)
  
  return(name_dt)
}

# ---------------------------------------------------------
# Full descriptive names for disease, vaccine, or vaccine type
# ---------------------------------------------------------
append_d_v_t_name = function(name_dt) {
  
  # All columns to update
  d_v_t = intersect(
    x = names(name_dt), 
    y = qc(disease, vaccine, type))
  
  # Convert column name type -> vaccine_type
  if ("type" %in% d_v_t)
    name_dt %<>% rename(vaccine_type = type)
  
  # Iterate through columns to update
  for (x in d_v_t) {
    
    # Vaccine type is a special case
    if (x == "type")
      x = "vaccine_type"
    
    # Column name of full description
    x_name = paste1(x, "name")
    
    # Plotting order
    x_ord = table(x)[[x_name]]
    
    # Append full name description
    name_dt %<>%
      left_join(table(x), by = x) %>%
      select(-all_of(x)) %>%
      rename(.x := all_of(x_name)) %>%
      mutate(.x = factor(.x, x_ord)) %>%
      rename(!!x := .x) %>%
      select(all_of(names(name_dt)))
  }
  
  # Convert back vaccine_type -> type
  if ("type" %in% d_v_t)
    name_dt %<>% rename(type = vaccine_type)
  
  return(name_dt)
}

# ---------------------------------------------------------
# Set natural order of d_v_a names - appending if necessary
# ---------------------------------------------------------
format_d_v_a_name = function(name_dt) {
  
  # Natural order of d_v_a_name (not necessarily alphabetical)
  d_v_a_dt = table("d_v_a") %>%
    select(d_v_a_id, d_v_a_name)
  
  # Check if d_v_a_name is already defined
  if (!"d_v_a_name" %in% names(name_dt)) {
    
    # Append name if it's not already in the data
    name_dt %<>%
      left_join(y  = d_v_a_dt, 
                by = "d_v_a_id")
  }
  
  # Set order using factors according to d_v_a_dt
  order_dt = name_dt %>%
    mutate(d_v_a_name = factor(
      x      = d_v_a_name, 
      levels = d_v_a_dt$d_v_a_name)) %>%
    arrange(d_v_a_name)
  
  return(order_dt)
}

# ---------------------------------------------------------
# Get vector of colours for a given variable
# ---------------------------------------------------------
get_palette = function(variable) {
  
  n_colours = list(
    disease = n_unique(table("disease")$disease), 
    region  = n_unique(table("country")$region), 
    income  = n_unique(table("income_status")$income))
  
  if (!variable %in% names(o$palette))
    stop("Input '", variable, "' must be one of: ", 
         paste(names(o$palette), collapse = ", "))
  
  p = o$palette[[variable]]
  n = n_colours[[variable]]
  
  colours = colour_scheme(p, n = n)
  
  return(colours)
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

