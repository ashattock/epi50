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
    gbd    = "Static modelling", 
    impute = "Geographic imputation model", 
    extrap = "Temporal extrapolation model")
  
  # Associated colours
  impact_colours = c("#EB7D5B", "#FED23F", "#B5D33D", "#6CA2EA", "#442288")
  
  # ---- Number of FVPs by pathogen ----
  
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
  
  # GBD approach
  gbd_dt = table("gbd_estimates") %>%
    select(disease, country, year) %>%
    unique() %>%
    mutate(class = "gbd")
  
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
  all_dt = rbind(gbd_dt, vimc_dt, impute_dt, extern_dt) %>%
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
      values = colour_scheme(
        map = "brewer::paired", 
        n   = length(age_groups))) +
    # Prettify legend...
    guides(fill = guide_legend(title = "Age group")) +
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
  save_fig(g, "Vaccine efficacy profiles", dir = "static")
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
    save_fig(g, save_name, by, dir = "static")
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
  save_fig(g, "Deaths averted by disease", dir = "static")
  
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
  save_fig(g, "Deaths averted by vaccine", dir = "static")
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
    impact_cum = "Cumuative deaths averted (in millions)", 
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

