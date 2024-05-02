###########################################################
# PLOTTING
#
# All plotting functionality in one place.
#
###########################################################

# **** Data visualisation ****

# ---------------------------------------------------------
# Plot methodology figure to be used in paper
# ---------------------------------------------------------
plot_scope = function() {
  
  message("  - Plotting country-disease scope")
  
  # Manually set tidy y axis limit
  y_max = 8  # In billions
  
  # Linear time interpolation for less pixilated figures
  smoothness = 10  # Higher value for smoother plot
  
  # Dictionary for full impact source descriptions
  impact_dict = c(
    extern = "Transmission modelling (Form 1)",
    vimc   = "Transmission modelling (Form 2)", 
    static = "Static modelling (Form 3)", 
    impute = "Geographic imputation model", 
    extrap = "Temporal extrapolation model")
  
  # Associated colours
  impact_colours = c("#EB7D5B", "#FED23F", "#B5D33D", "#6CA2EA", "#442288")
  
  # Alter figure dimensions
  save_height = 18
  
  # ---- Number of FVPs by pathogen ----
  
  # Number of FVPs by country
  fvps_dt = table("coverage") %>%
    filter(coverage > 0) %>%
    # Append disease details...
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_id") %>%
    # Concept here is FVP, so remove birth dose & boosters...
    filter(!grepl("_(bd|bx|px)$", vaccine)) %>% 
    # Summarise over age...
    group_by(disease, country, year) %>%
    summarise(fvps = sum(fvps)) %>%
    ungroup() %>%
    as.data.table()
  
  # Total cumulative FVPs per disease
  total_dt = fvps_dt %>%
    # Total by disease...
    group_by(disease, year) %>%
    summarise(total = sum(fvps)) %>%
    ungroup() %>%
    # Cumulative over time...
    group_by(disease) %>%
    mutate(cum_total = cumsum(total)) %>%
    ungroup() %>%
    as.data.table()
  
  # ---- Source of impact estimates ----
  
  # Static model settings
  static_dt = table("gbd_estimates") %>%
    select(disease, country, year) %>%
    unique() %>%
    filter(year %in% o$gbd_estimate_years) %>%
    mutate(class = "static")
  
  # External model settings
  extern_dt = table("extern_estimates") %>%
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_id") %>%
    select(disease, country, year) %>%
    unique() %>%
    mutate(class = "extern")
  
  # VIMC settings
  vimc_dt = table("vimc_estimates") %>%
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_id") %>%
    select(disease, country, year) %>%
    unique() %>%
    mutate(class = "vimc")
  
  # VIMC country imputation settings
  impute_dt = table("d_v_a") %>%
    filter(source == "vimc") %>%
    select(disease) %>%
    unique() %>%
    expand_grid(
      country = all_countries(), 
      year    = unique(vimc_dt$year)) %>%
    left_join(y  = vimc_dt, 
              by = c("disease", "country", "year")) %>%
    filter(is.na(class)) %>%
    mutate(class = "impute") %>%
    as.data.table()
  
  # ---- Construct plotting datatable ----
  
  # Combine all impact sources
  plot_dt = rbind(extern_dt, static_dt, vimc_dt, impute_dt) %>%
    # Append FVPs...
    right_join(y  = fvps_dt, 
               by = c("disease", "country", "year")) %>%
    # Anything not yet specified is time-exrapolated...
    mutate(class = ifelse(is.na(class) & year <= 2000, "extrap.1", class), 
           class = ifelse(is.na(class) & year >  2000, "extrap.2", class)) %>%
    # Summarise over all countries...
    group_by(disease, class, year) %>%
    summarise(fvps = sum(fvps)) %>%
    ungroup() %>%
    # Append cumulative total...
    left_join(y  = total_dt, 
              by = c("disease", "year")) %>%
    # Share of cumulative FVPs by class over time...
    mutate(share = fvps / total, 
           value = share * cum_total / 1e9) %>%
    # Expand for a more granular timescale...
    expand_grid(time = seq(
      from = 0 , 
      to   = 1 - 1 / smoothness, 
      by   = 1 / smoothness)) %>%
    mutate(time = year + time) %>%
    filter(time <= max(o$years)) %>%
    # Interpolate annual values for smoother plot...
    mutate(value = ifelse(time == year, value, NA)) %>%
    group_by(disease, class) %>%
    mutate(value = na_interpolation(value)) %>%
    ungroup() %>%
    # Use full disease names...
    left_join(y  = table("disease_name"), 
              by = "disease") %>%
    mutate(class = str_remove(class, "\\.[1-9]+$")) %>%
    select(disease = disease_name, class, time, value) %>%
    # Use impact source descriptions...
    mutate(class = recode(class, !!!impact_dict), 
           class = factor(class, rev(impact_dict))) %>%
    # Set disease order...
    arrange(desc(value)) %>%
    mutate(disease = fct_inorder(disease)) %>%
    # Tidy up...
    arrange(disease, class, time) %>%
    as.data.table()
  
  # ---- Construct label datatable ----
  
  # Year range of analysis
  year1 = min(o$years)
  year2 = max(o$years)
  
  # Label description string
  label = paste0("Total (", year1, "-", year2, "):")
  
  # Construct labels: total FVPs over analysis timeframe
  label_dt = plot_dt %>%
    filter(time == year2) %>%
    # Total cumulative FVPs by disease...
    group_by(disease) %>%
    summarise(total = sum(value)) %>%
    ungroup() %>%
    # Construct labels...
    mutate(total = round(total, 2),
           label = paste(label, total, "billion")) %>%
    # Set coordinates...
    mutate(time  = year1 + 0.01 * (year2 - year1),
           value = y_max * 0.9) %>%
    as.data.table()
  
  # ---- Produce plot ----
  
  # Plot FVP over time per pathogen and impact source
  g = ggplot(plot_dt) +
    aes(x = time, 
        y = value) +
    geom_bar(
      mapping  = aes(fill = class), 
      stat     = "identity", 
      position = "stack",  
      width    = 1 / smoothness) +
    # Add total labels...
    geom_text(
      data    = label_dt, 
      mapping = aes(label = label), 
      size    = 4.5,
      hjust   = 0, 
      vjust   = 1) + 
    # Some intricate faceting...
    facet_rep_wrap(
      facets = vars(disease), 
      ncol   = 1, 
      labeller = label_wrap_gen(width = 20), 
      strip.position = "right", 
      repeat.tick.labels = FALSE) + 
    # Set colours and legend title...
    scale_fill_manual(
      values = impact_colours, 
      name   = "Source of impact estimates") +
    guides(fill = guide_legend(
      reverse = TRUE, 
      byrow   = TRUE, 
      nrow    = 2)) +
    # Prettify x axis...
    scale_x_continuous(
      limits = c(year1 - 1 / smoothness, 
                 year2 + 1 / smoothness), 
      expand = expansion(mult = c(0, 0)), 
      breaks = seq(year1, year2, by = 5)) +  
    # Prettify y axis...
    scale_y_continuous(
      name   = paste("Cumulative number of people receiving",
                     "final primary dose (in billions)"), 
      limits = c(0, y_max),
      breaks = seq(2, y_max, by = 2),
      labels = comma,
      expand = expansion(mult = c(0, 0)))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(strip.text    = element_text(size = 14),
          strip.text.y  = element_text(angle = 0, hjust = 0),
          strip.background = element_blank(), 
          axis.title.x  = element_blank(),
          axis.title.y  = element_text(
            size = 18, margin = margin(l = 10, r = 20)),
          axis.text.x   = element_text(size = 12),
          axis.text.y   = element_text(size = 10),
          axis.ticks    = element_blank(), 
          axis.line     = element_line(linewidth = 0.25),
          panel.spacing = unit(-0.6, "lines"), 
          panel.grid.major.y = element_line(linewidth = 0.25),
          legend.title  = element_text(size = 14),
          legend.text   = element_text(size = 12),
          legend.key    = element_blank(),
          legend.position = "bottom", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save to file
  save_fig(g, "S06", height = save_height)
}

# ---------------------------------------------------------
# Plot total number of FVP over time for each d_v_a
# ---------------------------------------------------------
plot_total_fvps = function() {
  
  message("  - Plotting total number of FVP")
  
  # Flag for whether to plot FVPs cumulatively over time
  cumulative = TRUE
  
  # String to define total FVPs
  total = "All source total"
  
  # Total FVPs (sum of all sources)
  #
  # NOTE: Not necessarily equal to sum of all sources as
  #       SIA are assumed to be only partially targeted
  total_dt = table("coverage") %>%
    # Summarise over countries and age...
    group_by(d_v_a_id, year) %>%
    summarise(fvps = sum(fvps) / 1e9) %>%
    ungroup() %>%
    # Cumulative FVPs...
    group_by(d_v_a_id) %>%
    mutate(fvps_cum = cumsum(fvps)) %>%
    ungroup() %>%
    # Append d-v-a details...
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_id") %>%
    mutate(type = total) %>%
    # Tidy up...
    select(d_v_a_name, disease, type, 
           year, fvps, fvps_cum) %>%
    as.data.table()
  
  # Polio is a special case
  polio_dt = total_dt %>%
    filter(disease == "polio") %>%
    mutate(source = "POLIS") %>%
    select(d_v_a_name, source, year, fvps, fvps_cum)
  
  # Map external model vaccines to d-v-a ID
  extern_map = table("d_v_a_extern") %>%
    filter(disease != "polio") %>%
    select(disease, id = d_v_a_id) %>%
    left_join(y  = table("d_v_a"), 
              by = "disease") %>%
    select(disease, id, d_v_a_id)
  
  # Number of FVPs by source of data
  source_dt = table("coverage_source") %>%
    mutate(source = toupper(source)) %>%
    # Map external vaccines to d-v-a ID...
    rename(id = d_v_a_id) %>%
    left_join(y  = extern_map, 
              by = "id") %>%
    mutate(d_v_a_id = ifelse(
      test = is.na(d_v_a_id), 
      yes  = id, 
      no   = d_v_a_id)) %>%
    filter(d_v_a_id %in% table("d_v_a")$d_v_a_id) %>%
    # Summarise over countries and age...
    group_by(d_v_a_id, source, year) %>%
    summarise(fvps = sum(fvps) / 1e9) %>%
    ungroup() %>%
    # Cumulative FVPs...
    group_by(d_v_a_id, source) %>%
    mutate(fvps_cum = cumsum(fvps)) %>%
    ungroup() %>%
    # Append polio...
    format_d_v_a_name() %>%
    select(all_names(polio_dt)) %>%
    rbind(polio_dt) %>%
    as.data.table()
  
  # All sources being plotted
  sources = unique(source_dt$source)
  
  # Concatenate into plotting datatable
  plot_dt = total_dt %>%
    select(-disease) %>%
    bind_rows(source_dt) %>%
    replace_na(list(
      source = "NA", 
      type   = "NA")) %>%
    mutate(source = factor(source, c(sources, "NA")), 
           type   = factor(type,   c(total,   "NA")))
  
  # Metric to use for y axis
  y = ifelse(cumulative, "fvps_cum", "fvps")
  
  # Colours: named vector
  colours = colour_scheme(
    map = "brewer::set1", 
    n   = length(sources)) %>%
    c("black") %>%
    setNames(c(sources, "NA"))
  
  # Line types: named vector
  types = c("dashed", "solid") %>%
    setNames(c(total, "NA"))
  
  # Plot FVPs over time for each d_v_a
  g = ggplot(plot_dt) + 
    aes(x = year, 
        y = !!sym(y), 
        colour   = source, 
        linetype = type) + 
    geom_line(
      linewidth = 1.5) + 
    # Facet with strip text wrapping...
    facet_wrap(
      facets   = vars(d_v_a_name), 
      labeller = label_wrap_gen(width = 24), 
      scales   = "free_y") + 
    # Set colour scheme...
    scale_color_manual(
      breaks = sources,
      values = colours) +
    # Set line types...
    scale_linetype_manual(
      breaks = total,
      values = types) +
    # Prettify x axis...
    scale_x_continuous(
      limits = range(o$years), 
      expand = expansion(mult = c(0, 0)), 
      breaks = seq(
        from = min(o$years), 
        to   = max(o$years), 
        by   = 10)) +  
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
          legend.position = "bottom", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(3, "lines"))
  
  # Save to file
  save_fig(g, "S07")
}

# ---------------------------------------------------------
# Plot age targets as defined by WIISE and VIMC coverage data
# ---------------------------------------------------------
plot_coverage_age_density = function() {
  
  message("  - Plotting coverage data density by age")
  
  # Plot upto 2^x age
  log2_max = 6
  
  # Construct plotting datatable
  plot_dt = table("coverage_source") %>%
    mutate(trans_age = pmax(age, 1), .after = age) %>%
    filter(trans_age <= 2 ^ log2_max) %>%
    format_d_v_a_name() %>%
    filter(!is.na(d_v_a_name))
  
  # Plot age density of coverage data by source
  g = ggplot(plot_dt) +
    aes(x = trans_age, 
        y = after_stat(scaled), 
        colour = source,
        fill   = source) +
    geom_density(alpha = 0.2) +
    # Facet with strip text wrapping...
    facet_wrap(
      facets   = vars(d_v_a_name), 
      labeller = label_wrap_gen(width = 24)) + 
    # Prettify x axis...
    scale_x_continuous(
      name   = "Age (log2 scale)",
      trans  = "log2", 
      limits = c(1, 2 ^ log2_max), 
      expand = c(0, 0), 
      breaks = 2 ^ (0 : log2_max)) +  
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
  save_fig(g, "S08")
}

# **** Static models ****

# ---------------------------------------------------------
# Plot Global Burden of Disease burden estimates by age
# ---------------------------------------------------------
plot_gbd_estimates = function() {
  
  message("  - Plotting GBD burden estimates by age")
  
  # Flag for plotting recent extrapolation
  plot_extrap = FALSE
  
  # Define meaningful metric names
  metric_dict = c(
    deaths = "Deaths", 
    dalys  = "Disability-adjusted life years")
  
  # Define age groups and associated upper bounds
  age_groups = c(
    "Neonates"      = -1,
    "Other infants" = 0,
    "1-4 years"     = 5,
    "5-14"          = 15,
    "15-49"         = 50,
    "50-69"         = 70,
    "70+ years"     = max(o$ages))
  
  # Map each age bin to respective age group
  age_group_dt = data.table(age = c(-1, o$ages)) %>%
    mutate(group_idx = match(age, age_groups),
           group = names(age_groups[group_idx])) %>%
    fill(group, .direction = "up") %>%
    select(age, age_group = group)
  
  # Plot up to this year
  plot_to = ifelse(plot_extrap, max(o$years), max(o$gbd_estimate_years))
  
  # Load GBD estimates and categorise into age groups
  gbd_dt = table("gbd_estimates") %>%
    filter(year <= plot_to) %>%
    append_d_v_t_name() %>%
    # Append age group details...
    left_join(y  = age_group_dt,
              by = "age") %>%
    # Summarise for these broad age groups...
    group_by(disease, year, age_group) %>%
    summarise(deaths = sum(deaths_disease), 
              dalys  = sum(dalys_disease)) %>%
    ungroup() %>%
    # Melt to long format...
    pivot_longer(cols = c(deaths, dalys), 
                 names_to = "metric") %>%
    replace_na(list(value = 0)) %>%
    # Append metric names...
    mutate(metric = recode(metric, !!!metric_dict), 
           metric = factor(metric, metric_dict)) %>%
    # Set age group factors for meaningful plotting order...
    mutate(age_group = factor(age_group, names(age_groups))) %>%
    select(disease, metric, year, age_group, value) %>%
    arrange(disease, metric, year, age_group) %>%
    as.data.table()
  
  # Plot deaths over time by age group
  g = ggplot(gbd_dt) +
    aes(x = year, 
        y = value, 
        fill = age_group) +
    geom_bar(stat = "identity") +
    # Facet by disease...
    facet_grid2(
      cols   = vars(disease), 
      rows   = vars(metric), 
      scales = "free_y", 
      independent = "y") + 
    # Set colour scheme...
    scale_fill_manual(
      name   = "Age group",
      values = colour_scheme(
        map = "brewer::paired", 
        n   = length(age_groups))) +
    # Prettify y axis...
    scale_y_continuous(
      labels = comma,
      expand = expansion(mult = c(0, 0.05)), 
      breaks = pretty_breaks())
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title    = element_blank(),
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
          legend.position = "bottom", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save to file
  save_fig(g, "S12")
}

# ---------------------------------------------------------
# Plot static model pathogen vaccine efficacy with immunity decay
# ---------------------------------------------------------
plot_vaccine_efficacy = function() {
  
  message("  - Plotting vaccine efficacy profiles")
  
  # Dictionary for each vaccine schedule
  schedule_dict = c(
    x  = "Primary schedule",
    bx = "Booster schedule", 
    px = "Pregnancy schedule")
  
  # Function to extract vaccine efficacy data
  extract_data_fn = function(dt) {
    
    # Extract and format into datatable
    data = dt$data %>%
      unlist() %>%
      as_named_dt("value") %>%
      mutate(name    = c("y", "x"), 
             vaccine = dt$vaccine) %>%
      pivot_wider(id_cols = vaccine) %>%
      as.data.table()
    
    return(data)
  }
  
  # Function to group similar vaccines but split by dose
  schedule_fn = function(dt) {
    
    # Append descriptive columns
    shedule_dt = dt %>%
      # Primary schedule or booster dose...
      separate(col  = vaccine, 
               into = c("type", "schedule"), 
               sep  = "_", 
               fill = "right", 
               remove = FALSE) %>%
      replace_na(list(schedule = "x")) %>%
      mutate(schedule = recode(schedule, !!!schedule_dict), 
             schedule = factor(schedule, schedule_dict)) %>%
      # Append vaccine type name
      mutate(type = str_remove(type, "[0-9]")) %>%
      left_join(y  = table("type_name"), 
                by = "type") %>%
      select(-type) %>%
      rename(type = type_name) %>%
      mutate(type = fct_inorder(type)) %>%
      as.data.table()
    
    return(shedule_dt)
  }
  
  # Extract all vaccine efficacy data points
  data_dt = table("vaccine_efficacy") %>%
    dtapply(extract_data_fn) %>%
    rbindlist() %>%
    inner_join(y  = table("d_v_a"), 
               by = "vaccine") %>%
    schedule_fn() %>%
    select(type, schedule, time = x, value = y)
  
  # Load vaccine efficacy profiles
  profile_dt = table("vaccine_efficacy_profiles") %>%
    schedule_fn() %>%
    select(type, schedule, time, value = profile)
  
  # Plot vaccine efficacy with waning immunity (if any)
  g = ggplot(profile_dt) + 
    aes(x = time, 
        y = value) + 
    # Plot immunity profile...
    geom_line(
      colour    = "black",
      linewidth = 1.5) + 
    # Plot coloured data points...
    geom_point(
      data    = data_dt, 
      mapping = aes(colour = type),
      size    = 2.5) + 
    # Facet by vaccine and schedule...
    facet_grid(
      rows     = vars(type), 
      cols     = vars(schedule),
      labeller = label_wrap_gen(width = 20)) + 
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
    theme(axis.text     = element_text(size = 9),
          axis.title.x  = element_text(
            size = 16, margin = margin(b = 10, t = 20)),
          axis.title.y  = element_text(
            size = 16, margin = margin(l = 10, r = 20)),
          axis.line     = element_blank(),
          strip.text.x  = element_text(size = 13),
          strip.text.y  = element_text(size = 11),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(0.5, "lines"),
          panel.grid.major.y = element_line(linewidth = 0.25),
          legend.position = "none")
  
  # Save figure to file
  save_fig(g, "S09")
}

# ---------------------------------------------------------
# Plot effective coverage with waning immunity for static model pathogens
# ---------------------------------------------------------
plot_effective_coverage = function() {
  
  message("  - Plotting effective coverage by year and age")
  
  # Plot only up to a certain age
  age_max = 50
  
  # Manually define appropriate number of colours 
  colours = list(
    disease = "blues", 
    type    = "greens")
  
  # Number of breaks in continuous colour bar
  col_breaks = 40
  
  # Repeat for disease and vaccine type
  for (by in c("disease", "type")) {
    
    # Load previously calculated total coverage file
    effective_dt = read_rds("static", "effective_coverage", by)
    
    # Population weight over all countries
    plot_dt = effective_dt %>%
      append_d_v_t_name() %>%
      select(country, by = !!by, year, age, coverage) %>%
      filter(age >= 0,  # For simplicity, do not plot neonate immunity  
             age <= age_max) %>%
      # Append population size...
      left_join(y  = table("wpp_pop"),
                by = c("country", "year", "age")) %>%
      mutate(n = pop * coverage) %>%
      # Population weighted coverage...
      lazy_dt() %>%
      group_by(by, year, age) %>%
      summarise(effective_coverage = sum(n / sum(pop))) %>%
      ungroup() %>%
      as.data.table()
    
    # Construct continuous colour vector
    col_map = paste0("pals::brewer.", colours[by])
    col_vec = colour_scheme(col_map, n = col_breaks - 1)
    
    # Plot each pathogen by year ana age
    g = ggplot(plot_dt) + 
      aes(x = year, y = age, fill = effective_coverage) + 
      geom_tile() +
      facet_wrap(~by) + 
      # Set continuous colour bar...
      scale_fill_gradientn(
        colours = c("#FFFFFF", col_vec), 
        limits  = c(0, 1), 
        breaks  = pretty_breaks(), 
        label   = percent,
        guide   = guide_colourbar(
          title = "Effective coverage")) + 
      # Prettify x axis...
      scale_x_continuous(
        expand = c(0, 0), 
        breaks = seq(
          from = min(o$years), 
          to   = max(o$years), 
          by   = 5)) +
      # Prettify y axis...
      scale_y_continuous(
        name   = "Age (in years)", 
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
    if (by == "type")    save_fig(g, "S10")
    if (by == "disease") save_fig(g, "S11")
  }
}

# ---------------------------------------------------------
# Plot deaths and DALYs averted for static model pathogens
# ---------------------------------------------------------
plot_static = function() {
  
  message("  - Plotting static model impact results")
  
  # Flag for plotting recent extrapolation
  plot_extrap = FALSE
  
  # Disease burden / burden averted dictionary
  metric_dict = c(
    burden  = "Estimated disease-specific burden (GBD 2021)", 
    averted = "Estimated burden averted from static model")
  
  # Associated colours
  metric_colours = c("darkred", "navyblue")
  
  # ---- Plot by disease ----
  
  # Ensure consistent years of plotting
  plot_years = table("gbd_estimates") %>%
    lazy_dt() %>%
    pivot_longer(cols = -c(disease, country, year, age), 
                 names_to = "metric") %>%
    group_by(metric, disease, year) %>%
    summarise(value = sum(value)) %>%
    ungroup() %>%
    # Only years for which we have all data...
    filter(value > 0) %>%
    count(year) %>%
    filter(n == max(n)) %>%
    pull(year) %>%
    intersect(o$gbd_estimate_years)
  
  # Repeat for each metric
  for (metric in o$metrics) {
    
    # Load previously calculated total coverage file
    averted_dt = read_rds("static", metric, "averted_disease")
    
    # Summarise results over country and age
    disease_dt = averted_dt %>%
      lazy_dt() %>%
      filter(year %in% plot_years) %>%
      pivot_longer(cols = c(burden, averted), 
                   names_to = "metric") %>%
      # Summarise over countries...
      group_by(disease, year, metric) %>%
      summarise(value = sum(value) / 1e6) %>%
      ungroup() %>%
      # Recode deaths disease/averted...
      left_join(y  = table("disease_name"), 
                by = "disease") %>%
      mutate(metric = recode(metric, !!!metric_dict), 
             metric = factor(metric, metric_dict)) %>%
      # Tidy up...
      select(metric, disease = disease_name, year, value) %>%
      arrange(metric, disease, year) %>%
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
    
    # Set a figure title explaining metric
    g = g + ggtitle(
      label = table("metric_dict") %>%
        filter(metric == !!metric) %>%
        pull(metric_impact))
    
    # Prettify theme
    g = g + theme_classic() + 
      theme(plot.title    = element_text(size = 20, hjust = 0.5),
            axis.title.x  = element_blank(),
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
    
    # Save figures to file
    if (metric == "deaths") save_fig(g, "S13")
    if (metric == "dalys")  save_fig(g, "S14")
  }
}

# **** Regression (impute and infer) ****

# ---------------------------------------------------------
# Plot truth vs predicted for imputation training data
# ---------------------------------------------------------
plot_impute_quality = function(metric) {
  
  message("  - Plotting imputation quality of fit: ", metric)
  
  # ---- Load results from fitting ----
  
  # Function to load imputation results
  load_results_fn = function(id) {
    
    # Load file and extract model details
    result = try_load(
      pth  = o$pth$impute, 
      file = paste1("impute", metric, id), 
      throw_error = FALSE) %>%
      pluck("result")
    
    return(result)
  }
  
  # Load imputation results for all d-v-a
  results_dt = table("d_v_a") %>%
    filter(source == "vimc") %>%
    pull(d_v_a_id) %>%
    lapply(load_results_fn) %>%
    rbindlist()
  
  # ---- Construct plotting datatables ----
  
  # Prepare datatable for plotting
  plot_dt = results_dt %>%
    filter(target > 0, prediction > 0) %>%
    format_d_v_a_name() %>%
    select(d_v_a_name, target, prediction)
  
  # Maximum value in each facet (target or prediction)
  blank_dt = plot_dt %>%
    mutate(max_value = pmax(target, prediction)) %>%
    group_by(d_v_a_name) %>%
    summarise(max_value = max(max_value)) %>%
    ungroup() %>%
    expand_grid(type = c("target", "prediction")) %>%
    pivot_wider(names_from  = type, 
                values_from = max_value) %>%
    as.data.table()
  
  # Colour scheme
  colours = colour_scheme(
    map = "brewer::set1", 
    n   = n_unique(plot_dt$d_v_a_name))
  
  # ---- Produce plot ----
  
  # Single plot with multiple facets
  g = ggplot(plot_dt) +
    aes(x = target, 
        y = prediction, 
        color = d_v_a_name) +
    # Plot truth vs predicted...
    geom_point(
      alpha = 0.5, 
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
    # Set colour scheme...
    scale_colour_manual(
      values = colours) + 
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
    theme(axis.text     = element_text(size = 10),
          axis.text.x   = element_text(hjust = 1, angle = 50),
          axis.title.x  = element_text(
            size = 24, margin = margin(b = 10, t = 20)),
          axis.title.y  = element_text(
            size = 24, margin = margin(l = 10, r = 20)),
          axis.line     = element_blank(),
          strip.text    = element_text(size = 14),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"))
  
  # Save figure to file
  save_sub = letters[which(o$metrics %in% metric)]
  save_fig(g, paste0("S15", save_sub))
}

# --------------------------------------------------------
# Plot predictive performance for each country
# --------------------------------------------------------
plot_impute_perform = function(metric) {
  
  message("  - Plotting predictive performance by country: ", metric)
  
  # Function to load imputation results
  load_results_fn = function(id, pull) {
    
    # Load file and extract model details
    result = try_load(
      pth  = o$pth$impute, 
      file = paste1("impute", metric, id), 
      throw_error = FALSE) %>%
      pluck(pull)
    
    return(result)
  }
  
  # Data used to train regression models with associated fit
  train_dt = table("d_v_a") %>%
    filter(source == "vimc") %>%
    pull(d_v_a_id) %>%
    lapply(load_results_fn, pull = "model") %>%
    lapply(augment) %>%
    rbindlist() %>%
    rename(prediction = .fitted) %>%
    format_d_v_a_name() %>%
    select(d_v_a_name, region, country, 
           year, target, prediction)
  
  # Idnetify outliers for more meaningful plot
  outlier_dt = train_dt %>%
    group_by(d_v_a_name, region) %>%
    slice_max(prediction, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(outlier = target * 2) %>%
    select(d_v_a_name, region, outlier) %>%
    as.data.table()
  
  # Imputed predictions for countries without data
  impute_dt = table("d_v_a") %>%
    filter(source == "vimc") %>%
    pull(d_v_a_id) %>%
    lapply(load_results_fn, pull = "result") %>%
    rbindlist() %>%
    filter(is.na(target) & !is.na(prediction)) %>%
    format_d_v_a_name() %>%
    append_region_name() %>%
    left_join(y  = outlier_dt, 
              by = c("d_v_a_name", "region")) %>%
    filter(prediction <= outlier, 
           year >= 2000) %>%
    select(d_v_a_name, region, country, 
           year, target, prediction)
  
  # Construct colour scheme
  colours = colour_scheme(
    map = "brewer::set1", 
    n   = n_unique(train_dt$d_v_a_name))
  
  # Plot de-identified countries for all diseases
  g = ggplot(train_dt) + 
    aes(x = year, 
        y = prediction,
        group  = country,
        colour = d_v_a_name) + 
    # Plot training data and fit...
    geom_point(
      mapping = aes(y = target), 
      alpha   = 0.5, 
      shape   = 16) + 
    geom_line() + 
    # Plot imputed on top in black...
    geom_line(
      data   = impute_dt, 
      colour = "black") + 
    # Facet by disease and region...
    facet_grid(
      cols   = vars(region), 
      rows   = vars(d_v_a_name), 
      scales = "free_y", 
      labeller = label_wrap_gen(width = 17)) + 
    # Set colour scheme...
    scale_colour_manual(
      values = colours) + 
    # Prettify y axis...
    scale_y_continuous(
      name   = paste0(
        "Imputation prediction\n ", 
        "(ratio of cumulative impact and cumulative FVP)"), 
      labels = scientific,
      limits = c(0, NA),
      expand = c(0, 0), 
      breaks = pretty_breaks())
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.text     = element_text(size = 8),
          axis.text.x   = element_text(hjust = 1, angle = 50),
          axis.title.y  = element_text(
            size = 16, margin = margin(l = 10, r = 20)),
          axis.line     = element_blank(),
          strip.text    = element_text(size = 10),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(0.4, "lines"), 
          legend.position = "none")
  
  # Save figure to file
  save_sub = letters[which(o$metrics %in% metric)]
  save_fig(g, paste0("S16", save_sub))
}

# **** Impact functions ****

# ---------------------------------------------------------
# Plot impact function evaluation
# ---------------------------------------------------------
plot_model_fits = function(metric) {
  
  message("  - Plotting impact function fits: ", metric)
  
  # Load data used for impact function fitting
  data_dt = read_rds("impact", "impact", metric, "data") %>%
    format_d_v_a_name()
  
  # Largest cumulative FVP we need to plot up to
  max_dt = data_dt %>%
    group_by(d_v_a_id, country) %>%
    slice_max(fvps, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(d_v_a_id, country, x_max = fvps) %>%
    as.data.table()
  
  # Expand for a full series of evaluation points
  eval_dt = max_dt %>%
    expand_grid(scale = seq(0, 1, by = 0.01)) %>%
    mutate(fvps = x_max * scale) %>%
    select(d_v_a_id, country, fvps) %>%
    as.data.table()
  
  # Evaluate all points
  fit_dt = evaluate_impact_function(
    data   = eval_dt, 
    metric = metric, 
    uncert = FALSE)
  
  # Remove outliers for more meaningful plot
  outlier_dt = fit_dt %>%
    # Cumulative impact by country...
    group_by(d_v_a_id, country) %>%
    slice_max(impact, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    # Determine outlier threshold...
    group_by(d_v_a_id) %>%
    summarise(mean = mean(impact), 
              sd   = sd(impact)) %>%
    ungroup() %>%
    mutate(outlier = mean + sd * 3) %>%
    select(d_v_a_id, outlier) %>%
    # Determine outlier cases...
    left_join(y  = fit_dt, 
              by = "d_v_a_id") %>%
    filter(impact > outlier) %>%
    select(d_v_a_id, country) %>%
    unique() %>%
    mutate(outlier = TRUE) %>%
    as.data.table()
  
  # Remove outliers from data
  points_dt = data_dt %>%
    left_join(y  = outlier_dt, 
              by = c("d_v_a_id", "country")) %>%
    filter(is.na(outlier)) %>%
    select(d_v_a_name, country, fvps, impact)
  
  # Evaluate selected impact function
  lines_dt = fit_dt %>%
    left_join(y  = outlier_dt, 
              by = c("d_v_a_id", "country")) %>%
    filter(is.na(outlier)) %>%
    format_d_v_a_name() %>%
    select(d_v_a_name, country, fvps, impact)
  
  # Plot function evaluation against the data
  g = ggplot(lines_dt) +
    aes(x = fvps, 
        y = impact, 
        colour = country) +
    # Plot data, then fit on top...
    geom_point(
      data = points_dt,
      size = 0.75,
      alpha = 0.5,
      show.legend = FALSE) +
    geom_line(show.legend = FALSE) +
    # Faceting with wrap labelling...
    facet_wrap(
      facets   = vars(d_v_a_name), 
      labeller = label_wrap_gen(width = 24), 
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
      name   = "Cumulative impact per capita", 
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
  
  # Save figure to file
  save_sub = letters[which(o$metrics %in% metric)]
  save_fig(g, paste0("S17", save_sub))
}

# ---------------------------------------------------------
# Plot function selection statistics
# ---------------------------------------------------------
plot_model_selection = function(metric) {
  
  message("  - Plotting impact model selection: ", metric)
  
  # Define colour scheme
  colour_map = "brewer::set1"
  
  # Load stuff: best fit functions and associtaed coefficients
  choice_dt = read_rds("impact", "model_choice", metric) %>%
    format_d_v_a_name()
  
  # ---- Plot function count ----
  
  # Simple plotting function with a few features
  plot_selection = function(var, type = "count", stat = "number") {
    
    # Determine order - with 'focus' function first
    fn_dict = fn_set(dict = TRUE)
    
    # Number of times each model is optimal
    selection_dt = choice_dt %>% 
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
      arrange_at(intersect(
        x = names(fn_dict), 
        y = names(.))) %>%
      # Final formatting...
      pivot_longer(cols = -var, 
                   names_to  = "fn", 
                   values_to = "val") %>%
      mutate(fn  = recode(fn, !!!fn_dict), 
             fn  = factor(fn, fn_dict),
             var = fct_rev(fct_inorder(var))) %>%
      as.data.table()
    
    # Check figure type flag
    if (type == "count") {
      
      # Number of occurances
      g = ggplot(selection_dt) + 
        aes(x    = var, 
            y    = val, 
            fill = fct_rev(fn)) + 
        geom_col() + 
        coord_flip() + 
        # Set colour scheme...
        scale_fill_manual(
          values = colour_scheme(
            map = colour_map, 
            n   = n_unique(selection_dt$fn))) +
        # Prettify x axis (noting coord_flip)...
        scale_y_continuous(
          name   = paste(first_cap(stat), "of countries"), 
          expand = expansion(mult = c(0, 0.05)), 
          breaks = pretty_breaks()) + 
        # Prettify legend...
        guides(fill = guide_legend(
          reverse = TRUE))
      
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
  
  # Create a variety of plots
  g = list(
    
    # Plot by disease-vaccine-activity
    a = plot_selection("d_v_a_name", stat = "proportion"), 
    b = plot_selection("d_v_a_name", stat = "number"),
    
    # Plot by country
    c = plot_selection("country", type = "density"))
  
  # Save figures of interest
  save_sub = letters[which(o$metrics %in% metric)]
  save_fig(g$a, paste0("S18", save_sub))
}

# **** Historical impact ****

# ---------------------------------------------------------
# Main results plot - historical impact over time
# ---------------------------------------------------------
plot_historical_impact = function(by = NULL) {
  
  # ---- Figure properties ----
  
  # Dictionary for temporal and cumulative subplots
  metrics = qc(deaths, yll, dalys)
  
  # Set scale for set of metrics
  scale = list(
    deaths = list(val = 1e6, str = "million"),
    yll    = list(val = 1e9, str = "billion"), 
    dalys  = list(val = 1e9, str = "billion"))
  
  # Diseases to combine into one colour (set to NULL to turn off)
  grouping = qc(dip, hepb, je, mena, pcv, rota, rubella, yf) 
  
  # Descriptive name for this grouping
  other = "Other pathogens"
  
  # Custom y axes limits for set of metrics
  y_lims = c(200, 11, 11)
  
  # Year range string
  range = paste(range(o$years), collapse = "-")
  
  # ---- Regional or income scope ----
  
  # Grouping datatable: income and region of each country
  country_scope = table("income_status") %>%
    filter(year == max(year)) %>%
    left_join(y  = table("income_dict"), 
              by = "income") %>%
    select(country, income = income_name) %>%
    mutate(income = factor(
      x      = income, 
      levels = table("income_dict")$income_name)) %>%
    append_region_name()
  
  # Global scope
  if (is.null(by)) {
    scope_dt = country_scope %>%
      select(country) %>%
      mutate(scope = "global")
  }
  
  # Group by region or income
  if (!is.null(by)) {
    scope_dt = country_scope %>%
      select(country, scope = !!by)
  }
  
  # Display what we're plotting to user
  msg = "  - Plotting historical impact"
  message(paste(c(msg, by), collapse = ": "))
  
  # ---- Load results ----
  
  # Function to load results for a given metric
  load_fn = function(metric) {
    
    # Load result and append metric ID
    impact_dt = read_rds("history", "burden_averted", metric) %>%
      mutate(metric = !!metric)
    
    return(impact_dt)
  }
  
  # Load results and apply initial formatting
  impact_dt = lapply(metrics, load_fn) %>%
    rbindlist() %>%
    # Remove 'other' pathogens - these go last...
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_id") %>%
    left_join(y  = table("disease_name"), 
              by = "disease") %>%
    # Country grouping...
    left_join(y  = scope_dt, 
              by = "country") %>%
    # Tidy up...
    select(metric, disease, disease_name, 
           scope, year, impact)
  
  # ---- Summarize results ----
  
  # Cumulative impact by metric and disease group
  results_dt = impact_dt %>%
    # Combine subset of diseases...
    mutate(group = ifelse(
      test = disease %in% grouping, 
      yes  = other, 
      no   = disease_name)) %>%
    # Cumulative results for each disease...
    group_by(metric, group, scope, year) %>%
    summarise(impact = sum(impact)) %>%
    mutate(impact = cumsum(impact)) %>%
    ungroup() %>%
    rename(disease = group) %>%
    # All metrics bounded above by DALYs...
    group_by(disease, year) %>%
    mutate(impact = pmin(
      impact, impact[metric == "dalys"])) %>%
    ungroup() %>%
    as.data.table()
  
  # ---- Disease totals ----
  
  # Format scaling details into joinable datatable
  scale_dt = list2dt(scale) %>%
    mutate(metric = names(scale), .before = 1)
  
  # Construct labels: total FVPs over analysis timeframe
  label_dt = results_dt %>%
    # Total burden averted over all years...
    group_by(disease, metric) %>%
    summarise(total = max(impact)) %>%
    ungroup() %>%
    # Apply metric scaling...
    left_join(y  = scale_dt, 
              by = "metric") %>%
    mutate(total = paste(round(total / val, 1), str)) %>%
    select(-val, -str) %>%
    # Append metric and total strings...
    left_join(y  = table("metric_dict"), 
              by = "metric") %>%
    mutate(total = paste0(metric_short, ": ", total)) %>%
    select(-metric_short, -metric_long, -metric_impact) %>%
    # Format into single string per disease...
    pivot_wider(names_from  = metric,
                values_from = total) %>%
    mutate(name = paste0("**", disease, "**"), 
           .after = disease) %>%
    select(disease, name, all_of(metrics)) %>%
    # Collapse into single string with line breaks...
    unite(col = "label", names(.)[-1], 
          sep = "<br>", remove = FALSE) %>%
    select(disease, label) %>%
    as.data.table()
  
  # Simplify labels for non-global plot
  if (!is.null(by)) {
    label_dt %<>% 
      select(disease) %>%
      mutate(label = disease)
  }
  
  # ---- Disease colours ----
  
  # Common plotting order of diseases based on global deaths averted
  diseases = impact_dt %>%
    filter(metric == metrics[1], 
           !disease %in% grouping) %>%
    # Order by decreasing global impact...
    group_by(disease_name) %>%
    summarise(order = sum(impact)) %>%
    ungroup() %>%
    arrange(-order) %>%
    # Covert to vector
    pull(disease_name) %>%
    c(other)
  
  # Set colours as named vector for consistent colouring
  colours = colours_who("category", length(diseases))
  
  # ---- Construct plotting datatables ----
  
  # Set order and append total labels to plotting data
  plot_dt = results_dt %>%
    # Apply scaler...
    left_join(y  = scale_dt,
              by = "metric") %>%
    mutate(value = impact / val) %>%
    # Full metric description...
    left_join(y  = table("metric_dict"), 
              by = "metric") %>%
    mutate(metric_id = factor(metric, metrics), 
           metric = paste0(
             metric_impact, "\n(cumulative ", 
             range, ", in ", str, "s)")) %>%
    # Append labels...
    left_join(y  = label_dt, 
              by = "disease") %>%
    # Plotting order...
    mutate(disease = factor(disease, diseases)) %>%
    arrange(metric_id, disease) %>%
    mutate(metric = fct_inorder(metric), 
           label  = fct_inorder(label)) %>%
    select(label, metric_id, metric, scope, year, value)
  
  # Blank plot to get custom y axes limits
  blank_dt = plot_dt %>%
    select(metric) %>%
    unique() %>%
    mutate(year  = min(o$years), 
           value = y_lims)
  
  # We'll also plot metric totals in each facet
  text_dt = plot_dt %>%
    # Total impact by metric...
    filter(year %in% max(o$years)) %>%
    group_by(metric_id, metric) %>%
    summarise(value = round(sum(value), 1)) %>%
    ungroup() %>%
    # Construct total impact text...
    left_join(y  = scale_dt, 
              by = c("metric_id" = "metric")) %>%
    mutate(text = paste("Total:", value, str)) %>%
    select(metric, text) %>%
    # Set facet positioning...
    left_join(blank_dt, by = "metric") %>%
    mutate(value = 0.95 * value, 
           year  = 0.1 * diff(range(o$years)) + year) %>%
    select(metric, year, value, text)
  
  # ---- Produce plot ----
  
  # Stacked yearly bar plot
  g = ggplot(plot_dt) +
    aes(x = year, 
        y = value) + 
    # Plot bars...
    geom_col(
      mapping = aes(
        fill = label)) +
    # Set colours...
    scale_fill_manual(values = colours) + 
    # Prettiy x axis...
    scale_x_continuous(
      limits = c(min(o$years) - 1, 
                 max(o$years) + 1), 
      expand = expansion(mult = c(0, 0)), 
      breaks = seq(
        from = min(o$years), 
        to   = max(o$years), 
        by   = 5)) +
    # Prettify legend (needed for y spacing to take effect)...
    guides(fill = guide_legend(
      byrow = TRUE))
  
  # Prettify global plot
  if (is.null(by)) {
    g = g + 
      # Facet by metric only...
      facet_wrap(
        facets = vars(metric), 
        scales = "free_y", 
        nrow   = 1) + 
      # Add total text...
      geom_text(
        data    = text_dt, 
        mapping = aes(
          label = text), 
        hjust   = "left", 
        size    = 4) +
      # Prettify y axis...
      geom_blank(data = blank_dt) + 
      scale_y_continuous(
        labels = comma, 
        breaks = pretty_breaks(),
        expand = expansion(mult = c(0, 0)))
  }
  
  # Prettify disaggregated plot
  if (!is.null(by)) {
    g = g + 
      # Facet by metric and country-grouping...
      facet_grid(
        rows   = vars(metric),
        cols   = vars(scope), 
        scales = "free_y", 
        labeller = labeller(
          scope  = label_wrap_gen(20))) + 
      # Prettify y axis...
      scale_y_continuous(
        labels = comma, 
        breaks = pretty_breaks(),
        expand = expansion(mult = c(0, 0.05)))
  }
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title  = element_blank(),
          axis.text   = element_text(size = 9),
          axis.text.x = element_text(hjust = 1, angle = 50), 
          axis.line   = element_blank(),
          plot.title  = element_text(
            margin = margin(t = 10, b = 20), 
            size   = 18,
            hjust  = 0.5), 
          strip.text = element_text(size = 12),
          strip.background = element_blank(), 
          panel.border = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(0.5, "lines"),
          panel.grid.major.y = element_line(linewidth = 0.25),
          legend.title = element_blank(),
          legend.text  = element_markdown(size = 11),
          legend.key   = element_blank(),
          legend.position = "right", 
          legend.spacing.y  = unit(2, "lines"),
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save global results as main manuscript figure
  if (is.null(by)) save_fig(g, "1")
  
  # ... and disaggregated results as supplementary figures
  else if (by == "region") save_fig(g, "S01")
  else if (by == "income") save_fig(g, "S02")
}

# ---------------------------------------------------------
# Main results plot - historical impact over time
# ---------------------------------------------------------
plot_temporal_impact = function(metric) {
  
  message("  - Plotting temporal impact: ", metric)
  
  # Results by disease and region
  plot_dt = read_rds("history", "all_samples", metric) %>%
    # Summarising uncertainty...
    summarise_uncertainty() %>%
    # Append full region and disease names...
    append_region_name() %>%
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_id") %>%
    left_join(y  = table("disease_name"), 
              by = "disease") %>%
    # Results for each disease and each region...
    group_by(disease_name, region, year) %>%
    summarise(impact = sum(impact), 
              lower  = sum(lower), 
              upper  = sum(upper)) %>%
    ungroup() %>%
    rename(disease = disease_name) %>%
    # Remove only single data points...
    add_count(disease, region) %>%
    filter(n > 1) %>%
    select(-n) %>%
    as.data.table()
  
  # Colour scheme
  colours = colours_who_region() %>%
    enframe() %>%
    left_join(y  = table("region_dict"), 
              by = c("name" = "region")) %>%
    select(region_name, value) %>%
    deframe()
  
  # Y axis label describing metric
  axis_label = table("metric_dict") %>%
    filter(metric == !!metric) %>%
    pull(metric_impact)
  
  # Stacked yearly bar plot
  g = ggplot(plot_dt) +
    aes(x = year, 
        y = impact) + 
    # Plot uncertainty bands...
    geom_ribbon(
      mapping = aes(
        ymin = lower,
        ymax = upper,
        fill = region),
      colour = NA,
      alpha  = 0.5) +
    # Plot best estimate line...
    geom_line(
      mapping = aes(colour = region), 
      linewidth = 1.2) + 
    # Facet by temporal-cumulative metric...
    facet_wrap(
      facets = vars(disease), 
      scales = "free_y") + 
    # Set colours...
    scale_colour_manual(values = colours) + 
    scale_fill_manual(values = colours) + 
    # Prettify y axis...
    scale_y_continuous(
      name   = axis_label,
      labels = comma, 
      breaks = pretty_breaks(),
      expand = expansion(mult = c(0, 0.05))) +
    # Prettiy x axis...
    scale_x_continuous(
      limits = c(min(o$years), 
                 max(o$years)), 
      expand = expansion(mult = c(0, 0)), 
      breaks = seq(
        from = min(o$years), 
        to   = max(o$years), 
        by   = 5))
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title.x  = element_blank(),
          axis.title.y  = element_text(size = 20),
          axis.text     = element_text(size = 9),
          axis.text.x   = element_text(hjust = 1, angle = 50), 
          axis.line     = element_blank(),
          strip.text    = element_text(size = 12),
          strip.background = element_blank(), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(0.5, "lines"),
          panel.grid.major.y = element_line(linewidth = 0.25),
          legend.title  = element_blank(),
          legend.text   = element_text(size = 12),
          legend.key    = element_blank(),
          legend.position = "bottom", 
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  if (metric == "deaths") save_fig(g, "S04")
  if (metric == "dalys")  save_fig(g, "S05")
}

# ---------------------------------------------------------
# Plot change in infant mortality rates over time
# ---------------------------------------------------------
plot_infant_mortality = function() {
  
  message("  - Plotting infant mortality rates")
  
  # ---- Figure properties ----
  
  # Ages of interest
  age_bound = 0
  
  # Shorthand for whole analysis period
  year_range = paste(range(o$years), collapse = "-")
  
  # Construct subplot dictionary
  metric_dict = list(
    rate   = paste("Infant mortality rate", year_range), 
    deaths = paste("Cumulative infant deaths", year_range), 
    cov    = "Global vaccine coverage")
  
  # Metric formatting
  metric_lab = list(
    rate   = "percent", 
    deaths = "comma") 
  
  # Dictionary for each scenario
  line_dict = list(
    vaccine    = "Vaccination as obserevd",
    no_vaccine = "Hypothetical scenario:<br>no historical vaccination",
    no_other   = "Hypothetical scenario:<br>no improvement in infant survival")
  
  # Dictionary for each scenario
  fill_dict = list( 
    vaccine = "Reduction in infant mortality<br>attributable to vaccination", 
    other   = "Reduction in infant mortality<br>attributable to other causes")
  
  # Colours of vaccine-other contributions
  fill_colours = c(colours_who("logo", 1), "grey30")
  
  # Set colour scheme for vaccine coverage
  coverage_colours = "viridis::viridis"
  
  # ---- Compile results ----
  
  # Load mortality rates in vaccine and no vaccine scenarios
  mortality = mortality_rates(
    age_bound = age_bound, 
    grouping  = "none")
  
  # Mortality rate over time in each scenario
  rate_dt = mortality$rate %>%
    mutate(metric = "rate") %>%
    select(metric, scenario, year, value) %>%
    as.data.table()
  
  # Cumulative number of deaths in each scenario
  deaths_dt = mortality$value %>%
    group_by(scenario) %>%
    mutate(value = cumsum(value)) %>%
    ungroup() %>%
    mutate(metric = "deaths") %>%
    select(metric, scenario, year, value) %>%
    as.data.table()
  
  # Vaccine coverage over time - global average
  coverage_dt = table("coverage_global") %>%
    filter(coverage > 0,
           !grepl("_BD$", vaccine)) %>%
    # Combine DTP3 coverage...
    mutate(vaccine = ifelse(
      test = vaccine %in% qc(dip, tet, aper),
      yes  = "DTP3",
      no   = vaccine)) %>%
    group_by(vaccine, year) %>%
    summarise(coverage = mean(coverage)) %>%
    ungroup() %>%
    # Append descriptive vaccine names...
    left_join(y  = table("vaccine_name"),
              by = "vaccine") %>%
    select(vaccine = vaccine_name, year, coverage) %>%
    replace_na(list(vaccine = "DTP third dose")) %>%
    # Metric subplot order...
    mutate(metric = metric_dict$cov,
           metric = factor(metric, metric_dict)) %>%
    # Tidy up...
    arrange(vaccine, year) %>%
    as.data.table()
  
  # ---- Impact attributable to vaccination ----
  
  # Percentage of mortality decline attributable to vaccination
  attributable_dt = rate_dt %>%
    mutate(value = value * 100) %>%
    pivot_wider(names_from = scenario) %>%
    mutate(relative = 100 * (no_vaccine - vaccine) / 
             (no_vaccine[1] - vaccine)) %>%
    mutate(n_years = 1 : n(), 
           relative_mean = cumsum(relative) / n_years) %>%
    as.data.table()
  
  # Round vaccine attributable result
  attributable = round(mean(attributable_dt$relative))
  
  # Append this very important outcome to legend labels
  fill_dict$vaccine %<>% paste0(": **", attributable, "%**")
  fill_dict$other   %<>% paste0(": **", 100 - attributable, "%**")
  
  # ---- Construct plotting datatables ----
  
  # Combine into single plotting datatable
  lines_dt = rate_dt %>%
    bind_rows(deaths_dt) %>%
    # Smaller linewidth for hypothetical scenarios...
    mutate(size = ifelse(scenario == "vaccine", 2, 1), 
           .after = scenario) %>%
    # Scenario order...
    rename(line = scenario) %>%
    mutate(line = recode(line, !!!line_dict), 
           line = factor(line, line_dict)) %>%
    # Metric subplot order...
    mutate(metric = recode(metric, !!!metric_dict),
           metric = factor(metric, metric_dict)) %>%
    arrange(metric, line, year)
  
  # Second plotting datatable to highlight vaccine contribution
  area_dt = rate_dt %>%
    rbind(deaths_dt) %>%
    pivot_wider(names_from = scenario) %>%
    # Fill impact attributable to vaccination and other factors...
    mutate(y0.vaccine = pmin(vaccine, no_vaccine),
           y1.vaccine = pmax(vaccine, no_vaccine), 
           y0.other   = pmin(no_vaccine, no_other),
           y1.other   = pmax(no_vaccine, no_other)) %>%
    select(-vaccine, -no_vaccine, -no_other) %>%
    # Format into plottable datatable...
    pivot_longer(cols = -c(metric, year)) %>%
    separate(col  = "name", 
             into = c("var", "fill"), 
             sep  = "\\.") %>%
    pivot_wider(names_from = var) %>%
    # Scenario order...
    mutate(fill = recode(fill, !!!fill_dict), 
           fill = factor(fill, fill_dict)) %>%
    # Metric subplot order...
    mutate(metric = recode(metric, !!!metric_dict),
           metric = factor(metric, metric_dict)) %>%
    # Tidy up...
    select(metric, fill, year, y0, y1) %>%
    arrange(metric, fill, year) %>%
    as.data.table()
  
  # ---- Produce plot ----
  
  # Function to prettify theme
  theme_fn = function(g) {
    
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
            legend.text   = element_markdown(size = 11),
            legend.position   = "right",
            legend.spacing.y  = unit(0.5, "lines"),
            legend.key.height = unit(1.1, "lines"),
            legend.key.width  = unit(1.8, "lines"))
    
    return(g)
  }
  
  # Function to create subplot
  plot_fn = function(plot_metric) {
    
    # Subplot metric
    sub = metric_dict[[plot_metric]]
    
    # Plot infant mortality
    g = ggplot(area_dt[metric == sub]) +
      aes(x = year) +
      # Plot area...
      geom_ribbon(
        mapping = aes(
          ymin = y0,
          ymax = y1, 
          fill = fill),
        colour = NA, 
        alpha  = 0.5) +
      # Plot lines...
      geom_line(
        data    = lines_dt[metric == sub],
        mapping = aes(
          y = value,
          linetype  = line,
          linewidth = size), 
        colour = "black") +
      # Facet by metric...
      facet_wrap(
        facets = vars(metric), 
        scales = "free_y", 
        ncol   = 1) +
      # Set colour scheme...
      scale_fill_manual(
        values = fill_colours) +
      # Set linewidth values
      scale_linewidth(
        range = c(0.5, 1.2)) +
      # Prettify y axis...
      scale_y_continuous(
        labels = get(metric_lab[[plot_metric]]), 
        limits = c(0, NA),
        expand = expansion(mult = c(0, 0.05)), 
        breaks = pretty_breaks()) +
      # Prettiy x axis...
      scale_x_continuous(
        limits = c(min(o$years), max(o$years)), 
        expand = expansion(mult = c(0, 0)), 
        breaks = seq(min(o$years), max(o$years), by = 5)) +
      # Prettify legends...
      guides(
        linewidth = "none", 
        linetype = guide_legend(
          order     = 1, 
          keyheight = 1.8, 
          reverse   = TRUE),
        fill = guide_legend(
          order     = 2, 
          keyheight = 1.8, 
          reverse   = TRUE))
    
    # Apply theme to subplot
    g = theme_fn(g)
    
    return(g)
  }
  
  # Plot coverage over time
  g_cov = ggplot(coverage_dt) +
    aes(x = year, 
        y = coverage, 
        colour = vaccine) +
    # Plot lines...
    geom_line(linewidth = 1.2) +
    # Trivial facet...
    facet_wrap(
      facets = vars(metric)) +  
    # Set colour scheme...
    scale_color_manual(
      values = rev(colour_scheme(
        map = coverage_colours, 
        n   = n_unique(coverage_dt$vaccine)))) +
    # Set linewidth values
    scale_linewidth(
      range = c(0.5, 1.2)) +
    # Prettify y axis...
    scale_y_continuous(
      labels = percent,
      limits = c(0, 1),
      expand = c(0, 0),
      breaks = pretty_breaks()) + 
    # Prettiy x axis...
    scale_x_continuous(
      limits = c(min(o$years), max(o$years)), 
      expand = expansion(mult = c(0, 0)), 
      breaks = seq(min(o$years), max(o$years), by = 5))
  
  # Create subplots
  g1 = plot_fn("rate") 
  g2 = plot_fn("deaths") + theme(legend.position = "none")
  g3 = theme_fn(g_cov)
  
  # Patch subplots together
  g = g1 / g2 / g3
  
  # Save figure to file
  save_fig(g, "2")
}

# ---------------------------------------------------------
# Plot regional differences in infant mortality changes
# ---------------------------------------------------------
plot_mortality_change = function() {
  
  message("  - Plotting regional changes in infant mortality")
  
  # Description of metric and year scope
  metric_str = "infant mortality rate"
  years_str  = paste0("(", paste(range(o$years), collapse = "-"), ")")
  
  # Dictionary for metric type
  type_dict = list(
    abs = paste("Absolute decrease in", metric_str, years_str),
    rel = paste("Relative decrease in", metric_str, years_str), 
    att = paste("Contribution of vaccination to decrease in", 
                metric_str, years_str))
  
  # Metric type colour scheme (+1 for global average)
  colours = c("#EBAA2D", "#DF721F", "#71C2A9", "#808080")
  
  # Mortality rates - by region and global average
  region_dt = mortality_rates(grouping = "region")$rate
  world_dt  = mortality_rates(grouping = "none")$rate %>%
    mutate(group = "World")
  
  # Absolute decrease in infant mortality rate
  abs_dt = rbind(region_dt, world_dt) %>%
    filter(scenario == "vaccine", 
           year %in% range(o$years)) %>%
    group_by(group) %>%
    summarise(start = value[1], 
              end   = value[2], 
              abs   = start - end) %>%
    ungroup() %>%
    arrange(-abs) %>%
    as.data.table()
  
  # Relative decrease in infant mortality rate
  rel_dt = rbind(region_dt, world_dt) %>%
    filter(scenario == "vaccine", 
           year %in% range(o$years)) %>%
    group_by(group) %>%
    summarise(rel = 1 - min(value) / max(value)) %>%
    ungroup() %>%
    as.data.table()
  
  # Percentage of mortality decline attributable to vaccination
  att_dt = rbind(region_dt, world_dt) %>%
    pivot_wider(names_from = scenario) %>%
    select(group, year, vaccine, no_vaccine) %>%
    group_by(group) %>%
    mutate(att = (no_vaccine - vaccine) / 
             (no_vaccine[1] - vaccine)) %>%
    summarise(att = mean(att)) %>%
    ungroup() %>%
    as.data.table()
  
  # Order regions by absolute decrease in mortality rates
  order = abs_dt %>%
    filter(group != "World") %>%
    arrange(desc(abs)) %>%
    pull(group) %>%
    c("World")
  
  # Labels bars for clarity
  label_dt = abs_dt %>%
    # Convert to percentage...
    mutate(across(
      .cols = c(abs, start, end), 
      .fns  = ~ . * 100)) %>%
    # Construct string explaining start and end values...
    mutate(label = sprintf(
      "%.1f%%\n(%.1f%% to %.1f%%)", 
      abs, start, end)) %>%
    # Only required for absolute bars...
    mutate(type = "abs") %>%
    select(group, type, label)
  
  # Construct plotting datatable
  plot_dt = abs_dt %>%
    select(-start, -end) %>%
    # Append all results...
    left_join(rel_dt, by = "group") %>%
    left_join(att_dt, by = "group") %>%
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
  save_fig(g, "S03")
}

# ---------------------------------------------------------
# Plot absolute and relative probability of death in 2024
# ---------------------------------------------------------
plot_survival_increase = function(log_age = FALSE) {
  
  message("  - Plotting increase in childhood survival")
  
  # Year to create snapshot for - makes sence to be final year of analysis
  snapshot_year = max(o$years)
  
  # Age bound (nice round value or a square)
  age_bound = ifelse(log_age, 64, 50)
  
  # Construct pplot title
  title = "Historical vaccination compared to hypothetical no vaccination"
  
  # Construct a somewhat elaborate y axis label
  y_label = paste0(
    "Marginal increase in survival probability in n<sup>th</sup> ", 
    "year of life (", max(o$years), " snapshot)")
  
  # Dictionary for absolute and marginal metrics
  metric_dict = c(
    relative = "Relative", 
    absolute = "Absolute")
  
  # Estimated child deaths averted by vaccination
  averted_dt = read_rds("history", "burden_averted_deaths") %>%
    # Results in year of interest...
    filter(year == snapshot_year) %>%
    # Apply age effect...
    left_join(y  = table("impact_age_multiplier"), 
              by = "d_v_a_id", 
              relationship = "many-to-many") %>%
    mutate(impact = impact * scaler) %>%
    # Summarise results over region...
    append_region_name() %>%
    group_by(region, age) %>%
    summarise(averted = sum(impact)) %>%
    ungroup() %>%
    as.data.table()
  
  # Calculate absolute and relative effect in year of interest
  survival_dt = table("wpp_pop") %>%
    left_join(y  = table("wpp_deaths"), 
              by = c("country", "year", "age")) %>%
    # Retain only years and ages of interest...
    filter(year == snapshot_year, 
           age  <= age_bound - log_age * 1) %>%
    # Append region...
    append_region_name() %>%
    # Average prob of death by region and year...
    group_by(region, age) %>%
    summarise(pop    = sum(pop), 
              deaths = sum(deaths)) %>%
    ungroup() %>%
    # Append deaths averted results...
    left_join(y  = averted_dt, 
              by = c("region", "age")) %>%
    mutate(age = age + log_age * 1) %>%
    # Cumulative results...
    group_by(region) %>%
    mutate(c_pop   = cumsum(pop), 
           c_death = cumsum(deaths), 
           c_avert = cumsum(averted)) %>%
    ungroup() %>%
    # Compute absolute and relative results...
    mutate(vaccine    = c_death / c_pop, 
           no_vaccine = (c_death + c_avert) / c_pop, 
           absolute   = no_vaccine - vaccine,
           relative   = 1 - vaccine / no_vaccine) %>%
    # Melt to long format...
    select(region, age, pop, absolute, relative) %>%
    pivot_longer(cols = c(absolute, relative), 
                 names_to = "metric") %>%
    # Tidy up...
    select(metric, region, age, pop, value) %>%
    arrange(metric, region, age) %>%
    as.data.table()
  
  # Combine with global values
  plot_dt = survival_dt %>%
    # Calculate global average...
    group_by(metric, age) %>%
    mutate(weight = pop / sum(pop)) %>%
    summarise(value = sum(value * weight)) %>%
    ungroup() %>%
    # Append regional results...
    mutate(region = "World", 
           size   = 2) %>%
    bind_rows(survival_dt) %>%
    select(-pop) %>%
    replace_na(list(size = 1)) %>%
    # Metric order...
    mutate(metric = recode(metric, !!!metric_dict), 
           metric = factor(metric, metric_dict)) %>%
    # Tidy up...
    arrange(size, region) %>%
    mutate(region = fct_inorder(region)) %>%
    as.data.table()
  
  # ---- Produce plot ----
  
  # WHO regional colours, with global average in black
  colours = unname(c(colours_who_region(), "#000000"))
  
  # Plot absolute and relative difference in snapshot year
  g = ggplot(plot_dt) +
    aes(x = age, 
        y = value, 
        colour = fct_rev(region), 
        linewidth = size) +
    geom_line() +
    # Facet by absolute/relative...
    facet_wrap(
      facets = vars(metric), 
      scales = "free_y", 
      ncol   = 1) + 
    # Set colour scheme
    scale_colour_manual(
      values = rev(colours)) +
    scale_linewidth(
      range = c(1.2, 3)) +
    # Add plot title...
    ggtitle(title) +
    # Prettify y axis...
    scale_y_continuous(
      name   = y_label, 
      labels = percent, 
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05)), 
      breaks = pretty_breaks()) +  
    # Prettify legend...
    guides(
      linetype  = "none", 
      linewidth = "none", 
      colour = guide_legend(
        reverse = TRUE, 
        override.aes = list(
          linewidth = 2)))
  
  # Transform x axis to log2 scale
  if (log_age == TRUE) {
    g = g + scale_x_continuous(
      name   = "Age (log2 scale)",
      trans  = "log2", 
      breaks = 2 ^ (1 : (sqrt(age_bound) - 1) - 1))
  }
  
  # Keep x axis on real number line
  if (log_age == FALSE) {
    g = g + scale_x_continuous(
      name   = "Age in years", 
      expand = expansion(add = 1),
      breaks = seq(0, age_bound, by = 5))
  }
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(axis.title   = element_text(size = 16),
          axis.title.x = element_text(
            margin = margin(t = 20, b = 10)),
          axis.title.y = element_markdown(
            margin = margin(l = 10, r = 20)),
          axis.text = element_text(size = 11),
          axis.line  = element_blank(),
          strip.text = element_text(size = 14),
          strip.background = element_blank(), 
          plot.title    = element_text(
            margin = margin(t = 10, b = 20), 
            size   = 18,
            hjust  = 0.5), 
          panel.border  = element_rect(
            linewidth = 0.5, fill = NA),
          panel.spacing = unit(1, "lines"),
          legend.title  = element_blank(),
          legend.text   = element_text(size = 13),
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines"))
  
  # Save figure to file
  save_fig(g, "3")
}

# **** Helper functions ****

# ---------------------------------------------------------
# Append metric (deaths and DALYs) descriptive names
# ---------------------------------------------------------
append_metric_name = function(dt) {
  
  # Append metric names
  name_dt = dt %>%
    left_join(y  = table("metric_dict"), 
              by = "metric") %>%
    select(-metric, -metric_short) %>%
    rename(metric = metric_long) %>%
    mutate(metric = fct_inorder(metric))
  
  return(name_dt)
}

# ---------------------------------------------------------
# Append full region names
# ---------------------------------------------------------
append_region_name = function(dt) {
  
  # Throw error if region already appended - too confusing
  if ("region" %in% names(dt))
    stop("Provide plotting datatable without 'region'")
  
  # Country - region name mapping
  region_dt = table("country") %>%
    left_join(y  = table("region_dict"), 
              by = "region") %>%
    select(region = region_name, country) %>%
    # Convert to factors...
    arrange(region, country) %>%
    mutate(region = fct_inorder(region))
  
  # Append region names
  name_dt = dt %>%
    left_join(y  = region_dt, 
              by = "country")
  
  return(name_dt)
}

# ---------------------------------------------------------
# Full descriptive names for disease, vaccine, or vaccine type
# ---------------------------------------------------------
append_d_v_t_name = function(dt) {
  
  # All columns to update
  d_v_t = intersect(
    x = names(dt), 
    y = qc(disease, vaccine, type))
  
  # Iterate through columns to update
  for (x in d_v_t) {
    
    # Column name of full description
    x_name = paste1(x, "name")
    
    # Plotting order
    x_ord = table(x_name)[[x_name]]
    
    # Append full name description
    dt %<>%
      lazy_dt() %>%
      left_join(y  = table(x_name), 
                by = x) %>%
      select(-all_of(x)) %>%
      # Rename and reorder columns...
      rename(.x := all_of(x_name)) %>%
      mutate(.x = factor(.x, x_ord)) %>%
      rename(!!x := .x) %>%
      select(all_names(dt)) %>%
      as.data.table()
  }
  
  return(dt)
}

# ---------------------------------------------------------
# Set natural order of d_v_a names - appending if necessary
# ---------------------------------------------------------
format_d_v_a_name = function(dt) {
  
  # Natural order of d_v_a_name (not necessarily alphabetical)
  d_v_a_dt = table("d_v_a") %>%
    select(d_v_a_id, d_v_a_name)
  
  # Check if d_v_a_name is already defined
  if (!"d_v_a_name" %in% names(dt)) {
    
    # Append name if it's not already in the data
    dt %<>%
      left_join(y  = d_v_a_dt, 
                by = "d_v_a_id")
  }
  
  # Set order using factors according to d_v_a_dt
  order_dt = dt %>%
    mutate(d_v_a_name = factor(
      x      = d_v_a_name, 
      levels = d_v_a_dt$d_v_a_name)) %>%
    arrange(d_v_a_name)
  
  return(order_dt)
}

# ---------------------------------------------------------
# Save a ggplot figure to file with default settings
# ---------------------------------------------------------
save_fig = function(g, fig_num, width = NULL, height = NULL) {
  
  # Set default height and width if undefined
  if (is.null(width))  width  = o$save_width
  if (is.null(height)) height = o$save_height
  
  # Figure formats to save with
  formats = o$figure_format
  if (grepl("^[0-9]+", fig_num)) 
    formats = c(formats, o$manuscript_format)
  
  # Repeat the saving process for each image format in figure_format
  for (fig_format in formats) {
    
    # File name and extension
    file_name = paste("Figure", fig_num)
    file_ext  = paste0(".", fig_format)
    
    # Concatenate with path
    full_path = paste0(o$pth$figures, file_name, file_ext)
    
    # Save figure (size specified in options.R)
    ggsave(full_path, 
           plot   = g, 
           width  = width, 
           height = height, 
           device = fig_format, 
           dpi    = o$save_resolution, 
           units  = o$save_units)
  }
}

