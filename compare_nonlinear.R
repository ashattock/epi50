
# ---------------------------------------------------------
# Parent function for comparing linear to non-linear approach
# ---------------------------------------------------------
compare_nonlinear = function() {
  
  # Load stuff up front
  load_tables("wpp_input", "coverage")
  
  # We'll calculate DALYs for all viable countries
  countries = loc_table$location_iso3
  
  # An example country for diagnostic plots
  plot_country = "ETH"  # countries[1]
  plot_disease = c("Hib", "Measles", "PCV", "Rota")  # Childhood diseases
  
  # ---- Load coverage ----
  
  # Load population size of each country over time
  pop_dt = wpp_input %>%
    filter(country %in% countries) %>%
    group_by(country, year) %>%
    summarise(pop = sum(nx)) %>%
    ungroup() %>%
    as.data.table()
  
  # Load FVPs over time
  fvps_dt = get_temporal_fvps() %>%
    filter(country %in% countries) %>%
    # Coverage over lifetime?...
    # total_coverage() %>%
    # Cumulative sum of FVPs...
    group_by(country, d_v_a) %>%
    mutate(fvps_cum = cumsum(fvps)) %>%
    ungroup() %>%
    # ... relative to 100k people...
    left_join(y  = pop_dt, 
              by = c("country", "year")) %>%
    mutate(fvps_rel = o$per_person * fvps_cum / pop) %>%
    select(-pop) %>%
    as.data.table()
  
  # y = coverage %>%
  #   filter(country %in% countries, 
  #          v_at_id == 3) %>%
  #   total_coverage()
  
  # x = get_scenario_fvps()
  
  # ---- Method 1: impact factors ----
  
  # NOTE: Here we are using impact factors as per Austin's first iteration
  
  results_2021 = readRDS(file.path(o$pth$main, "impact_factors_2019.rds"))
  # results_2021 = try_load(o$pth$impact_factors, "impact_dt")
  
  # TEMP: Recode GBD disease names to old format until it's changed everywhere
  temp_dict = c("D" = "Dip", "T" = "Tet", "P" = "Per", "TB" = "TB")
  
  # Load impact factors calculated in step 2
  m1_dt = results_2021 %>%
    # Append country 'names'...
    left_join(y  = loc_table[, .(location_id, country = location_iso3)], 
              by = "location_id") %>%
    filter(country %in% countries) %>%
    # TEMP: Recode GBD disease names to old format until it's changed everywhere
    mutate(disease = recode(disease, !!!temp_dict)) %>%
    # Apply decent D-V-A names...
    left_join(y  = d_v_a_name(), 
              by = c("disease", "vaccine", "activity_type")) %>%
    # Take only what we need
    select(country, d_v_a, impact_factor)
  
  # Calculate impact using impact factors
  m1_impact = fvps_dt %>%
    select(country, d_v_a, year, fvps, fvps_cum) %>%
    inner_join(y  = m1_dt, 
               by = c("country", "d_v_a")) %>%
    mutate(impact     = impact_factor * fvps, 
           impact_cum = impact_factor * fvps_cum) %>%
    # Tidy up...
    select(country, d_v_a, year, 
           fvps,   fvps_cum, 
           impact, impact_cum) %>%
    mutate(method = "Impact factor method")
  
  # ---- Method 2: model fitting ----
  
  # NOTE: Here we are fitting to cumulative FVPs (could be non-linear)
  
  # Evaluate best fitting model at cumulative FVPs per population-person
  m2_dt = evaluate_best_model(country = countries, x = fvps_dt)
  
  m2_impact = fvps_dt %>%
    inner_join(y  = m2_dt, 
               by = c("country", "d_v_a", "fvps_rel")) %>%
    # Population weight for population-level impact...
    left_join(y  = pop_dt, 
              by = c("country", "year")) %>%
    mutate(impact_cum = impact_rel * pop / o$per_person) %>%
    # Revere cumsum to derive annual impact...
    group_by(country, d_v_a) %>%
    mutate(impact = rev_cumsum(impact_cum)) %>%
    ungroup() %>%
    # Tidy up...
    select(country, d_v_a, year, 
           fvps,   fvps_cum, 
           impact, impact_cum) %>%
    mutate(method = "Non-linear method")
  
  # ---- Data ----
  
  # Combine both methods for plotting
  impact_dt = rbind(m1_impact, m2_impact) %>%
    left_join(y  = d_v_a_name(), 
              by = "d_v_a") %>%
    filter(#!grepl(".* campaign$", d_v_a), 
           disease %in% plot_disease) %>%
    arrange(method, country, d_v_a, year)
  
  # Also load the data - we'll plot fits against this
  vimc_dt = readRDS(paste0(o$pth$testing, "vimc_dt.rds"))
  
  # VIMC 'data' (ie impact estimates)
  data_dt = vimc_dt %>%
    left_join(y  = d_v_a_name(), 
              by = "d_v_a") %>%
    filter(#!grepl(".* campaign$", d_v_a),
           disease %in% plot_disease, 
           country %in% countries)
  
  # ---- Plot impact comparison ----
  
  plot_impact = impact_dt[country == plot_country]
  plot_data   = data_dt[country   == plot_country]
  
  # Is this the primary figure? Cumulative impact over time?
  g = (ggplot(plot_impact) +
         aes(x = year, y = impact_cum, colour = method) +
         geom_line(size = 1.5) +
         geom_point(data = plot_data, colour = "black") +
         facet_wrap(~d_v_a, scales = "free_y")) %>%
    prettify3(save = c("Compare", "year v impact cumulative"))
  
  g1 = (ggplot(plot_impact) +
          aes(x = year, y = impact, colour = method) +
          geom_line(size = 1.5) +
          geom_point(data = plot_data, colour = "black") +
          facet_wrap(~d_v_a, scales = "free_y")) %>%
    prettify3(save = c("Compare", "year v impact"))
  
  g2 = (ggplot(plot_impact) +
          aes(x = fvps, y = impact, colour = method) +
          geom_line(size = 1.5) +
          geom_point(data = plot_data, colour = "black") +
          facet_wrap(~d_v_a, scales = "free")) %>%
    prettify3(save = c("Compare", "FVP v impact"))
  
  g3 = (ggplot(plot_impact) +
          aes(x = fvps_cum, y = impact, colour = method) +
          geom_line(size = 1.5) +
          geom_point(data = plot_data, colour = "black") +
          facet_wrap(~d_v_a, scales = "free")) %>%
    prettify3(save = c("Compare", "FVP cumulative v impact"))
  
  g4 = (ggplot(plot_impact) +
          aes(x = fvps_cum, y = impact_cum, colour = method) +
          geom_line(size = 1.5) +
          geom_point(data = plot_data, colour = "black") +
          facet_wrap(~d_v_a, scales = "free")) %>%
    prettify3(save = c("Compare", "FVP cumulative v impact cumulative"))
  
  # ---- Total impact comparison ----
  
  annual_dt = impact_dt %>%
    filter(!grepl(".* campaign$", d_v_a)) %>%
    group_by(method, d_v_a, year) %>%
    summarise(impact = sum(impact)) %>%
    ungroup() %>%
    as.data.table()
  
  plot_years = data.table(
    from = c(2000, 2000, 2020), 
    to   = c(2022, 2039, 2030)) 
  
  # Funcion to prettify the following special plots
  prettify_fill = function(g, type, y0, y1) {
    
    # Axes label dictionary
    lab_dict = c(
      area = "Deaths averted",
      bar  = "Cumulative deaths averted")
    
    # COnstruct y label for this plot
    y_lab = paste0(lab_dict[type], " (", y0, "-", y1, ")")
    
    # Prettyify axes
    g = g + scale_y_continuous(
      name   = y_lab,
      limits = c(0, NA), 
      expand = expansion(mult = c(0, 0.05)),
      breaks = pretty_breaks(), 
      labels = comma)
    
    # Prettify x axis for area plots only
    if (type == "area")
      g = g + scale_x_continuous(expand = c(0, 0), breaks = pretty_breaks())
    
    # Apply colours
    cols = colour_scheme("viridis::viridis", n = n_unique(annual_dt$d_v_a))
    g = g + scale_fill_manual(values = cols)
    
    # Prettify theme
    g = g + theme_classic() + 
      theme(strip.text    = element_text(size = 12),
            axis.title.y  = element_text(size = 16),
            axis.title.x  = element_blank(),
            axis.text     = element_text(size = 10),
            axis.line     = element_blank(),
            panel.border  = element_rect(size = 1, colour = "black", fill = NA),
            panel.spacing = unit(1, "lines"),
            strip.background = element_blank(), 
            legend.title  = element_blank(),
            legend.text   = element_text(size = 11),
            legend.key    = element_blank(),
            legend.key.height = unit(2, "lines"),
            legend.key.width  = unit(2, "lines")) 
    
    return(g)
  }
  
  # Iterate through plots to produce
  for (i in seq_len(nrow(plot_years))) {
    y0 = plot_years$from[i]
    y1 = plot_years$to[i]
    
    years_dt = annual_dt %>%
      filter(year >= y0, 
             year <= y1)
    
    total_dt = years_dt %>%
      group_by(method, d_v_a) %>%
      summarise(impact_cum = sum(impact)) %>%
      ungroup() %>%
      as.data.table()
    
    ga = (ggplot(years_dt) +
            aes(x = year, y = impact, fill = d_v_a) + 
            geom_area() + 
            facet_wrap(~method)) %>%
      prettify_fill("area", y0, y1)
    
    gb = (ggplot(total_dt) +
            aes(x = method, y = impact_cum, fill = d_v_a) + 
            geom_bar(stat = "identity")) %>%
      prettify_fill("bar", y0, y1)
    
    g = ggarrange(ga, gb, ncol = 1, 
                  common.legend = TRUE, 
                  legend = "right")
    
    fig_save(g, dir = "testing", "Non-linear effect", y0, y1)
  }
}

# ---------------------------------------------------------
# Load coverage data and extract FVPs over time
# ---------------------------------------------------------
get_temporal_fvps = function() {
  
  # Load and format FVPs over time
  fvps_dt = coverage %>%
    # Summarise over age...
    group_by(country, v_at_id, year, sex_id) %>%
    summarise(fvps = sum(fvps)) %>%
    ungroup() %>%
    # Deal with gender groupings...
    mutate(sex_id = paste0("s", sex_id)) %>%
    pivot_wider(names_from  = sex_id, 
                values_from = fvps, 
                values_fill = 0) %>%
    mutate(fvps = ifelse(s1 == 0, s2 + s3, s1)) %>%
    select(country, v_at_id, year, fvps) %>%
    # Apply decent D-V-A names...
    left_join(y  = v_at_table, 
              by = "v_at_id") %>%
    left_join(y  = d_v_at_table, 
              by = c("vaccine", "activity_type")) %>%
    left_join(y  = d_v_a_name(), 
              by = c("disease", "vaccine", "activity_type")) %>%
    # Reduce down to what we're interested in...
    filter(year >= 2000) %>%
    select(country, d_v_a, year, fvps) %>%
    arrange(country, d_v_a, year) %>%
    as.data.table()
  
  return(fvps_dt)
}

# ---------------------------------------------------------
# Apply colour scheme and tidy up axes - impact plots
# ---------------------------------------------------------
prettify3 = function(g, save = NULL) {
  
  # Construct manual colour scheme
  cols = c("gold", "forestgreen")
  
  # Axes label dictionary
  lab_dict = c(
    year       = "Year",
    fvps       = "Fully vaccinated persons (FVPs) in one year",
    fvps_cum   = "Cumulative fully vaccinated persons (FVPs)",
    impact     = "Deaths averted (per FVP)",
    impact_cum = "Cumulative deaths averted (per FVP)")
  
  # Extract info from the plot
  g_info = ggplot_build(g)
  
  # Prettyify axes
  g = g + 
    scale_x_continuous(
      name   = lab_dict[g_info$plot$labels$x], 
      expand = expansion(mult = c(0, 0.05)),
      labels = comma) +
    scale_y_continuous(
      name   = lab_dict[g_info$plot$labels$y], 
      limits = c(0, NA), 
      expand = expansion(mult = c(0, 0.05)),
      labels = comma)
  
  # Apply colours
  g = g + scale_colour_manual(values = cols)
  
  # Prettify theme
  g = g + theme_classic() + 
    theme(strip.text    = element_text(size = 11),
          axis.title    = element_text(size = 15),
          axis.text     = element_text(size = 8, angle = 30, hjust = 1),
          axis.line     = element_blank(),
          panel.border  = element_rect(size = 1, colour = "black", fill = NA),
          panel.spacing = unit(1, "lines"),
          strip.background = element_blank(), 
          legend.title  = element_blank(),
          legend.text   = element_text(size = 12),
          legend.key    = element_blank(),
          legend.key.height = unit(2, "lines"),
          legend.key.width  = unit(2, "lines")) 
  
  # Save plots to file
  if (!is.null(save))
    fig_save(g, dir = "testing", save)
  
  return(g)
}

