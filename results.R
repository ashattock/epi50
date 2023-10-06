###########################################################
# RESULTS
#
# xxxxxxx
#
###########################################################

# ---------------------------------------------------------
# xxxxxxxx
# Called by: launch.R, main.R (and other launch-style scripts)
# ---------------------------------------------------------
run_results = function() {
  
  # Only continue if specified by do_step
  if (!is.element(4, o$do_step)) return()
  
  message("* Producing results")
  
  # Load tables up front (see db_utils.R)
  load_tables("wpp_input")
  
  # Load impact factors calculated in step 2
  # impact_dt       = try_load(o$pth$impact_factors, "impact_dt")
  scenario_impact = try_load(o$pth$impact_factors, "scenario_impact")
  
  # Load uncertainty generated in step 3
  draws_dt = try_load(o$pth$uncertainty, "draws")
  
  # ---- Reference results ----
  
  message(" - Saving results")
  
  browser() # Tidy up...
  
  country_dt = country_table %>%
    mutate(gavi = gavi73 == 1) %>%
    select(country, country_name, region, economy, gavi)
  
  # Number of births per year
  cohort_dt = wpp_input %>%
    filter(age == 0) %>% # Why??
    group_by(country, year) %>%
    summarise(cohort_size = sum(nx)) %>%
    ungroup() %>%
    as.data.table()
  
  g = ggplot(cohort_dt, aes(x = year, y = cohort_size, group = country, colour = country)) + 
    geom_line()
  
  results_dt = draws_dt %>%
    left_join(y  = country_table, 
              by = "country") %>%
    select()

  ## Save for sharing
  out_dt <- merge(
    draws_dt,
    country_table[, .(country, country_name, region, economy, gavi73)],
    by = "country"
  )
  out_dt <- out_dt[order(country_name, disease, vaccine, activity_type)]
  out_dt[is.na(deaths_averted), deaths_averted := 0]
  
  out_dt <- merge(out_dt, cohort_dt, by = c("country", "year"))
  
  total_dt <- wpp_input[, .(nx = sum(nx)), by = .(country, year)]
  setnames(total_dt, "nx", "total_pop")
  out_dt <- merge(out_dt, total_dt, by = c("country", "year"))
  setcolorder(
    out_dt, 
    c(
      "country_name", "country", "year", "disease", 
      "vaccine", "activity_type", "age", "impact_factor", "fvps", "deaths_averted",
      "region", "economy", "gavi73", "cohort_size", "total_pop", paste0("draw_", 1:200)
    )
  )
  
  # Save results to file
  save_file(out_dt, o$pth$results, "reference_results")
  
  # Upload results if desired
  if (o$results_upload)
    upload_object(out_dt, "ia2030_reference_results")
  
  # ---- Impact by region ----
  
  no_lin_range_cov <- gen_ia2030_goals(linear = F, no_covid_effect = 2022, 
                                       intro_year = 2025, intro_range = T)
  scen_dt <- merge(
    no_lin_range_cov,
    country_table[, .(country, country_name)],
    by = "country"
  )
  
  browser() # v_at_table is now v_a_table
  
  scen_dt <- merge(scen_dt, v_at_table, by = "v_at_id")
  write.csv(scen_dt, "outputs/ia2030_cov_trajectories.csv", row.names = F)
  
  
  dt <- scenario_impact
  dt <- merge(dt, country_table[, .(country, region)], by = "country")
  baseline_dt <- dt[year == 2019]
  setnames(baseline_dt, "deaths_averted", "baseline_deaths_averted")
  baseline_dt[, year := NULL]
  table_dt <- merge(dt[year %in% 2021:2030], baseline_dt,
                    by = c("country", "region", "disease", "vaccine", "activity_type"),
                    all.x = T)
  table_dt[is.na(baseline_deaths_averted), baseline_deaths_averted := 0]
  table_dt[is.na(deaths_averted), deaths_averted := 0]
  table_dt[activity_type %in% c("routine", "combined"), 
           incremental := deaths_averted - baseline_deaths_averted]
  table_totals <- table_dt[, .(total = sum(deaths_averted, na.rm = T),
                               incremental = sum(incremental, na.rm = T)), by = .(region, disease, year)][
                                 order(region, disease, year)
                               ]
  
  # NOTE: Used to produce markdown document
  # write.csv(table_totals, "outputs/detailed_results.csv", row.names = F)
  save_file(table_totals, o$pth$results, "detailed_results")
  
  # ---- xxxxxx ----
  
  dt <- scenario_impact
  dt <- merge(dt, country_table[, .(country, economy)], by = "country")
  global_dt <- dt[, .(deaths_averted = sum(deaths_averted, na.rm = T)), by = .(year, disease, vaccine, activity_type)]
  global_dt[, economy := "Global"]
  economy_dt <- dt[, .(deaths_averted = sum(deaths_averted, na.rm = T)), by = .(economy, year, disease, vaccine, activity_type)]
  all_dt <- rbind(global_dt, economy_dt)
  baseline_dt <- all_dt[year == 2019]
  setnames(baseline_dt, "deaths_averted", "baseline_deaths_averted")
  baseline_dt[, year := NULL]
  table_dt <- merge(all_dt[year %in% 2021:2030], baseline_dt,
                    by = c("economy", "disease", "vaccine", "activity_type"), all.x = T)
  table_dt[, incremental := deaths_averted - baseline_deaths_averted]
  table_totals <- table_dt[, .(total = sum(deaths_averted, na.rm = T) / 1e5,
                               incremental = sum(incremental, na.rm = T) / 1e5), by = economy]
  
  # write.csv(table_totals, "outputs/results_table_income.csv", row.names = F)
  
  # ---- Impact by income ----
  
  dt <- scenario_impact
  dt <- merge(dt, country_table[, .(country, economy)], by = "country")
  baseline_dt <- dt[year == 2019]
  setnames(baseline_dt, "deaths_averted", "baseline_deaths_averted")
  baseline_dt[, year := NULL]
  table_dt <- merge(dt[year %in% 2021:2030], baseline_dt,
                    by = c("country", "economy", "disease", "vaccine", "activity_type"),
                    all.x = T)
  table_dt[is.na(baseline_deaths_averted), baseline_deaths_averted := 0]
  table_dt[is.na(deaths_averted), deaths_averted := 0]
  table_dt[activity_type %in% c("routine", "combined"), 
           incremental := deaths_averted - baseline_deaths_averted]
  table_totals <- table_dt[, .(total = sum(deaths_averted, na.rm = T),
                               incremental = sum(incremental, na.rm = T)), by = .(economy, disease, year)][
                                 order(economy, disease, year)]
  
  # NOTE: Used to produce markdown document
  # write.csv(table_totals, "outputs/detailed_results_income.csv", row.names = F)
  save_file(table_totals, o$pth$results, "detailed_results_income")
  
  # ---- Produce markdown results document ----
  
  # Check flag
  if (o$results_markdown) {
    
    # Render document
    #
    # NOTE: Uses file outputs/detailed_results.csv
    rmarkdown::render(input             = paste0(o$pth$code, "results.Rmd"), 
                      # envir             = globalenv(), # globalenv
                      output_dir        = o$pth$figures, 
                      knit_root_dir     = o$pth$figures,
                      intermediates_dir = o$pth$figures, 
                      # clean             = FALSE,  # Save intermediates
                      quiet             = TRUE)
  }
  
  # ---- All diagnostic plots ----
  
  # Check flag
  if (o$plot_diagnostics) {
    
    message(" - Plotting diagnostics")
    
    # Plot parameters of fitted beta distribution to vaccine efficacy (GBD diseases only)
    fig_name = "Uncertainty distributions - GBD diseases"
    plot_gbd_uncertainty_dist(fig_name)  # See plotting.R
    
    # Plot uncertainty draws for all diseases
    fig_name = "Uncertainty draws - All diseases"
    plot_draws(fig_name)  # See plotting.R
    
    # Plot annual totals to check alignment of means
    fig_name = "Uncertainty bounds - Annual total"
    plot_annual_total(fig_name)  # See plotting.R
  }
}

