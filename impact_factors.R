###########################################################
# IMPACT FACTORS
#
# xxxxxxx
#
###########################################################

# ---------------------------------------------------------
# Parent function to calculate impact factors from relative risk
# Called by: launch.R, main.R (and other launch-style scripts)
# ---------------------------------------------------------
run_impact_factors = function() {
  
  # Only continue if specified by do_step
  if (!is.element(2, o$do_step)) return()
  
  message("* Calculating impact factors")
  
  # Load tables up front (see db_utils.R)
  load_tables("coverage")  # Needed by calc_impact_factors()
  
  # ---- Calculate impact factors ----
  
  # Load relative risk calculations and predictions from step 1
  rr_dt = try_load(o$pth$relative_risk, "relative_risk")
  
  # Calculate impact factors and rake to VIMC
  #
  # TODO: What is meant by 'rake' here?
  impact_factors = calc_impact_factors(rr_dt)
  impact_dt      = rake_impact(impact_factors)
  
  # Save to file
  save_file(impact_dt, o$pth$impact_factors, "impact_dt")
  
  # ---- Use impact factors to determine deaths averted ----
  
  # ?? What is the different between impact_dt and scenario_impact? - Age! I think
  
  # Get vaccine coverage and FVPs for all years up to 2030
  scenario_dt = get_scenario_fvps()  # See sceanrios.R
  
  # Calculate deaths averted by multiplying FVPs by impact factor
  scenario_impact = calc_scenario_impact(scenario_dt, impact_dt)
  
  # Summarise deaths averted to annual totals across all diseases and countries
  scenario_total = calc_scenario_total(scenario_impact = scenario_impact) 
  
  # Save both of these outputs to file
  save_file(scenario_impact, o$pth$impact_factors, "scenario_impact")
  save_file(scenario_total,  o$pth$impact_factors, "scenario_total")
}

# ---------------------------------------------------------
# xxxxxxxxx
# Called by: xxxxxx
# ---------------------------------------------------------
calc_impact_factors <- function(dt) {
  
  dt <- merge(dt, d_v_at_table, by = "d_v_at_id")
  
  # Routine deaths averted
  routine_averted <- dt[
    (year - age) %in% 2000:2030 & activity_type == "routine",
    .(total_averted = sum(averted, na.rm = T)), by = .(country, d_v_at_id)
  ]
  
  # Campaign deaths averted
  campaign_averted <- dt[
    year %in% 2000:2030 & activity_type == "campaign",
    .(total_averted = sum(averted, na.rm = T)), by = .(country, d_v_at_id)
  ]
  
  total_averted <- rbind(routine_averted, campaign_averted)
  total_averted <- merge(total_averted, d_v_at_table, by = "d_v_at_id")
  
  browser() # v_at_table is now v_a_table
  
  # Totals FVPs
  total_fvps <- coverage[year %in% 2000:2030, .(total_fvps = sum(fvps,  na.rm = T)),
                         by = .(country, v_at_id)]
  total_fvps <- merge(total_fvps, v_at_table)
  
  # Combine
  total_dt <- merge(
    total_averted[total_averted > 0], total_fvps,
    by = c("vaccine", "activity_type", "country"),
    all.x = T
  )
  
  browser() # v_at_table is now v_a_table
  
  # Collapse Rubella to combined
  rub_dt <- total_dt[vaccine == "Rubella", .(total_fvps = sum(total_fvps, na.rm = T), total_averted = sum(total_averted, na.rm = T)), by = .(vaccine, country, disease)]
  rub_dt[, activity_type := "combined"]
  rub_dt <- merge(rub_dt, v_at_table, by = c("vaccine", "activity_type"))
  rub_dt <- merge(rub_dt, d_v_at_table, by = c("disease", "vaccine", "activity_type"))
  total_dt <- rbind(total_dt[vaccine != "Rubella"], rub_dt, fill = T)
  
  total_dt[, pred_deaths_averted_rate := total_averted / total_fvps]
  
  return(total_dt[])
}

# ---------------------------------------------------------
# xxxxxxxxx
# Called by: xxxxxx
# ---------------------------------------------------------
rake_impact <- function(impact_factors) {
  
  load_tables("vimc_yov_impact") # Now just vimc_yov
  
  browser()
  
  yov_dt <- unique(vimc_yov_impact[, .(country, vaccine, activity_type, deaths_averted_rate)])
  
  # Take the mean here for the Rubella mishap
  yov_dt <- yov_dt[, lapply(.SD, mean), by = .(country, vaccine, activity_type)]
  
  impact_dt <- merge(impact_factors, yov_dt,
                     by = c("country", "vaccine", "activity_type"), all.x = T)
  impact_dt[, raking_factor := deaths_averted_rate / pred_deaths_averted_rate]
  raking_summary <- summarize_raking(impact_dt)
  impact_dt <- merge(impact_dt, raking_summary[, .(vaccine, activity_type, median)],
                     by = c("vaccine", "activity_type"), all.x = T)
  overall_routine_median <- median(impact_dt[activity_type == "routine"]$raking_factor, na.rm = T)
  impact_dt[is.na(median), median := overall_routine_median]
  impact_dt[is.na(deaths_averted_rate), raked_deaths_averted_rate := pred_deaths_averted_rate * median]
  impact_dt[, impact_factor := ifelse(is.na(deaths_averted_rate), raked_deaths_averted_rate, deaths_averted_rate)]
  
}

# ---------------------------------------------------------
# Calculate deaths averted by multiplying FVPs by impact factor
# Called by: run_impact_factors(), calc_scenario_total()
# ---------------------------------------------------------
calc_scenario_impact <- function(scenario_dt, impact_dt) {
  
  message(" - Calculating scenario impact")

  browser() # v_at_table is now v_a_table
  
  # Calculate deaths averted 
  scenario_impact = scenario_dt %>%
    filter(fvps > 0) %>%
    # Get vaccine-activity details...
    # TODO: Feel this wouldn't be necessary without the 'Rubella hiccup'...
    left_join(y  = v_at_table,
              by = "v_at_id") %>%
    mutate(activity_type = ifelse(vaccine == "Rubella", "combined", activity_type)) %>%  # TODO: This *shouldn't* be needed
    # Join with pre-calculated impact factors...
    inner_join(y  = impact_dt, 
               by = c("vaccine", "activity_type", "country")) %>%
    select(disease, vaccine, activity_type, country, impact_factor, year, age, fvps) %>%
    # Now a simple calculation for deaths averted...
    mutate(deaths_averted = fvps * impact_factor) %>%
    arrange(disease, vaccine, country, year, age)
  
  # Plotting impact_dt...
  plot_dt1 = impact_dt %>%
    group_by(disease, vaccine, activity_type) %>%
    summarise(deaths_averted = sum(total_averted)) %>%
    ungroup() %>%
    unite("d_v_at", disease, vaccine, activity_type) %>%
    as.data.table()

  g1 = ggplot(plot_dt1, aes(x = d_v_at, y = deaths_averted, fill = d_v_at)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    theme(axis.text.x = element_text(angle = 50, hjust = 1))

  # Plotting scenario_impact...
  plot_dt2 = scenario_impact %>%
    group_by(disease, vaccine, activity_type) %>%
    summarise(deaths_averted = sum(deaths_averted)) %>%
    ungroup() %>%
    unite("d_v_at", disease, vaccine, activity_type) %>%
    as.data.table()

  g2 = ggplot(plot_dt2, aes(x = d_v_at, y = deaths_averted, fill = d_v_at)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    theme(axis.text.x = element_text(angle = 50, hjust = 1))

  # Sanity check: no NA deaths averted
  if (any(is.na(scenario_impact$deaths_averted)))
    stop("NA deaths averted calculated")
  
  # Sanity check: all values are positive
  if (any(scenario_impact$deaths_averted <= 0))
    stop("Zero or negative deaths averted calculated")

  return(scenario_impact)
}

# ---------------------------------------------------------
# Summarise total deaths averted per year
# Called by: run_impact_factors()
# ---------------------------------------------------------
calc_scenario_total = function(scenario_impact = NULL, ...) {
  
  # If scenario_impact is not provided, calculated it using mandatory additional arguments
  if (is.null(scenario_impact))
    scenario_impact = calc_scenario_impact(...)
  
  # Total number of deaths averted per year
  scenario_total = scenario_impact %>%
    group_by(year) %>%
    summarise(deaths_averted = sum(deaths_averted)) %>%
    ungroup() %>%
    arrange(year) %>%
    as.data.table()
  
  return(scenario_total)
}

