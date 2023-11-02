###########################################################
# SCENARIOS
#
# xxxxxxx
#
###########################################################

# ---------------------------------------------------------
# Use 2030 coverage targets for DTP vaccine to project future coverage by country
# ---------------------------------------------------------
gen_ia2030_goals <- function(linear = T, no_covid_effect = 2022, intro_year = 2025, intro_range = T) {
  
  # Load 2019 coverage
  load_tables("coverage")
  cov_dt <- coverage[year == 2019]
  
  browser() # v_at_table is now v_a_table
  
  # Iterate through each vaccine (except HPV which is special)
  vaccs <- setdiff(unique(v_at_table[activity_type == "routine"]$v_at_id), 2)
  dt <- rbindlist(lapply(vaccs, function(v) {
    v_dt <- cov_dt[v_at_id == v][order(country)]
    if (length(unique(v_dt$gender)) > 1) {
      v_dt <- v_dt[gender == 2]
    }
    missing_locs <- setdiff(country_table$country, v_dt$country)
    missing_dt <- data.table(country = missing_locs, coverage = 0,
                             v_at_id = v, year = 2019, age = 0, gender = unique(v_dt$gender),
                             fvps = 0)
    v_dt <- rbind(v_dt, missing_dt, use.names = T)
    v_dt <- v_dt[order(country)]
    # Repeat out to 'no covid effect' year
    n_covid <- no_covid_effect - 2019
    covid_mat <- matrix(rep(v_dt$coverage, n_covid), ncol = n_covid)
    
    # Increase to goal with zeros delayed to intro_year
    n_increase <- 2030 - no_covid_effect + 1
    zero_idx <- which(v_dt$coverage == 0)
    # Change intro_year to 2024 for YF
    if (v == 19) {
      temp_intro_year <- intro_year - 1
    } else {
      temp_intro_year <- intro_year
    }
    if (intro_range) {
      # Range of intro years split up by quintile of coverage goal
      zero_locs <- v_dt[zero_idx]$country
      
      browser() # ia2030_dtp_goal now part of country_table, but needs year = 2030 to be appended
      
      # 2030 coverage targets for DTP vaccine...
      ordered_locs <- ia2030_dtp_goal[country %in% zero_locs][rev(order(value))]$country
      split_locs <- split(ordered_locs, floor(5 * seq.int(0, length(ordered_locs) - 1) / length(ordered_locs)))
      names(split_locs) <- (-2:2 + temp_intro_year)[1:length(split_locs)]
    } else {
      zero_n <- 2030 - temp_intro_year + 1
    }
    
    browser() # ia2030_dtp_goal now part of country_table, but needs year = 2030 to be appended
    
    setnames(v_dt, "coverage", "current")
    roc_dt <- merge(
      v_dt[, .(country, current)],
      ia2030_dtp_goal[, .(country, value)],
      by = "country"
    )
    roc_dt[, n := n_increase]
    t_mat <- matrix(
      1:n_increase, byrow = T, ncol = n_increase, 
      nrow = nrow(roc_dt)
    )
    if (length(zero_idx > 0)) {
      if (intro_range) {
        for (i in zero_idx) {
          i_year <- as.integer(names(split_locs)[unlist(
            lapply(split_locs, function(s) {
              v_dt[i,]$country %in% s
            })
          )])
          zero_n <- 2030 - i_year + 1
          t_mat[i, ] <- c(rep(0, n_increase - zero_n), 1:zero_n)
          roc_dt[i, n := zero_n]
        }
      } else {
        t_mat[zero_idx,] <- matrix(
          c(rep(0, n_increase - zero_n), 1:zero_n),
          nrow = length(zero_idx), ncol = n_increase, byrow = T)
        roc_dt[zero_idx, n := zero_n]
      }
    }
    # Handle regionally-specific vaccines
    if (v %in% c(19, 12, 7)) {
      if (v == 19) {
        reg_locs <- country_table[yf == 1]$country
        # Remove Argentina and Kenya
        reg_locs <- setdiff(reg_locs, c(7, 90))
      } else if (v == 12){
        reg_locs <- country_table[mena == 1]$country
      } else if (v == 7) {
        reg_locs <- country_table[je == 1]$country
        # Remove Russia and Pakistan
        reg_locs <- setdiff(reg_locs, c(131, 144))
      }
      reg_idx <- which(v_dt$country %in% reg_locs)
      roc_dt[!reg_idx, value := current]
    }
    # Do no introduce HepB or BCG in any countries
    if (v %in% c(3, 21)) {
      roc_dt[zero_idx, value := current]
    }
    # Hold medium/low BCG coverage constant in France, Ireland, Sweden
    if (v == 21) {
      hold_locs <- c(63, 83, 168)
      hold_idx <- which(v_dt$country %in% hold_locs)
      roc_dt[hold_idx, value := current]
    }
    # Keep current JE levels in countries with only subnational endemicity
    # Indonesia, India, and Malaysia
    if (v == 7) {
      hold_locs <- c(79, 80, 104)
      hold_idx <- which(v_dt$country %in% hold_locs)
      roc_dt[hold_idx, value := current]
    }
    # Linear vs. non-linear
    if (linear) {
      roc_dt[, roc := ((value - current) / n)]
      roc_dt[roc < 0, roc := 0]
      inc_mat <- v_dt$current + roc_dt$roc * t_mat
    } else {
      roc_dt[, roc := log((1 - value) /(1 - current)) /  n]
      roc_dt[roc > 0, roc := 0]
      inc_mat <- 1 - (1 - v_dt$current) * exp(roc_dt$roc * t_mat)
    }
    # Combine and convert to data.table
    c_mat <- cbind(covid_mat, inc_mat)
    colnames(c_mat) <- 2019:2030
    # matplot(t(c_mat), type = "l")
    dt <- cbind(v_dt[, -c("year", "current", "fvps"), with = F], c_mat)
    melt_dt <- melt(dt,
                    id.vars = c("country", "v_at_id", "age", "gender"),
                    variable.name = "year"
    )
    return(melt_dt)
  }))
  
  ## Add HPV
  dt <- rbind(dt, hpv_target, use.names = T)
  
  dt[, year := as.integer(as.character(year))]
  dt <- dt[year > 2019]
  return(dt)
}

# ---------------------------------------------------------
# Get vaccine coverage and FVPs for all years up to 2030
# ---------------------------------------------------------
get_scenario_fvps = function() {
  
  # Past coverage and FVPs comes directly from the data
  past_fvps = coverage %>%
    filter(year >= 2000, 
           year <  2020)  # TODO: Should this be more flexible? Depends on gen_ia2030_goals I guess
  
  # Generate target coverage using seperate function
  #
  # TODO: Move this coverage renaming to occur inside of gen_ia2030_goals()
  future_coverage = gen_ia2030_goals(linear = FALSE) %>%
    rename(coverage = value)
  
  # Convert target coverage to FVPs using pop size projections
  future_fvps = cov2fvp(future_coverage)
  
  # Bind together
  scenario_dt = rbind(past_fvps, future_fvps)
  
  return(scenario_dt)
}


