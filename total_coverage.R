###########################################################
# TOTAL COVERAGE
#
# Determine effective vaccine coverage for each cohort per year, 
# possible over different vaccine activities (eg routine & campaign).
#
###########################################################

# ---------------------------------------------------------
# Calculate total coverage for every vaccine, country, gender
# ---------------------------------------------------------
total_coverage2 <- function(coverage_dt) {
  
  # TODO: Allow each d_v_a to be 'targeted' or 'non-targeted'
  
  # Create full combination table
  #
  # NOTE: Final result is sparse => most age-year values will be zero
  full_dt = expand_grid(
    country = unique(coverage_dt$country),
    year    = o$data_years, 
    age     = o$data_ages) %>%
    as.data.table()
  
  # Function to extract total coverage for each data point
  total_coverage_fn = function(i) {
    
    # Data relating to this row of coverage_dt
    data = coverage_dt[i, ]
    
    # Indicies for years and ages
    year_idx = match(data$year, o$data_years) : length(o$data_years)
    age_idx  = match(data$age,  o$data_ages)  : length(o$data_ages)
    
    # Index upto only the smallest of these two vectors
    vec_idx = 1 : min(length(year_idx), length(age_idx))
    
    # These form the only non-trivial entries
    total_dt = data.table(
      country = data$country, 
      year    = o$data_years[year_idx[vec_idx]],
      age     = o$data_ages[age_idx[vec_idx]], 
      value   = data$coverage)
    
    return(total_dt)
  }

  # Coverage data values to work through  
  total_idx = seq_len(nrow(coverage_dt))
  
  # Apply total coverage function to each row
  total_dt = lapply(total_idx, total_coverage_fn) %>%
    rbindlist() %>%
    # First summarise for cumulative total coverage...
    group_by(country, year, age) %>%
    # summarise(value = 1 - prod(1 - value)) %>%  # Assumes non-targeted vaccination
    summarise(value = min(sum(value), 1)) %>%   # Assumes targeted vaccination
    ungroup() %>%
    # Then join with full grid...
    full_join(y  = full_dt, 
              by = names(full_dt)) %>%
    replace_na(list(value = 0)) %>%
    arrange(country, year, age) %>%
    # Append d_v_a info...
    mutate(d_v_a_id = unique(coverage_dt$v_a_id), 
           .after = 1) %>%
    as.data.table()
  
  # TODO: Set a cap on BCG effect at age 15
  
  return(total_dt)
}

# ---------------------------------------------------------
# Determine non-trivial year-age vaccine effects by cohort
# ---------------------------------------------------------
calc_total_cov2 <- function(dt) {
  

  
  return(cov_dt)
}

# ---------------------------------------------------------
# v1 function
# ---------------------------------------------------------
total_coverage <- function(coverage) {
  
  browser() # v_at_table is now v_a_table
  
  cov_dt <- merge(coverage, v_at_table)
  total_coverage <- rbindlist(lapply(unique(cov_dt$vaccine), function(v) {
    v_dt <- cov_dt[vaccine == v]
    vacc_dt <- rbindlist(lapply(unique(v_dt$activity_type), function(a) {
      a_dt <- v_dt[activity_type == a]
      act_dt <- rbindlist(lapply(unique(a_dt$country), function(l) {
        l_dt <- v_dt[country == l]
        loc_dt <- rbindlist(lapply(unique(l_dt$gender), function(s) {
          dt <- l_dt[gender == s]
          total_dt <- calc_total_cov(dt)
          total_dt[, gender := s]
        }))
        loc_dt[, country := l]
        return(loc_dt)
      }))
      act_dt[, activity_type := a]
      return(act_dt)
    }))
    vacc_dt <- vacc_dt[, vaccine := v]
    # Set a cap on BCG effect at age 15
    if (v == "BCG") {
      vacc_dt[age >= 15, value := 0]
    }
    return(vacc_dt[])
  }))
  return(total_coverage[])
}

# ---------------------------------------------------------
# v1 function
# ---------------------------------------------------------
calc_total_cov <- function(dt) {
  ages <- 0:95
  years <- 2000:2039
  activity_types <- unique(dt$activity_type)
  if (length(activity_types) != 1)
    browser()
  full_dt <- CJ(age = ages, year = years, activity_type = activity_types)
  merge_dt <- merge(
    full_dt,
    dt[, .(activity_type, age, year, coverage)],
    by = c("activity_type", "age", "year"), all.x = T
  )
  merge_dt[is.na(coverage), coverage := 0]
  nrow <- length(ages)
  ncol <- length(years)
  mat <- matrix(0, nrow, ncol)
  rownames(mat) <- ages
  colnames(mat) <- years
  if ("routine" %in% activity_types){
    r_mat <- as.matrix(
      dcast(
        merge_dt[activity_type == "routine"],
        age ~ year, value.var = "coverage"
      )
    )[, -1, drop = F]
    # Routine -- take the max
    for (i in (ncol - 1):1) {
      vec <- r_mat[, i]
      for (j in (i + 1):min(nrow, ncol)) {
        mat[max(1, (j - i + 1)):nrow, j] <- pmax(
          mat[max(1, (j - i + 1)):nrow, j],
          vec[1:min(nrow, (nrow - j + i))]
        )
      }
    }
  }
  if ("campaign" %in% activity_types) {
    c_mat <- as.matrix(
      dcast(
        merge_dt[activity_type == "campaign"],
        age ~ year, value.var = "coverage"
      )
    )[, -1, drop = F]
    # Campaign - Assume independence
    for (i in (ncol - 1):1) {
      vec <- c_mat[, i]
      for (j in (i + 1):min(nrow, ncol)) {
        mat[max(1, (j - i + 1)):nrow, j] <- 1 -
          (1 - mat[max(1, (j - i + 1)):nrow, j]) *
          (1 - vec[1:min(nrow, (nrow - j + i))])
      }
    }
  }
  total_dt <- melt(
    as.data.table(cbind(age = 0:95, mat)),
    id.vars = "age", variable.name = "year"
  )
  total_dt[, year := as.integer(as.character(year))]
  
  return(total_dt[])
}

# ---------------------------------------------------------
# Calculate FVPs using coverage and demogrpahic data
# Called by: get_scenario_fvps() and similar FVPs-generating functions
# ---------------------------------------------------------
cov2fvp = function(coverage_dt) {
  
  # Load demographic data
  wpp_input = table("wpp_input")
  
  # Total number of people per country (both genders combined)
  #
  # NOTE: nx := number of people
  both_dt = wpp_input[, .(nx = sum(nx)), .(country, age, year)]
  both_dt[, gender := 3]
  
  # Combine so we have both genders seperate and combined
  pop_dt = rbind(wpp_input, both_dt, fill = T)
  
  # Join with coverage details
  fvp_dt = merge(coverage_dt, pop_dt[, .(country, age, gender, year, nx)])
  
  # Then just a simple calculation for FVPs
  fvp_dt[, fvps := coverage * nx]
  fvp_dt[, nx := NULL]
  
  return(fvp_dt)
}

