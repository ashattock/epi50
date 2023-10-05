###########################################################
# TOTAL COVERAGE
#
# Determine effective vaccine coverage for each cohort per year, 
# possible over different vaccine activities (eg routine & campaign).
#
###########################################################

# ---------------------------------------------------------
#' Calculate total coverage for every vaccine, country, sex
# Called by: prep_rr()
# ---------------------------------------------------------
total_coverage2 <- function(coverage) {
  
  # TODO: Apply with dtplyr to achieve full effect...
  
  # Coverage data of all activities for this vaccine
  all_data = coverage %>%
    left_join(v_at_table, by = "v_at_id") %>%
    select(-v_at_id) %>%
    arrange(country, sex_id, year, age)
  
  # Sanity check that only one vaccine type has been provided
  if (length(unique(all_data$vaccine)) > 1)
    stop("This function handles only one vaccine type at a time")
  
  # All unique countries and sex combinations for this vaccine
  unique_data = all_data %>%
    select(country, sex_id, vaccine) %>%
    unique()
  
  # Initiate list of total coverage results
  cov_list = list()
  
  # Loop through unique country-sex combinations
  for (i in seq_len(nrow(unique_data))) {
    
    # Shorthand for this activity
    this_data = unique_data[i, ]
    
    # Determine effective vaccine coverage for each cohort per year
    #
    # NOTE: Vaccine has already been filtered by parent function
    cov_list[[i]] = all_data %>%
      filter(country %in% this_data$country, 
             sex_id  %in% this_data$sex_id) %>%  # Data for this country and sex
      group_by(year, age, activity_type) %>%
      summarise(coverage = max(coverage)) %>%  # Max of any duplicate campaigns
      ungroup() %>%
      as.data.table() %>%
      calc_total_cov2() %>%  # Determines only non-trivial values in year x age matrix
      cbind(this_data)
  }
  
  # All possible combinations including year and age
  #
  # NOTE: Sparse result, so most of these age-year pairs will have zero value
  full_dt = unique_data %>%
    expand_grid(year = o$data_years, 
                age  = o$data_ages)
  
  # Join total coverage results to full combination table
  total_dt = rbindlist(cov_list) %>%
    full_join(y  = full_dt, 
              by = names(full_dt)) %>%
    mutate(value = ifelse(is.na(value), 0, value)) %>%
    arrange(country, sex_id, year, age)
  
  # TODO: Set a cap on BCG effect at age 15
  
  return(total_dt)
}

# ---------------------------------------------------------
# Determine non-trivial year-age vaccine effects by cohort
# Called by: total_coverage()
# ---------------------------------------------------------
calc_total_cov2 <- function(dt) {
  
  # Initiate list
  cov_list = list()
  
  # Loop through the different activities
  for (i in seq_len(nrow(dt))) {
    
    # Indicies for years and ages
    year_idx = (1 + match(dt[i, year], o$data_years)) : length(o$data_years)
    age_idx  = (1 + match(dt[i, age],  o$data_ages))  : length(o$data_ages)
    
    # Index upto only the smallest of these two vectors
    vec_idx = 1 : min(length(year_idx), length(age_idx))
    
    # These form the only non-trivial entries
    cov_list[[i]] = 
      data.table(year = o$data_years[year_idx[vec_idx]],
                 age  = o$data_ages[age_idx[vec_idx]], 
                 value   = dt[i, coverage], 
                 actvity = dt[i, activity_type])
  }
  
  # Calculate coverage for any identicial year and age activities from different activities
  #
  # NOTE: Assuming independence of 'routine' and 'campaign' hence the 1 - prod(1 - value)
  cov_dt = rbindlist(cov_list) %>%
    group_by(year, age) %>%
    summarise(value = 1 - prod(1 - value)) %>%
    ungroup() %>%
    as.data.table()
  
  return(cov_dt)
}

# ---------------------------------------------------------
# v1 function
# ---------------------------------------------------------
total_coverage <- function(coverage) {
  cov_dt <- merge(coverage, v_at_table)
  total_coverage <- rbindlist(lapply(unique(cov_dt$vaccine), function(v) {
    v_dt <- cov_dt[vaccine == v]
    vacc_dt <- rbindlist(lapply(unique(v_dt$activity_type), function(a) {
      a_dt <- v_dt[activity_type == a]
      act_dt <- rbindlist(lapply(unique(a_dt$country), function(l) {
        l_dt <- v_dt[country == l]
        loc_dt <- rbindlist(lapply(unique(l_dt$sex_id), function(s) {
          dt <- l_dt[sex_id == s]
          total_dt <- calc_total_cov(dt)
          total_dt[, sex_id := s]
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
  load_tables("wpp_input")
  
  # Total number of people per country (both sexes combined)
  #
  # NOTE: nx := number of people
  both_dt = wpp_input[, .(nx = sum(nx)), .(country, age, year)]
  both_dt[, sex_id := 3]
  
  # Combine so we have both genders seperate and combined
  pop_dt = rbind(wpp_input, both_dt, fill = T)
  
  # Join with coverage details
  fvp_dt = merge(coverage_dt, pop_dt[, .(country, age, sex_id, year, nx)])
  
  # Then just a simple calculation for FVPs
  fvp_dt[, fvps := coverage * nx]
  fvp_dt[, nx := NULL]
  
  return(fvp_dt)
}

