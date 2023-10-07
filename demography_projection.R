###########################################################
# POPULATION PROJECTION
#
# xxx
#
###########################################################

#' Execute cohort-component projection model
#' @param nx Population matrix
#' @param sx Survival matrix
#' @param fx Fertility matrix
#' @param mig Migration matrix
#' @return A list with population and births matrices
#' @export
get_ccpm <- function(nx, sx, fx, mig, n) {
  z <- dim(sx)[1] / 2
  n <- dim(sx)[2]
  srb <- 1.05
  bx <- array(dim = n)
  p <- array(dim = c(2 * z, n + 1))
  p[, 1] <- nx[, 1]
  for (i in 1:n) {
    next_pop <- update_pop(p[, i], sx[, i], fx[, i], mig[, i], z, srb)
    p[, i + 1] <- next_pop[[1]]
    bx[i] <- next_pop[[2]]
  }
  
  ccpm <- list(population = p[, 2:(n + 1)], births = bx)
  
  return(ccpm)
}

# ---------------------------------------------------------
# xxx
# Called by: xxx
# ---------------------------------------------------------
update_pop <- function(p_in, sx, fx, mig, z, srb) {
  p_out <- vector()
  bx_out <- 0
  for (sex in c("female", "male")) {
    if (sex == "male") {
      first_idx <- 1
      mid_idx <- 2:(z - 1)
      final_idx <- z
    } else {
      first_idx <- z + 1
      mid_idx <- z + 2:(z - 1)
      final_idx <- z + z
    }
    
    # Middle age groups
    p_out[mid_idx] <- p_in[mid_idx - 1] * (
      sx[mid_idx - 1] * (1 + .5 * mig[mid_idx]) +
        .5 * mig[mid_idx]
    )
    
    # Final age group
    p_out[final_idx] <- p_in[final_idx - 1] * (
      sx[final_idx - 1] * (1 + .5 * mig[final_idx]) +
        .5 * mig[final_idx]
    ) +
      p_in[final_idx] * (
        sx[final_idx] * (1 + .5 * mig[final_idx]) +
          .5 * mig[final_idx]
      )
    
    # First age group
    bx <- calc_births(srb, fx, sx, p_in, mig, z, sex)
    bx_out <- bx_out + bx
    p_out[first_idx] <- bx * (sx[first_idx] * (1 + .5 * mig[first_idx]) +
                                .5 * mig[first_idx])
  }
  
  
  return(list(p_out, bx_out))
}

# ---------------------------------------------------------
# xxx
# Called by: xxx
# ---------------------------------------------------------
calc_births <- function(srb, fx, sx, p_in, mig, z, sex = "female") {
  # Calculate births
  fxb <- ((1 + srb) ^ (-1) * (
    fx[z + 10:54] +
      fx[z + 11:55] * sx[z + 11:55]
  ) * 0.5)
  
  if (sex == "male") {
    fxb <- fxb * srb
  }
  
  bx <- sum(
    fxb * p_in[z + 10:54] * (1 + .5 * mig[z + 10:54])
  )
  
  return(bx)
}


# TODO: This could be more clear
#' Not sure why this is needed
get_person <- function(x) {
  x[1:96] + x[-c(1:96)]
}

#' Calculate life-table values from mortality rates
#' @param y Vector of mortality rates
#' @return  A vector of life-expectancy, 1q0, and 5q0
lt_est <- function(y) {
  px <- exp(-y)
  lx <- c(1, cumprod(px))
  ex <- round(sum(head(lx, -1) + tail(lx, -1)) / 2, 3)
  q5 <- 1000 * (1 - lx[6])
  q1 <- 1000 * (1 - lx[2])
  lt <- c(ex, q1, q5)
  
  return(lt)
}

#' Run the population projection for a particular country and year range
#' @param iso Country ISO3
#' @param y0 Start year of projection
#' @param y1 End year of projection
#' @param wpp_input Input WPP data
#' @param scen Name of the scenario being run
#' @return A list of tables with population projection results
project_pop <- function(iso, y0, y1, wpp_input, scen = "Default") {
  
  wpp_ina <- wpp_input  %>%
    filter(country == iso, 
           year %in% (y0 - 1) : y1) %>%
    arrange(sex_id, age, year)
  
  nx <- wpp_ina %>%
    select(sex_id, age, year, nx) %>%
    mutate(nx = ifelse(nx == 0, 1e-09, nx)) %>%
    spread(year, nx) %>%
    select(-sex_id, -age) %>%
    as.matrix()
  
  wpp_ina <- wpp_ina %>%
    filter(year %in% y0:y1)
  
  fx <- wpp_ina %>%
    select(sex_id, age, year, fx) %>%
    spread(year, fx) %>%
    select(-sex_id, -age) %>%
    as.matrix()
  
  mig <- wpp_ina %>%
    select(sex_id, age, year, mig) %>%
    spread(year, mig) %>%
    select(-sex_id, -age) %>%
    as.matrix()
  
  mx <- wpp_ina %>%
    select(sex_id, age, year, mx) %>%
    mutate(mx = ifelse(mx == 0, 1e-09, mx)) %>%
    spread(year, mx) %>%
    select(-sex_id, -age) %>%
    as.matrix()
  
  n <- y1 - y0 + 1
  
  if (scen != "Default")
    mx <- get_mx_scen(iso, y0, y1, scen, mx, nx[, 2 : n])
  
  sx <- exp(-mx)
  
  ccpm <- get_ccpm(nx, sx, fx, mig)
  
  deaths = -1 * log(sx) * ccpm$population[, 1 : n]
  
  output = c(ccpm, list(deaths = deaths, mx = mx))
  
  return(output)
}

# ---------------------------------------------------------
# xxx
# Called by: xxx
# ---------------------------------------------------------
add_lt <- function(projected_pop, is, y0, y1) {
  population <- projected_pop$population
  deaths <- projected_pop$deaths
  births <- projected_pop$births
  mx <- projected_pop$mx
  
  pop_tot <- apply(population, 2, get_person)
  deaths_tot <- apply(deaths, 2, get_person)
  mx_both  <- deaths_tot / pop_tot
  
  lt_out_both <- apply(mx_both, 2, lt_est)
  lt_out_fmle <- apply(mx[1:96, ], 2, lt_est)
  lt_out_mle <- apply(mx[97:192, ], 2, lt_est)
  
  deaths_both <- apply(deaths, 2, sum)
  deaths_male <- apply(deaths[97:192, ], 2, sum)
  deaths_female <- apply(deaths[1:96, ], 2, sum)
  
  e0 <- lt_out_both[1, ]
  e0_male <- lt_out_mle[1, ]
  e0_female <- lt_out_fmle[1, ]
  
  imr <- lt_out_both[2, ]
  u5mr <- lt_out_both[3, ]
  
  locs <- country_table %>%
    filter(country_name == is) %>%
    select(country, country_name)
  out_df <- data.table(
    year = y0:y1,
    deaths_both = deaths_both,
    deaths_male = deaths_male,
    deaths_female = deaths_female,
    e0 = e0,
    e0_male = e0_male,
    e0_female = e0_female,
    imr = imr, u5mr = u5mr, births = births, locs
  )
  
  return(out_df)
}

# ---------------------------------------------------------
# Calculate all single-year deaths
# ---------------------------------------------------------
get_all_deaths = function(y0, y1, wpp_input) {
  
  # INPUTS:
  #  y0 := Start year of projection
  #  y1 := End year of projection
  
  # OUTPUT:
  #  A datatable with single-year deaths
  
  # Initiate list to store deaths for each country
  deaths_list = list()
  
  # All countries to calculate this for
  countries = unique(wpp_input$country)
  
  # Iterate through countries
  for (country in countries) {
    
    # Use cohort-component projection model to project deaths
    ccpm = project_pop(country, y0, y1, wpp_input)
    
    # Melt down death projections to tidy format
    deaths_list[[country]] = ccpm$deaths %>%
      as_named_dt(paste0(y0 : y1)) %>%
      mutate(age    = rep(0 : 95, times = 2), 
             sex_id = rep(1 : 2,  each  = 96)) %>%
      gather(year, deaths, -age, -sex_id) %>%
      mutate(year = as.numeric(year),
             country = !!country) %>%
      arrange(sex_id, age, year) %>%
      as.data.table()
  }
  
  # Bind all country results into single datatable
  deaths_dt = rbindlist(deaths_list) %>%
    mutate(deaths = ifelse(is.finite(deaths), deaths, NA), 
           deaths = pmax(deaths, 0))
  
  return(deaths_dt)
}

