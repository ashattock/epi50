###########################################################
# DEMOGRAPHY
#
# xxxxxx
#
###########################################################

# ---------------------------------------------------------
# Prepare demography-related estimates from WPP
# ---------------------------------------------------------
prepare_demography = function() {
  
  message(" - Demography data")
  
  # TEMP: Load table from IA2030 database
  readRDS("temp/wpp_input.rds")$db_dt %>%
    mutate(gender = c("m", "f", "b")[sex_id]) %>%
    select(country, year, gender, age, nx, mx, fx, mig) %>%
    save_table("wpp_input")

  # TEMP: Load table from IA2030 database
  readRDS("temp/all_deaths.rds")$db_dt %>%
    # Deal with 'NA' deaths being stored as character...
    mutate(deaths = ifelse(deaths == "NA", 0, deaths), 
           deaths = as.numeric(deaths)) %>%
    # Summarise over gender...
    group_by(country, year, age) %>%
    summarise(deaths_allcause = sum(deaths)) %>%
    ungroup() %>%
    arrange(country, year, age) %>%
    as.data.table() %>%
    save_table("deaths_allcause")

  return()
  
  # Load historical and projected pop sizes from WPP for both genders
  wpp_dt = load_wpp_data() %>%
    append_mx() %>%  # Append mortality
    append_fx()      # Append fertility
  
  browser()
  
  # ---- Append migration details ----
  
  # loop through pulling migration
  isc <- sort(unique(wpp_dt$country))
  isn <- length(isc)
  wpp_in_list <- list(isn)
  
  for (c in 1:isn) {
    is <- isc[c]
    wpp_ina <- wpp_dt  %>%
      filter(country == is) %>%
      arrange(gender, age, year)
    fx <- wpp_ina %>%
      select(gender, age, year, fx) %>%
      tidyr::spread(year, fx) %>%
      select(-c(gender, age)) %>%
      as.matrix()
    nx <- wpp_ina %>%
      select(gender, age, year, nx) %>%
      tidyr::spread(year, nx) %>%
      select(-c(gender, age)) %>%
      as.matrix()
    nx[nx == 0] <- 1e-09
    mx <- wpp_ina %>%
      select(gender, age, year, mx) %>%
      tidyr::spread(year, mx) %>%
      select(-c(gender, age)) %>%
      as.matrix()
    mx[mx == 0] <- 1e-09
    sx <- exp(-mx)
    z <- length(0:95)
    n <- ncol(nx) - 1
    
    migs <- get_mig(is, nx, sx, fx, z) %>% as.data.table()
    
    yrv  <- paste0(1979 + 1:n)
    sxv  <- rep(c(2, 1), each = 96)
    agv  <- c(0:95, 0:95)
    
    colnames(migs) <- yrv
    
    migs <- migs %>%
      mutate(age = agv, gender = sxv) %>%
      tidyr::gather(year, mig, -age, -gender) %>%
      mutate(year = as.numeric(year)) %>%
      arrange(gender, age, year)
    
    wpp_in_list[[c]] <- right_join(
      wpp_ina,
      migs,
      by = c("year", "gender", "age")
    ) %>%
      select(country, gender, age, year, nx, mx, fx, mig)
  }
  
  # ---- Prepare outputs ----
  
  browser()
  
  # Bind into single datatable
  wpp_input = rbindlist(wpp_in_list) %>%
    mutate(across(.cols = c(age, year, nx), 
                  .fns  = as.integer)) %>%
    arrange(country, year, age, gender)
  
  # Through an error if any NA entries for number of people
  if (any(is.na(wpp_input$nx)))
    stop("NA nx entries")
  
  # Project population estimates and predict deaths
  all_deaths = get_all_deaths(2000, 2095, wpp_input) %>%
    select(country, year, age, gender, deaths) %>%
    arrange(country, year, age, gender)
  
  browser()
  
  # Save objects as tables in cache
  save_table(wpp_input,  "wpp_input")
  save_table(all_deaths, "deaths_allcause")
}

# ---------------------------------------------------------
# Load population data from large WPP data files
# ---------------------------------------------------------
load_wpp_data = function() {
  
  message("  > Load WWP data")
  
  # Initiate list to store results
  wpp_list = list()
  
  # Iterate through gender and timeframe
  for (time in c("past", "future")) {
    for (gender in c("male", "female")) {
      
      # Path to streamlined WPP pop estimates - by gender and time frame
      wpp_file = paste0(o$pth$input, "wpp19_", gender, "_", time, ".csv")
      
      # Load data and melt to tidy format
      wpp_list[[paste1(time, gender)]] = fread(wpp_file) %>%
        pivot_longer(cols     = -c(country, year),
                     names_to = "age") %>%
        mutate(age    = as.integer(age),
               gender = substr(gender, 1, 1), 
               nx     = value * 1000) %>%
        select(country, year, gender, age, nx) %>%
        arrange(country, year, age) %>%
        as.data.table()
    }
  }
  
  # Squash into single datatable
  wpp_dt = rbindlist(wpp_list)
  
  return(wpp_dt)
}

# ---------------------------------------------------------
# Append mortality rates to population data
# ---------------------------------------------------------
append_mx = function(wpp_dt) {
  
  message("  > Append mx")
  
  data(mxF, package = "wpp2019")
  data(mxM, package = "wpp2019")
  
  mx_f <- mxF %>%
    gather(year, mx, -country_code, -name, -age) %>%
    mutate(year = as.numeric(substr(year, 1, 4)) + 2.5) %>%
    filter(year > 1970) %>%
    spread(year, mx) %>%
    select(-name) %>%
    mutate(gender = "f")
  
  mx_m <- mxM %>%
    distinct() %>%
    gather(year, mx, -country_code, -name, -age) %>%
    mutate(year = as.numeric(substr(year, 1, 4)) + 2.5) %>%
    filter(year > 1970) %>%
    spread(year, mx) %>%
    select(-name) %>%
    mutate(gender = "m")
  
  wppmx42   <- rbind(mx_f, mx_m)
  
  for (j in seq(1977.5, 2092.5, 5)) {
    for (i in (j + 0.5):(j + 4.5)) {
      k <- j + 5
      eval(
        parse(
          text = paste(
            paste(
              "wppmx42 <- wppmx42 %>% mutate(`",
              i,
              "` = `",
              j,
              "`*exp((i - j)*1/5*log(`",
              k,
              "`/`",
              j,
              "`)))",
              sep = ""
            ),
            collapse = ";"
          )
        )
      )
    }
  }
  
  wppmx42 <- wppmx42 %>%
    select(country_code, gender, age, paste0(1980:2096))
  
  browser()
  
  wppmx_list = list()
  
  for (code in unique(wppmx42$country_code)) {
    
    dpred <- wppmx42 %>%
      filter(country_code == code) %>%
      arrange(gender, age) %>%
      select(-country_code, -gender, -age) %>%
      as.matrix()
    
    dpred[is.nan(dpred) | is.na(dpred)] <- 0.5
    dpred_mat <- apply(dpred, 2, split_rate)
    
    wppmx_list[[code]] = dpred_mat %>%
      as.data.table() %>%
      mutate(country_code = code,
             gender = rep(c(2, 1), each = 101),
             age    = rep(0 : 100, times = 2)) %>%
      pivot_longer(cols = -c(country_code, gender, age), 
                   names_to  = "year", 
                   values_to = "mx") %>%
      mutate(year = as.numeric(year))
  }
  
  wppmx = rbindlist(wppmx_list) %>% 
    rename(wpp_country_code = country_code)
  
  browser()
  
  # g1 = mx_f %>%
  #   filter(country_code == 4) %>%
  #   select(-country_code, -gender) %>%
  #   pivot_longer(cols = -age, 
  #                names_to = "year") %>%
  #   mutate(year = as.numeric(year), 
  #          age  = as.factor(age)) %>%
  #   ggplot(aes(x = year, y = value, colour = age)) +
  #   geom_line()
  
  # g2 = wppmx %>%
  #   filter(wpp_country_code == 4, 
  #          gender == 2) %>%
  #   select(-wpp_country_code, -gender) %>%
  #   mutate(age  = as.factor(age)) %>%
  #   ggplot(aes(x = year, y = mx, colour = age)) +
  #   geom_line()
  
  
  
  
  
  
  wpp_dt %<>% 
    left_join(y  = wppmx,
              by = c("wpp_country_code", "gender", "year", "age")) %>%
    mutate(dx  = mx * nx, 
           age = ifelse(age > 95, 95, age)) %>%  # Set upper age bound
    group_by(wpp_country_code, country, gender, year, age, ) %>%
    summarise(nx = sum(nx, na.rm = T), 
              dx = sum(dx, na.rm = T)) %>%
    ungroup() %>%
    mutate(mx = ifelse(nx == 0, 0, dx / nx))
  
  return(wpp_dt)
}

# ---------------------------------------------------------
# Append fertility rate
# ---------------------------------------------------------
append_fx = function(wpp_dt) {
  
  data("tfr", package = "wpp2019")
  data("tfrprojMed", package = "wpp2019")
  
  tfra <- tfr %>%
    select(-c(last.observed)) %>%
    tidyr::gather(year, tfr, -country_code, -name) %>%
    mutate(year = as.numeric(substr(year, 1, 4)) + 2.5) %>%
    filter(year > 1977) %>%
    tidyr::spread(year, tfr) %>%
    select(-c(name))
  tfrb <- tfrprojMed %>%
    tidyr::gather(year, tfr, -country_code, -name) %>%
    mutate(year = as.numeric(substr(year, 1, 4)) + 2.5) %>%
    filter(year > 1977) %>%
    tidyr::spread(year, tfr) %>%
    select(-c(name))
  tfrall   <- left_join(tfra, tfrb, by = "country_code")
  for (j in seq(1977.5, 2092.5, 5)) {
    for (i in (j + 0.5):(j + 4.5)) {
      k <- j + 5
      eval(
        parse(
          text = paste(
            paste(
              "tfrall <- tfrall %>% mutate(`",
              i,
              "` = `",
              j,
              "`*exp((i - j)*1/5*log(`",
              k,
              "`/`",
              j,
              "`)))",
              sep = ""
            ),
            collapse = ";"
          )
        )
      )
    }
  }
  tfrall <- tfrall %>%
    select(country_code, paste0(1980:2096)) %>%
    tidyr::gather(year, tfr, -country_code) %>%
    mutate(year = as.numeric(year))
  
  data("percentASFR", package = "wpp2019")
  
  agef <- tibble(
    agec = paste0(seq(15, 45, 5), "-", seq(15, 45, 5) + 4),
    age = seq(15, 45, 5)
  )
  asfr <- percentASFR %>%
    tidyr::gather(year, afr, -country_code, -name, -age) %>%
    mutate(year = as.numeric(substr(year, 1, 4)) + 2.5) %>%
    filter(year > 1977) %>%
    tidyr::spread(year, afr)  %>%
    rename(agec = age) %>%
    mutate(agec = paste0(agec)) %>%
    left_join(agef, by = "agec") %>%
    select(-c(name, agec)) %>%
    arrange(country_code, age)
  
  age_asfr <- tibble(agen = rep(seq(15, 45, 5), each = 5), age = 15:49)
  
  for (j in seq(1977.5, 2092.5, 5)) {
    for (i in (j + 0.5):(j + 4.5)) {
      k <- j + 5
      eval(
        parse(
          text = paste(
            paste(
              "asfr <- asfr %>% mutate(`",
              i,
              "` = `",
              j,
              "`*exp((i - j)*1/5*log(`",
              k,
              "`/`",
              j,
              "`)))",
              sep = ""
            ),
            collapse = ";"
          )
        )
      )
    }
  }
  
  asfrall <- asfr %>%
    select(country_code, age, paste0(1980:2096)) %>%
    rename(agen = age) %>%
    left_join(age_asfr, by = "agen") %>%
    select(-c(agen)) %>%
    gather(year, fx, -country_code, -age) %>%
    group_by(country_code, year) %>%
    mutate(sfr = sum(fx, na.rm = T)) %>%
    ungroup() %>%
    mutate(asfr = fx / sfr)  %>%
    select(country_code, age, year, asfr) %>%
    mutate(year = as.numeric(year))
  
  wppfx = tfrall %>%
    left_join(y  = asfrall, 
              by = c("country_code", "year")) %>%
    mutate(fx = asfr * tfr, 
           gender = 2) %>%
    select(wpp_country_code = country_code, 
           gender, age, year, fx) %>%
    arrange(wpp_country_code, age) %>%
    setDT()
  
  browser()
  
  
  
  
  wpp_dt %<>%
    left_join(y  = wppfx,
              by = c("wpp_country_code", "gender", "year", "age")) %>%
    mutate(fx = ifelse(is.na(fx), 0, fx)) %>%
    select(wpp_country_code, country, gender, age, year, nx, mx, fx) %>%
    filter(!is.na(gender) & !is.na(country)) %>%
    arrange(country, gender, age)
  
  return(wpp_dt)
}

# ---------------------------------------------------------
# xxxxxxx
# ---------------------------------------------------------
get_mig = function(is, nx, sx, fx, z) {
  nxf <- nx[1:z, ]
  nxm <- nx[(z + 1):(2 * z), ]
  sxf <- sx[1:z, ]
  sxm <- sx[(z + 1):(2 * z), ]
  
  n <- ncol(nx) - 1
  migm <- migf <- array(dim = c(z, n))
  srb <- 1.05
  # Females
  for (i in 1:n) {
    num <- 2 * (nxf[2:(z - 1), i + 1] - nxf[1:(z - 2), i] * sxf[1:(z - 2), i])
    denom <- (nxf[1:(z - 2), i] * (1 + sxf[1:(z - 2), i]))
    migf[2:(z - 1), i] <- num / denom
    migf[z, i] <- migf[z - 1, i]
    fxbf <- (
      (1 + srb) ^ (-1) *
        (fx[10:54, i] + fx[11:55, i] * sxf[11:55, i]) *
        0.5
    )
    bxfs <- sum(
      fxbf *
        nxf[10:54, i] *
        (1 + .5 * migf[10:54, i])
    )
    migf[1, i] <- 2 * (nxf[1, i + 1] / bxfs - sxf[1, i]) / (1 + sxf[1, i])
  }
  
  # Males
  for (i in 1:n) {
    num <- 2 * (nxm[2:(z - 1), i + 1] - nxm[1:(z - 2), i] * sxm[1:(z - 2), i])
    denom <- (nxm[1:(z - 2), i] * (1 + sxm[1:(z - 2), i]))
    migm[2:(z - 1), i] <- num / denom
    migm[z, i] <- migm[z - 1, i]
    fxbm <- (
      srb * (1 + srb) ^ (-1) *
        (fx[10:54, i] + fx[11:55, i] * sxf[11:55, i]) *
        0.5
    )
    bxms <- sum(
      fxbm *
        nxf[10:54, i] *
        (1 + .5 * migf[10:54, i])
    )
    migm[1, i] <- 2 * (nxm[1, i + 1] / bxms - sxm[1, i]) / (1 + sxm[1, i])
  }
  
  mig <- rbind(migf, migm)
  
  return(mig)
}

# ---------------------------------------------------------
# xxxxxxx
# ---------------------------------------------------------
split_rate = function(mx) {
  pop <- log(mx)
  pop[pop < -13] <- -13
  pop[pop > -0.0001] <- -0.0001
  m1 <- predict(
    pspline::smooth.Pspline(
      c(0, 2, seq(7, 97, 5), 100),
      pop[1:22],
      spar = 0.1
    ),
    0:100
  )
  m2 <- predict(
    pspline::smooth.Pspline(
      c(0, 2, seq(7, 97, 5), 100),
      pop[23:44], spar = 0.1
    ),
    0:100
  )
  pop2 <- c(m1, m2)
  pop2[pop2 > -0.0001] <- -0.0001
  exp(pop2)
}

