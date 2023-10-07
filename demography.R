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
  
  message(" - Demoraphy data")
  
  # ---- Load population from WPP ----
  
  # Load historical and projected pop sizes for both genders
  pm1 = get_wpp_pop("male",   "past")
  pm2 = get_wpp_pop("male",   "future")
  pf1 = get_wpp_pop("female", "past")
  pf2 = get_wpp_pop("female", "future")
  
  browser()
  
  # Concatenate into single datatable
  wpp_in = rbind(pm1, pm2, pf1, pf2) # %>% 
  # left_join(y  = country_table[, .(country, wpp_country_code)], 
  #           by = "wpp_country_code")
  
  browser()
  
  # ---- Append mx ----
  
  # Append mx - what is mx??
  wppmx <- get_mx()
  wpp_in %<>% 
    left_join(y  = wppmx,
              by = c("wpp_country_code", "sex_id", "year", "age")) %>%
    mutate(dx  = mx * nx, 
           age = ifelse(age > 95, 95, age)) %>%  # Set upper age bound
    group_by(wpp_country_code, country, sex_id, year, age, ) %>%
    summarise(nx = sum(nx, na.rm = T), 
              dx = sum(dx, na.rm = T)) %>%
    ungroup() %>%
    mutate(mx = ifelse(nx == 0, 0, dx / nx))
  
  # ---- Append fx ----
  
  # Append fx - what is fx??
  wppfx <- get_fx()
  wpp_in %<>%
    left_join(y  = wppfx,
              by = c("wpp_country_code", "sex_id", "year", "age")) %>%
    mutate(fx = ifelse(is.na(fx), 0, fx)) %>%
    select(wpp_country_code, country, sex_id, age, year, nx, mx, fx) %>%
    filter(!is.na(sex_id) & !is.na(country)) %>%
    arrange(country, sex_id, age)
  
  # ---- Append migration details ----
  
  # loop through pulling migration
  isc <- sort(unique(wpp_in$country))
  isn <- length(isc)
  wpp_in_list <- list(isn)
  
  for (c in 1:isn) {
    is <- isc[c]
    wpp_ina <- wpp_in  %>%
      filter(country == is) %>%
      arrange(sex_id, age, year)
    fx <- wpp_ina %>%
      select(sex_id, age, year, fx) %>%
      tidyr::spread(year, fx) %>%
      select(-c(sex_id, age)) %>%
      as.matrix()
    nx <- wpp_ina %>%
      select(sex_id, age, year, nx) %>%
      tidyr::spread(year, nx) %>%
      select(-c(sex_id, age)) %>%
      as.matrix()
    nx[nx == 0] <- 1e-09
    mx <- wpp_ina %>%
      select(sex_id, age, year, mx) %>%
      tidyr::spread(year, mx) %>%
      select(-c(sex_id, age)) %>%
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
      mutate(age = agv, sex_id = sxv) %>%
      tidyr::gather(year, mig, -age, -sex_id) %>%
      mutate(year = as.numeric(year)) %>%
      arrange(sex_id, age, year)
    
    wpp_in_list[[c]] <- right_join(
      wpp_ina,
      migs,
      by = c("year", "sex_id", "age")
    ) %>%
      select(country, sex_id, age, year, nx, mx, fx, mig)
  }
  
  # ---- Prepare outputs ----
  
  # Bind into single datatable
  wpp_input = rbindlist(wpp_in_list) %>%
    mutate(across(.cols = c(age, year, nx), 
                  .fns  = as.integer)) %>%
    arrange(country, year, age, sex_id)
  
  # Through an error if any NA entries for number of people
  if (any(is.na(wpp_input$nx)))
    stop("NA nx entries")
  
  # Project population estimates and predict deaths
  all_deaths = get_all_deaths(2000, 2095, wpp_input) %>%
    select(country, year, age, sex_id, deaths) %>%
    arrange(country, year, age, sex_id)
  
  # Upload both objects to database
  upload_object(wpp_input,  "wpp_input")
  upload_object(all_deaths, "all_deaths")
}

# ---------------------------------------------------------
# Load population data from large WPP data files
# ---------------------------------------------------------
get_wpp_pop = function(sex, time_frame) {
  
  # # Details for selecting historical data
  # if (time_frame == "past") {
  #   sheet_name = "ESTIMATES"
  #   cell_range = "A17:DE18122"
  #   
  #   # Years to select from
  #   year_from  = 1950
  #   year_until = 2020  # This is the latest year in the ESTIMATES sheet
  # }
  # 
  # # Details for selecting future data
  # if (time_frame == "future") {
  #   sheet_name = "MEDIUM VARIANT"
  #   cell_range = "A17:DE20672"
  #   
  #   # Years to select from
  #   year_from  = 2020  # This is the earliest year in the MEDIUM VARIANT sheet
  #   year_until = 2097
  # }
  # 
  # browser()
  # 
  # # Construct path to specific input file
  # xls_file = paste0(o$pth$input, "wpp19_", sex, ".xlsx")
  
  csv_file = paste0(o$pth$input, "wpp19_", sex, "_", time_frame, ".csv")
  
  # browser() # Append country table and filter
  
  # Load data and melt to tidy format
  # pop_dt = readxl::read_excel(path  = xls_file, 
  #                             sheet = sheet_name, 
  #                             range = cell_range) %>%
  # pop_dt = 
    fread(csv_file) %>%
    select(wpp_country_code = "Country code", 
           year = "Reference date (as of 1 July)",
           all_of(as.character(0 : 100))) %>%
    mutate(across(.cols = all_of(as.character(0 : 100)), 
                  .fns  = function(x) str_remove(as.character(x), " "))) %>%
    inner_join(y  = table("country")[, .(country, wpp_country_code)], 
               by = "wpp_country_code") %>%
    select(country, year, all_of(as.character(0 : 100))) %>%
    write_delim(csv_file, delim = ",")
    
    # pivot_longer(cols     = -c(wpp_country_code, year), 
    #              names_to = "age") %>%
    # filter(year >= year_from, 
    #        year <  year_until) %>%
  #   mutate(age = as.numeric(age), 
  #          sex_id = which(c("male", "female") == sex)) %>%
  #   group_by(wpp_country_code, year, sex_id, age) %>%
  #   summarise(nx = sum(value * 1000)) %>%
  #   ungroup() %>%
  #   setDT()
  # 
  # return(pop_dt)
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

# ---------------------------------------------------------
# xxxxxxx
# ---------------------------------------------------------
get_mx = function() {
  
  data(mxF, package = "wpp2019")
  data(mxM, package = "wpp2019")

  mx_f <- mxF %>%
    gather(year, mx, -country_code, -name, -age) %>%
    mutate(year = as.numeric(substr(year, 1, 4)) + 2.5) %>%
    filter(year > 1977) %>%
    spread(year, mx) %>%
    select(-name) %>%
    mutate(sex_id = 2)
  
  mx_m <- mxM %>%
    distinct() %>%
    gather(year, mx, -country_code, -name, -age) %>%
    mutate(year = as.numeric(substr(year, 1, 4)) + 2.5) %>%
    filter(year > 1977) %>%
    spread(year, mx) %>%
    select(-name) %>%
    mutate(sex_id = 1)
  
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
    select(country_code, sex_id, age, paste0(1980:2096))
  
  wppmx_list = list()
  
  for (code in unique(wppmx42$country_code)) {
    
    dpred <- wppmx42 %>%
      filter(country_code == code) %>%
      arrange(sex_id, age) %>%
      select(-country_code, -sex_id, -age) %>%
      as.matrix()
    
    dpred[is.nan(dpred) | is.na(dpred)] <- 0.5
    dpred_mat <- apply(dpred, 2, split_rate)
    
    wppmx_list[[code]] = dpred_mat %>%
      as.data.table() %>%
      mutate(country_code = code,
             sex_id = rep(c(2, 1), each = 101),
             age    = rep(0 : 100, times = 2)) %>%
      pivot_longer(cols = -c(country_code, sex_id, age), 
                   names_to  = "year", 
                   values_to = "mx") %>%
      mutate(year = as.numeric(year))
  }
  
  wppmx = rbindlist(wppmx_list) %>% 
    rename(wpp_country_code = country_code)

  return(wppmx)
}

# ---------------------------------------------------------
# xxxxxxx
# ---------------------------------------------------------
get_fx = function() {
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
           sex_id = 2) %>%
    select(wpp_country_code = country_code, 
           sex_id, age, year, fx) %>%
    arrange(wpp_country_code, age) %>%
    setDT()

  return(wppfx)
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

