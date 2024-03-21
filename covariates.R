###########################################################
# COVARIATES
#
# Prepare covariates fron various sources for use in regression
# modelling. That is, for geographical imputation, and also
# for inferring drivers of impact.
#
###########################################################

# ---------------------------------------------------------
# Parent functio for prepare all covariate data
# ---------------------------------------------------------
prepare_covariates = function() {
  
  message(" > Regression covariates")
  
  download_gapminder()
  
  c1 = covariates_gapminder()
  
  browser()
  
  c2 = covariates_gapminder()
  c3 = covariates_gbd()
  
  browser()
  
  interpolate_covariates()
  
  # # Expand to single year age bins...
  # complete(country, year = min(year) : max(o$years)) %>%
  # # Interpolate timeseries trends...
  # as_tsibble(index = year, key = country) %>%
  # interp_ts_trend() %>%
  
  
  
  
  
  # ---- Combine and store ----
  
  # TODO: Combine into single datatable
  
  browser()
  
  # Save in tables cache
  save_table(c1_dt, "gapminder_covariates")
  save_table(c2_dt, "unicef_covariates")
  save_table(c3_dt, "gbd_covariates")
}

# ---------------------------------------------------------
# Download all gapminer data from github
# ---------------------------------------------------------
download_gapminder = function() {
  
  gapminder_dict = table("gapminder_dict")
  
  files = paste0(o$pth$gapminder, gapminder_dict$var, ".rds")
  
  # Skip process if all files have already been downloaded
  if (all(file.exists(files)) && !o$force_download_gapminder)
    return()
  
  # ---- Download all files ----
  
  gapminder_url = paste0(
    "https://raw.githubusercontent.com/open-numbers/", 
    "ddf--gapminder--systema_globalis/", 
    "master/countries-etc-datapoints/", 
    "ddf--datapoints--[name]--by--geo--time.csv")
  
  for (i in seq_row(gapminder_dict)) {
    
    data_url = str_replace(
      string  = gapminder_url, 
      pattern = "\\[name\\]", 
      replacement = gapminder_dict$file[i])
    
    data = read_csv(data_url)
    
    # Save rds file locally for easy re-loading
    saveRDS(data, file = files[i])
  }
}

# ---------------------------------------------------------
# xxxxxxx
# ---------------------------------------------------------
covariates_gapminder = function() {
  
  message("  - Gapminder")
  
  # Initiate covariates list
  covariates = list()
  
  # ---- Load raw data ----
  
  gapminder_dict = table("gapminder_dict")
  
  raw_files = paste0(o$pth$gapminder, gapminder_dict$var, ".rds")
  raw_data  = rbindlist(lapply(raw_files, readRDS))
  
  browser()
  
  # Prepare Gapminder data for use as predictors
  
  # Gini coefficient
  gini_dt = fread(paste0(o$pth$input, "ddf--datapoints--gapminder_gini--by--geo--time.csv")) %>%
    rename(gini = gapminder_gini)
  
  # GPD per capita US$ inflation-adjusted
  gdp_dt = fread(paste0(o$pth$input, "ddf--datapoints--gdppercapita_us_inflation_adjusted--by--geo--time.csv")) %>%
    rename(gdp = gdppercapita_us_inflation_adjusted)
  
  # Literacy rate (female) aged 15+ years
  literacy_female_15_dt = fread(paste0(o$pth$input, "ddf--datapoints--literacy_rate_adult_female_percent_of_females_ages_15_above--by--geo--time.csv")) %>%
    rename(lit_female = literacy_rate_adult_female_percent_of_females_ages_15_above)
  
  # Literacy rate (male) aged 15+ years
  literacy_male_15_dt = fread(paste0(o$pth$input, "ddf--datapoints--literacy_rate_adult_male_percent_of_males_ages_15_and_above--by--geo--time.csv")) %>%
    rename(lit_male = literacy_rate_adult_male_percent_of_males_ages_15_and_above)
  
  # Doctors per 1000 population
  doctors_per_1000_dt = fread(paste0(o$pth$input, "ddf--datapoints--medical_doctors_per_1000_people--by--geo--time.csv")) %>%
    rename(doctors_per_1000 = medical_doctors_per_1000_people)
  
  # Underweight children
  underweight_children_dt = fread(paste0(o$pth$input, "ddf--datapoints--underweight_children--by--geo--time.csv")) %>%
    rename(underwt_children = underweight_children)
  
  # Population aged 0 to 14
  pop_0_to_14_dt = fread(paste0(o$pth$input, "ddf--datapoints--population_aged_0_14_years_both_sexes_percent--by--geo--time.csv")) %>%
    rename(pop_0to14 = population_aged_0_14_years_both_sexes_percent)
  
  # Population density
  pop_density_dt = fread(paste0(o$pth$input, "ddf--datapoints--population_density_per_square_km--by--geo--time.csv")) %>%
    rename(pop_density = population_density_per_square_km)
  
  # Urban population (%)
  urban_dt = fread(paste0(o$pth$input, "ddf--datapoints--urban_population_percent_of_total--by--geo--time.csv")) %>%
    rename(urban_percent = urban_population_percent_of_total)
  
  # Rural poverty (% rural people below poverty line)
  rural_poverty_dt = fread(paste0(o$pth$input, "ddf--datapoints--rural_poverty_percent_rural_people_below_national_rural--by--geo--time.csv")) %>%
    rename(rural_poverty = rural_poverty_percent_rural_people_below_national_rural)
  
  # HIV prevalence
  hiv_prev_dt = fread(paste0(o$pth$input, "ddf--datapoints--adults_with_hiv_percent_age_15_49--by--geo--time.csv")) %>%
    rename(hiv_prev = adults_with_hiv_percent_age_15_49)
  
  # HIV mortality
  hiv_mortality_dt = fread(paste0(o$pth$input, "ddf--datapoints--annual_hiv_deaths_number_all_ages--by--geo--time.csv")) %>%
    rename(hiv_mortality = annual_hiv_deaths_number_all_ages)
  
  # Maternal mortality
  maternal_mortality_dt = fread(paste0(o$pth$input, "ddf--datapoints--maternal_mortality_ratio_who--by--geo--time.csv")) %>%
    rename(maternal_mortality = maternal_mortality_ratio_who)
  
  # Health spending ($)
  health_spending_dt = fread(paste0(o$pth$input, "ddf--datapoints--total_health_spending_per_person_us--by--geo--time.csv")) %>%
    rename(health_spending = total_health_spending_per_person_us)
  
  # Private health spending (%)
  private_health_spending_dt = fread(paste0(o$pth$input, "ddf--datapoints--private_share_of_total_health_spending_percent--by--geo--time.csv")) %>%
    rename(private_health = private_share_of_total_health_spending_percent)
  
  # At least basic sanitation (%)
  sanitation_dt = fread(paste0(o$pth$input, "ddf--datapoints--at_least_basic_sanitation_overall_access_percent--by--geo--time.csv")) %>%
    rename(basic_sanitation = at_least_basic_sanitation_overall_access_percent)
  
  # At least basic water source (%)
  water_dt = fread(paste0(o$pth$input, "ddf--datapoints--at_least_basic_water_source_overall_access_percent--by--geo--time.csv")) %>%
    rename(basic_water = at_least_basic_water_source_overall_access_percent)
  
  # Human development index (life expectancy, education, per-person income)
  hdi_dt = fread(paste0(o$pth$input, "ddf--datapoints--hdi_human_development_index--by--geo--time.csv")) %>%
    rename(HDI = hdi_human_development_index)
  
  # Attended births
  attended_births_dt = fread(paste0(o$pth$input, "ddf--datapoints--births_attended_by_skilled_health_staff_percent_of_total--by--geo--time.csv")) %>%
    rename(attended_births = births_attended_by_skilled_health_staff_percent_of_total)
  
  # Internet users
  internet_users_dt = fread(paste0(o$pth$input, "ddf--datapoints--internet_users--by--geo--time.csv")) 
  
  
  # Create table of Gapminder covariates
  gapminder_dt = gini_dt %>%
    full_join(gdp_dt, by=c("geo", "time")) %>%
    full_join(literacy_male_15_dt, by=c("geo", "time")) %>%
    full_join(literacy_female_15_dt, by=c("geo", "time")) %>%
    full_join(doctors_per_1000_dt, by=c("geo", "time")) %>%
    full_join(health_spending_dt, by=c("geo", "time")) %>%
    full_join(private_health_spending_dt, by=c("geo", "time")) %>%
    full_join(pop_0_to_14_dt, by=c("geo", "time")) %>%
    full_join(pop_density_dt, by=c("geo", "time")) %>%
    full_join(urban_dt, by=c("geo", "time")) %>%
    full_join(hiv_prev_dt, by=c("geo", "time")) %>%
    full_join(hiv_mortality_dt, by=c("geo", "time")) %>%
    full_join(maternal_mortality_dt, by=c("geo", "time")) %>%
    full_join(internet_users_dt, by=c("geo", "time")) %>%
    full_join(rural_poverty_dt, by=c("geo", "time")) %>%
    full_join(underweight_children_dt, by=c("geo", "time")) %>%
    full_join(attended_births_dt, by=c("geo", "time")) %>%
    full_join(water_dt, by=c("geo", "time")) %>%
    full_join(sanitation_dt, by=c("geo", "time")) %>%
    full_join(hdi_dt, by=c("geo", "time")) %>%
    mutate(country_code = toupper(geo)) %>%
    select(-geo) %>%
    relocate(country_code) %>%
    rename(year = time) %>%
    full_join(WHO_regions_dt, by="country_code", relationship = "many-to-many") %>%
    select(-country) %>%
    rename(country = country_code) %>%
    arrange(country, year) %>%
    filter(year >= 1964 & year <= 2024)  %>%# include ten years before EPI to allow for historical effect of covariates
    as.data.table()              
  
  # Check for Gapminder countries not linked to WHO regions
  gapminder_dt %>% filter(is.na(region_short)) %>%
    select(country, region_short) %>%
    unique()
  
  gapminder_dt = gapminder_dt %>%
    filter(!is.na(region_short))
  
}

# ---------------------------------------------------------
# xxxxxxx
# ---------------------------------------------------------
covariates_unicef = function() {
  
  message("  - UNICEF")
  
  # Initiate covariates list
  covariates = list()
  
  # ---- Stunting ----
  
  # Read in UNICEF stunting data
  covariates$stunting = 
    fread(paste0(o$pth$input, "unicef_stunting.csv")) %>%
    # Format column names...
    setnames(names(.), tolower(names(.))) %>% 
    rename(country = "iso code") %>%
    # Reduce down to values of interest...
    filter(country %in% all_countries(), 
           indicator == "Stunting", 
           estimate  == "Point Estimate") %>%
    select(country, matches("^[0-9]+")) %>% 
    # Melt to long format...
    pivot_longer(cols = -country,
                 names_to = "year") %>% 
    # Format values...
    replace_with_na(list(value = "-")) %>% 
    mutate(value  = as.numeric(value), 
           year   = as.integer(year), 
           metric = "stunting") %>%
    as.data.table()
  
  # ---- Maternal mortality ----
  
  # Read in UNICEF maternal mortality data
  covariates$maternal_mortality = 
    fread(paste0(o$pth$input, "unicef_maternal_mortality.csv")) %>%
    # Reduce down to values of interest...
    select(country = iso, matches("^y[0-9]+")) %>%
    filter(country %in% all_countries()) %>%
    # Melt to long format...
    pivot_longer(cols = -country,
                 names_to = "year") %>%
    mutate(year   = as.integer(str_remove(year, "y")), 
           metric = "maternal_mortality") %>%
    as.data.table()
  
  return(covariates)
}

# ---------------------------------------------------------
# xxxxxxx
# ---------------------------------------------------------
covariates_gbd = function() {
  
  message("  - GBD")
  
  # Initiate covariates list
  covariates = list()
  
  # ---- Health and quality index ----
  
  # Prepare GBD 2019 HAQI for use as a covariate
  covariates$haqi = 
    fread(paste0(o$pth$input, "gbd19_haqi.csv")) %>%
    # Countries of interest...
    rename(country_name = location_name) %>%
    inner_join(y  = table("country"),
               by = "country_name") %>%
    # Tidy up...
    select(country, year = year_id, value = val) %>%
    mutate(value  = value / 100, 
           metric = "haqi") %>%
    arrange(country, year)
  
  # ---- Socio-demographic index ----
  
  # Function to load each SDI file
  load_sdi = function(name) {
    
    # Full file name
    file = paste0(o$pth$input, name, ".csv")
    
    # Load and convert to long form to allow binding
    sdi_dt = fread(file, header = TRUE) %>%
      rename(gbd_alt_name = Location) %>%
      pivot_longer(cols = -gbd_alt_name,
                   names_to  = "year") %>%
      as.data.table()
    
    return(sdi_dt)
  }
  
  # GBD country names
  gbd_name_dt = table("country") %>%
    mutate(gbd_alt_name = ifelse(
      test = is.na(gbd_alt_name), 
      yes  = country_name, 
      no   = gbd_alt_name)) %>%
    select(country, gbd_alt_name)
  
  # Prepare GBD 2019 SDI for use as a covariate
  covariates$sdi = 
    rbind(load_sdi("gbd19_sdi_1970"), 
          load_sdi("gbd19_sdi_1990")) %>%
    # Countries of interest...
    inner_join(y  = gbd_name_dt,
               by = "gbd_alt_name") %>%
    # Years of interest...
    mutate(year = as.integer(year)) %>%
    filter(year %in% o$years) %>%
    # Tidy up...
    select(country, year, value) %>%
    mutate(metric = "sdi") %>%
    arrange(country, year) %>%
    as.data.table()
  
  return(covariates)
}

