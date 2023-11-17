###########################################################
# PREPARE
#
# Prepare various data sources for use throughout the analysis.
# The idea is that this process needs to be done only once, 
# and that the prepared inputs are streamlined for quick loading.
#
###########################################################

# ---------------------------------------------------------
# Parent function for preparing model inputs from various data sources
# ---------------------------------------------------------
run_prepare = function() {
  
  # Only continue if specified by do_step
  if (!is.element(1, o$do_step)) return()
  
  message("* Preparing input data")
  
  # Convert config yaml files to datatables
  prepare_config_tables()
  
  # Streamline VIMC impact estimates for quick loading
  prepare_vimc_estimates()
  
  # Parse vaccine efficacy profile for non-VIMC pathogens
  prepare_vaccine_efficacy()

  # Prepare GBD estimates of deaths for non-VIMC pathogens
  prepare_gbd_estimates()
  
  # Prepare GBD covariates for extrapolating to non-VIMC countries
  prepare_gbd_covariates()
  
  # Prepare demography-related estimates from WPP
  prepare_demography()
  
  # Prepare historical vaccine coverage
  prepare_coverage()  # See coverage.R for coverage-related functions
}

# ---------------------------------------------------------
# Convert config yaml files to datatables 
# ---------------------------------------------------------
prepare_config_tables = function() {
  
  message(" - Config files")
  
  # NOTE: Convert from yaml (/config) to rds (/tables) for fast loading
  
  # List of config yaml files to convert
  config_files = o$pth$config %>%
    list.files(pattern = ".+\\.yaml$") %>%
    str_remove(".yaml$")
  
  # Iterate through these files
  for (file in config_files) {
    
    # Load the yaml file
    yaml_file = paste0("config/", file, ".yaml")
    yaml_data = read_yaml(yaml_file)$table
    
    # Convert to datatable
    config_dt = yaml_data %>%
      lapply(as.data.table) %>%
      rbindlist(fill=TRUE)
    
    # Save in tables cache
    save_table(config_dt, file)
  }
  
  # Other config-related tables...
  
  # Disease-vaccine table
  table("d_v_a") %>%
    select(disease, vaccine) %>%
    unique() %>%
    mutate(d_v_id = 1 : n(), 
           .before = 1) %>%
    save_table("d_v")
  
  # Vaccine-activity table
  table("d_v_a") %>%
    select(vaccine, activity) %>%
    unique() %>%
    mutate(v_a_id = 1 : n(), 
           .before = 1) %>%
    save_table("v_a")
}

# ---------------------------------------------------------
# Streamline VIMC impact estimates for quick loading
# ---------------------------------------------------------
prepare_vimc_estimates = function() {
  
  message(" - VIMC estimates")
  
  # TODO: In raw form, we could instead use vimc_estimates.csv
  
  # Prepare VIMC vaccine impact estimates
  read_rds("input", "vimc_estimates") %>%
    left_join(y  = table("d_v_a"), 
              by = c("disease", "vaccine", "activity")) %>%
    filter(!is.na(d_v_a_id), 
           year %in% o$analysis_years) %>%
    select(country, d_v_a_id, year, age, deaths_averted) %>%
    arrange(country, d_v_a_id, age, year) %>%
    save_table("vimc_estimates")
  
  # Simply store VIMC in it's current form
  read_rds("input", "vimc_uncertainty") %>%
    save_table("vimc_uncertainty")
}

# ---------------------------------------------------------
# Parse vaccine efficacy profile for non-modelled pathogens
# ---------------------------------------------------------
prepare_vaccine_efficacy = function() {
  
  message(" - Vaccine efficacy")
  
  # Vaccine efficacy details
  pars_dt = table("vaccine_efficacy") %>%
    left_join(y  = table("d_v"), 
              by = "vaccine") %>%
    select(disease, vaccine, 
           a = efficacy, 
           x = decay_x, 
           y = decay_y)
  
  # We'll use these to represent waning immunity profile
  #
  # NOTE: We represent immunity decay using an exponential y = ae^(-rx)
  profile_list = list()
  
  # Points at which to evaluate exponential function
  #
  # NOTE: These represent immunity in the years following vaccination
  x_eval = seq_along(o$analysis_years) - 1
  
  # Iterate through diseases
  for (i in seq_row(pars_dt)) {
    pars = pars_dt[i, ]
    
    # If values are NA, interpret as no immunity decay
    has_decay = !any(is.na(pars[, .(x, y)]))
    
    # In this case, constant efficacy for n years
    if (!has_decay)
      profile = rep(pars$a, length(x_eval))
    
    # Otherwise represent this decay
    if (has_decay) {
      
      # Solve exponential function for r
      r = -log(pars$y / pars$a) / pars$x
      
      # Evaluate exponential using this rate
      profile = pars$a * exp(-r * x_eval)
    }
    
    # Form profile into a datatable
    profile_list[[i]] = pars %>%
      select(disease, vaccine) %>%
      expand_grid(time = x_eval) %>%
      mutate(profile = profile) %>%
      as.data.table()
  }
  
  # Bind into single datatable and save
  rbindlist(profile_list) %>%
    save_table("vaccine_efficacy_profiles")
  
  # Plot these profiles
  plot_vaccine_efficacy()
}

# ---------------------------------------------------------
# Prepare GBD estimates of deaths for non-modelled pathogens
# ---------------------------------------------------------
prepare_gbd_estimates = function() {
  
  message(" - GBD estimates")
  
  # NOTE: We are using GBD's 'number' metric here, where IA2030
  #       pipeline uses 'rate' (which is number per 100k people)
  
  # Parse specific age strings
  age_dict = c(
    "<1 year" = "0 to 1", 
    "80 plus" = "80 to 95")
  
  # Parse diseases
  disease_dict = c(
    "Diphtheria"     = "Dip",
    "Tetanus"        = "Tet",
    "Whooping cough" = "Per", 
    "Tuberculosis"   = "TB")
  
  # Load GBD estimates of deaths for relevant diseases
  gbd_dt = fread(paste0(o$pth$input, "gbd19_deaths.csv")) %>%
    # Parse disease and countries...
    mutate(disease = recode(cause, !!!disease_dict), 
           country = countrycode(
             sourcevar   = location,
             origin      = "country.name", 
             destination = "iso3c")) %>%
    filter(country %in% all_countries()) %>%
    # Parse age groups...
    mutate(age     = recode(age, !!!age_dict), 
           age_bin = str_extract(age, "^[0-9]+"), 
           age_bin = as.numeric(age_bin)) %>%
    # Tidy up...
    select(country, disease, year, age_bin, val) %>%
    arrange(country, disease, year, age_bin)
  
  # Construct age datatable to expand age bins to single years
  age_bins = sort(unique(gbd_dt$age_bin))
  age_dt   = data.table(age = o$data_ages) %>%
    mutate(age_bin = ifelse(age %in% age_bins, age, NA)) %>%
    fill(age_bin, .direction = "down") %>%
    group_by(age_bin) %>%
    add_count(age_bin) %>%
    ungroup() %>%
    as.data.table()
  
  # Expand to all ages and store
  gbd_dt %>%
    full_join(y  = age_dt, 
              by = "age_bin", 
              relationship = "many-to-many") %>%
    arrange(country, disease, year, age) %>%
    mutate(deaths_disease = val / n) %>%
    select(country, disease, year, age, deaths_disease) %>%
    save_table("gbd_estimates")
}

# ---------------------------------------------------------
# Prepare GBD covariates for extrapolating to non-modelled countries
# ---------------------------------------------------------
prepare_gbd_covariates = function() {
  
  message(" - GBD covariates")
  
  # Prepare GBD 2019 HAQI for use as a covariate
  haqi_dt = fread(paste0(o$pth$input, "gbd19_haqi.csv")) %>%
    # Countries of interest...
    rename(country_name = location_name) %>%
    inner_join(y  = table("country"),
               by = "country_name") %>%
    # Tidy up...
    select(country, year = year_id, haqi = val) %>%
    mutate(haqi = haqi / 100) %>%
    arrange(country, year)

  # Function to load each SDI file
  load_sdi = function(name) {
    
    # Full file name
    file = paste0(o$pth$input, name, ".csv")
    
    # Load and convert to long form to allow binding
    sdi_dt = fread(file, header = TRUE) %>%
      rename(gbd_alt_name = Location) %>%
      pivot_longer(cols = -gbd_alt_name,
                   names_to  = "year",
                   values_to = "sdi") %>%
      as.data.table()
    
    return(sdi_dt)
  }
  
  # Prepare GBD 2019 SDI for use as a covariate
  sdi_dt = rbind(load_sdi("gbd19_sdi_1"), 
        load_sdi("gbd19_sdi_2")) %>%
    # Countries of interest...
    inner_join(y  = table("country"),
               by = "gbd_alt_name") %>%
    select(country, year, sdi) %>%
    # Years of interest...
    mutate(year = as.integer(year)) %>%
    filter(year %in% o$analysis_years) %>%
    # Tidy up...
    arrange(country, year) %>%
    as.data.table()
  
  # Join metrics into single datatable
  gbd_covariates = 
    full_join(x  = sdi_dt, 
              y  = haqi_dt, 
              by = c("country", "year")) %>%
    arrange(country, year)
  
  # Save in tables cache
  save_table(gbd_covariates, "gbd_covariates")
  
  # Plot SDI - HAQi relationship
  plot_sdi_haqi()
}

# ---------------------------------------------------------
# Prepare demography-related estimates from WPP
# ---------------------------------------------------------
prepare_demography = function() {
  
  # TODO: Could we instead use the 'wpp2022' package?
  
  message(" - Demography data")
  
  # SOURCE: https://population.un.org/wpp/Download/Standard
  
  # File names parts for WPP data
  file_names = list(
    pop   = "Population1January",
    death = "Deaths")
  
  # Files name years - past and future
  file_years = c("1950-2021", "2022-2100")
  
  # Loop through data types
  for (type in names(file_names)) {
    
    # Filename part of this datataype
    metric = file_names[[type]]
    
    # Initiate list to store data
    data_list = list()
    
    # Loop through past and future data
    for (year in file_years) {
      
      # Construct full file name
      name = paste0("WPP2022_", metric, "BySingleAgeSex_Medium_", year, ".csv")
      file = paste0(o$pth$input, file.path("wpp", name))
      
      # Stop here if file missing - ask user to download raw data
      if (!file.exists(file))
        stop("Please first download the file '", name, "' from",
             " https://population.un.org/wpp/Download/Standard",  
             " and copy to the /input/wpp/ directory")
      
      # Construct name of key data column
      data_name = paste0(first_cap(type), "Total")
      
      # Load the file and wrangle what we need
      data_list[[name]] = fread(file) %>%
        select(country = ISO3_code,
               year    = Time,
               age     = AgeGrp,
               metric  = !!data_name) %>%
        # Scale metrics by factor of 1k...
        mutate(metric = metric * 1e3) %>%
        rename_with(~type, metric) %>%
        # Only countries and years of interest...
        filter(country %in% all_countries(),
               year    %in% o$analysis_years) %>%
        mutate(age = ifelse(age == "100+", 100, age),
               age = as.integer(age))
    }
    
    # Combine past and future data and save to file
    rbindlist(data_list) %>%
      arrange(country, year, age) %>%
      save_table(paste1("wpp", type))
  }
}

# ---------------------------------------------------------
# Simple wrapper to load all countries
# ---------------------------------------------------------
all_countries = function() {
  
  # Pull all countries defined in config file
  countries = table("country")$country
  
  return(countries)
}

# ---------------------------------------------------------
# Save table in cache directory for quick loading
# ---------------------------------------------------------
save_table = function(x, table) {
  
  # Save table in tables cache directory
  save_rds(x, "tables", table, "table")
}

# ---------------------------------------------------------
# Load and return cached datatable
# ---------------------------------------------------------
table = function(table) {
  
  # Construct file path
  file = paste0(o$pth$tables, table, "_table.rds")
  
  # Throw an error if this file doesn't exist
  if (!file.exists(file))
    stop("Table ", table, " has not been cached - have you run step 0?")
  
  # Load rds file
  y = read_rds(file)
  
  return(y)
}

