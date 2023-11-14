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
# Parse vaccine efficacy profile for non-VIMC pathogens
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
# Prepare GBD estimates of deaths for non-VIMC pathogens
# ---------------------------------------------------------
prepare_gbd_estimates = function() {
  
  message(" - GBD estimates")
  
  # Load GBD estimates of deaths for relevant diseases
  gbd_dt = read_rds("input", "gbd19_estimates")
  
  # Construct age datatable to expand age bins to single years
  age_bins = sort(unique(gbd_dt$age))
  age_dt   = data.table(age = o$data_ages) %>%
    mutate(age_bin = ifelse(age %in% age_bins, age, NA)) %>%
    fill(age_bin, .direction = "down") %>%
    group_by(age_bin) %>%
    add_count(age_bin) %>%
    ungroup() %>%
    as.data.table()
  
  # Expand to all ages and store
  gbd_dt %>%
    rename(age_bin = age) %>%
    full_join(y  = age_dt, 
              by = "age_bin", 
              relationship = "many-to-many") %>%
    arrange(country, disease, year, age) %>%
    mutate(deaths_disease = value / n) %>%
    # NOTE: OK to join only on disease as d_v_a is unique for GBD pathogens...
    left_join(y  = table("d_v_a"), 
              by = "disease", 
              relationship = "many-to-many") %>%
    select(country, d_v_a_id, year, age, deaths_disease) %>%
    save_table("gbd_estimates")
}

# ---------------------------------------------------------
# Prepare GBD covariates for extrapolating to non-VIMC countries
# ---------------------------------------------------------
prepare_gbd_covariates = function() {
  
  # TODO: Covariates stop at 2019... can we get an updated version?
  
  message(" - GBD covariates")
  
  # Prep GBD 2019 SDI for use as a covariate
  
  # fread(paste0(o$pth$input, "gbd19_sdi.csv"), header = TRUE) %>%
  #   mutate(n = 1 : n()) %>%
  #   mutate(Location = ifelse(n == 654, "Côte d'Ivoire", Location), 
  #          Location = ifelse(n == 664, "São Tomé and PrÍncipe", Location)) %>%
  #   filter(n != 105) %>%
  #   select(-n) %>%
  #   rename(gbd_alt_name = Location) %>%
  #   inner_join(y  = table("country")[, .(country, gbd_alt_name)],
  #              by = "gbd_alt_name") %>%
  #   select(country, all_of(as.character(1990 : 2019))) %>%
  #   pivot_longer(cols = -country,
  #                names_to = "year",
  #                values_to = "sdi") %>%
  #   mutate(year = as.integer(year)) %>%
  #   arrange(country, year) %>%
  #   as.data.table() %>%
  #   save_rds("input", "gbd19_sdi")
  
  # Prep GBD 2019 HAQI for use as a covariate
  
  # fread(paste0(o$pth$input, "gbd19_haqi.csv")) %>%
  #   mutate(n = 1 : n()) %>%
  #   filter(!n %in% (8641 : 8680)) %>%
  #   select(-n) %>%
  #   rename(country_name = location_name) %>%
  #   inner_join(y  = table("country")[, .(country, country_name)],
  #              by = "country_name") %>%
  #   select(country, year = year_id, haqi = val) %>%
  #   mutate(haqi = haqi / 100) %>%
  #   arrange(country, year) %>%
  #   save_rds("input", "gbd19_haqi")
  
  # NOTE: We're missing SDI for two countries: 
  gbd_sdi  = read_rds("input", "gbd19_sdi")
  gbd_haqi = read_rds("input", "gbd19_haqi")
  
  # Join metrics into single datatable
  gbd_covariates = 
    inner_join(x  = gbd_sdi, 
               y  = gbd_haqi, 
               by = c("country", "year")) %>%
    arrange(country, year)
  
  # Save in tables cache
  save_table(gbd_covariates, "gbd_covariates")
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
    name = file_names[[type]]
    
    # Initiate list to store data
    data_list = list()
    
    # Loop through past and future data
    for (year in file_years) {
      
      # Construct full file name
      file = paste0("WPP2022_", name, "BySingleAgeSex_Medium_", year, ".csv")
      
      # Stop here if file missing - ask user to download raw data
      if (!file.exists(paste0(o$pth$data, file)))
        stop("Please first download the file '", file, "' from",
             " https://population.un.org/wpp/Download/Standard",  
             " and copy to the /data/ directory")
      
      # Construct name of key data column
      data_name = paste0(first_cap(type), "Total")
      
      # Load the file and wrangle what we need
      data_list[[file]] = fread(paste0(o$pth$data, file)) %>%
        select(country = ISO3_code,
               year    = Time,
               age     = AgeGrp,
               metric  = !!data_name) %>%
        # Scale metrics by factor of 1k...
        mutate(metric = metric * 1e3) %>%
        rename_with(~type, metric) %>%
        # Only countries and years of interest...
        filter(country %in% table("country")$country,
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

