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
  
  # create_yaml(table_name, group_by, type = "rds")
  
  # Convert config yaml files to datatables
  prepare_config_tables()
  
  # Streamline VIMC impact estimates for quick loading
  prepare_vimc_impact()
  
  # Prepare GBD estimates of deaths for non-VIMC pathogens
  prepare_gbd_estimates()
  
  # Prepare GBD covariates for extrapolating to non-VIMC countries
  prepare_gbd_covariates()
  
  # Prepare demography-related estimates from WPP
  prepare_demography()
  
  # Prepare historical vaccine coverage
  prepare_coverage()  # See coverage.R for coverage-related functions
  
  # TODO: Is this needed? Can we just use WUENIC coverage instead?
  # prepare_hpv_target()
}

# ---------------------------------------------------------
# Convert config yaml files to datatables 
# ---------------------------------------------------------
prepare_config_tables = function() {
  
  message(" - Config files")
  
  # NOTE: Here we load in yaml files (/config/) and create
  #       rds files (/cache/) for fast loading
  
  # List of config yaml files to convert
  config_files = o$pth$config %>%
    list.files(pattern = ".+\\.yaml$") %>%
    str_remove(".yaml$")
  
  # Iterate through these files
  for (file in config_files) {
    
    # Load the yaml file
    y_file = paste0("config/", file, ".yaml")
    y_data = read_yaml(y_file)$table
    
    # Convert to datatable
    dt = rbindlist(lapply(y_data, as.data.table))
    
    # Save in cache
    save_table(dt, file)
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
prepare_vimc_impact = function() {
  
  message(" - VIMC estimates")
  
  # TODO: In raw form, we could instead use vimc_estimates.csv
  
  # All diseases of interest, not necessarily everything produced by VIMC
  diseases = table("disease")$disease
  
  # Prepare VIMC vaccine impact estimates
  readRDS(paste0(o$pth$input, "vimc_estimates.rds")) %>%
    filter(disease %in% diseases) %>%
    left_join(y  = table("d_v_a"), 
              by = c("disease", "vaccine", "activity")) %>%
    select(country, d_v_a_id, year, age, deaths_averted) %>%
    arrange(d_v_a_id, country, age, year) %>%
    save_table("vimc_impact")
  
  # Prepare VIMC year-of-vaccination results - take the mean across models
  readRDS(paste0(o$pth$input, "vimc_yov.rds")) %>%
    filter(disease %in% diseases) %>%
    left_join(y  = table("d_v_a"), 
              by = c("disease", "vaccine", "activity")) %>%
    select(country, d_v_a_id, model, year, deaths_averted, deaths_averted_rate) %>%
    group_by(country, d_v_a_id, year) %>%
    summarise(deaths_averted      = mean(deaths_averted,      na.rm = TRUE),
              deaths_averted_rate = mean(deaths_averted_rate, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(!is.na(deaths_averted_rate)) %>%
    arrange(d_v_a_id, country, year) %>%
    as.data.table() %>%
    save_table("vimc_yov")
  
  # Simply store VIMC in it's current form
  readRDS(paste0(o$pth$input, "vimc_uncertainty.rds")) %>%
    # select() %>%
    save_table("vimc_uncertainty")
}

# ---------------------------------------------------------
# Prepare GBD estimates of deaths for non-VIMC pathogens
# ---------------------------------------------------------
prepare_gbd_estimates = function() {
  
  message(" - GBD estimates")
  
  # Load GBD estimates of deaths for relevant diseases
  gbd_dt = readRDS(paste0(o$pth$input, "gbd19_estimates.rds"))
  
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
    full_join(age_dt, by = "age_bin", 
              relationship = "many-to-many") %>%
    arrange(country, disease, year, age) %>%
    mutate(deaths_disease = value / n) %>%
    # NOTE: OK to join only on disease as d_v_a is unique for GBD pathogens...
    left_join(y  = table("d_v_a"), 
              by = "disease") %>%
    select(country, d_v_a_id, year, age, deaths_disease) %>%
    save_table("gbd_estimates")
}

# ---------------------------------------------------------
# Prepare GBD covariates for extrapolating to non-VIMC countries
# ---------------------------------------------------------
prepare_gbd_covariates = function() {
  
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
  #   saveRDS("input/gbd19_sdi.rds")
  
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
  #   saveRDS("input/gbd19_haqi.rds")
  
  # NOTE: We're missing SDI for two countries: 
  gbd_sdi  = readRDS(paste0(o$pth$input, "gbd19_sdi.rds"))
  gbd_haqi = readRDS(paste0(o$pth$input, "gbd19_haqi.rds"))
  
  # Join metrics into single datatable
  gbd_covariates = 
    inner_join(x  = gbd_sdi, 
               y  = gbd_haqi, 
               by = c("country", "year")) %>%
    arrange(country, year)
  
  # Save in cache
  # save_table(gbd_covariates, "gbd_covariates")
  
  # ---- Forecast out to 2100 ----
  
  # TODO: Is this necessary for EPI50?? - probably not!
  
  forecast_gbd_cov <- function(gbd_covariates) {
    
    years_back = 5
    
    max_year = max(gbd_covariates$year)
    
    covariates = setdiff(names(gbd_covariates), qc(country, year))
    
    fcast_list = list()
    
    for (covariate in covariates) {
      
      gbd_mat = gbd_covariates %>%
        select(country, year, all_of(covariate)) %>%
        pivot_wider(names_from  = year, 
                    values_from =  all_of(covariate)) %>%
        select(-country) %>%
        as.matrix()
      
      start_vec = 1 - gbd_mat[, dim(gbd_mat)[2] - years_back]
      end_vec   = 1 - gbd_mat[, dim(gbd_mat)[2]]
      
      t_mat = matrix(data = 1 : (2100 - max_year), 
                     ncol = (2100 - max_year),
                     nrow = length(end_vec), 
                     byrow = TRUE)
      
      aroc = log(end_vec / start_vec) / years_back
      
      pred_mat = 1 - end_vec * exp(aroc * t_mat)
      colnames(pred_mat) = (max_year + 1):2100
      
      if (max(pred_mat) >= 1)
        warning("Forecast greater than or equal to 1: consider forecasting in logit space")
      
      out_dt = cbind(gbd_mat, pred_mat) %>%
        as.data.table() %>%
        mutate(country = unique(gbd_covariates$country)) %>%
        pivot_longer(cols = -country, 
                     names_to = "year") %>%
        mutate(metric = !!covariate, 
               year = as.integer(year)) %>%
        as.data.table()
      
      g = ggplot(out_dt) +
        aes(x = year, y = value, colour = country) +
        geom_line(show.legend = FALSE) +
        geom_vline(xintercept = max_year)
      
      fcast_list[[covariate]] <- out_dt
    }
    
    fcast_dt = rbindlist(fcast_list) %>%
      pivot_wider(names_from = metric) %>%
      as.data.table()
    
    return(fcast_dt)
  }
  
  gbd_covariates_fcast = forecast_gbd_cov(gbd_covariates)
  
  # Save in cache
  save_table(gbd_covariates_fcast, "gbd_covariates")
}

# ---------------------------------------------------------
# Prepare demography-related estimates from WPP
# ---------------------------------------------------------
prepare_demography = function() {
  
  message(" - Demography data")
  
  # SOURCE: https://population.un.org/wpp/Download/Standard
  
  # File names parts for WPP data
  file_names = list(
    pop   = "Population1January",
    death = "Deaths")
  
  # Files name years - past and future
  file_years = c("1950-2021", "2022-2100")
  
  # Years to extract (100 years worth)
  extract_years = 0 : 100 + min(o$analysis_years)
  
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
      if (!file.exists(paste0(o$pth$input, file)))
        stop("Please first download the file '", file, "' from",
             " https://population.un.org/wpp/Download/Standard",  
             " and copy to the /data/ directory")
      
      # Construct name of key data column
      data_name = paste0(first_cap(type), "Total")
      
      # Load the file and wrangle what we need
      data_list[[file]] = fread(paste0(o$pth$input, file)) %>%
        select(country = ISO3_code,
               year    = Time,
               age     = AgeGrp,
               metric  = !!data_name) %>%
        # Scale metrics by factor of 1k...
        mutate(metric = metric * 1e3) %>%
        rename_with(~type, metric) %>%
        # Only countries and years of interest...
        filter(country %in% table("country")$country,
               year    %in% extract_years) %>%
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
  
  # Save table in cache directory
  saveRDS(x, paste0(o$pth$cache, table, "_table.rds"))
}

# ---------------------------------------------------------
# Load and return cached datatable
# ---------------------------------------------------------
table = function(table) {
  
  # Construct file path
  file = paste0(o$pth$cache, table, "_table.rds")
  
  # Throw an error if this file doesn't exist
  if (!file.exists(file))
    stop("Table ", table, " has not been cached - have you run step 0?")
  
  # Load rds file
  y = readRDS(file)
  
  return(y)
}

# ---------------------------------------------------------
# TEMP: Convert input to yaml format
# ---------------------------------------------------------
create_yaml = function(name, id, type = "rds") {
  
  read_file = paste0("input/", name, "_table.", type) # Or data/
  yaml_file = paste0("config/", name, ".yaml")
  
  if (type == "csv") load_fn = "fread"
  if (type == "rds") load_fn = "readRDS"
  
  get(load_fn)(read_file) %>%
    split(.[[id]]) %>%
    unname() %>%
    as.yaml() %>%
    write_yaml(yaml_file)
}

