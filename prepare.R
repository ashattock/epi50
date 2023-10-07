
# ---------------------------------------------------------
# xxxxxxxxxxxxxxxxx
# ---------------------------------------------------------
run_prepare = function() {
  
  # Only continue if specified by do_step
  if (!is.element(0, o$do_step)) return()
  
  # Convert config yaml files to datatables 
  prepare_config_tables()
  
  # # Streamline VIMC impact estimates for quick loading
  prepare_vimc_impact()

  # TODO: Is this needed? Can we just use WUENIC coverage instead?
  prepare_hpv_target()

  # Prepare GBD estimates of deaths for non-VIMC pathogens
  prepare_gbd_estimates()
  
  # 
  prepare_gbd_covariates()
  
  browser()
  
  # 
  create_coverage()
  
}

# ---------------------------------------------------------
# Convert config yaml files to datatables 
# ---------------------------------------------------------
prepare_config_tables = function() {
  
  # create_yaml("country", "country", type = "rds")
  
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
  load_table("d_v_a") %>%
    select(disease, vaccine) %>%
    unique() %>%
    mutate(d_v_id = 1 : n(), 
           .before = 1) %>%
    save_table("d_v")
  
  # Vaccine-activity table
  load_table("d_v_a") %>%
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
  
  # TODO: In raw form, we could instead use vimc_estimates.csv
  
  # Prepare VIMC vaccine impact estimates
  readRDS(paste0(o$pth$input, "vimc_estimates.rds")) %>%
    left_join(y  = load_table("d_v_a"), 
              by = c("disease", "vaccine", "activity")) %>%
    select(country, d_v_a_id, year, age, deaths_averted) %>%
    arrange(d_v_a_id, country, age, year) %>%
    save_table("vimc_impact")
  
  # Prepare VIMC year-of-vaccination results - take the mean across models
  readRDS(paste0(o$pth$input, "vimc_yov.rds")) %>%
    left_join(y  = load_table("d_v_a"), 
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
# xxxxxxxx
# ---------------------------------------------------------
prepare_hpv_target = function() {
  
  # TODO: De we need these targets for EPI50 - probably not
  
  # Prep HPV target coverage for IA2030 coverage scenario
  fread(paste0(o$pth$input, "hpv_target_coverage.csv")) %>%
    pivot_longer(cols = -c(country, sex, age, doses), 
                 names_to = "year") %>%
    mutate(year = as.integer(year)) %>%
    replace_na(list(value = 0)) %>%
    # Keep the max of the one and two does coverage levels...
    group_by(country, year, sex, age) %>%
    summarise(value = max(value)) %>%
    ungroup() %>%
    # Append country ISO...
    rename(hpv_name = country) %>%
    left_join(y  = load_table("country")[, .(country, hpv_name)], 
              by = c("hpv_name")) %>%
    select(-hpv_name) %>%
    # Breakdown aggregate into distinct gender...
    mutate(sex = tolower(substr(sex, 1, 1))) %>%
    pivot_wider(names_from = sex) %>%
    mutate(f = ifelse(is.na(f), b, f), 
           m = ifelse(is.na(m), b, m)) %>%
    select(-b) %>%
    pivot_longer(cols = -c(country, year, age), 
                 names_to = "sex") %>%
    filter(!is.na(value)) %>%
    # Final formatting...
    select(country, year, age, sex, hpv_target = value) %>%
    arrange(country, year, age, sex) %>%
    as.data.table() %>%
    save_table("hpv_target")
}

# ---------------------------------------------------------
# Prepare GBD estimates of deaths for non-VIMC pathogens
# ---------------------------------------------------------
prepare_gbd_estimates = function() {
  
  # Load GBD estimates of deaths from relevant diseases
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
    arrange(country, disease, year, sex, age) %>%
    filter(!is.na(value)) %>%
    mutate(value := value / n) %>%
    select(-age_bin, n) %>%
    left_join(y  = load_table("d_v_a"), 
              by = "disease") %>%
    select(d_v_a_id, country, year, sex, age, strata_deaths = value) %>%
    save_table("gbd_estimates")
}

# ---------------------------------------------------------
# xxxxxxxx
# ---------------------------------------------------------
prepare_gbd_covariates = function() {
  
  # Prep GBD 2019 SDI for use as a covariate
  
  # fread(paste0(o$pth$input, "gbd19_sdi.csv"), header = TRUE) %>%
  #   mutate(n = 1 : n()) %>%
  #   mutate(Location = ifelse(n == 654, "Côte d'Ivoire", Location), 
  #          Location = ifelse(n == 664, "São Tomé and PrÍncipe", Location)) %>%
  #   filter(n != 105) %>%
  #   select(-n) %>%
  #   rename(gbd_alt_name = Location) %>%
  #   inner_join(y  = load_table("country")[, .(country, gbd_alt_name)],
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
  #   inner_join(y  = load_table("country")[, .(country, country_name)],
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
load_table = function(table) {
  
  # Construct file path
  file = paste0(o$pth$cache, table, "_table.rds")
  
  # Throw an error if this file doesn't exist
  if (!file.exists(file))
    stop("Table ", table, " has not been cached - have you run step 0?")
  
  # Load rds file
  y = readRDS(file)
  
  return(y)
}

