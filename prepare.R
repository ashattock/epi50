
# ---------------------------------------------------------
# xxxxxxxxxxxxxxxxx
# ---------------------------------------------------------
run_prepare = function() {
  
  # Convert config yaml files to datatables 
  prepare_config_tables()
  
  # Streamline VIMC impact estimates for quick loading
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
  
  ## Read data
  dt <- fread(
    system.file("extdata", "gbd19_sdi.csv", package = "vieIA2030"),
    header = T)
  
  browser()
}

# ---------------------------------------------------------
# TEMP: Convert input to yaml format
# ---------------------------------------------------------
create_yaml = function(name, id, type = "rds") {
  
  read_file = paste0("input/", name, "_table.", type)
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

