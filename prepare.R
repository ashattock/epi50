
# ---------------------------------------------------------
# xxxxxxxxxxxxxxxxx
# ---------------------------------------------------------
run_prepare = function() {
  
  # Convert config yaml files to datatables 
  prepare_config_tables()
  
  # Streamline VIMC impact estimates for quick loading
  prepare_vimc_impact()
  
  prepare_hpv_target()
  
  browser()
  
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
  
  # ---- Config-related datatables ----
  
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
# xxxxxxxxxxxxxxxxx
# ---------------------------------------------------------
load_table = function(table) {
  
  y = readRDS(paste0(o$pth$cache, table, "_table.rds"))
  
  return(y)
}

