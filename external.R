###########################################################
# EXTERNAL
#
# Prepare, simulate (if appropriate), and extract results
# for external polio and measles models.
#
###########################################################

# ---------------------------------------------------------
# Parent function for extracting results for all external models
# ---------------------------------------------------------
run_external = function() {
  
  # Only continue if specified by do_step
  if (!is.element(2, o$do_step)) return()
  
  message("* Preparing external models")
  
  browser()
  
  # Simulate and extract results from DynaMICE
  simulate_dynamice()
  extract_dynamice()
  
  extract_measles()
  
  extract_polio()
}

# ---------------------------------------------------------
# Prepare and simulate DynaMICE measles model
# ---------------------------------------------------------
simulate_dynamice = function() {
  
  browser()
  
  # TODO: How important are 'mid_day' and 'coverage_subnat' 
  #       and how should they be derived?
  
  # Extract path of local DynaMICE repo
  repo_path = repo_exists("dynamice")
  
  # Return out now if repo doesn't exist locally
  if (is.null(repo_path))
    return()
  
  # ---- Dynamice coverage inputs ----
  
  message(" - Preparing inputs for DynaMICE measles model")
  
  # Convert EPI50-DynaMICE vaccine references
  dynamice_dict = c(
    MCV1    = "MCV1", 
    MCV2    = "MCV2", 
    Measles = "SIA")
  
  # Load EPI50 coverage details
  coverage_dt = table("coverage") %>%
    left_join(y  = table("v_a"), 
              by = "v_a_id") %>%
    filter(vaccine %in% names(dynamice_dict)) %>%
    mutate(vaccine = recode(vaccine, !!!dynamice_dict)) %>%
    select(vaccine, country, year, age, coverage)
  
  # Routine coverage - no age disaggregation required
  routine_dt = coverage_dt %>%
    filter(vaccine != "SIA") %>%
    select(vaccine, country, year, coverage) %>%
    # Summarise over age groups (noting that coverage the same)...
    group_by(vaccine, country, year) %>%
    summarise(coverage = mean(coverage)) %>%
    ungroup() %>%
    as.data.table()
  
  # Non-routine coverage - requires additonal details
  sia_dt = coverage_dt %>%
    filter(vaccine == "SIA") %>%
    # Convert coverage to char to enable pivot...
    mutate(coverage = round(coverage, 8), 
           coverage = as.character(coverage)) %>%
    # Group by all but age to find age bounds per campaign...
    group_by(vaccine, country, year, coverage) %>%
    mutate(age_first = min(age), 
           age_last  = max(age)) %>%
    ungroup() %>%
    # Reduce down to individual campaigns...
    select(-age) %>%
    unique() %>%
    # Convert coverage back to numeric...
    mutate(coverage = as.numeric(coverage)) %>%
    # Append any other additional variables needed...
    mutate(mid_day = 180) %>%             # NOTE: A placeholder assumption
    # Assume trivial subntaional coverage...
    mutate(coverage_subnat = coverage,    # NOTE: A placeholder assumption
           .after = coverage) %>%  
    as.data.table()
  
  # ---- Save input files to DynaMICE repo ----
  
  # Concatenate routine and non-routine coverage data
  data_dt = bind_rows(routine_dt, sia_dt) %>%
    arrange(vaccine, country, year, 
            age_first, age_last)
  
  # Inputs for 'all vaccine' and 'no vaccine' scenarios
  data_list = list(
    mcv1_mcv2_sia = data_dt, 
    nomcv         = data_dt[coverage == 0])
  
  # Iterate through scenarios
  for (scenario in names(data_list)) {
    
    # Construct file path to save to (in DynaMICE repo)
    save_path = file.path(repo_path, "input", "coverage", "coverage")
    save_file = paste0(paste1(save_path, scenario), ".csv")
    
    # Save data as a csv
    fwrite(data_list[[scenario]], file = save_file)
  }
  
  # ---- Other config files ----
  
  # Save full EPI50 country list
  country_dt   = data.table(country = all_countries())
  country_file = file.path(repo_path, "config", "countries.csv")
  
  # Also save associated regions
  region_dt   = table("country")[, .(country, region)]
  region_file = file.path(repo_path, "config", "regions.csv")
  
  # Write to config folder
  fwrite(country_dt, file = country_file)
  fwrite(region_dt,  file = region_file)
}

# ---------------------------------------------------------
# Extract outputs from DynaMICE measles model
# ---------------------------------------------------------
extract_dynamice = function() {
  
  browser()
}

# ---------------------------------------------------------
# Determine if specific repo exists locally
# ---------------------------------------------------------
repo_exists = function(repo) {
  
  # Path the parent directory on this EPI50 repository
  parent_path = str_remove(o$pth$code, "[a-z,A-Z,0-9]+/$")
  
  # Path to the repo in question
  repo_path = paste0(parent_path, repo)
  
  # If repo exists, return path
  if (dir.exists(repo_path))
    return(repo_path)
  
  # If it doesn't exist, return trivial
  if (!dir.exists(repo_path))
    return(NULL)
}

