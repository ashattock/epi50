###########################################################
# INTERFACE
#
# Set of transmission-model specific helper functions to interface
# between EPI50 pipeline and the specific model itself.
#
###########################################################

# ---------------------------------------------------------
# Prepare inputs and outputs for DynaMICE measles model
# ---------------------------------------------------------
prepare_dynamice = function() {
  
  # TODO: How important are 'mid_day' and 'coverage_subnat' 
  #       and how should they be derived?
  
  # Extract path of local DynaMICE repo
  repo_path = repo_exists("dynamice")
  
  # Return out now if repo doesn't exist locally
  if (is.null(repo_path))
    return()
  
  # ---- Dynamice inputs ----
  
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
    mutate(coverage = mean(coverage)) %>%
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
           age_last  = max(age), 
           .after = age) %>%
    ungroup() %>%
    # Reduce down to individual campaigns...
    select(-age) %>%
    unique() %>%
    # Convert coverage back to numeric...
    mutate(coverage = as.numeric(coverage)) %>%
    # Append any other additional variables needed...
    mutate(mid_day = 180,                   # NOTE: An assumption
           .before = coverage) %>%
    # Assume trivial subntaional coverage...
    mutate(coverage_subnat = coverage) %>%  # NOTE: An assumption
    as.data.table()
  
  # ---- Save input files to DynaMICE repo ----
  
  # Construct list of all coverage data
  coverage_data = list(
    # All vaccine scenario...
    mcv1_mcv2_sia = list(
      routine = routine_dt, 
      sia     = sia_dt), 
    # No vaccine scenario...
    nomcv = list(
      routine = routine_dt[coverage == 0], 
      sia     = sia_dt[coverage == 0]))
  
  # Iterate through scenarios
  for (scenario in names(coverage_data)) {
    
    # Iterate through activity types
    for (activity in names(coverage_data[[scenario]])) {
      
      # Extract the relevant data
      save_data = coverage_data %>%
        pluck(scenario) %>%
        pluck(activity)
      
      # Construct file path to save to (in DynaMICE repo)
      save_path = file.path(repo_path, "input", "coverage", activity)
      save_file = paste0(paste1(save_path, scenario), ".csv")
      
      # Save data as a csv
      fwrite(save_data, file = save_file)
    }
  }
}

# ---------------------------------------------------------
# Determine if specific repo exists locally
# ---------------------------------------------------------
repo_exists = function(repo) {
  
  parent_path = str_remove(o$pth$code, "[a-z,A-Z,0-9]+/$")
  
  repo_path = paste0(parent_path, repo)
  
  if (dir.exists(repo_path))
    return(repo_path)
  
  if (!dir.exists(repo_path))
    return(NULL)
}

