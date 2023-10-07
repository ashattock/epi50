###########################################################
# PREPARE WPP INPUT POPULATION
#
# Load historical and projected future population values.
#
###########################################################

# ---------------------------------------------------------
# Load population data from large WPP data files
# ---------------------------------------------------------
get_pop = function(data, sex, time_frame) {

  # Details for selecting historical data
  if (time_frame == "past") {
    sheet_name = "ESTIMATES"
    cell_range = "A17:DE18122"
    
    # Years to select from
    year_from  = 1980
    year_until = 2020  # This is the latest year in the ESTIMATES sheet
  }
  
  # Details for selecting future data
  if (time_frame == "future") {
    sheet_name = "MEDIUM VARIANT"
    cell_range = "A17:DE20672"
    
    # Years to select from
    year_from  = 2020  # This is the earliest year in the MEDIUM VARIANT sheet
    year_until = 2097
  }
  
  # Load data and melt to tidy format
  pop_dt = readxl::read_excel(path  = data[[sex]], 
                              sheet = sheet_name, 
                              range = cell_range) %>%
    select(wpp_country_code = "Country code", 
           year = "Reference date (as of 1 July)",
           all_of(as.character(0 : 100))) %>%
    pivot_longer(cols     = -c(wpp_country_code, year), 
                 names_to = "age") %>%
    filter(year >= year_from, 
           year <  year_until) %>%
    mutate(age = as.numeric(age), 
           sex_id = which(c("male", "female") == sex)) %>%
    group_by(wpp_country_code, year, sex_id, age) %>%
    summarise(nx = sum(value * 1000)) %>%
    ungroup() %>%
    setDT()
  
  return(pop_dt)
}

# ---- Load and format population data ----

# Path to large WPP poulation data files
data_pth = "supp_data/wpp_input_pop/WPP2019_INT_"

# Specify path for large WPP population data files
data = list(
  male   = paste0(data_pth, "F03_2_POPULATION_BY_AGE_ANNUAL_MALE.xlsx"), 
  female = paste0(data_pth, "F03_3_POPULATION_BY_AGE_ANNUAL_FEMALE.xlsx"))

# Load historical and projected pop sizes for both genders
pm1 = get_pop(data, "male",   "past")
pm2 = get_pop(data, "male",   "future")
pf1 = get_pop(data, "female", "past")
pf2 = get_pop(data, "female", "future")

# Concatenate into single datatable
wpp_input_pop = rbind(pm1, pm2, pf1, pf2)

# Save this datatable so it can be used directly within the package
save(wpp_input_pop, file = "inst/extdata/wpp_input_pop.RData")

