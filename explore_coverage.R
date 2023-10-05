


explore_coverage = function() {
  
  # Load raw data
  #
  # SOURCE: xxxx
  file_raw = "coverage-data.xlsx"
  path_raw = system.file("sia", file_raw, package = "vieIA2030")
  data_raw = read_excel(path_raw)
  
  # Load data dictionary
  # file_dict = "SIA_data_dictionary.csv"
  # data_dict = fread(system.file("sia", file_dict, package = "vieIA2030"))
  
  # ----
  
  data = data_raw %>%
    select(group    = GROUP, 
           country  = CODE, 
           year     = YEAR, 
           vaccine  = ANTIGEN,
           type     = COVERAGE_CATEGORY, 
           coverage = COVERAGE) %>%
    filter(group == "COUNTRIES", 
           type  == "OFFICIAL", 
           !is.na(coverage)) %>%
    select(-group, -type)
  
  # ----
  
  # # Load and format data dictionary
  # dict_dt = data_dict %>%
  #   mutate(intervention = tolower(intervention)) %>%
  #   fill(intervention) %>%
  #   select(-notes)
  # 
  # no_int = setdiff(data$intervention, dict_dt$intervention)
  
  # ----
  
  g = ggplot(data) +
    aes(x = year, 
        y = coverage, 
        colour = country) +
    geom_point(show.legend = FALSE) +
    facet_wrap(~vaccine)
  
  g = g + ylim(0, 100)
  
  browser()
  
  
}
