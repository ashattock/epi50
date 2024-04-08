
# NOTE: These raw results are not git committed to the repository due to sheer file size

age_dict = c(
  "Early Neonatal"  = "<28 days", 
  "Late Neonatal"   = "<28 days", 
  "1-5 months"      = "28-364 days", 
  "6-11 months"     = "28-364 days", 
  "12 to 23 months" = "1-4", 
  "2 to 4"          = "1-4", 
  "95 plus"         = "95+")

gbd_path = paste0(o$pth$input, "gbd")

gbd_data = fread(file.path(gbd_path, "gbd21_raw.csv")) %>%
  lazy_dt() %>%
  rename(cause    = cause_name,
         location = location_name, 
         age      = age_group_name, 
         year     = year_id) %>%
  mutate(cause = ifelse(
    test = grepl("Tuberculosis", cause), 
    yes  = "Tuberculosis", 
    no   = cause)) %>%
  mutate(country = countrycode(
    sourcevar   = location,
    origin      = "country.name", 
    destination = "iso3c")) %>%
  filter(country %in% all_countries()) %>%
  mutate(age = recode(age, !!!age_dict), 
         age = str_replace(age, " to ", "-")) %>%
  group_by(measure, cause, country, age, year) %>%
  summarise(value = sum(value),
            lower = sum(lower),
            upper = sum(upper)) %>%
  ungroup() %>%
  as.data.table() %>%
  setkey(NULL) %>%
  split(.$measure) %>%
  lapply(function(x) select(x, -measure))

save_rds(gbd_data$death, file.path(gbd_path, "gbd21_deaths.rds"))
save_rds(gbd_data$daly,  file.path(gbd_path, "gbd21_dalys.rds"))

