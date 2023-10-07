file <- system.file("extdata", "obs_wpp_demo.csv", package = "vieIA2030")
obs_wpp <- fread(file) %>%
  filter(!(variant == "Medium variant" & year == 2020)) %>%
  select(-variant) %>%
  filter(year %in% 2000:2095) %>%
  rename(wpp_country_code = location_code) %>%
  right_join(country_table[, .(country, wpp_country_code)], by = "wpp_country_code") %>%
  select(-wpp_country_code) %>%
  filter(!is.na(year)) %>%
  mutate(
    deaths_both = deaths_both * 1e03,
    deaths_male = deaths_male * 1e03,
    deaths_female = deaths_female * 1e03,
    births = births * 1e03)

upload_object(obs_wpp, "obs_wpp")