###########################################################
# SIA DATA EXPLORATION
#
# Sources of SIA data:
#   1) https://extranet.who.int/xmart4/WIISE/data/SIA_MAIN
#   2) https://www.who-immunization-repository.org/
#
###########################################################

# ---------------------------------------------------------
# Extract SIA coverage data
# ---------------------------------------------------------
coverage_sia = function() {
  
  message("  > SIA coverage")
  
  # ---- Load and filter raw data ----
  
  # TODO: There is nothing we can do (in terms of impact estimation)
  #       for SIA events not modelled by VIMC, so filter out these.
  
  data_dict = table("sia_dictionary") %>%
    select(-notes) %>%
    mutate(activity = "campaign") %>%
    inner_join(y  = table("d_v_a"), 
               by = c("disease", "vaccine", "activity"))
  
  # Entries to set as NA
  na_var = c("unknown", "undefined", "")
  
  # Load raw data
  sia_raw = fread(paste0(o$pth$input, "sia_data.csv"))
  
  # browser()
  
  # Select data of interest
  sia_dt = sia_raw %>%
    select(country      = ISO3_CODE,          # Country ISO3 codes
           intervention = INTERVENTION_CODE,
           sia_type     = ACTIVITY_TYPE_CODE, 
           status       = ACTIVITY_STATUS,    # Done, Ongoing, Planned, etc
           age_group    = ACTIVITY_AGE_GROUP,
           cohort       = TARGET, 
           doses        = DOSES, 
           coverage     = COVERAGE, 
           matches("YEAR$|MONTH$|DAY$")) %>%
    # Format data types...
    mutate_if(is.character, tolower) %>%
    mutate(across(.cols = c(cohort, coverage, doses), 
                  .fns  = as.numeric)) %>%
    # Attempt to impute doses and then filter out the unworkable...
    mutate(doses = ifelse(is.na(doses), cohort * coverage, doses)) %>%
    filter(doses > 0) %>%
    # Remove any activities not of interest...
    filter(!status   %in% c("cancelled", "postponed"), 
           !sia_type %in% c("routine")) %>%
    select(-status, -sia_type) %>%
    # Remove any unknown countries...
    mutate(country = toupper(country)) %>%
    filter(country %in% table("country")$country) %>%
    # Remove any unknown interventions...
    filter(intervention %in% unique(data_dict$intervention)) %>%
    arrange(country, intervention) %>%
    # Deal with other unknown entries...
    mutate(across(.cols = where(is.character), 
                  .fns  = ~if_else(. %in% na_var, NA, .)), 
           across(.cols = where(is.character), 
                  .fns  = ~if_else(is.na(.), "unknown", .))) %>%
    # Deal with dates...
    format_sia_dates() %>%
    impute_sia_dates() %>%
    # Parse age groups into age ranges...
    parse_age_groups() %>%
    # Then covert doses to FVP...
    convert_fvps() %>%
    # Convert to d_v_a...
    left_join(y  = data_dict, 
              by = "intervention", 
              relationship = "many-to-many") %>%
    select(country, d_v_a_id, year, age, fvps, cohort, coverage) %>%
    mutate(source = "sia")
  
  browser()
  
  # ---- Fully Vaccinated Persons (FVP) ----
  
  # Load doses per FVP table
  doses_per_fvp = table("sia_doses") %>%
    rename(doses_fvps = doses) %>%
    select(-notes)
  
  fvps_dt = sia_dt %>%
    left_join(doses_per_fvp, 
              by = "disease") %>%
    mutate(fvps = (doses / doses_fvps) / 1e6) %>%
    group_by(region, disease) %>%
    summarise(fvps = sum(fvps)) %>%
    ungroup() %>%
    complete(region, disease, fill = list(fvps = 0)) %>%
    as.data.table()
  
  g1 = ggplot(fvps_dt) + 
    aes(x = 1, y = fvps, fill = region) + 
    geom_bar(stat = "identity", colour = "black", size = 0.05) + 
    facet_wrap(~disease)
  
  # Save figures to file
  save_fig(g1, dir = "sia_data", "SIA FVPs")
  
  # ---- Doses per month ----
  
  # Plot durations for all entries
  plot_durations(sia_dt)
  
  # Then same plots grouped by different variables
  plot_durations(sia_dt, "extent")
  plot_durations(sia_dt, "sia_type")
  plot_durations(sia_dt, "disease")
  
  # ---- Plot by country ----
  
  # Calculate doses and cumulative doses for all possible months
  #
  # NOTE: we only keep trivial months for nicer plotting of cumulative doses
  country_dt1 = sia_month_dt %>%
    # unite("d_v", disease, vaccine) %>%
    # group_by(country, d_v, month) %>%
    group_by(country, disease, month) %>%
    summarise(doses = sum(doses)) %>%  # Doses over time country-disease-vaccine
    mutate(cum_doses = cumsum(doses)) %>%  # Cumulative doses country-disease-vaccine
    ungroup() %>%
    as.data.table()
  
  # Remove (most of) the zeros for nicer line plots
  country_dt2 = country_dt1 %>%
    group_by(country, disease) %>%
    filter(lead(cum_doses) > 0) %>%  # Remove all but most recent trailing zeros (for pretty plotting)
    ungroup() %>%
    filter(!(doses == 0 & cum_doses > 0)) %>%  # Remove leading zeros
    as.data.table()
  
  # Get colours - one per country
  cols = get_colours(length(levels(country_dt2$country)))
  
  # Plot area of cumulative doses per country over time
  g1 = ggplot(country_dt1) + 
    aes(x = month, y = cum_doses / 1e6, fill = country) +
    geom_area() +
    facet_wrap(~disease, scales = "free_y")
  
  # Plot cumulative doses per country over time
  g2 = ggplot(country_dt2) + 
    aes(x = month, y = cum_doses / 1e6, colour = country) +
    geom_line() +
    facet_wrap(~disease, scales = "free_y") 
  
  # Save figures to file
  save_fig(g1, dir = "sia_data", "SIA doses by country", "area")
  save_fig(g2, dir = "sia_data", "SIA doses by country", "line")
  
  # ---- Plot by disease ----
  
  # Total doses for each d_v - sum over countries
  disease_dt1 = sia_month_dt %>%
    # unite("d_v", disease, vaccine) %>%
    # group_by(d_v, month) %>%
    group_by(disease, month) %>%
    summarise(doses = sum(doses)) %>%
    mutate(cum_doses = cumsum(doses)) %>%
    ungroup() %>%
    # Number of FVPs...
    left_join(doses_per_fvp, 
              by = "disease") %>%
    mutate(fvps     = (doses     / doses_fvps), 
           cum_fvps = (cum_doses / doses_fvps)) %>%
    replace_na(list(fvps = 0, cum_fvps = 0)) %>%
    as.data.table()
  
  # Remove trivial entires for nicer line plot
  disease_dt2 = disease_dt1 %>%
    filter(doses > 0, cum_doses > 0)
  
  # Get colours - one per disease
  n_cols = length(unique(disease_dt2$disease))
  cols   = get_colours(n_cols)
  
  # Plot doses as areas using full datatable
  g1 = ggplot(disease_dt1) +
    aes(x = month, y = cum_doses / 1e9, fill = disease) +
    geom_area() +
    scale_y_continuous(name   = "Cumulative doses (billions)",
                       expand = expansion(mult = c(0, 0.05)),
                       labels = comma)
  
  # Plot FVPs as areas using full datatable
  g2 = ggplot(disease_dt1) +
    aes(x = month, y = cum_fvps / 1e9, fill = disease) +
    geom_area() +
    scale_y_continuous(name   = "Cumulative FVPs (billions)",
                       expand = expansion(mult = c(0, 0.05)),
                       labels = comma)
  
  # Plot lines using reduced datatable
  g3 = ggplot(disease_dt2) +
    aes(x = month, y = cum_doses, colour = disease) +
    geom_line(size = 2) +
    scale_y_continuous(name   = "Cumulative doses (log10 scale)", 
                       trans  = "log10", 
                       limits = c(1, NA),
                       expand = expansion(mult = c(0, 0.05)),
                       labels = comma)
  
  # Save figures to file
  save_fig(g1, dir = "sia_data", "SIA doses by disease", "area")
  save_fig(g2, dir = "sia_data", "SIA FVPs by disease",  "area")
  save_fig(g3, dir = "sia_data", "SIA doses by disease", "line")
  
  # ---- Plot by data source ----
  
  # First determine countries reported by VIMC
  vimc_countries = table("vimc_estimates") %>%
    select(country, d_v_a_id) %>%
    unique() %>%
    left_join(y  = table("d_v_a"), 
              by = "d_v_a_id") %>%
    select(country, disease) %>%
    unique() %>%
    mutate(country_source = "vimc")
  
  # Join VIMC disease and countries to SIA data
  vimc_dt = sia_month_dt %>%
    # Number of FVPs...
    left_join(doses_per_fvp, 
              by = "disease") %>%
    mutate(fvps = (doses / doses_fvps)) %>%
    replace_na(list(fvps = 0)) %>%
    select(-doses_fvps) %>%
    # Information source...
    left_join(vimc_countries[, .(country, disease, country_source)], 
              by = c("country", "disease")) %>%
    left_join(table("disease")[, .(disease, source)], 
              by = "disease") %>%
    rename(impact_source = source) %>%
    mutate(country_source = ifelse(is.na(country_source), "nosource", country_source),
           impact_source  = ifelse(is.na(impact_source),  "nosource", impact_source), 
           country_source = paste0("country_", country_source), 
           impact_source  = paste0("impact_",  impact_source))
  
  # Group by source of data for disease and country
  source_dt = vimc_dt %>%
    group_by(disease, country_source, impact_source, month) %>%
    summarise(doses = sum(doses), 
              fvps  = sum(fvps)) %>%
    mutate(cum_doses = cumsum(doses), 
           cum_fvps  = cumsum(fvps)) %>%
    ungroup() %>%
    as.data.table()
  
  # Plot by data source - filled by d_v
  g1 = ggplot(source_dt) + 
    aes(x = month, y = cum_fvps / 1e9, fill = disease) + 
    geom_area() + 
    facet_grid(country_source ~ impact_source)
  
  # Save figures to file
  save_fig(g1, dir = "sia_data", "SIA doses by source")
  
  
  
  
  
  
  return(sia_dt)
}

# ---------------------------------------------------------
# SIA database has numerous date columns - combine into useable format
# ---------------------------------------------------------
format_sia_dates = function(sia_dt) {
  
  # TODO... do any of the entries with missing start dates have end dates??
  # If so, we could impute the start date as we do for missing end dates.
  
  # Define date columns in raw data set
  date_cols = c("plan", "postponed", "done")
  date_type = c("start", "end")
  
  # Function to create date strings from multiple columns
  date_fn = function(col) {
    
    # Combine all columns to create y-m-d string
    date_str = paste(
      sia_dt[[toupper(paste0(col, "_year"))]], 
      sia_dt[[toupper(paste0(col, "_month"))]], 
      sia_dt[[toupper(paste0(col, "_day"))]], 
      sep = ".")
    
    # Trivialise any dates containing NA
    date_str[grepl("NA", date_str)] = ""
    
    return(date_str)
  }
  
  # Loop through date columns to create a single variable
  for (i in date_cols) {
    for (j in date_type) {
      
      # Construct new column name
      col = paste(i, j, sep = "_")
      
      # Apend this newly compiled column
      sia_dt[[col]] = date_fn(col)
    }
  }
  
  # First and last dates we are interested in
  data_from = paste(min(o$analysis_years),     1, 1, sep = ".")
  data_to   = paste(max(o$analysis_years + 1), 1, 1, sep = ".")
  
  # Create single start and end columns
  date_dt = sia_dt %>%
    mutate(start_date = paste0(plan_start, postponed_start, done_start), 
           end_date   = paste0(plan_end,   postponed_end,   done_end)) %>%
    # Convert to date format...
    mutate(start_date = format_date(start_date), 
           end_date   = format_date(end_date)) %>%
    # Remove the numerous now-redundant date columns...
    select(-matches("YEAR$|MONTH$|DAY$"), 
           -ends_with("start"), 
           -ends_with("end")) %>%
    # Remove ineligible dates...
    filter(!is.na(start_date),
           start_date >= format_date(data_from), 
           start_date <  format_date(data_to)) %>%
    arrange(country, intervention, start_date)
  
  return(date_dt)
}

# ---------------------------------------------------------
# Impute missing end dates, distribute over time, and sum over years
# ---------------------------------------------------------
impute_sia_dates = function(sia_dt) {
  
  average_fn = "mean"  # OPTION: "mean" or "median"
  
  # Impute missing end dates
  impute_dt = sia_dt %>%
    # Calculate average duration...
    mutate(days    = as.numeric(end_date - start_date), 
           average = get(average_fn)(days, na.rm = TRUE)) %>%
    # Fill in any missing end dates with duration average...
    mutate(fix_date = start_date + average, 
           end_date = if_else(is.na(end_date), fix_date, end_date)) %>%
    select(-days, -average, -fix_date)
  
  # Melt monthly dates to tidy format...
  
  # Single datatable column of all possible months
  all_months_dt = 
    seq(from = floor_date(min(impute_dt$start_date), "month"), 
        to   = floor_date(max(impute_dt$end_date),   "month"), 
        by   = "month") %>%
    as.character() %>%
    as_named_dt("month")
  
  # All months to distibute doses across (for which campaign has been 'run')
  run_months_dt = impute_dt %>%
    mutate(start_date = floor_date(start_date, "month"),  # Beginning of month
           end_date   = floor_date(end_date,   "month"),  # Beginning of month
           end_date   = pmax(start_date, end_date)) %>%   # In case of end_date < start_date
    rowwise() %>%
    mutate(run_months = seq(start_date, end_date, by = "month") %>%  # All months to distibute across
             paste(collapse = " & ") %>%
             as.character(), 
           n_months = str_count(run_months, "&") + 1) %>%  # Number of months to distibute across
    ungroup() %>%
    as.data.table()
  
  # Expand for all possible and distrubte doses across months
  #
  # NOTES: 
  #  - We'll use this 'all months' dt for pretty plotting
  #  - Whilst this works, there could well be a more efficient way
  sia_month_dt = run_months_dt %>%
    expand_grid(all_months_dt) %>%                     # Full factorial for all possible months
    mutate(value = str_detect(run_months, month)) %>%  # Successful matches
    mutate(month = format_date(month), 
           doses = (doses / n_months) * value) %>%     # Divide total doses across the months
    select(-start_date, -end_date, -run_months, -n_months) %>%
    arrange(country, intervention, month) %>%
    as.data.table()
  
  # Remove these trivial dose entries and sum over year
  sia_year_dt = sia_month_dt %>%
    filter(doses > 0) %>%
    mutate(year = year(month)) %>%
    group_by(country, intervention, year, age_group) %>%
    summarise(doses = sum(doses)) %>%
    ungroup() %>%
    as.data.table()
  
  # Sanity check that we haven't changed number of doses
  # if (abs(sum(sia_year_dt$doses) - check_doses) > 1e-3)
  #   stop("We seem to have gained/lost doses here")
  
  return(sia_year_dt)
}

# ---------------------------------------------------------
# Parse age groups into age ranges
# ---------------------------------------------------------
parse_age_groups = function(sia_dt) {
  
  # All unique age group strings to parse
  age_groups = sort(unique(sia_dt$age_group))
  
  # Initialise working datatable
  group_dt = data.table(
    group = age_groups, 
    age   = age_groups)
  
  # ---- Parse age group stings ----
  
  # Regular expression to remove
  exp_rm = c(",.+", "\n.+", "\\+", "=", "[a-z, ]")
  
  # Regular expression to substitute
  exp_sub = c(
    "&"  = "and",
    "^<" = "0&", 
    "^>" = "100&", 
    "-<" = "&", 
    "->" = "&", 
    "-"  = "&", 
    " y" = "#", 
    " m" = "@") 
  
  # Several named special cases
  exp_txt = c(
    "all ages"    = "0-100",
    "adolescents" = "10-19",
    "adults"      = "18-60",
    "children"    = "1-12",
    "elderly"     = "60-100",
    "school"      = "5-16",
    "women"       = "18-60")
  
  # Parse key characters
  for (exp in names(exp_sub))
    group_dt[grepl(exp, age), age := gsub(exp, exp_sub[[exp]], age)]
  
  # Remove certain characters
  for (exp in exp_rm)
    group_dt[grepl(exp, age), age := gsub(exp, "", age)]
  
  # Parse strings that we don't yet have info for
  for (exp in names(exp_txt))
    group_dt[grepl(exp, group) & age == "", 
           age := gsub("-", "&", exp_txt[[exp]])]
  
  # ---- Apply parsed ages to data ----
  
  # Construct age group - age range dictionary
  age_dict = group_dt %>%
    # Split at denominator if it exists...
    separate_wider_delim(
      cols  = age, 
      delim = "&", 
      names = c("a1", "a2"), 
      too_few = "align_start", 
      cols_remove = FALSE) %>%
    # Extract units - years (#) or months (@)...
    mutate(u1 = str_extract(a1, "#|@"), 
           u2 = str_extract(a2, "#|@"), 
           u1 = ifelse(is.na(u1), u2, u1), 
           u2 = ifelse(is.na(u2), u1, u2)) %>%
    # Assume units are years if still trivial...
    mutate(u1 = ifelse(is.na(u1), "#", u1), 
           u2 = ifelse(is.na(u2), "#", u2)) %>%
    # Convert ages to numeric...
    mutate(a1 = str_remove(a1, "(#|@).*"), 
           a2 = str_remove(a2, "(#|@).*"), 
           a1 = suppressWarnings(as.numeric(a1)), 
           a2 = suppressWarnings(as.numeric(a2))) %>%
    # Convert all ages to year format...
    mutate(y1 = ifelse(u1 == "@", a1 / 12, a1), 
           y2 = ifelse(u2 == "@", a2 / 12, a2)) %>%
    # Deal with missing and infeasible values...
    mutate(y1 = ifelse(y1 > 100, NA, y1), 
           y2 = ifelse(is.na(y2), y1, y2)) %>%
    replace_na(list(y1 = 0, y2 = 0)) %>%
    # Ensure correct order and round...
    mutate(age_min = round(pmin(y1, y2)), 
           age_max = round(pmax(y1, y2))) %>%
    # Expand out for single age bins...
    select(age_group = group, age_min, age_max) %>%
    expand_grid(age = 0 : 100) %>%
    filter(age >= age_min, 
           age <= age_max) %>%
    # Tidy up...
    select(age_group, age) %>%
    as.data.table()
  
  # Convert to long form and distribute across age bins
  age_dt = sia_dt %>%
    rename(total_doses = doses) %>%
    # Expand with parsed age bins...
    left_join(y  = age_dict, 
              by = "age_group", 
              relationship = "many-to-many") %>%
    # Append population size...
    left_join(y  = table("wpp_pop"), 
              by = c("country", "year", "age")) %>%
    # Distribute across age bins...
    group_by(country, intervention, year, age_group) %>%
    mutate(doses = total_doses * pop / sum(pop)) %>%
    ungroup() %>%
    # Tidy up...
    select(country, intervention, year, 
           age, doses, cohort = pop) %>%
    as.data.table()
  
  # Sanity check that we haven't lost/gained doses
  dose_diff = sum(sia_dt$doses) - sum(age_dt$doses)
  if (abs(dose_diff) > 1e-6)
    stop("Age disaggregation failed")
    
  return(age_dt)
}

# ---------------------------------------------------------
# Convert number of doses to FVPs
# ---------------------------------------------------------
doses2fvps = function() {
  
  browser()
}

# ---------------------------------------------------------
# Simple wrapper around colour_scheme()
# ---------------------------------------------------------
get_colours = function(n) {
  
  # Generate colour palettes (see auxiliary.R)
  cols = colour_scheme(o$palette_sia, n = n)
  
  return(cols)
}

# ---------------------------------------------------------
# Plot duration (in days) of campaigns - can be grouped
# ---------------------------------------------------------
plot_durations = function(dt, by = NULL, zoom = TRUE) {
  
  # Calculate duration in days
  plot_dt = dt %>%
    mutate(duration = as.numeric(end_date - start_date)) %>%
    filter(duration > 0)
  
  # Check if plotting all data together
  if (is.null(by)) {
    
    # Produce single density - looks better than boxplot
    g = ggplot(plot_dt, aes(x = duration)) + 
      stat_density(adjust = 5, alpha = 0.5)
    
  } else {  # Otherwise plot by group
    
    # Produce plot by group - boxplot looks good for this
    g = ggplot(plot_dt) + 
      aes_string(x = by, y = "duration", fill = by) + 
      geom_boxplot(show.legend = FALSE)
    
    if (by == "disease") {
      
      # Get colours - one per disease
      n_cols = length(unique(plot_dt$disease))
      cols   = get_colours(n_cols)
      
      g = g + scale_fill_manual(values = cols)
    }
  }
  
  # Prettify plot
  g = g + theme_classic() + 
    theme(axis.title    = element_text(size = 20),
          axis.text     = element_text(size = 12), 
          axis.line     = element_blank(),
          panel.border  = element_rect(linewidth = 1, colour = "black", fill = NA),
          panel.spacing = unit(1, "lines"))
  
  # Check if we want a zoomed in version too
  if (zoom && !is.null(by)) {
    
    # Determine IQR for a second 'zoomed in' plot
    iqr_dt = plot_dt %>%
      group_by_at(by) %>%
      summarise(p = list(boxplot.stats(duration)$stats)) %>%
      ungroup() %>%
      unnest_wider(p, names_sep = "") %>%
      select(lb = p2, ub = p4)
    
    # Min and max IQR - across groups if need be
    iqr_lim = c(min(iqr_dt$lb), max(iqr_dt$ub))
    
    # Zoom in for the second plot
    g_zoom = g + coord_cartesian(ylim = c(0, max(iqr_lim)))  # Actually, just reset the upper limit
    
    # Place one on top of the other
    g = ggpubr::ggarrange(g, g_zoom, ncol = 1, align = "v")
  }
  
  # Save figure to file
  save_fig(g, dir = "sia_data", "SIA durations", by)
}

