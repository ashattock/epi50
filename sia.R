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
  
  # ---- User options ----
  
  # Ignore campaigns starting (or ending??) before this date
  from_date = "1980-01-01"
  
  # Flag for throwing data warnings
  throw_warnings = FALSE
  
  # Entries to remove from data
  rm_var = list(
    status   = c("cancelled", "postponed"),
    activity = c("routine"))
  
  # Entries to set as NA
  na_var = c("unknown", "undefined", "")
  
  # ---- Load stuff ----
  
  # Load raw data
  data_raw = fread(paste0(o$pth$input, "sia_data.csv"))
  
  # Load data dictionary table
  data_dict = table("sia_dictionary") %>%
    select(-notes)
  
  # Load doses per FVP table
  data_fvps = table("sia_doses") %>%
    rename(doses_fvps = doses) %>%
    select(-notes)
  
  # ---- Format dates ----
  
  # Define date columns in raw data set
  date_cols = c("plan", "postponed", "done")
  date_type = c("start", "end")
  
  # Format y-m-d columns to date class
  date_fn = function(dt, col) {
    paste(dt[[paste0(toupper(col), "_YEAR")]], 
          dt[[paste0(toupper(col), "_MONTH")]], 
          dt[[paste0(toupper(col), "_DAY")]], 
          sep = ".")
  }
  
  # Loop through date columns to create a single variable
  dates_dt = copy(data_raw)
  for (i in date_cols) {
    for (j in date_type) {
      
      # Construct new column name
      col = paste(i, j, sep = "_")
      
      # Extract and format date, if any
      date_col = date_fn(data_raw, col)
      date_col[grepl("NA", date_col)] = ""
      
      # Apend this newly compiled column
      dates_dt[[col]] = date_col
    }
  }
  
  # TODO... do any of the entries with missing start dates have end dates??
  # If so, we could impute the start date as we do for missing end dates.
  
  # Create single start and end columns
  dates_dt %<>%
    mutate(start_date = paste0(plan_start, postponed_start, done_start), 
           end_date   = paste0(plan_end,   postponed_end,   done_end)) %>%
    mutate(start_date = format_date(start_date), 
           end_date   = format_date(end_date)) %>%
    # Remove ineligible dates...
    filter(!is.na(start_date),
           start_date >= format_date(from_date))
  
  # ---- Select data of interest ----
  
  # Select data of interest
  sia_dt = dates_dt %>%
    select(country      = ISO3_CODE,          # Country ISO3 codes
           status       = ACTIVITY_STATUS,    # Done, Ongoing, Planned, Postponed, Cancelled, Unknown
           activity     = ACTIVITY_TYPE_CODE, 
           extent       = ACTIVITY_EXTENT,
           intervention = INTERVENTION_CODE,
           target       = TARGET, 
           doses        = DOSES, 
           coverage     = COVERAGE, 
           start_date, end_date) %>%
    # Format columns...
    mutate(across(c(target, coverage, doses), as.numeric)) %>%
    mutate_if(is.character, tolower) %>%
    mutate(country = toupper(country))
  
  # Attempt to fill missing doses data with target-coverage info
  sia_dt %<>%
    mutate(doses = ifelse(is.na(doses), target * coverage, doses)) %>%
    select(-target, -coverage) %>%
    filter(!is.na(doses), doses > 0)
  
  # ---- Recode key variables ----
  
  # Remove any strings defined in rm_var
  for (var in names(rm_var))
    sia_dt = sia_dt[!get(var) %in% rm_var[[var]]]
  
  # Recode extent, and replace unknown entries with NA
  sia_dt %<>%
    inner_join(table("country")[, .(country, region, economy)],
               by = "country") %>%  
    mutate(extent = str_extract(extent, get_words("national", "subnational")), 
           across(where(is.character), ~if_else(. %in% na_var, NA, .)), 
           across(where(is.character), ~if_else(is.na(.), "unknown", .)))
  
  # ---- Recode interventions ----
  
  # TODO: Deal with partial dosing
  
  # Throw a warning if any unidentified interventions
  no_int = setdiff(sia_dt$intervention, data_dict$intervention)
  if (throw_warnings && length(no_int) > 0)
    warning("Undefined interventions: ", paste(no_int, collapse = ", "))
  
  # Join the dict to the data for disease-vaccine activities
  sia_dt = data_dict %>%
    inner_join(y  = table("d_v"), 
               by = c("disease", "vaccine")) %>%
    inner_join(y  = sia_dt, 
               by = "intervention", 
               relationship = "many-to-many") %>%
    select(-intervention, -d_v_id) %>%  # Experimenting with removing vaccine variable
    mutate_if(is.character, as.factor)
  
  # Number of doses of interest - used for sanity checking
  #
  # NOTE: These are disease-specific doses, and NOT number of injections
  n_doses = sum(sia_dt$doses)
  
  # summary(sia_dt)
  
  # ---- Colour scheme ----
  
  # Function for constructing colour vector
  colour_fn = function(var) {
    
    # All variable states
    lv = levels(sia_dt[[var]])
    na = lv == "unknown"

    # Extract colour paletter string
    pal = o[[paste0("palette_", var)]]
    
    # Construct colours (plus possible trivial grey)
    col = rep(NA, length(lv))
    col[!na] = colour_scheme(pal, n = length(lv) - sum(na))
    col[na]  = "grey80"
    
    return(col)
  }
  
  # Construct colour vectors for plotting variables
  vars = c("disease", "country", "region", "economy")
  pals = lapply(vars, colour_fn) %>% setNames(vars)
  
  # ---- Plot frequency of variables ----
  
  # Activity frequency by disease-vaccine
  freq_dt1 = sia_dt %>%
    group_by(disease, activity) %>%
    summarise(doses = sum(doses) / 1e9) %>%
    ungroup() %>%
    complete(disease, activity, 
             fill = list(doses = 0)) %>%
    as.data.table()
  
  # Activity frequency by region
  freq_dt2 = sia_dt %>%
    group_by(region, activity, extent) %>%
    summarise(doses = sum(doses) / 1e9) %>%
    ungroup() %>%
    complete(region, activity, extent, 
             fill = list(doses = 0)) %>%
    as.data.table()
  
  # Activity frequency by economy
  freq_dt3 = sia_dt %>%
    group_by(economy, activity, extent) %>%
    summarise(doses = sum(doses) / 1e9) %>%
    ungroup() %>%
    complete(economy, activity, extent, 
             fill = list(doses = 0)) %>%
    as.data.table()
  
  # Plot frequency of disease - stacked
  g1 = ggplot(freq_dt1) + 
    aes(x = activity, y = doses, fill = disease) + 
    geom_bar(stat = "identity", colour = "black", size = 0.05)
  
  # Save figures to file
  save_fig(g1, dir = "sia_data", "SIA activity by disease")
  
  # Plot activity by region
  g2 = ggplot(freq_dt2) + 
    aes(x = 1, y = doses, fill = region) + 
    geom_bar(stat = "identity", position = "dodge", 
             colour = "black", size = 0.05) + 
    facet_grid(extent ~ activity)
  
  # Plot activity by economy
  g3 = ggplot(freq_dt3) + 
    aes(x = 1, y = doses, fill = economy) + 
    geom_bar(stat = "identity", position = "dodge", 
             colour = "black", size = 0.05) + 
    facet_grid(extent ~ activity)
  
  # Save figures to file
  save_fig(g2, dir = "sia_data", "SIA activity by region")
  save_fig(g3, dir = "sia_data", "SIA activity by economy")
  
  # ---- Doses per month ----
  
  # Plot durations for all entries
  plot_durations(sia_dt)
  
  # Then same plots grouped by different variables
  plot_durations(sia_dt, "extent")
  plot_durations(sia_dt, "activity")
  plot_durations(sia_dt, "disease")
  
  # Impute missing end dates...
  
  # First compute average durations - grouped by extent
  duration_dt = sia_dt %>%
    mutate(days = as.numeric(end_date - start_date)) %>%
    group_by(extent) %>%
    summarise(dur_mean   = mean(days,   na.rm = TRUE), 
              dur_median = median(days, na.rm = TRUE)) %>%
    ungroup() %>%
    as.data.table()
  
  # Fill in any missing end dates with extent average
  sia_dt %<>%
    left_join(duration_dt, by = "extent") %>%
    mutate(fix_date = start_date + dur_mean, 
           end_date = if_else(is.na(end_date), fix_date, end_date)) %>%
    select(-dur_mean, -dur_median, -fix_date) %>%
    arrange(country, disease, start_date)
  
  # Melt monthly dates to tidy format...
  
  # Single datatable column of all possible months
  all_months_dt = seq(from = floor_date(min(sia_dt$start_date), "month"), 
                      to   = floor_date(max(sia_dt$end_date),   "month"), 
                      by   = "month") %>%
    as.character() %>%
    as_named_dt("month")
  
  # All months to distibute doses across (for which campaign has been 'run')
  run_months_dt = sia_dt %>%
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
    arrange(country, disease, month) %>%
    as.data.table()
  
  # Remove these trivial dose entries from the main datatable
  #
  # NOTE: This reduced dt would be the desired output of this exploration
  sia_dt = sia_month_dt %>%
    filter(doses > 0) %>%
    select(country, disease, activity, extent, 
           month, doses, status, region, economy)
  
  # Sanity check that we haven't changed number of doses
  if (abs(sum(sia_dt$doses) - n_doses) > 1e-3)
    stop("We seem to have gained/lost doses here")
  
  # ---- Fully Vaccinated Persons (FVP) ----
  
  fvps_dt = sia_dt %>%
    left_join(data_fvps, 
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
    left_join(data_fvps, 
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
  
  browser() # disease_table
  
  # Join VIMC disease and countries to SIA data
  vimc_dt = sia_month_dt %>%
    # Number of FVPs...
    left_join(data_fvps, 
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
}

# ---------------------------------------------------------
# Generate regular expression for finding one of several strings
# ---------------------------------------------------------
get_words = function(...) {
  
  # Concatenate word bound char and collapse to single string
  regexp = paste0("\\b", list(...), "\\b", collapse = "|")
  
  return(regexp)
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

