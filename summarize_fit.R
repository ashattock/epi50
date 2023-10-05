###########################################################
# SUMMARIZE FIT
#
# xxx
#
###########################################################

# ---------------------------------------------------------
# xxx
# Called by: xxx
# ---------------------------------------------------------
summarize_fit <- function(pred_all) {
  
  summary <- pred_all[strata_deaths_averted > 0,
                      .(
                        mse = mean((averted - strata_deaths_averted)**2) / mean(strata_deaths_averted**2),
                        min_rr = min(pred_rr)
                      ),
                      by = .(d_v_at_id)]
  summary <- merge(summary, d_v_at_table)
  
  return(summary)
}

# ---------------------------------------------------------
# Summarise raking factor for each vaccine-activity
# Called by: rake_impact(), run_uncertainty()
# ---------------------------------------------------------
summarize_raking <- function(impact_dt) {
  
  # Summarise raking factor for each vaccine-activity
  raking_dt = impact_dt %>%
    filter(!is.na(raking_factor), 
           pred_deaths_averted_rate != 0) %>%
    group_by(vaccine, activity_type) %>% 
    summarise(median = median(raking_factor),
              mean   = mean(raking_factor),
              min    = min(raking_factor),
              max    = max(raking_factor),
              sd     = sd(raking_factor),
              cv     = sd / mean) %>%
    ungroup() %>%
    as.data.table()
  
  return(raking_dt)
}

