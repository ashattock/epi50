###########################################################
# DATABASE UTILITIES
#
# Functions required to connect to SQLite databases, pull
# the relevant data, and load to the global environment.
#
###########################################################

# ---------------------------------------------------------
# Parent function to load (possibly several) database tables to global environment
# Called by: numerous funcions
# ---------------------------------------------------------
load_tables <- function(...) {
  
  # All tables we want to load
  tables_to_load = list(...) %>%
    unlist() %>%
    setdiff(ls(envir = .GlobalEnv))  # Check if any have already been loaded
  
  # Loop through tables to load
  for (table in tables_to_load) {
    
    # Either load from cache...
    db_dt = cache_load(o, table)
    
    # ... or pull from database
    if (is.null(db_dt))
      db_dt = db_pull(table)
    
    # Assign to global environment
    assign(table, db_dt, envir = .GlobalEnv)
  }
}

# ---------------------------------------------------------
# Pull data from database for all or specific countries
# Called by: load_tables(), numerous other funcions
# ---------------------------------------------------------
db_pull <- function(table, iso3_list = NULL, append_names = F) {
  
  message("  > Loading table from database: ", table)
  
  # Append file name to differentiate in database
  db_file = paste0(table, "_", o$db_version)
  
  # Open database connection
  db_con = open_connection()
  
  # Pull the desired table from the database
  db_dt = tbl(db_con, db_file) %>%
    collect() %>%
    quiet() %>%
    as.data.table()
  
  # Close database connection
  DBI::dbDisconnect(db_con)
  
  # Save the datatable in cache
  cache_save(o, table, db_dt)
  
  # Check if we want to append country and/or vaccine details
  if (append_names) {
    
    browser()
    
    # Details we may wish to append
    country_dt = country_table[, .(country, country_name)]
    vaccine_dt = vaccine_table[, .(vaccine_id, vaccine, vaccine_name)]
    
    # Append country details
    if ("country" %in% names(db_dt))
      db_dt %<>% left_join(country_dt, by = "country")
    
    # Append vaccine details
    if ("vaccine_id" %in% names(db_dt))
      db_dt %<>% left_join(vaccine_dt, by = "vaccine_id")
  }
  
  return(db_dt)
}

# ---------------------------------------------------------
# Open database connection
# Called by: db_pull(), list_db_tables()
# ---------------------------------------------------------
open_connection <- function() {
  
  # Fetch token defined in gcs_creds.json file to authorise database query
  bigrquery::bq_auth(path = o$pth$credentials)
  
  # Create a database connecion
  db_con = DBI::dbConnect(bigrquery::bigquery(),
                          project = "vaccine-impact",
                          dataset = "data")
  
  return(db_con)
}

# ---------------------------------------------------------
#' Generates the SQLite database with input data
# ---------------------------------------------------------
gen_db <- function() {
  response <- menu(c("Yes", "No"), title = "Delete and rebuild database?")
  if (response == 1) {
    message("   - Prepping vaccine coverage...")
    source("data-raw/coverage.R")
    
    message("   - Prepping WPP inputs and all-cause deaths...")
    source("data-raw/wpp_input.R")
    
    message("   - Prepping WPP observed...")
    source("data-raw/obs_wpp.R")
    
    if (file.exists("inst/extdata/vimc_estimates.csv")) {
      message("   - Prepping VIMC impact estimates...")
      source("data-raw/vimc_impact.R")
    } else {
      warning("Add VIMC estimates to inst/extdata")
    }
    
    message("Done!")
  }
}

# ---------------------------------------------------------
# Upload table to datbase
# Called by: various functions in data-raw directory - often via gen_db()
# ---------------------------------------------------------
upload_object = function(object, name) {
  
  message("  > Uploading '", name, "' to database")
  
  # Append version number to differentiate in database
  file = paste0(name, "_", o$db_version)
  temp = paste0(file, ".csv")
  
  # Temporarily save object as a csv file
  write.csv(object, temp, row.names = FALSE)
  
  # The uploading is done with Python
  #
  # NOTE: I'm having package issues so going via bash as a workaround
  py_file = system.file("upload_file.py", package = "vieIA2030")
  sh_file = system.file("upload_file.sh", package = "vieIA2030")
  
  # Call the bash file, which itself calles the upload_file.py file
  system(paste("bash", sh_file, py_file, file, o$pth$credentials))

  # We're down with the temporary file - delete it
  file.remove(temp)
}

# ---------------------------------------------------------
# Helper function to list available database tables
# Called by: none
# ---------------------------------------------------------
list_db_tables = function() {
  
  # Open database connection
  db_con = open_connection()
  
  # Return tables accessible through this connection
  tables = DBI::dbListTables(db_con)
  
  # Close the connection
  DBI::dbDisconnect(db_con)
  
  return(tables)
}

# ---------------------------------------------------------
# Store recently pulled database table in data cache
# Called by: db_pull()
# ---------------------------------------------------------
cache_save = function(o, table, db_dt) {
  
  # Store datatable along with timestamp so we know when it was pulled
  cache_info = list(db_dt = db_dt, timestamp = Sys.time())
  
  # Save list to file in cache dir
  save_file(cache_info, o$pth$cache, table)
}

# ---------------------------------------------------------
# Attempt to load database table from cache (conditions must be met)
# Called by: load_tables()
# ---------------------------------------------------------
cache_load = function(o, table) {
  
  # NOTE: See force_db_pull and cache_hour_limit limit in options.R
  
  # Initiate trivial output
  db_dt = NULL
  
  # Path to (possibly) cached file
  cache_pth = paste0(o$pth$cache, table, ".rds")
  
  # Only continue if not forcing and cache exists
  if (!o$force_db_pull && file.exists(cache_pth)) {
    
    message("  > Loading table from cache: ", table)
    
    # Load cached table from file
    cache_info = readRDS(paste0(o$pth$cache, table, ".rds"))
    
    # Time limit for this cached table to remain valid for loading
    cache_limit = cache_info$timestamp + lubridate::hours(o$cache_hour_limit)
    
    # Assign datatable if within the time limit
    if (Sys.time() < cache_limit) {
      db_dt = cache_info$db_dt
      
      # Otherwise inform user that we'll need a new database pull
    } else {
      message("  ! Cache time limit exceeded - fresh database pull required")
    }
  }
  
  return(db_dt)  
}

