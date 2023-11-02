###########################################################
# HISTORY
#
# xxxxxxxxxxxxxxxxx
#
###########################################################

# ---------------------------------------------------------
# Parent function for xxxxxxxxxxx
# ---------------------------------------------------------
run_history = function() {
  
  # Only continue if specified by do_step
  if (!is.element(5, o$do_step)) return()
  
  message("* Calculating impact of historical coverage")
  
  # ---- Load stuff ----

  coef_dt = read_rds("impact", "coef") 

  best_dt = read_rds("impact", "best_model")
  
  browser()
  

}

