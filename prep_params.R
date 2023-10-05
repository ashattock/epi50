###########################################################
# PREPARE PARAMETERS
#
# @ Austin: why are alpha and beta always 1? Starting point for
# some sort of fitting?
# 
###########################################################

prep_params <- function() {
  params_list <- lapply(d_v_at_table$d_v_at_id, function(id) {
    list(
      disease = d_v_at_table[d_v_at_id == id]$disease,
      vaccine = d_v_at_table[d_v_at_id == id]$vaccine,
      activity_type = d_v_at_table[d_v_at_id == id]$activity_type,
      alpha = 1, beta = 1, age_knots = c(2, 5, 10, 25)
    )
  })
  names(params_list) <- d_v_at_table$d_v_at_id
  jsonlite::write_json(params_list, o$pth$parameters, pretty = T)
}

