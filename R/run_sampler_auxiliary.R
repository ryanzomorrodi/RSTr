#' Append new values to plots
#' @noRd
update_plots <- function(plots, output, batch) {
  start <- min(batch * 100 / 2, 2000) + 10
  plots <- rbind(plots, sapply(output, extract_last_margin))
  if (start < 2000) plots <- plots[-(1:5), ]
  ts(plots, start = start, frequency = 0.1)
}

#' Append new values to output
#' @noRd
append_to_output <- function(output, RSTr_obj) {
  current_sample <- RSTr_obj$current_sample
  along <- sapply(current_sample, \(par) length(dim(par)) + 1) 
  mapply(abind::abind, output, current_sample, along = along)
}

#' Update params list inside of RSTr_obj
#' @noRd
update_params <- function(RSTr_obj, current_batch) {
  params <- RSTr_obj$params
  params$total <- params$total + 100
  params$batch <- current_batch
  params
}

#' Extract all samples for first instance of parameter
#' @noRd
extract_last_margin <- function(arr) {
  idx <- c(rep(list(1), length(dim(arr)) - 1), list(seq_len(rev(dim(arr))[1])))
  do.call(`[`, c(list(arr), idx))
}