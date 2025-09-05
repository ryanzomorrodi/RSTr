#' Calculate Medians
#' @param output Array of output to generate estimates for.
#'
#' @export
get_medians <- function(output) {
  ndims <- length(dim(output)) - 1
  apply(output, 1:ndims, stats::median)
}

#' Calculate Relative Precisions
#' @param output Array of output to generate relative precisions for.
#' @param medians Array of medians generated from output
#' @param perc_ci Number from 0 to 1. Determines width of credible interval.
#'
#' @export
get_relative_precision <- function(output, medians, perc_ci = 0.95) {
  ndims <- length(dim(output)) - 1
  alpha <- (1 - perc_ci) / 2
  hi    <- apply(output, 1:ndims, stats::quantile, 1 - alpha)
  lo    <- apply(output, 1:ndims, stats::quantile, alpha)
  medians / (hi - lo)
}

#' Calculate Credible Intervals
#' @param output Array of output to generate relative precisions for.
#' @param perc_ci Number from 0 to 1. Determines width of credible interval.
#'
#' @export
get_credible_interval <- function(output, perc_ci = 0.95) {
  ndims <- length(dim(output)) - 1
  alpha <- (1 - perc_ci) / 2
  apply(output, 1:ndims, stats::quantile, c(alpha, 1 - alpha))
}
