#' @noRd
get_medians <- function(samples) {
  apply(samples, 1:(length(dim(samples)) - 1), stats::median)
}

#' @noRd
get_credible_interval <- function(sample, perc_ci = 0.95) {
  ndims <- length(dim(sample)) - 1
  alpha <- (1 - perc_ci) / 2
  list(
    lower = apply(sample, 1:ndims, stats::quantile, alpha),
    upper = apply(sample, 1:ndims, stats::quantile, 1 - alpha)
  )
}


#' @noRd
get_relative_precision <- function(medians, credible_interval) {
  medians / (credible_interval$upper - credible_interval$lower)
}
