#' Logit transformation
#'
#' @param x the value to be logit-transformed
#' @noRd
logit <- function(x) {
  log(x / (1 - x))
}

#' Expit transformation
#'
#' @param x the value to be expit-transformed
#' @noRd
expit <- function(x) {
  1 / (1 + exp(-x))
}
