#' Logit transformation
#' @noRd
logit <- function(x) {
  log(x / (1 - x))
}

#' Expit transformation
#' @noRd
expit <- function(x) {
  1 / (1 + exp(-x))
}

#' Consolidated transformation down function
#' @noRd
log_logit <- function(x, method) {
  if (method == "binomial") logit(x) else if (method == "poisson") log(x)
}

#' Consolidated transformation up function
#' @noRd
exp_expit <- function(x, method) {
  if (method == "binomial") expit(x) else if (method == "poisson") exp(x)
}