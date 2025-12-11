get_priors <- function(RSTr_obj, priors) {
  UseMethod("get_priors")
}

#' @export
get_priors.ucar <- function(RSTr_obj, priors) {
  if (is.null(priors$lambda_sd)) priors$lambda_sd <- array(0.025, dim = dim(RSTr_obj$data$Y))
  priors$lambda_accept <- array(0, dim = dim(priors$lambda_sd))
  if (is.null(priors$tau_a)) priors$tau_a <- 0.001
  if (is.null(priors$tau_b)) priors$tau_b <- 0.001
  if (is.null(priors$sig_a)) priors$sig_a <- 0.001
  if (is.null(priors$sig_b)) priors$sig_b <- 0.001
  RSTr_obj$priors <- priors
  RSTr_obj
}

#' @export
get_priors.mcar <- function(RSTr_obj, priors) {
  num_group <- dim(RSTr_obj$data$Y)[2]
  if (is.null(priors$lambda_sd)) priors$lambda_sd <- array(0.025, dim = dim(RSTr_obj$data$Y))
  priors$lambda_accept <- array(0, dim = dim(priors$lambda_sd))
  if (is.null(priors$tau_a)) priors$tau_a <- 0.001
  if (is.null(priors$tau_b)) priors$tau_b <- 0.001
  if (is.null(priors$G_df)) priors$G_df <- num_group + 2
  if (is.null(priors$G_scale)) priors$G_scale <- diag(1 / 7, num_group)
  RSTr_obj$priors <- priors
  RSTr_obj
}

#' @export
get_priors.mstcar <- function(RSTr_obj, priors) {
  num_group <- dim(RSTr_obj$data$Y)[2]
  if (is.null(priors$lambda_sd)) priors$lambda_sd <- array(0.025, dim = dim(RSTr_obj$data$Y))
  priors$lambda_accept <- array(0, dim = dim(priors$lambda_sd))
  if (is.null(priors$tau_a)) priors$tau_a <- 0.001
  if (is.null(priors$tau_b)) priors$tau_b <- 0.001
  if (is.null(priors$G_df)) priors$G_df <- num_group + 2
  if (is.null(priors$G_scale)) priors$G_scale <- diag(1 / 7, num_group)
  if (is.null(priors$Ag_scale)) priors$Ag_scale <- diag(1 / 7, num_group)
  if (is.null(priors$Ag_df)) priors$Ag_df <- num_group + 2
  RSTr_obj$priors <- priors
  RSTr_obj
}

#' @export
get_priors.mstcar_update_rho <- function(RSTr_obj, priors) {
  num_group <- dim(RSTr_obj$data$Y)[2]
  if (is.null(priors$lambda_sd)) priors$lambda_sd <- array(0.025, dim = dim(RSTr_obj$data$Y))
  priors$lambda_accept <- array(0, dim = dim(priors$lambda_sd))
  if (is.null(priors$tau_a)) priors$tau_a <- 0.001
  if (is.null(priors$tau_b)) priors$tau_b <- 0.001
  if (is.null(priors$G_df)) priors$G_df <- num_group + 2
  if (is.null(priors$G_scale)) priors$G_scale <- diag(1 / 7, num_group)
  if (is.null(priors$Ag_scale)) priors$Ag_scale <- diag(1 / 7, num_group)
  if (is.null(priors$Ag_df)) priors$Ag_df <- num_group + 2
  if (is.null(priors$rho_a)) priors$rho_a <- 95
  if (is.null(priors$rho_b)) priors$rho_b <- 5
  if (is.null(priors$rho_sd)) priors$rho_sd <- matrix(0.05, num_group, 1)
  priors$rho_accept <- array(0, dim = dim(priors$rho_sd))
  RSTr_obj$priors <- priors
  RSTr_obj
}
