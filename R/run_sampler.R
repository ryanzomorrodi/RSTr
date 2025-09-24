#' Run Gibbs sampler
#'
#' \code{run_sampler()} generates samples for model \code{name} in \code{dir}. The model used to generate samples (e.g., MSTCAR, MCAR, UCAR) along with the model's other parameters are specified in \code{initialize_model()}.
#' @param name Name of model and corresponding folder
#' @param dir Directory where model lives
#' @param iterations Specifies number of iterations to run
#' @param .show_plots If set to \code{FALSE}, hides traceplots
#' @param .show_progress If set to \code{FALSE}, hides progress bar
#' @param .discard_burnin If set to \code{TRUE}, won't save burn-in samples
#' @returns No output, saves sampler output to \code{dir}
#' @examples
#' data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
#' adj_min <- list(2, 1)
#' initialize_model("test", tempdir(), data_min, adj_min)
#' run_sampler("test", .show_plots = FALSE, .show_progress = FALSE)
#' \dontshow{
#' unlink(paste0(tempdir(), "\\test"), recursive = TRUE)
#' }
#' @export
run_sampler <- function(name, dir = tempdir(), iterations = 6000, .show_plots = TRUE, .show_progress = TRUE, .discard_burnin = FALSE) {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  if (substr(dir, nchar(dir), nchar(dir)) != "/") {
    dir <- paste0(dir, "/")
  }
  iterations <- iterations - iterations %% 100
  model <- readRDS(paste0(dir, name, "/params.Rds"))$model
  if (model == "ucar") {
    run_sampler_u(name, dir, iterations, .show_plots, .show_progress, .discard_burnin)
  }
  if (model == "mcar") {
    run_sampler_m(name, dir, iterations, .show_plots, .show_progress, .discard_burnin)
  }
  if (model == "mstcar") {
    run_sampler_mst(name, dir, iterations, .show_plots, .show_progress, .discard_burnin)
  }
}

#' Impute event values
#' @noRd
impute_missing_events <- function(data, inits, params, miss) {
  if (params$method == "binomial") {
    rate <- expit(inits$theta[miss])
    rp <- stats::runif(
      length(miss),
      stats::pbinom(params$impute_lb - 0.1, round(data$n[miss]), rate),
      stats::pbinom(params$impute_ub + 0.1, round(data$n[miss]), rate)
    )
    data$Y[miss] <- stats::qbinom(rp, round(data$n[miss]), rate)
  }
  if (params$method == "poisson") {
    rate <- data$n[miss] * exp(inits$theta[miss])
    rp <- stats::runif(
      length(miss),
      stats::ppois(params$impute_lb - 0.1, rate),
      stats::ppois(params$impute_ub + 0.1, rate)
    )
    data$Y[miss] <- stats::qpois(rp, rate)
  }
  data$Y
}

#' Tune metropolis standard deviation
#' @noRd
tune_metropolis_sd <- function(sd, accept) {
  accept <- pmin(pmax(accept, 1 / 6), 0.75)
  sd <- ifelse(accept > 0.5, sd * accept / 0.5, sd)
  sd <- ifelse(accept < 0.35, sd * accept / 0.35, sd)
  sd
}

#' Append new values to output
#' @noRd
append_to_output <- function(output, inits, output_mar) {
  onames = names(output)
  output <- lapply(names(output), \(par) abind::abind(output[[par]], inits[[par]], along = output_mar[par]))
  names(output) <- onames
  output
}

#' Append new values to plots
#' @noRd
append_to_plots <- function(plots, inits) {
  pnames = names(plots)
  plots <- lapply(names(plots), \(par) c(plots[[par]], inits[[par]][1]))
  names(plots) <- pnames
  plots
}

#' MSTCAR parameter updates
#' @noRd
update_inits_mst <- function(data, inits, spatial_data, priors, params, miss, t_accept, r_accept) {
  if (length(miss)) {
    data$Y <- impute_missing_events(data, inits, params, miss)
  }
  inits$beta <- update_beta_mst(inits, spatial_data)
  inits$Z <- update_Z_mst(inits, spatial_data)
  inits$G <- update_G_mst(inits, priors, spatial_data)
  inits$Ag <- update_Ag_mst(inits, priors)
  inits$tau2 <- update_tau2_mst(inits, priors, spatial_data)
  inits$theta <- update_theta_mst(inits, data, priors, spatial_data, params, t_accept)
  if (rho_up) {
    inits$rho <- update_rho_mst(inits, priors, spatial_data, r_accept)
  }
  inits
}