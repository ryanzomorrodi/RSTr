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
impute_missing_events <- function(Y, n, theta, miss, method, impute_lb, impute_ub) {
  if (method == "binomial") {
    rate <- expit(theta[miss])
    rp <- stats::runif(
      length(miss),
      stats::pbinom(impute_lb - 0.1, round(n[miss]), rate),
      stats::pbinom(impute_ub + 0.1, round(n[miss]), rate)
    )
    Y[miss] <- stats::qbinom(rp, round(n[miss]), rate)
  }
  if (method == "poisson") {
    rate <- n[miss] * exp(theta[miss])
    rp <- stats::runif(
      length(miss),
      stats::ppois(impute_lb - 0.1, rate),
      stats::ppois(impute_ub + 0.1, rate)
    )
    Y[miss] <- stats::qpois(rp, rate)
  }
  Y
}