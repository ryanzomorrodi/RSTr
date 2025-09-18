#' Run Gibbs sampler
#' @param name Name of model and corresponding folder
#' @param dir Directory where model lives
#' @param iterations Specifies number of iterations to run
#' @param .show_plots Show or hide traceplots as 
#' @param .discard_burnin If set to \code{TRUE}, won't save burn-in samples
#'
#' @export
run_sampler <- function(name, dir = tempdir(), iterations = 6000, .show_plots = TRUE, .show_progress = TRUE, .discard_burnin = FALSE) {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  if (substr(dir, nchar(dir), nchar(dir)) != "/") {
    dir <- paste0(dir, "/")
  }
  iterations = iterations - iterations %% 100
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
