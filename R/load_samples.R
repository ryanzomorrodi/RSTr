#' Load MCMC samples
#'
#' \code{load_samples()} gathers samples saved for model \code{name} in directory \code{dir}. By default, loads the rate estimate samples \code{lambda}, but any model parameters can be loaded. Users can also specify a burn-in period.
#' @param name  Name of model
#' @param dir   Directory where model lives
#' @param param Which parameter samples to load
#' @param burn  Numer of burn-in samples to discard
#'
#' @returns An \code{array} of samples from model \code{name}
#' @examples
#' # prepare truncated dataset
#' data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
#' adj_min <- list(2, 1)
#' initialize_model("test", tempdir(), data_min, adj_min, show_plots = FALSE)
#' run_sampler("test", show_plots = FALSE, show_progress = FALSE)
#' samples <- load_samples("test", tempdir()) * 1e5
#' \dontshow{
#' unlink(paste0(tempdir(), "\\test"), recursive = TRUE)
#' }
#' @export
load_samples <- function(name, dir = tempdir(), param, burn = 2000) {
  if (substr(dir, nchar(dir), nchar(dir)) != "/") {
    dir <- paste0(dir, "/")
  }
  params <- readRDS(paste0(dir, name, "/params.Rds"))
  mar <- c("lambda" = 4, "beta" = 4, "Z" = 4, "G" = 4, "Ag" = 3, "tau2" = 3, "sig2" = 3, "rho" = 2)
  if (params$model == "mstcar") mar["tau2"] = 2
  batch <- which(1:params$batch * 100 > burn)
  files <- paste0(dir, name, "/", param, "/", param, "_out_", batch, ".Rds")
  output <- abind::abind(lapply(files, readRDS), along = mar[param])
  dims <- params$dimnames
  its <- seq(burn + 10, max(batch) * 100, by = 10)
  if (param == "beta") {
    num_island <- readRDS(paste0(dir, name, "/spatial_data.Rds"))$num_island
    dimnames(output) <- list(island = 1:num_island, group = dims[[2]], time = dims[[3]], its = its)
  } else if (param %in% c("Z", "lambda")) {
    dimnames(output) <- c(dims, list(its = its))
  } else if (param == "rho") {
    dimnames(output) <- list(group = dims[[2]], its = its)
  } else if (param == "tau2") {
    if (params$model == "mstcar") {
      dimnames(output) <- list(group = dims[[2]], its = its)
    } else if (params$model %in% c("ucar", "mcar")) {
      dimnames(output) <- list(group = dims[[2]], time = dims[[3]], its = its)
    }
  } else if (param == "Ag") {
    dimnames(output) <- list(group1 = dims[[2]], group2 = dims[[2]], its = its)
  } else if (param == "G") {
    dimnames(output) <- list(group1 = dims[[2]], group2 = dims[[2]], time = dims[[3]], its = its)
  } else if (param == "sig2") {
    dimnames(output) <- list(group = dims[[2]], time = dims[[3]], its = its)
  }
  output
}