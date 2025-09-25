#' Load MCMC samples
#'
#' \code{load_samples()} gathers samples saved for model \code{name} in directory \code{dir}. By default, loads the rate estimate samples \code{theta}, but any model parameters can be loaded. Users can also specify a burn-in period.
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
#' run_sampler("test", show_plots = FALSE, .show_progress = FALSE)
#' theta <- load_samples("test", tempdir()) * 1e5
#' \dontshow{
#' unlink(paste0(tempdir(), "\\test"), recursive = TRUE)
#' }
#' @export
load_samples <- function(name, dir = tempdir(), param = "theta", burn = 2000) {
  if (substr(dir, nchar(dir), nchar(dir)) != "/") {
    dir <- paste0(dir, "/")
  }
  model <- readRDS(paste0(dir, name, "/params.Rds"))$model
  if (model == "ucar") {
    samples <- load_samples_u(name, dir, param, burn)
  }
  if (model == "mcar") {
    samples <- load_samples_m(name, dir, param, burn)
  }
  if (model == "mstcar") {
    samples <- load_samples_mst(name, dir, param, burn)
  }
  samples
}

#' Load MCMC samples, UCAR
#'
#' @noRd
load_samples_u <- function(name, dir, param, burn) {
  mar <- c("theta" = 2, "beta" = 2, "Z" = 2, "sig2" = 1, "tau2" = 1)
  params <- readRDS(paste0(dir, name, "/params.Rds"))
  batch <- which(1:params$batch * 100 > burn)
  if (substr(dir, nchar(dir), nchar(dir)) != "/") {
    dir <- paste0(dir, "/")
  }
  files <- paste0(dir, name, "/", param, "/", param, "_out_", batch, ".Rds")
  output <- abind::abind(lapply(files, readRDS), along = mar[param])
  if (param %in% c("theta", "beta")) {
    if (params$method == "binomial") output <- expit(output)
    if (params$method == "poisson") output <- exp(output)
  }
  dims <- params$dimnames
  if (!is.null(dims)) {
    its <- seq(burn + 10, max(batch) * 100, by = 10)
    if (param == "beta") {
      num_island <- readRDS(paste0(dir, name, "/spatial_data.Rds"))$num_island
      dimnames(output) <- list(island = 1:num_island, its = its)
    }
    if (param %in% c("Z", "theta")) {
      dimnames(output) <- list(region = dims[[1]], its = its)
    }
  }
  output
}

#' Load MCMC samples, MCAR
#'
#' @noRd
load_samples_m <- function(name, dir, param, burn) {
  mar <- c("theta" = 3, "beta" = 3, "Z" = 3, "G" = 3, "tau2" = 2)
  params <- readRDS(paste0(dir, name, "/params.Rds"))
  batch <- which(1:params$batch * 100 > burn)
  if (substr(dir, nchar(dir), nchar(dir)) != "/") {
    dir <- paste0(dir, "/")
  }
  files <- paste0(dir, name, "/", param, "/", param, "_out_", batch, ".Rds")
  output <- abind::abind(lapply(files, readRDS), along = mar[param])
  if (param %in% c("theta", "beta")) {
    if (params$method == "binomial") output <- expit(output)
    if (params$method == "poisson") output <- exp(output)
  }
  dims <- params$dimnames
  if (!is.null(dims)) {
    its <- seq(burn + 10, max(batch) * 100, by = 10)
    if (param == "beta") {
      num_island <- readRDS(paste0(dir, name, "/spatial_data.Rds"))$num_island
      dimnames(output) <- list(island = 1:num_island, group = dims[[2]], time = dims[[3]], its = its)
    }
    if (param %in% c("Z", "theta")) {
      dimnames(output) <- c(dims, list(its = its))
    }
    if (param %in% c("tau2")) {
      dimnames(output) <- list(group = dims[[2]], its = its)
    }
    if (param == "G") {
      dimnames(output) <- list(group1 = dims[[2]], group2 = dims[[2]], its = its)
    }
  }
  output
}

#' Load MCMC samples, MSTCAR
#'
#' @noRd
load_samples_mst <- function(name, dir, param, burn) {
  mar <- c("theta" = 4, "beta" = 4, "Z" = 4, "G" = 4, "Ag" = 3, "tau2" = 2, "rho" = 2)
  params <- readRDS(paste0(dir, name, "/params.Rds"))
  batch <- which(1:params$batch * 100 > burn)
  if (substr(dir, nchar(dir), nchar(dir)) != "/") {
    dir <- paste0(dir, "/")
  }
  files <- paste0(dir, name, "/", param, "/", param, "_out_", batch, ".Rds")
  output <- abind::abind(lapply(files, readRDS), along = mar[param])
  if (param %in% c("theta", "beta")) {
    if (params$method == "binomial") output <- expit(output)
    if (params$method == "poisson") output <- exp(output)
  }
  dims <- params$dimnames
  if (!is.null(dims)) {
    its <- seq(burn + 10, max(batch) * 100, by = 10)
    if (param == "beta") {
      num_island <- readRDS(paste0(dir, name, "/spatial_data.Rds"))$num_island
      dimnames(output) <- list(island = 1:num_island, group = dims[[2]], time = dims[[3]], its = its)
    }
    if (param %in% c("Z", "theta")) {
      dimnames(output) <- c(dims, list(its = its))
    }
    if (param %in% c("tau2", "rho")) {
      dimnames(output) <- list(group = dims[[2]], its = its)
    }
    if (param == "Ag") {
      dimnames(output) <- list(group1 = dims[[2]], group2 = dims[[2]], its = its)
    }
    if (param == "G") {
      dimnames(output) <- list(group1 = dims[[2]], group2 = dims[[2]], time = dims[[3]], its = its)
    }
  }
  output
}
