#' @noRd
validate_model <- function(RSTr_obj) {
  errout <- NULL
  errout <- check_data(RSTr_obj, errout)
  errout <- check_initial_values(RSTr_obj, errout)
  errout <- check_priors(RSTr_obj, errout)
  display_errors(errout)
}

#' \code{update_model()} generates further samples for model \code{name} in \code{dir}. The model used to generate samples (e.g., MSTCAR, MCAR, UCAR) along with the model's other parameters are specified in \code{*car()}.
#' @param RSTr_obj The \code{RSTr} model object to generate samples for
#' @param iterations Specifies number of iterations to run
#' @param show_plots If set to \code{FALSE}, hides traceplots
#' @param verbose If set to \code{FALSE}, hides progress bar and other messages
#' @returns No output, saves sampler output to \code{dir}
#' @examples
#' data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
#' adj_min <- list(2, 1)
#' mod_mst <- mstcar("mod_mst", data_min, adj_min, tempdir())
#' mod_mst <- update_model(mod_mst, show_plots = FALSE, verbose = FALSE)
#' \dontshow{
#' unlink(paste0(tempdir(), "\\mod_mst"), recursive = TRUE)
#' }
#' @export
update_model <- function(RSTr_obj, iterations = 6000, show_plots = TRUE, verbose = TRUE) {
  RSTr_obj <- run_sampler(RSTr_obj, iterations, show_plots, verbose)
  if (verbose) message("Generating estimates...")
  RSTr_obj <- get_estimates(RSTr_obj)
  save_model(RSTr_obj)
  if (verbose) message("Model finished at ", format(Sys.time(), "%a %b %d %X"))
  RSTr_obj
}

#' @noRd
prepare_data <- function(data) {
  if (is.null(dim(data$Y))) {
    data <- lapply(data, \(x) array(x, dim = c(length(x), 1, 1), dimnames = list(names(x))))
  } else if (length(dim(data$Y)) == 2) {
    data <- lapply(data, \(x) array(x, dim = c(dim(x), 1), dimnames = dimnames(x)))
  }
  data
}

#' @noRd
make_data_plots <- function(data) {
  data_series <- sapply(data, apply, 3, sum, na.rm = TRUE)
  data_series <- ts(data_series, start = rownames(data_series)[1])
  plot(data_series, type = "p", main = NA, nc = 2)
}

#' @noRd
get_estimates <- function(RSTr_obj) {
  name <- RSTr_obj$params$name
  dir <- RSTr_obj$params$dir
  samples <- load_samples(name, dir)
  medians <- get_medians(samples)
  credible_interval <- get_credible_interval(samples)
  relative_precision <- get_relative_precision(medians, credible_interval)
  RSTr_obj$medians <- medians
  RSTr_obj$credible_interval <- credible_interval
  RSTr_obj$relative_precision <- relative_precision
  RSTr_obj
}