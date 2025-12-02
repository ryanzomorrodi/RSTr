#' Run Gibbs sampler
#'
#' \code{run_sampler()} generates samples for model \code{name} in \code{dir}. The model used to generate samples (e.g., MSTCAR, MCAR, UCAR) along with the model's other parameters are specified in \code{*car()}.
#' @param name Name of model and corresponding folder
#' @param dir Directory where model lives
#' @param iterations Specifies number of iterations to run
#' @param show_plots If set to \code{FALSE}, hides traceplots
#' @param verbose If set to \code{FALSE}, hides progress bar
#' @returns No output, saves sampler output to \code{dir}
#' @examples
#' data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
#' adj_min <- list(2, 1)
#' mstcar("test", data_min, adj_min, tempdir())
#' run_sampler("test", show_plots = FALSE, verbose = FALSE)
#' \dontshow{
#' unlink(paste0(tempdir(), "\\test"), recursive = TRUE)
#' }
#' @useDynLib RSTr, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom RcppDist bayeslm
#' @importFrom RcppArmadillo fastLm
#' @export
run_sampler <- function(RSTr_obj, iterations = 6000, show_plots = TRUE, verbose = TRUE) {
  iterations <- iterations - iterations %% 100
  sampler_start <- Sys.time()
  missing_Y <- RSTr_obj$params$missing_Y
  start_batch <- RSTr_obj$params$batch
  total <- RSTr_obj$params$total
  method <- RSTr_obj$params$method
  batches <- seq(start_batch + 1, start_batch + iterations / 100)

  if (verbose) message("Starting sampler on Batch ", start_batch + 1, " at ", format(Sys.time(), "%a %b %d %X"))
  for (batch in batches) {
    if (verbose) display_progress(batch, max(batches), total, 0, sampler_start)
    output <- setNames(vector("list", length(RSTr_obj$current_sample)), names(RSTr_obj$current_sample))
    RSTr_obj$current_sample$lambda <- log_logit(RSTr_obj$current_sample$lambda, method)
    for (it in 1:100) {
      if (missing_Y) RSTr_obj <- impute_missing_data(RSTr_obj)
      RSTr_obj <- update_current_sample(RSTr_obj)
      if (it %% 10 == 0) output <- append_to_output(output, RSTr_obj)
      if (verbose) display_progress(batch, max(batches), total, it, sampler_start)
    }
    output <- prepare_output(output, method)
    RSTr_obj$current_sample$lambda <- exp_expit(RSTr_obj$current_sample$lambda, method)
    RSTr_obj <- update_priors_sd(RSTr_obj)
    RSTr_obj <- update_params(RSTr_obj, batch)
    save_model(RSTr_obj)
    save_output(output, batch, RSTr_obj$params$dir, RSTr_obj$params$name)

    if (show_plots) {
      if (!exists("plots")) plots <- NULL
      plots <- update_plots(plots, output, batch)
      plot(plots, xlab = "Iterations", main = "Traceplots")
    } 
  }
  RSTr_obj
}