#' Run Gibbs sampler
#' @useDynLib RSTr, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom RcppDist bayeslm
#' @importFrom RcppArmadillo fastLm
#' @noRd
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