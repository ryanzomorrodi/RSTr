#' Update model
#' 
#' \code{update_model()} generates additional samples for model \code{RSTr_obj}.
#' 
#' @param RSTr_obj The \code{RSTr} model object to generate samples for.
#' @param iterations Number of iterations to run.
#' @param show_plots If set to \code{FALSE}, hides traceplots.
#' @param verbose If set to \code{FALSE}, hides progress bar and other messages.
#' @returns An \code{RSTr} model object.
#' @examples
#' data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
#' adj_min <- list(2, 1)
#' mod_mst <- mstcar("test", data_min, adj_min, tempdir())
#' mod_mst <- update_model(mod_mst, iterations = 1000, show_plots = FALSE, verbose = FALSE)
#' \dontshow{
#' unlink(paste0(tempdir(), "\\test"), recursive = TRUE)
#' }
#' @export
update_model <- function(RSTr_obj, iterations = 6000, show_plots = TRUE, verbose = TRUE) {
  RSTr_obj <- run_sampler(RSTr_obj, iterations, show_plots, verbose)
  if (verbose) {
    if (RSTr_obj$params$age_standardized) {
      message("Age-standardized estimates detected. Updating estimates...")
    } else if (RSTr_obj$params$suppressed) {
      message("Suppressed estimates detected. Updating estimates...")
    } else {
      message("Generating estimates...")
    }
  }
  RSTr_obj <- post_sampler_output(RSTr_obj)
  save_model(RSTr_obj)
  if (RSTr_obj$params$age_standardized) {
    old_class <- class(RSTr_obj)
    RSTr_obj <- RSTr_obj[-grep(".*_as", names(RSTr_obj))]
    class(RSTr_obj) <- old_class
    for (agegroup in RSTr_obj$age_metadata$names) {
      RSTr_obj <- age_standardize(RSTr_obj, RSTr_obj$age_metadata$std_pop[[agegroup]], agegroup, RSTr_obj$age_metadata$groups[[agegroup]])
    }
    if (RSTr_obj$params$suppressed) {
      RSTr_obj <- suppress_estimates(RSTr_obj, RSTr_obj$params$suppress_threshold)
    }
  } else if (RSTr_obj$params$suppressed) {
    old_class <- class(RSTr_obj)
    RSTr_obj <- RSTr_obj[-grep(".*_suppressed", names(RSTr_obj))]
    class(RSTr_obj) <- old_class
    RSTr_obj <- suppress_estimates(RSTr_obj, RSTr_obj$params$suppress_threshold)
  }
  save_model(RSTr_obj)
  if (verbose) message("Model finished at ", format(Sys.time(), "%a %b %d %X"))
  RSTr_obj
}