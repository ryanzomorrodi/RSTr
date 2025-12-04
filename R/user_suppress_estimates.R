#' Suppress estimates based on reliability criteria
#' 
#' Generates suppressed estimates for an \code{RSTr} model object with a given relative precision and population threshold.
#' 
#' While the \code{pop_threshold} argument is optional, population thresholds are necessary for non-enhanced models. Population thresholds should only be omitted for enhanced CAR models, such as the EUCAR.
#' 
#' @param RSTr_obj An \code{RSTr} model object.
#' @param rates_per The desired scaling for estimate rates.
#' @param standardized If \code{RSTr_obj} contains age-standardized rates, shows the age-standardized rates. If set to \code{FALSE}, always shows the non-age-standardized rates.
#' @returns A long \code{table} containing region/group/time period names, estimates, credible intervals, relative precisions, and the associated event/population counts.
#' @examples
#' std_pop <- c(113154, 100640, 95799)
#' data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
#' adj_min <- list(2, 1)
#' mod_mst <- mstcar("test", data_min, adj_min, tempdir(), show_plots = FALSE)
#' estimates_table <- get_estimates(mod_mst)
#' mod_mst <- age_standardize(mod_mst, std_pop, "35-64")
#' estimates_table_as <- get_estimates(mod_mst)
#' \dontshow{
#' unlink(paste0(tempdir(), "\\test"), recursive = TRUE)
#' }
#' @export
suppress_estimates <- function(RSTr_obj, pop_threshold = Inf) {
  RSTr_obj$params$suppressed <- TRUE
  RSTr_obj$params$suppress_threshold <- pop_threshold
  if (is.infinite(pop_threshold) & !(RSTr_obj$params$model %in% c("eucar"))) {
    warning("Suppressing estimates without a population threshold is not recommended for non-enhanced models. Specify a `pop_threshold` or re-run with enhanced model")
  }
  if (is.finite(pop_threshold) & (RSTr_obj$params$model %in% c("eucar"))) {
    warning("Suppressing estimates with a population threshold not necessary for EUCAR models")
  }
  medians_suppressed <- RSTr_obj$medians
  medians_suppressed[(RSTr_obj$relative_precision < 1) | (RSTr_obj$data$n < pop_threshold)] <- NA
  RSTr_obj$medians_suppressed <- medians_suppressed
  if (RSTr_obj$params$age_standardized) {
    medians_suppressed_as <- RSTr_obj$medians_as
    medians_suppressed_as[(RSTr_obj$relative_precision_as < 1) | (RSTr_obj$data_as$n < pop_threshold)] <- NA
    RSTr_obj$medians_suppressed_as <- medians_suppressed_as
  }
  RSTr_obj
}