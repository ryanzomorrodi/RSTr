#' Age-standardize samples
#' 
#' Age-standardizes samples using a standard population for either an \code{RSTr} model object or an \code{array} of samples.
#' 
#' @usage
#' ## S3 method for class 'RSTr':
#' age_standardize(object, std_pop, new_name, groups = NULL)
#' ## Default S3 method:
#' age_standardize(object, std_pop, margin, groups = NULL, bind_new = FALSE, new_name = NULL)
#' @param object An \code{RSTr} model object or an \code{array} of samples imported with \code{load_samples()}.
#' @param std_pop A vector of standard populations.
#' @param new_name The name to assign to the age-standardized group.
#' @param groups A vector of either indices for each group or a vector of strings for each group name. If set to \code{NULL}, will use all groups in the dataset.
#' @param margin The margin on which the groups of interest are stratified.
#' @param bind_new If set to \code{TRUE}, will bind to the original sample dataset. Otherwise, will generate a standalone array of samples.
#' @returns An \code{RSTr} object or an \code{array} of age-standardized samples.
#' @examples
#' std_pop <- c(113154, 100640, 95799)
#' data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
#' adj_min <- list(2, 1)
#' mod_mst <- mstcar("test", data_min, adj_min, tempdir(), show_plots = FALSE)
#' # age-standardize by all age groups
#' mod_mst <- age_standardize(mod_mst, std_pop, "35-64")
#' # Add onto age-standardized estimates. Age-standardize only by the first two age groups
#' mod_mst <- age_standardize(mod_mst, std_pop[1:2], "35-54", groups = 1:2)
#' 
#' std_pop <- c(113154, 100640, 95799)
#' age_margin <- 2
#' # age-standardize by all age groups
#' samples_3564 <- age_standardize(minsample, std_pop, margin = age_margin)
#' # age-standardize only by the first two age groups
#' samples_3554 <- age_standardize(minsample, std_pop[1:2], groups = 1:2, margin = age_margin)
#' # bind age-standardized samples to original samples
#' samples_as <- age_standardize(
#'   minsample,
#'   std_pop,
#'   new_name = "35-64",
#'   margin = age_margin,
#'   bind_new = TRUE
#' )
#' \dontshow{
#' unlink(paste0(tempdir(), "\\test"), recursive = TRUE)
#' }
#' @export
age_standardize <- function(object, std_pop, new_name, groups, ...) {
  UseMethod("age_standardize")
}

#' @export
age_standardize.default <- function(object, std_pop, new_name = NULL, groups = NULL, margin, bind_new = FALSE) {
  mar <- seq_along(dim(object))[-margin]
  wts <- std_pop / sum(std_pop)
  sub_sample <- object
  if (!is.null(groups)) {
    sub_sample <- subset_array(object, margin, groups)
  }
  if (bind_new == TRUE) {
    new_dim <- dim(object)
    new_dim[margin] <- 1
    agg_sample <- array(
      apply(sweep(sub_sample, margin, wts, "*"), mar, sum, na.rm = TRUE),
      dim <- new_dim
    )
    array_new <- abind::abind(object, agg_sample, along = margin)
    newnames <- c(dimnames(object)[[margin]], new_name)
    dimnames(array_new)[[margin]] <- newnames
  } else {
    array_new <- apply(sweep(sub_sample, margin, wts, "*"), mar, sum, na.rm = TRUE)
  }
  array_new
}

#' @export
age_standardize.RSTr <- function(object, std_pop, new_name, groups = NULL) {
  object$params$age_standardized <- TRUE
  samples <- load_samples(object)
  if (is.null(groups)) groups <- 1:dim(samples)[2]
  data <- lapply(object$data, aggregate_count, 2, groups, TRUE, new_name)
  data <- lapply(data, \(x) x[, new_name, , drop = FALSE])
  samples <- subset_array(samples, 2, groups)
  samples <- age_standardize(samples, std_pop, new_name, groups, 2, TRUE)[, new_name, , , drop = FALSE]
  medians <- get_medians(samples)
  credible_interval <- get_credible_interval(samples)
  relative_precision <- get_relative_precision(medians, credible_interval)
  object$medians_as <- erase_duplicates(abind::abind(object$medians_as, medians, along = 2))
  object$data_as$Y <- erase_duplicates(abind::abind(object$data_as$Y, data$Y, along = 2))
  object$data_as$n <- erase_duplicates(abind::abind(object$data_as$n, data$n, along = 2))
  object$credible_interval_as$lower <- erase_duplicates(abind::abind(object$credible_interval_as$lower, credible_interval$lower, along = 2))
  object$credible_interval_as$upper <- erase_duplicates(abind::abind(object$credible_interval_as$upper, credible_interval$upper, along = 2))
  object$relative_precision_as <- erase_duplicates(abind::abind(object$relative_precision_as, relative_precision, along = 2))
  object$age_metadata$names <- colnames(object$medians_as)
  object$age_metadata$std_pop[[new_name]] <- std_pop
  object$age_metadata$groups[[new_name]] <- groups
  if (object$params$suppressed) object <- suppress_estimates(object, object$params$suppress_threshold)
  object
}

#' @noRd
erase_duplicates <- function(arr) {
  arr_groups <- which(!duplicated(dimnames(arr)[[2]], fromLast = TRUE))
  arr[, arr_groups, , drop = FALSE]
}