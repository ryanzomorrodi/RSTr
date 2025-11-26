#' Check data
#' @noRd
check_data <- function(RSTr_obj, errout = NULL) {
  data <- RSTr_obj$data
  check_missing_data_objects(data)
  # Check for warnings
  check_unused_data(data)
  # Check for errors
  errout <- check_mismatched_dimensions(data, errout) 
  errout <- check_invalid_Y(data$Y, errout)
  errout <- check_invalid_n(data$n, errout)
  errout <- check_zero_events(RSTr_obj, errout)
  errout
}

#' Check for unused elements
#' @noRd
check_unused_data <- function(data) {
  chk_elem <- which(!(names(data) %in% c("Y", "n")))
  if (length(chk_elem)) {
    warning("Unused elements of list 'data':", paste(names(data)[chk_elem], collapse = ", "))
  }
}

#' Check for missing elements
#' @noRd
check_missing_data_objects <- function(data) {
  miss <- sapply(c("Y", "n"), \(x) !any(names(data) == x))
  if (any(miss)) {
    stop("One or more objects missing from list 'data': ", paste(chk[miss], collapse = ", "))
  }
}

#' Check mismatched Y/n dimensions
#' @noRd
check_mismatched_dimensions <- function(data, errout) {
  if (any(dim(data$Y) != dim(data$n))) {
    errtxt <- "Data not same dimensions. Ensure dim(Y) == dim(n)"
    errout <- c(errout, errtxt)
  }
  errout
}

#' Check for invalid Y values
#' @noRd
check_invalid_Y <- function(Y, errout) {
  Ychk <- Y[(!is.na(Y)) & (!is.null(Y))]
  if (any((Ychk < 0) | is.infinite(Ychk))) {
    errtxt <- ": Invalid Y values. Check that all Y's are at least 0 and finite"
    errout <- c(errout, errtxt)
  }
  errout
}

#' Check invalid n values
#' @noRd
check_invalid_n <- function(n, errout) {
  if (any((n < 0) | is.infinite(n))) {
    errtxt <- "Invalid n values. Check that all n's are at least 0 and finite"
    errout <- c(errout, errtxt)
  }
  errout
}

#' Check for zero events
#' @noRd
check_zero_events <- function(RSTr_obj, errout) {
  UseMethod("check_zero_events")
}

#' Check for zero events, ucar
#' @noRd
check_zero_events.ucar <- function(RSTr_obj, errout) {
  if (any(apply(RSTr_obj$data$Y, 2:3, sum) == 0)) {
    errtxt <- "At least one set of regions has no events. Ensure that Y has at least one event for each set of regions"
    errout <- c(errout, errtxt)
  }
  errout
}

#' Check for zero events, mcar
#' @noRd
check_zero_events.mcar <- function(RSTr_obj, errout) {
  if (any(apply(RSTr_obj$data$Y, 3, sum) == 0)) {
    errtxt <- "No events in Y for at least one time period. Ensure that Y has at least one event for each time period"
    errout <- c(errout, errtxt)
  }
  errout
}

#' Check for zero events, mstcar
#' @noRd
check_zero_events.mstcar <- function(RSTr_obj, errout) {
  if (sum(RSTr_obj$data$Y) == 0) {
    errtxt <- "No events in Y. Ensure that Y has at least one event"
    errout <- c(errout, errtxt)
  }
  errout
}