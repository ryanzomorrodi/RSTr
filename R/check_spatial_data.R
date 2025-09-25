#' Check spatial data
#'
#' @noRd
check_spatial_data <- function(adjacency) {
  message("Checking spatial data...")
  # Check for errors
  errout <- NULL
  errct <- 0
  # Regions have no neighbors
  if (any(sapply(adjacency, length) == 0)) {
    errct <- errct + 1
    errtxt <- paste(errct, ": Some regions have no neighbors. Ensure all regions have at least 1 neighbor. Check vignette('RSTr-adj') for more information")
    errout <- c(errout, errtxt)
  }
  if (errct) {
    stop(paste(errct, "error(s) found in list 'data':\n", paste(errout, collapse = "\n ")))
  }
}
