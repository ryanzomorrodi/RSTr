#' @export
print.RSTr <- function(x, ...) {
  cat("RSTr object:\n\n")
  cat("Model name:", x$params$name, "\n")
  cat("Model type:", toupper(x$params$model), "\n")
  cat("Data likelihood:", x$params$method, "\n")
  cat("Estimate Credible Interval:", paste0(round(x$params$perc_ci * 100, 2), "%"), "\n")
  cat("Number of geographic units:", dim(x$data$Y)[1], "\n")
  #cat("Model informativeness:", "TBD", "\n")
  cat("Number of samples:", x$params$total, "\n")
  cat("Estimates age-standardized:", ifelse(x$params$age_standardized, "Yes", "No"), "\n")
  if (x$params$age_standardized) cat("Age-standardized groups:", colnames(x$medians_as), "\n")
  cat("Estimates suppressed:", ifelse(x$params$suppressed, "Yes", "No"), "\n")
  if (x$params$suppressed) {
    if (x$params$age_standardized) {
      tot_reliable <- sum(!is.na(x$medians_suppressed_as))
      tot_estimate <- length(x$medians_suppressed_as)
      pct <- round(tot_reliable / tot_estimate * 100, 1)
      cat("Number of reliable age-standardized rates:", tot_reliable, "/", tot_estimate, paste0("(", pct, "%)\n"))
    } else {
      tot_reliable <- sum(!is.na(x$medians_suppressed))
      tot_estimate <- length(x$medians_suppressed)
      pct <- round(tot_reliable / tot_estimate * 100, 1)
      cat("Number of reliable rates:", tot_reliable, "/", tot_estimate, paste0("(", pct, "%)\n"))
    }
  }
}