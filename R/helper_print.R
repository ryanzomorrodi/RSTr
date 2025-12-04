#' @export
print.RSTr <- function(RSTr_obj) {
  cat("RSTr object:\n\n")
  cat("Model name:", RSTr_obj$params$name, "\n")
  cat("Model type:", toupper(RSTr_obj$params$model), "\n")
  cat("Data likelihood:", RSTr_obj$params$method, "\n")
  cat("Number of geographic units:", dim(RSTr_obj$data$Y)[1], "\n")
  #cat("Model informativeness:", "TBD", "\n")
  cat("Number of samples:", RSTr_obj$params$total, "\n")
  cat("Estimates age-standardized:", ifelse(RSTr_obj$params$age_standardized, "Yes", "No"), "\n")
  if (RSTr_obj$params$age_standardized) cat("Age-standardized groups:", colnames(RSTr_obj$medians_as), "\n")
  cat("Estimates suppressed:", ifelse(RSTr_obj$params$suppressed, "Yes", "No"), "\n")
  if (RSTr_obj$params$suppressed) {
    if (RSTr_obj$params$age_standardized) {
      tot_reliable <- sum(!is.na(RSTr_obj$medians_suppressed_as))
      tot_estimate <- length(RSTr_obj$medians_suppressed_as)
      pct <- round(tot_reliable / tot_estimate * 100, 1)
      cat("Number of reliable age-standardized rates:", tot_reliable, "/", tot_estimate, paste0("(", pct, "%)\n"))
    } else {
      tot_reliable <- sum(!is.na(RSTr_obj$medians_suppressed))
      tot_estimate <- length(RSTr_obj$medians_suppressed)
      pct <- round(tot_reliable / tot_estimate * 100, 1)
      cat("Number of reliable rates:", tot_reliable, "/", tot_estimate, paste0("(", pct, "%)\n"))
    }
  }
}