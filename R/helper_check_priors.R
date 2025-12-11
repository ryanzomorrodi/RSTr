#' @noRd
check_priors <- function(RSTr_obj, errout) {
  UseMethod("check_priors")
}

#' @noRd
check_priors.ucar <- function(RSTr_obj, errout) {
  priors <- RSTr_obj$priors
  chk <- c("tau_a", "tau_b", "lambda_sd", "lambda_accept", "sig_a", "sig_b")
  check_missing_priors(priors, chk)
  # Check for warnings
  check_unused_priors(priors, chk)
  # Check for errors
  errout <- check_tau_a(priors$tau_a, errout)
  errout <- check_tau_b(priors$tau_b, errout)
  errout <- check_lambda_sd(priors$lambda_sd, RSTr_obj$data$Y, errout)
  errout <- check_sig_a(priors$sig_a, errout)
  errout <- check_sig_b(priors$sig_b, errout)
  errout
}

#' @noRd
check_priors.mcar <- function(RSTr_obj, errout) {
  priors <- RSTr_obj$priors
  num_group <- dim(RSTr_obj$data$Y)[2]
  chk <- c("tau_a", "tau_b", "lambda_sd", "lambda_accept", "G_scale", "G_df")
  check_missing_priors(priors, chk)
  # Check for warnings
  check_unused_priors(priors, chk)
  # Check for errors
  errout <- check_tau_a(priors$tau_a, errout)
  errout <- check_tau_b(priors$tau_b, errout)
  errout <- check_lambda_sd(priors$lambda_sd, RSTr_obj$data$Y, errout)
  errout <- check_G_scale(priors$G_scale, num_group, errout)
  errout <- check_G_df(priors$G_df, num_group, errout)
  errout
}

#' @noRd
check_priors.mstcar <- function(RSTr_obj, errout) {
  priors <- RSTr_obj$priors
  num_group <- dim(RSTr_obj$data$Y)[2]
  chk <- c("tau_a", "tau_b", "lambda_sd", "lambda_accept", "G_scale", "G_df", "Ag_scale", "Ag_df")
  check_missing_priors(priors, chk)
  # Check for warnings
  check_unused_priors(priors, chk)
  # Check for errors
  errout <- check_tau_a(priors$tau_a, errout)
  errout <- check_tau_b(priors$tau_b, errout)
  errout <- check_lambda_sd(priors$lambda_sd, RSTr_obj$data$Y, errout)
  errout <- check_G_scale(priors$G_scale, num_group, errout)
  errout <- check_G_df(priors$G_df, num_group, errout)
  errout <- check_Ag_scale(priors$Ag_scale, num_group, errout)
  errout <- check_Ag_df(priors$Ag_df, num_group, errout)
  errout
}

#' @noRd
check_priors.mstcar_update_rho <- function(RSTr_obj, errout) {
  priors <- RSTr_obj$priors
  num_group <- dim(RSTr_obj$data$Y)[2]
  chk <- c("tau_a", "tau_b", "lambda_sd", "lambda_accept", "G_scale", "G_df", "Ag_scale", "Ag_df", "rho_a", "rho_b", "rho_sd", "rho_accept")
  check_missing_priors(priors, chk)
  # Check for warnings
  check_unused_priors(priors, chk)
  # Check for errors
  errout <- check_tau_a(priors$tau_a, errout)
  errout <- check_tau_b(priors$tau_b, errout)
  errout <- check_lambda_sd(priors$lambda_sd, RSTr_obj$data$Y, errout)
  errout <- check_G_scale(priors$G_scale, num_group, errout)
  errout <- check_G_df(priors$G_df, num_group, errout)
  errout <- check_Ag_scale(priors$Ag_scale, num_group, errout)
  errout <- check_Ag_df(priors$Ag_df, num_group, errout)
  errout <- check_rho_a(priors$rho_a, errout)
  errout <- check_rho_b(priors$rho_b, errout)
  errout <- check_rho_sd(priors$rho_sd, num_group, errout)
  errout
}

#' @noRd
check_missing_priors <- function(priors, chk) {
  miss <- sapply(chk, \(x) !any(names(priors) == x))
  if (any(miss)) {
    stop("One or more objects missing from list 'priors': ", toString(chk[miss]))
  }
}

#' @noRd
check_unused_priors <- function(priors, chk) {
  chk_elem <- !(names(priors) %in% chk)
  if (any(chk_elem)) {
    warning("Unused elements of list 'priors':", toString(names(priors)[chk_elem]))
  }
}

#' @noRd
check_tau_a <- function(tau_a, errout) {
  # is non-positive or infinite
  if ((tau_a <= 0) || !is.finite(tau_a)) {
    errtxt <- "tau_a is not positive. Ensure tau_a > 0 and not infinite or use default value"
    errout <- c(errout, errtxt)
  }
  errout
}

#' @noRd
check_tau_b <- function(tau_b, errout) {
  # is non-positive or infinite
  if ((tau_b <= 0) || !is.finite(tau_b)) {
    errtxt <- "tau_b is not positive. Ensure tau_b > 0 and not infinite or use default value"
    errout <- c(errout, errtxt)
  }
  errout
}

#' @noRd
check_lambda_sd <- function(lambda_sd, Y, errout) {
  # dim not num_time num_region
  if (!all(dim(lambda_sd) == dim(Y))) {
    errtxt <- "lambda_sd has different length than data. Ensure length(lambda_sd) == length(Y) or use default value"
    errout <- c(errout, errtxt)
  }
  # is non-positive or infinite
  if (any((lambda_sd <= 0) | !is.finite(lambda_sd))) {
    errtxt <- "lambda_sd contains non-positive values. Ensure all(lambda_sd > 0) and not infinite or use default value"
    errout <- c(errout, errtxt)
  }
  errout
}

#' @noRd
check_sig_a <- function(sig_a, errout) {
  # is non-positive or infinite
  if ((sig_a <= 0) || !is.finite(sig_a)) {
    errtxt <- "sig_a is not positive. Ensure sig_a > 0 and not infinite or use default value"
    errout <- c(errout, errtxt)
  }
  errout
}

#' @noRd
check_sig_b <- function(sig_b, errout) {
  # is non-positive or infinite
  if ((sig_b <= 0) || !is.finite(sig_b)) {
    errtxt <- "sig_b is not positive. Ensure sig_b > 0 and not infinite or use default value"
    errout <- c(errout, errtxt)
  }
  errout
}

#' @noRd
check_G_scale <- function(G_scale, num_group, errout) {
  # dimensions don't match num_group num_group
  if (!all(dim(G_scale) == c(num_group, num_group))) {
    errtxt <- "G_scale is not an num_group x num_group matrix. Ensure dim(G_scale) == num_group x num_group or use default value"
    errout <- c(errout, errtxt)
  }
  # matrix is not symmetric
  if (!isSymmetric(G_scale)) {
    errtxt <- "G_scale is not symmetric. Ensure G_scale is symmetric or use default value"
    errout <- c(errout, errtxt)
  }
  # values are infinite
  if (!all(is.finite(G_scale))) {
    errtxt <- "G_scale contains infinite values. Ensure G_scale is finite or use default value"
    errout <- c(errout, errtxt)
  }
  # diagonals are not positive
  if (any(diag(G_scale) <= 0)) {
    errtxt <- "diag(G_scale) contains non-positive values. Ensure diag(G_scale) is all positive or use default value"
    errout <- c(errout, errtxt)
  }
  errout
}

#' @noRd
check_G_df <- function(G_df, num_group, errout) {
  # is not a whole number
  if (floor(G_df) != G_df) {
    errtxt <- "G_df is not a whole number. Ensure G_df is whole number or use default value"
    errout <- c(errout, errtxt)
  }
  # is less than df
  if (G_df <= num_group - 1) {
    errtxt <- "G_df too small. Ensure G_df > num_group - 1 or use default value"
    errout <- c(errout, errtxt)
  }
  errout
}

#' @noRd
check_Ag_scale <- function(Ag_scale, num_group, errout) {
  # Ag_scale
  # dimensions don't match num_group num_group
  if (!all(dim(Ag_scale) == c(num_group, num_group))) {
    errtxt <- "Ag_scale is not an num_group x num_group matrix. Ensure dim(Ag_scale) == num_group x num_group or use default value"
    errout <- c(errout, errtxt)
  }
  # matrix is not symmetric
  if (!isSymmetric(Ag_scale)) {
    errtxt <- "Ag_scale is not symmetric. Ensure Ag_scale is symmetric or use default value"
    errout <- c(errout, errtxt)
  }
  # values are infinite
  if (!all(is.finite(Ag_scale))) {
    errtxt <- "Ag_scale contains infinite values. Ensure Ag_scale is finite or use default value"
    errout <- c(errout, errtxt)
  }
  # diagonals are not positive
  if (any(diag(Ag_scale) <= 0)) {
    errtxt <- "diag(Ag_scale) contains non-positive values. Ensure diag(Ag_scale) is positive or use default value"
    errout <- c(errout, errtxt)
  }
  errout
}

#' @noRd
check_Ag_df <- function(Ag_df, num_group, errout) {
  # is not a whole number
  if (floor(Ag_df) != Ag_df) {
    errtxt <- "Ag_df is not a whole number. Ensure Ag_df is whole number or use default value"
    errout <- c(errout, errtxt)
  }
  # is less than df
  if (Ag_df <= num_group - 1) {
    errtxt <- "Ag_df too small. Ensure Ag_df > num_group - 1 or use default value"
    errout <- c(errout, errtxt)
  }
  errout
}

#' @noRd
check_rho_a <- function(rho_a, errout) {
  # is non-positive or infinite
  if ((rho_a <= 0) || !is.finite(rho_a)) {
    errtxt <- "rho_a is not positive. Ensure rho_a > 0 and not infinite or use default value"
    errout <- c(errout, errtxt)
  }
  errout
}

#' @noRd
check_rho_b <- function(rho_b, errout) {
  # is non-positive or infinite
  if ((rho_b <= 0) || !is.finite(rho_b)) {
    errtxt <- "rho_b is not positive. Ensure rho_b > 0 and not infinite or use default value"
    errout <- c(errout, errtxt)
  }
  errout
}

#' @noRd
check_rho_sd <- function(rho_sd, num_group, errout) {
  # length not num_group
  if (length(rho_sd) != num_group) {
    errtxt <- "rho_sd is not length num_group. Ensure length(rho_sd) == num_group or use default value"
    errout <- c(errout, errtxt)
  }
  # is non-positive or infinite
  if (any((rho_sd <= 0) | !is.finite(rho_sd))) {
    errtxt <- "rho_sd contains non-positive values. Ensure all(rho_sd > 0) and not infinite or use default value"
    errout <- c(errout, errtxt)
  }
  errout
}