#' Check initial values
#' @noRd
check_initial_values <- function(RSTr_obj, errout = NULL) {
  UseMethod("check_initial_values")
}

#' Check initial values UCAR
#' @noRd
check_initial_values.ucar <- function(RSTr_obj, errout = NULL) {
  Y <- RSTr_obj$data$Y
  method <- RSTr_obj$params$method
  num_region <- dim(Y)[[1]]
  num_group <- dim(Y)[[2]]
  num_time <- dim(Y)[[3]]
  num_island <- RSTr_obj$spatial_data$num_island
  chk <- c("beta", "tau2", "lambda", "Z", "sig2")
  check_missing_initial_values(RSTr_obj, chk)
  # Check for warnings
  check_unused_initial_values(RSTr_obj, chk)
  # Check for errors
  errout <- check_beta(RSTr_obj$initial_values$beta, num_island, num_group, num_time, errout)
  errout <- check_lambda(RSTr_obj$initial_values$lambda, Y, method, errout)
  errout <- check_sig2(RSTr_obj$initial_values$sig2, errout)
  errout <- check_tau2(RSTr_obj$initial_values$tau2, errout)
  errout <- check_Z(RSTr_obj$initial_values$Z, Y, errout)
  errout
}

#' Check initial values MCAR
#' @noRd
check_initial_values.mcar <- function(RSTr_obj, errout = NULL) {
  Y <- RSTr_obj$data$Y
  method <- RSTr_obj$params$method
  num_region <- dim(Y)[[1]]
  num_group <- dim(Y)[[2]]
  num_time <- dim(Y)[[3]]
  num_island <- RSTr_obj$spatial_data$num_island
  chk <- c("beta", "tau2", "lambda", "Z", "G")
  check_missing_initial_values(RSTr_obj, chk)
  # Check for warnings
  check_unused_initial_values(RSTr_obj, chk)

  # Check for errors
  errout <- check_beta(RSTr_obj$initial_values$beta, num_island, num_group, num_time, errout)
  errout <- check_lambda(RSTr_obj$initial_values$lambda, Y, method, errout)
  errout <- check_G(RSTr_obj$initial_values$G, errout)
  errout <- check_tau2(RSTr_obj$initial_values$tau2, errout)
  errout <- check_Z(RSTr_obj$initial_values$Z, Y, errout)
  errout
}

#' Check initial values MSTCAR
#' @noRd
check_initial_values.mstcar <- function(RSTr_obj, errout = NULL) {
  Y <- RSTr_obj$data$Y
  method <- RSTr_obj$params$method
  num_region <- dim(Y)[[1]]
  num_group <- dim(Y)[[2]]
  num_time <- dim(Y)[[3]]
  num_island <- RSTr_obj$spatial_data$num_island
  chk <- c("beta", "tau2", "lambda", "Z", "G", "rho", "Ag")
  check_missing_initial_values(RSTr_obj, chk)
  # Check for warnings
  check_unused_initial_values(RSTr_obj, chk)
  # Check for errors
  errout <- check_beta(RSTr_obj$initial_values$beta, num_island, num_group, num_time, errout)
  errout <- check_lambda(RSTr_obj$initial_values$lambda, Y, method, errout)
  errout <- check_G(RSTr_obj$initial_values$G, errout)
  errout <- check_rho(RSTr_obj$initial_values$rho, errout)
  errout <- check_tau2(RSTr_obj$initial_values$tau2, errout)
  errout <- check_Z(RSTr_obj$initial_values$Z, Y, errout)
  errout <- check_Ag(RSTr_obj$initial_values$Ag, errout)
  errout
}

#' Check for missing elements
#' @noRd
check_missing_initial_values <- function(RSTr_obj, chk) {
  miss <- sapply(1:length(chk), \(x) !any(names(RSTr_obj$initial_values) == chk[x]))
  if (sum(miss)) {
    stop("One or more objects missing from list 'initial_values': ", paste(chk[miss], collapse = ", "))
  }
}

#' Check for unused elements
#' @noRd
check_unused_initial_values <- function(RSTr_obj, chk) {
  chk_elem <- which(!(names(RSTr_obj$initial_values) %in% chk))
  if (length(chk_elem)) {
    warning(paste("Unused elements of list 'initial_values':", paste(names(RSTr_obj$initial_values)[chk_elem], collapse = ", ")))
  }
}

#' Check beta
#' @noRd
check_beta <- function(beta, num_island, num_group, num_time, errout) {
  # dimensions don't match num_island num_group num_time
  if (!all(dim(beta) == c(num_island, num_group, num_time))) {
    errtxt <- "beta is not an num_island x num_group x num_time array. Ensure dim(beta) == num_island x num_group x num_time or use default value"
    errout <- c(errout, errtxt)
  }
  # values are infinite
  if (any(!is.finite(beta))) {
    errtxt <- "beta contains infinite values. Ensure all(is.finite(beta)) or use default value"
    errout <- c(errout, errtxt)
  }
}

#' Check lambda
#' @noRd
check_lambda <- function(lambda, Y, method, errout) {
  # dimensions don't match num_region num_group num_time
  if (!all(dim(lambda) == dim(Y))) {
    errtxt <- "lambda is not a num_region x num_group x num_time array. Ensure dim(lambda) == dim(Y) or use default value"
    errout <- c(errout, errtxt)
  }
  # values are unsupported
  lower_lim <- 0
  upper_lim <- ifelse(method == "binomial", 1, Inf)
  if (any((lambda <= lower_lim) | (lambda >= upper_lim))) {
    errtxt <- "lambda contains unsupported values. Ensure lambdas are within range (0, 1) for `method = binomial` or (0, Inf) for `method = poisson` or use default value"
    errout <- c(errout, errtxt)
  }
  errout
}

#' Check sig2
#' @noRd
check_sig2 <- function(sig2, errout) {
  # is non-positive or infinite
  if (any(sig2 <= 0) | any(!is.finite(sig2))) {
    errtxt <- "sig2 contains non-positive or infinite values. Ensure all sig2 > 0 and not infinite or use default value"
    errout <- c(errout, errtxt)
  }
  errout
}

#' Check tau2
#' @noRd
check_tau2 <- function(tau2, errout) {
  # is non-positive or infinite
  if (any(tau2 <= 0) | any(!is.finite(tau2))) {
    errtxt <- "Some or all tau2 are non-positive or infinite. Ensure all tau2 > 0 and not infinite or use default value"
    errout <- c(errout, errtxt)
  }
  errout
}

#' Check Z
#' @noRd
check_Z <- function(Z, Y, errout) {
  # dimensions don't match num_region num_group num_time
  if (!all(dim(Z) == dim(Y))) {
    errtxt <- "Z is not an num_region x num_group x num_time array. Ensure dim(Z) == dim(Y) or use default value"
    errout <- c(errout, errtxt)
  }
  # values are infinite
  if (any(!is.finite(Z))) {
    errtxt <- "Z contains infinite values. Ensure all(is.finite(Z)) or use default value"
    errout <- c(errout, errtxt)
  }
  errout
}

#' Check G
#' @noRd
check_G <- function(G, errout) {
  sig2 <- apply(G, 3, diag)
  gcor <- apply(G, 3, \(G) G[lower.tri(G)])
  # diagonals non-positive or infinite
  if (any((sig2 <= 0) | !is.finite(sig2))) {
    errtxt <- "Diagonals of G contain non-positive values. Ensure all diag(G) > 0 and not infinite or use default value"
    errout <- c(errout, errtxt)
  }
  # off-diagonal values are infinite
  if (any(!is.finite(gcor))) {
    errtxt <- "Off-diagonals of G contain infinite values. Ensure all(is.finite(G)) or use default value"
    errout <- c(errout, errtxt)
  }
  errout
}

#' Check rho
#' @noRd
check_rho <- function(rho, errout) {
  # is non-positive or infinite
  if (any((rho <= 0) | !is.finite(rho))) {

    errtxt <- "rho contains non-positive values. Ensure all(rho > 0) and not infinite or use default value"
    errout <- c(errout, errtxt)
  }
  errout
}

#' Check Ag
#' @noRd
check_Ag <- function(Ag, errout) {
  # matrix is not symmetric
  if (!isSymmetric(Ag)) {
    errtxt <- "Ag is not symmetric. Ensure Ag is symmetric or use default value"
    errout <- c(errout, errtxt)
  }
  # values are infinite
  if (any(!is.finite(Ag))) {
    errtxt <- "Ag contains infinite values. Ensure Ag is finite or use default value"
    errout <- c(errout, errtxt)
  }
  # diagonals are not positive
  if (any(diag(Ag) <= 0)) {
    errtxt <- "diag(Ag) contains non-positive values. Ensure diag(Ag) is positive or use default value"
    errout <- c(errout, errtxt)
  }
  errout
}

#' Display errors
#' @noRd
display_errors <- function(errout) {
  if (length(errout)) {
    stop(paste0(length(errout), "error(s) found:\n", paste(errout, collapse = "\n ")))
  }
}