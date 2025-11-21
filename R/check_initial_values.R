#' Check initial values
#'
#' @noRd
#'
check_initial_values <- function(initial_values, data, num_island, model, method) {
  message("Checking initial_values...")
  Y <- data$Y
  num_region <- dim(Y)[[1]]
  num_group <- dim(Y)[[2]]
  num_time <- dim(Y)[[3]]
  beta <- initial_values$beta
  tau2 <- initial_values$tau2
  lambda <- initial_values$lambda
  Z <- initial_values$Z
  chk <- list(
    "ucar" = c("beta", "tau2", "lambda", "Z", "sig2"),
    "mcar" = c("beta", "tau2", "lambda", "Z", "G"),
    "mstcar" = c("beta", "tau2", "lambda", "Z", "G", "Ag", "rho")
  )[[model]]
  miss <- sapply(1:length(chk), \(x) !any(names(initial_values) == chk[x]))
  if (sum(miss)) {
    stop("One or more objects missing from list 'initial_values': ", paste(chk[miss], collapse = ", "))
  }
  # Check for warnings
  warnout <- NULL
  warnct <- 0
  # Check for unused elements in 'initial_values'
  chk_elem <- which(!(names(initial_values) %in% chk))
  if (length(chk_elem)) {
    warnct <- warnct + 1
    warntxt <- paste(warnct, ": Unused elements of list 'initial_values':", paste(names(initial_values)[chk_elem], collapse = ", "))
    warnout <- c(warnout, warntxt)
  }
  if (warnct) {
    warning(paste(warnct, "warning(s) found in list 'initial_values':\n", paste(warnout, collapse = "\n ")))
  }

  # Check for errors
  errout <- NULL
  errct <- 0
  # beta
  # dimensions don't match num_island num_group num_time
  if (!all(dim(beta) == c(num_island, num_group, num_time))) {
    errct <- errct + 1
    errtxt <- paste(errct, ": beta is not an num_island x num_group x num_time array. Ensure dim(beta) == num_island x num_group x num_time or use default value")
    errout <- c(errout, errtxt)
  }
  # values are infinite
  if (any(!is.finite(beta))) {
    errct <- errct + 1
    errtxt <- paste(errct, ": beta contains infinite values. Ensure all(is.finite(beta)) or use default value")
    errout <- c(errout, errtxt)
  }
  # lambda
  # dimensions don't match num_region num_group num_time
  if (!all(dim(lambda) == dim(Y))) {
    errct <- errct + 1
    errtxt <- paste(errct, ": lambda is not a num_region x num_group x num_time array. Ensure dim(lambda) == dim(Y) or use default value")
    errout <- c(errout, errtxt)
  }
  # values are unsupported
  lower_lim <- 0
  upper_lim <- ifelse(method == "binomial", 1, Inf)
  if (any((lambda <= lower_lim) | (lambda >= upper_lim))) {
    errct <- errct + 1
    errtxt <- paste(errct, ": lambda contains unsupported values. Ensure lambdas are within range (0, 1) for `method = binomial` or (0, Inf) for `method = poisson` or use default value")
    errout <- c(errout, errtxt)
  }
  # sig2/G
  if (model == "ucar") {
    sig2 <- initial_values$sig2
    # is non-positive or infinite
    if (any(sig2 <= 0) | any(!is.finite(sig2))) {
      errct <- errct + 1
      errtxt <- paste(errct, ": sig2 contains non-positive or infinite values. Ensure all sig2 > 0 and not infinite or use default value")
      errout <- c(errout, errtxt)
    }
  } else if (model %in% c("mcar", "mstcar")) {
    G <- initial_values$G
    sig2 <- apply(G, 3, diag)
    gcor <- apply(G, 3, \(G) G[lower.tri(G)])
    # diagonals non-positive or infinite
    if (any((sig2 <= 0) | !is.finite(sig2))) {
      errct <- errct + 1
      errtxt <- paste(errct, ": Diagonals of G contain non-positive values. Ensure all diag(G) > 0 and not infinite or use default value")
      errout <- c(errout, errtxt)
    }
    # off-diagonal values are infinite
    if (any(!is.finite(gcor))) {
      errct <- errct + 1
      errtxt <- paste(errct, ": Off-diagonals of G contain infinite values. Ensure all(is.finite(G)) or use default value")
      errout <- c(errout, errtxt)
    }
  }
  # tau2
  # is non-positive or infinite
  if (any(tau2 <= 0) | any(!is.finite(tau2))) {
    errct <- errct + 1
    errtxt <- paste(errct, ": Some or all tau2 are non-positive or infinite. Ensure all tau2 > 0 and not infinite or use default value")
    errout <- c(errout, errtxt)
  }
  # Z
  # dimensions don't match num_region num_group num_time
  if (!all(dim(Z) == dim(Y))) {
    errct <- errct + 1
    errtxt <- paste(errct, ": Z is not an num_region x num_group x num_time array. Ensure dim(Z) == dim(Y) or use default value")
    errout <- c(errout, errtxt)
  }
  # values are infinite
  if (any(!is.finite(Z))) {
    errct <- errct + 1
    errtxt <- paste(errct, ": Z contains infinite values. Ensure all(is.finite(Z)) or use default value")
    errout <- c(errout, errtxt)
  }
  if (model == "mstcar") {
    rho <- initial_values$rho
    # rho
    # is non-positive or infinite
    if (any((rho <= 0) | !is.finite(rho))) {
      errct <- errct + 1
      errtxt <- paste(errct, ": rho contains non-positive values. Ensure all(rho > 0) and not infinite or use default value")
      errout <- c(errout, errtxt)
    }
    Ag <- initial_values$Ag
    # Ag
    # dimensions don't match num_group num_group
    if (!all(dim(Ag) == c(num_group, num_group))) {
      errct <- errct + 1
      errtxt <- paste(errct, ": Ag is not an num_group x num_group matrix. Ensure dim(Ag) == num_group x num_group or use default value")
      errout <- c(errout, errtxt)
    }
    # matrix is not symmetric
    if (!isSymmetric(Ag)) {
      errct <- errct + 1
      errtxt <- paste(errct, ": Ag is not symmetric. Ensure Ag is symmetric or use default value")
      errout <- c(errout, errtxt)
    }
    # values are infinite
    if (any(!is.finite(Ag))) {
      errct <- errct + 1
      errtxt <- paste(errct, ": Ag contains infinite values. Ensure Ag is finite or use default value")
      errout <- c(errout, errtxt)
    }
    # diagonals are not positive
    if (any(diag(Ag) <= 0)) {
      errct <- errct + 1
      errtxt <- paste(errct, ": diag(Ag) contains non-positive values. Ensure diag(Ag) is positive or use default value")
      errout <- c(errout, errtxt)
    }
  }
  if (errct) {
    stop(paste(errct, "error(s) found in list 'initial_values':\n", paste(errout, collapse = "\n ")))
  }
}