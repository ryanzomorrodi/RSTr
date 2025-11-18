#' Check initial values
#'
#' @noRd
#'
check_inits <- function(inits, data, num_island, model) {
  message("Checking inits...")
  Y <- data$Y
  num_region <- dim(Y)[[1]]
  num_group <- dim(Y)[[2]]
  num_time <- dim(Y)[[3]]
  beta <- inits$beta
  tau2 <- inits$tau2
  theta <- inits$theta
  Z <- inits$Z
  chk <- list(
    "ucar" = c("beta", "tau2", "theta", "Z", "sig2"),
    "mcar" = c("beta", "tau2", "theta", "Z", "G"),
    "mstcar" = c("beta", "tau2", "theta", "Z", "G", "Ag", "rho")
  )[[model]]
  miss <- sapply(1:length(chk), \(x) !any(names(inits) == chk[x]))
  if (sum(miss)) {
    stop("One or more objects missing from list 'inits': ", paste(chk[miss], collapse = ", "))
  }
  # Check for warnings
  warnout <- NULL
  warnct <- 0
  # Check for unused elements in 'inits'
  chk_elem <- which(!(names(inits) %in% chk))
  if (length(chk_elem)) {
    warnct <- warnct + 1
    warntxt <- paste(warnct, ": Unused elements of list 'inits':", paste(names(inits)[chk_elem], collapse = ", "))
    warnout <- c(warnout, warntxt)
  }
  if (warnct) {
    warning(paste(warnct, "warning(s) found in list 'inits':\n", paste(warnout, collapse = "\n ")))
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
  # sig2/G
  if (model == "ucar") {
    sig2 <- inits$sig2
    # is non-positive or infinite
    if (any(sig2 <= 0) | any(!is.finite(sig2))) {
      errct <- errct + 1
      errtxt <- paste(errct, ": sig2 contains non-positive or infinite values. Ensure all sig2 > 0 and not infinite or use default value")
      errout <- c(errout, errtxt)
    }
  } else if (model %in% c("mcar", "mstcar")) {
    G <- inits$G
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
  # theta
  # dimensions don't match num_region num_group num_time
  if (!all(dim(theta) == dim(Y))) {
    errct <- errct + 1
    errtxt <- paste(errct, ": theta is not a num_region x num_group x num_time array. Ensure dim(theta) == dim(Y) or use default value")
    errout <- c(errout, errtxt)
  }
  # values are infinite
  if (any(!is.finite(theta))) {
    errct <- errct + 1
    errtxt <- paste(errct, ": theta contains infinite values. Ensure all(is.finite(theta)) or use default value")
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
    rho <- inits$rho
    # rho
    # is non-positive or infinite
    if (any((rho <= 0) | !is.finite(rho))) {
      errct <- errct + 1
      errtxt <- paste(errct, ": rho contains non-positive values. Ensure all(rho > 0) and not infinite or use default value")
      errout <- c(errout, errtxt)
    }
    Ag <- inits$Ag
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
    stop(paste(errct, "error(s) found in list 'inits':\n", paste(errout, collapse = "\n ")))
  }
}