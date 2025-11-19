#' Get initial values
#'
#' @noRd
get_initial_values <- function(initial_values, data, spatial_data, model, method, ignore_checks) {
  Y <- data$Y
  n <- data$n
  island_id <- spatial_data$island_id
  num_region <- dim(Y)[[1]]
  num_group <- dim(Y)[[2]]
  num_time <- dim(Y)[[3]]
  num_island <- length(unique(island_id))
  # Prepare initial values
  initmiss <- NULL
  # beta
  if (is.null(initial_values$beta)) {
    beta <- apply(Y, 2:3, sum, na.rm = TRUE) / apply(n, 2:3, sum)
    if (method == "poisson") {
      beta <- array(log(beta), dim = c(num_group, num_time, num_island))
      beta[!is.finite(beta)] <- log(sum(Y, na.rm = TRUE) / sum(n))
    }
    if (method == "binomial") {
      beta <- array(logit(beta), dim = c(num_group, num_time, num_island))
      beta[!is.finite(beta)] <- logit(sum(Y, na.rm = TRUE) / sum(n))
    }
    beta <- aperm(beta, c(3, 1, 2))
    initial_values$beta <- beta
    initmiss <- c(initmiss, "beta")
  }
  # theta
  if (is.null(initial_values$theta)) {
    if (method == "poisson") theta <- log(Y / n)
    if (method == "binomial") theta <- logit(Y / n)
    theta[!is.finite(theta)] <- beta[island_id + 1, , ][which(!is.finite(theta))]
    initial_values$theta <- theta
    initmiss <- c(initmiss, "theta")
  }
  # Z
  if (is.null(initial_values$Z)) {
    initial_values$Z <- initial_values$theta - initial_values$beta[island_id + 1, , , drop = FALSE]
    initmiss <- c(initmiss, "Z")
  }
  # tau2
  if (is.null(initial_values$tau2)) {
    tau2_cols <- ifelse(model == "mstcar", 1, num_time)
    initial_values$tau2 <- matrix(1 / 100, num_group, tau2_cols)
    initmiss <- c(initmiss, "tau2")
  }
  # sig2/G
  if (model == "ucar") {
    if (is.null(initial_values$sig2)) {
      initial_values$sig2 <- matrix(1 / 100, num_group, num_time)
      initmiss <- c(initmiss, "sig2")
    }
  } else if (model %in% c("mcar", "mstcar")) {
    if (is.null(initial_values$G)) {
      initial_values$G <- array(diag(num_group) / 7, dim = c(num_group, num_group, num_time))
      initmiss <- c(initmiss, "G")
    }
  }
  if (model == "mstcar") {
    # rho
    if (is.null(initial_values$rho)) {
      initial_values$rho <- matrix(0.95, 1, num_group)
      initmiss <- c(initmiss, "rho")
    }
    # Ag
    if (is.null(initial_values$Ag)) {
      initial_values$Ag <- diag(1 / 7, num_group)
      initmiss <- c(initmiss, "Ag")
    }
  }
  if (!ignore_checks) {
    check_initial_values(initial_values, data, num_island, model)
  }
  if (!is.null(initmiss)) {
    message("The following objects were created using defaults in 'initial_values': ", paste(initmiss, collapse = " "))
  }
  initial_values
}
