#' Get priors
#' @noRd
get_priors <- function(priors, data, model, .ignore_checks) {
  # Prepare Hyperparameters
  primiss <- NULL
  if (is.null(priors$theta_sd)) {
    priors$theta_sd <- data$Y
    priors$theta_sd[] <- 0.025
    primiss <- c(primiss, "theta_sd")
  }
  if (is.null(priors$tau_a)) {
    priors$tau_a <- 0.001
    primiss <- c(primiss, "tau_a")
  }
  if (is.null(priors$tau_b)) {
    priors$tau_b <- 0.001
    primiss <- c(primiss, "tau_b")
  }
  if (model == "ucar") {
    num_region <- length(data$Y) # this will be phased out

    if (is.null(priors$sig_a)) {
      priors$sig_a <- 0.001
      primiss <- c(primiss, "sig_a")
    }
    if (is.null(priors$sig_b)) {
      priors$sig_b <- 0.001
      primiss <- c(primiss, "sig_b")
    }
    if (!.ignore_checks) {
      check_priors_u(priors, num_region)
    }
  } else if (model == "mcar") {
    num_region <- dim(data$Y)[1]
    num_group <- dim(data$Y)[2]

    if (is.null(priors$G_scale)) {
      priors$G_scale <- diag(1 / 7, num_group)
      primiss <- c(primiss, "G_scale")
    }
    if (is.null(priors$G_df)) {
      priors$G_df <- num_group + 2
      primiss <- c(primiss, "G_df")
    }
    if (!.ignore_checks) {
      check_priors_m(priors, num_region, num_group)
    }
  } else if (model == "mstcar") {
    num_region <- dim(data$Y)[1]
    num_group <- dim(data$Y)[2]
    num_time <- dim(data$Y)[3]
    
    if (is.null(priors$Ag_scale)) {
      priors$Ag_scale <- diag(1 / 7, num_group)
      primiss <- c(primiss, "Ag_scale")
    }
    if (is.null(priors$Ag_df)) {
      priors$Ag_df <- num_group + 2
      primiss <- c(primiss, "Ag_df")
    }
    if (is.null(priors$G_df)) {
      priors$G_df <- num_group + 2
      primiss <- c(primiss, "G_df")
    }
    if (is.null(priors$rho_a)) {
      priors$rho_a <- 95
      primiss <- c(primiss, "rho_a")
    }
    if (is.null(priors$rho_b)) {
      priors$rho_b <- 5
      primiss <- c(primiss, "rho_b")
    }
    if (is.null(priors$rho_sd)) {
      priors$rho_sd <- rep(0.05, num_group)
      primiss <- c(primiss, "rho_sd")
    }
    if (!.ignore_checks) {
      check_priors_mst(priors, num_region, num_group, num_time)
    }
  }
  if (!is.null(primiss)) {
    message("The following objects were created using defaults in 'priors': ", paste(primiss, collapse = " "))
  }
  priors
}
