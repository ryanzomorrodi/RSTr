#' Update Z
#'
#' @noRd
update_Z <- function(initial_values, spatial_data, params) {
  if (params$model == "ucar") {
    update_Z_ucar(initial_values, spatial_data)
  } else if (params$model == "mcar") {
    update_Z_mcar(initial_values, spatial_data)
  } else if (params$model == "mstcar") {
    update_Z_mstcar(initial_values, spatial_data)
  }
}

#' Update tau2
#'
#' @noRd
update_tau2 <- function(initial_values, spatial_data, priors, params) {
  if (params$model == "ucar") {
    if (params$restricted) {
      update_tau2_ucar_restricted(initial_values, spatial_data, priors, params)
    } else {
      update_tau2_ucar(initial_values, spatial_data, priors)
    }
  } else if (params$model == "mcar") {
    update_tau2_ucar(initial_values, spatial_data, priors)
  } else if (params$model == "mstcar") {
    update_tau2_mstcar(initial_values, spatial_data, priors)
  }
}

#' Update beta
#'
#' @noRd
update_beta <- function(initial_values, spatial_data, params) {
  if (params$model == "ucar") {
    if (params$restricted) {
      update_beta_ucar_restricted(initial_values, spatial_data, params)
    } else {
      update_beta_ucar(initial_values, spatial_data)
    }
  } else if (params$model == "mcar") {
    update_beta_ucar(initial_values, spatial_data)
  } else if (params$model == "mstcar") {
    update_beta_mstcar(initial_values, spatial_data)
  }
}

#' Update sig2
#'
#' @noRd
update_sig2 <- function(initial_values, spatial_data, priors, params) {
  if (params$restricted) {
    update_sig2_ucar_restricted(initial_values, spatial_data, priors, params)
  } else {
    update_sig2_ucar(initial_values, spatial_data, priors)
  }
}

#' Update G
#'
#' @noRd
update_G <- function(initial_values, spatial_data, priors, params) {
  if (params$model == "mcar") {
    update_G_mcar(initial_values, spatial_data, priors)
  } else if (params$model == "mstcar") {
    update_G_mstcar(initial_values, spatial_data, priors)
  }
}
