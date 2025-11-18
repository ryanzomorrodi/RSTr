#' Update Z
#'
#' @noRd
update_Z <- function(inits, spatial_data, params) {
  if (params$model == "ucar") {
    update_Z_ucar(inits, spatial_data)
  } else if (params$model == "mcar") {
    update_Z_mcar(inits, spatial_data)
  } else if (params$model == "mstcar") {
    update_Z_mstcar(inits, spatial_data)
  }
}

#' Update tau2
#'
#' @noRd
update_tau2 <- function(inits, spatial_data, priors, params) {
  if (params$model == "ucar") {
    if (params$restricted) {
      update_tau2_ucar_restricted(inits, spatial_data, priors, params)
    } else {
      update_tau2_ucar(inits, spatial_data, priors)
    }
  } else if (params$model == "mcar") {
    update_tau2_ucar(inits, spatial_data, priors)
  } else if (params$model == "mstcar") {
    update_tau2_mstcar(inits, spatial_data, priors)
  }
}

#' Update beta
#'
#' @noRd
update_beta <- function(inits, spatial_data, params) {
  if (params$model == "ucar") {
    if (params$restricted) {
      update_beta_ucar_restricted(inits, spatial_data, params)
    } else {
      update_beta_ucar(inits, spatial_data)
    }
  } else if (params$model == "mcar") {
    update_beta_ucar(inits, spatial_data)
  } else if (params$model == "mstcar") {
    update_beta_mstcar(inits, spatial_data)
  }
}

#' Update sig2
#'
#' @noRd
update_sig2 <- function(inits, spatial_data, priors, params) {
  if (params$restricted) {
    update_sig2_ucar_restricted(inits, spatial_data, priors, params)
  } else {
    update_sig2_ucar(inits, spatial_data, priors)
  }
}

#' Update G
#'
#' @noRd
update_G <- function(inits, spatial_data, priors, params) {
  if (params$model == "mcar") {
    update_G_mcar(inits, spatial_data, priors)
  } else if (params$model == "mstcar") {
    update_G_mstcar(inits, spatial_data, priors)
  }
}
