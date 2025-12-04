#' @noRd
create_new_model <- function(model, data, restricted = NULL, update_rho = NULL) {
  data <- prepare_data(data)
  if (model == "ucar") {
    if (restricted) new_ucar_restricted(data) else new_ucar(data)
  } else if (model == "mcar") {
    new_mcar(data)
  } else if (model == "mstcar") {
    if (update_rho) new_mstcar_update_rho(data) else new_mstcar(data)
  }
}

#' @noRd
new_RSTr <- function(data, subclass = character()) {
  structure(list(data = data), class = c(subclass, "RSTr"))
}

#' @noRd
new_ucar <- function(data, subclass = character()) {
  new_RSTr(data, subclass = c(subclass, "ucar"))
}

#' @noRd
new_ucar_restricted <- function(data) {
  new_ucar(data, subclass = "eucar")
}

#' @noRd
new_mcar <- function(data, subclass = character()) {
  new_RSTr(data, subclass = c(subclass, "mcar"))
}

#' @noRd
new_mstcar <- function(data, subclass = character()) {
  new_RSTr(data, subclass = c(subclass, "mstcar"))
}

#' @noRd
new_mstcar_update_rho <- function(data) {
  new_mstcar(data, subclass = "mstcar_update_rho")
}

