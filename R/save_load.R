#' @noRd
save_output <- function(output, batch, dir, name) {
  UseMethod("save_output")
}

#' @export
save_output.default <- function(output, batch, dir, name) {
  if (substr(dir, nchar(dir), nchar(dir)) != "/") dir <- paste0(dir, "/")
  for (par in names(output)) saveRDS(output[[par]], paste0(dir, name, "/", par, "/", par, "_out_", batch, ".Rds"))
}

#' @noRd
save_model <- function(RSTr_obj) {
  dir <- RSTr_obj$params$dir
  name <- RSTr_obj$params$name
  if (substr(dir, nchar(dir), nchar(dir)) != "/") dir <- paste0(dir, "/")
  saveRDS(RSTr_obj, paste0(dir, name, "/", name, ".Rds"))
}

#' @noRd
load_model <- function(name, dir = tempdir()) {
  if (substr(dir, nchar(dir), nchar(dir)) != "/") dir <- paste0(dir, "/")
  RSTr_obj <- readRDS(paste0(dir, name, "/", name, ".Rds"))
  RSTr_obj$params$name <- name
  RSTr_obj$params$dir <- dir
  RSTr_obj
}

#' @noRd
create_model_directory <- function(name, dir, pars) {
  if (substr(dir, nchar(dir), nchar(dir)) != "/") dir <- paste0(dir, "/")
  if (!dir.exists(paste0(dir, name))) dir.create(paste0(dir, name))
  for (par in pars) dir.create(paste0(dir, name, "/", par), showWarnings = FALSE)
}