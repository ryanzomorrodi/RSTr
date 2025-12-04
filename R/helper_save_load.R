#' @noRd
save_output <- function(output, batch, dir, name) {
  for (par in names(output)) saveRDS(output[[par]], paste0(dir, "/", name, "/", par, "/", par, "_out_", batch, ".Rds"))
}

#' @noRd
save_model <- function(RSTr_obj) {
  dir <- RSTr_obj$params$dir
  name <- RSTr_obj$params$name
  saveRDS(RSTr_obj, paste0(dir, "/", name, "/", name, ".Rds"))
}

#' @noRd
create_model_directory <- function(name, dir, pars) {
  if (!dir.exists(paste0(dir, name))) dir.create(paste0(dir, "/", name))
  for (par in pars) dir.create(paste0(dir, "/", name, "/", par), showWarnings = FALSE)
}