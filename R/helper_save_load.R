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
  new_dir <- paste0(dir, "/", name)
  if (!dir.exists(new_dir)) dir.create(new_dir)
  param_dirs <- paste0(dir, "/", name, "/", pars)
  for (par in param_dirs) if(!dir.exists(par)) dir.create(par)
}