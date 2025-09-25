#' Initialize CAR model
#'
#' This function performs checks and prepares data for use with either an MSTCAR, MCAR, or UCAR model. This function additionally specifies all of the model parameters, such as model type, event data type, intensity of smoothing in the UCAR model, and more.
#'
#' @param name Name of model and corresponding folder
#' @param dir Directory where model will live
#' @param data Dataset including mortality (Y) and population (n) information
#' @param adjacency Dataset including adjacency information
#' @param inits Optional list of initial conditions for each parameter
#' @param priors Optional list of priors for updates
#' @param model Run model as an MSTCAR/UCAR/MCAR model
#' @param method Run model with either Binomial data or Poisson data
#' @param m0 Baseline neighbor count by region
#' @param A Describes intensity of smoothing between regions
#' @param rho_up Controls whether rho update is performed for MSTCAR models
#' @param impute_lb If counts are suppressed for privacy reasons, \code{impute_lb} is lower bound of suppression, typically 0 or 1
#' @param impute_ub If counts are suppressed for privacy reasons, \code{impute_ub} is upper bound of suppression, typically 10
#' @param seed Set of random seeds to use for data replication
#' @param .show_plots If set to \code{FALSE}, suppresses check plots generated for MSTCAR models
#' @param .ignore_checks If set to \code{TRUE}, ignores data checks. Only use if you are certain that your input data is
#' correct and you are encountering bugs during setup
#' @returns No output, only sets up model and saves files to directory
#' @examples
#' # Initialize an MSTCAR model
#' initialize_model("test", dir = tempdir(), data = miheart, adjacency = miadj)
#' # Initialize an MCAR model
#' data_m <- lapply(miheart, \(x) x[, , "1979"])
#' initialize_model("test", tempdir(), data_m, miadj, model = "mcar")
#' # Initialize an MCAR model with Poisson-distributed event data
#' initialize_model("test", tempdir(), data_m, miadj, model = "mcar", method = "poisson")
#' # Initialize a restricted UCAR model
#' data_u <- lapply(miheart, \(x) x[, "65-74", "1979"])
#' initialize_model("test", tempdir(), data_u, miadj, model = "ucar", A = 6)
#' \dontshow{
#' unlink(paste0(tempdir(), "\\test"), recursive = TRUE)
#' }
#' @export
initialize_model <- function(
    name,
    dir = tempdir(),
    data,
    adjacency,
    inits = NULL,
    priors = NULL,
    model = c("mstcar", "ucar", "mcar"),
    method = c("binomial", "poisson"),
    m0 = NULL,
    A = NULL,
    rho_up = FALSE,
    impute_lb = 1,
    impute_ub = 9,
    seed = 1234,
    .show_plots = TRUE,
    .ignore_checks = FALSE) {
  method <- match.arg(method)
  model <- match.arg(model)
  miss <- which(!is.finite(data$Y))
  if (!.ignore_checks) {
    check_data(data)
  }
  if (model == "mstcar" & .show_plots) {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
    par(mfrow = c(1, 2))
    plot(dimnames(data$Y)[[3]], apply(data$Y, 3, sum, na.rm = TRUE), xlab = "Year", ylab = "Events")
    plot(dimnames(data$Y)[[3]], apply(data$n, 3, sum), xlab = "Year", ylab = "Population")
  }
  if (substr(dir, nchar(dir), nchar(dir)) != "/") {
    dir <- paste0(dir, "/")
  }
  if (!dir.exists(paste0(dir, name))) {
    dir.create(paste0(dir, name))
  }
  pars <- list(
    "ucar" = c("theta", "beta", "Z", "sig2", "tau2"),
    "mcar" = c("theta", "beta", "Z", "G", "tau2"),
    "mstcar" = c("theta", "beta", "Z", "G", "Ag", "tau2")
  )[[model]]
  if (model == "mstcar" & rho_up) pars <- c(pars, "rho")
  for (par in pars) {
    if (!dir.exists(paste0(dir, name, "/", par))) {
      dir.create(paste0(dir, name, "/", par))
    }
  }
  set.seed(seed)
  params <- list(
    seed = .Random.seed,
    batch = 0,
    total = 0,
    model = model,
    method = method
  )
  if (model == "ucar") { # the num_region/num_group/num_time will phase out from this statement
    if (is.null(m0)) m0 <- 3
    if (is.null(A)) A <- 6
    params$m0 <- m0
    params$A <- A
    params$dimnames <- names(data$Y)
  } else if (model == "mcar") {
    params$dimnames <- dimnames(data$Y)
  } else if (model == "mstcar") {
    params$rho_up = rho_up
    params$dimnames = dimnames(data$Y)
  }
  if (length(miss)) {
    params$impute_lb = impute_lb
    params$impute_ub = impute_ub
  }

  spatial_data <- get_spatial_data(adjacency, .ignore_checks)
  if (model == "ucar") {
    priors <- get_priors_u(priors, data, .ignore_checks) # unique to UCAR
    inits <- get_inits_u(inits, data, spatial_data$island_id, method, .ignore_checks) # unique to UCAR
  } else if (model == "mcar") {
    priors <- get_priors_m(priors, data, .ignore_checks) # unique to MCAR
    inits <- get_inits_m(inits, data, spatial_data$island_id, method, .ignore_checks) # unique to MCAR
  } else if (model == "mstcar") {
    priors <- get_priors_mst(priors, data, .ignore_checks) # unique to MSTCAR
    inits <- get_inits_mst(inits, data, spatial_data$island_id, method, .ignore_checks) # unique to MSTCAR
  }

  saveRDS(data, file = paste0(dir, name, "/data.Rds"))
  saveRDS(params, file = paste0(dir, name, "/params.Rds"))
  saveRDS(spatial_data, file = paste0(dir, name, "/spatial_data.Rds"))
  saveRDS(priors, file = paste0(dir, name, "/priors.Rds"))
  saveRDS(inits, file = paste0(dir, name, "/inits.Rds"))
  message("Model ready!")
}
