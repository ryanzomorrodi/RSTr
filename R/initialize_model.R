initialize_model <- function(
  name,
  data,
  adjacency,
  dir = tempdir(),
  show_plots = TRUE,
  ignore_checks = FALSE,
  method = "binomial",
  impute_lb = 0,
  impute_ub = 10,
  seed = 1234,
  initial_values = NULL,
  priors = NULL,
  model = c("mstcar", "ucar", "mcar"),
  pars,
  restricted = NULL,
  A = NULL,
  m0 = NULL,
  update_rho = NULL
) {
  RSTr_obj <- create_new_model(model, data, restricted, update_rho)
  if (show_plots & (dim(RSTr_obj$data$Y)[3] > 1)) make_data_plots(RSTr_obj$data)
  RSTr_obj$params <- get_params(RSTr_obj$data, seed, method, model, name, dir, restricted, A, m0, update_rho, impute_lb, impute_ub)
  RSTr_obj$spatial_data <- get_spatial_data(adjacency)
  RSTr_obj <- get_priors(RSTr_obj, priors)
  RSTr_obj$initial_values <- get_initial_values(RSTr_obj, initial_values, method)
  RSTr_obj$current_sample <- RSTr_obj$initial_values
  if (!ignore_checks) validate_model(RSTr_obj)
  create_model_directory(name, dir, pars)
  save_model(RSTr_obj)
  RSTr_obj
}

#' Initialize CAR model
#' This function performs checks and prepares data for use with either an MSTCAR, MCAR, or UCAR model. This function additionally specifies all of the model parameters, such as model type, event data type, intensity of smoothing in the UCAR model, and more.
#' @param name Name of model and corresponding folder
#' @param data Dataset including mortality (Y) and population (n) information
#' @param adjacency Dataset including adjacency information
#' @param dir Directory where model will live
#' @param show_plots If set to \code{FALSE}, suppresses check plots generated for MSTCAR models
#' @param ignore_checks If set to \code{TRUE}, ignores data checks. Only use if you are certain that your input data is
#' correct and you are encountering bugs during setup
#' @param method Run model with either Binomial data or Poisson data
#' @param impute_lb If counts are suppressed for privacy reasons, \code{impute_lb} is lower bound of suppression, typically 0 or 1
#' @param impute_ub If counts are suppressed for privacy reasons, \code{impute_ub} is upper bound of suppression, typically 10
#' @param seed Set of random seeds to use for data replication
#' @param initial_values Optional list of initial conditions for each parameter
#' @param priors Optional list of priors for updates
#' @param m0 For restricted UCAR models, baseline neighbor count by region
#' @param A For restricted UCAR models, describes intensity of smoothing between regions
#' @param update_rho For MSTCAR models, controls whether rho update is performed for MSTCAR models
#' @returns No output, only sets up model and saves files to directory
#' @examples
#' # Initialize an MSTCAR model
#' mstcar(name = "test", data = miheart, adjacency = miadj, dir = tempdir())
#' # Initialize an MCAR model
#' data_m <- lapply(miheart, \(x) x[, , "1979"])
#' mcar("test", data_m, miadj, tempdir())
#' # Initialize an MCAR model with Poisson-distributed event data
#' mcar("test", data_m, miadj, tempdir(), method = "poisson")
#' # Initialize a restricted UCAR model
#' data_u <- lapply(miheart, \(x) x[, "65-74", "1979"])
#' rucar("test", data_u, miadj, tempdir(), A = 6)
#' \dontshow{
#' unlink(paste0(tempdir(), "\\test"), recursive = TRUE)
#' }
#' @export
ucar <- function(
  name,
  data,
  adjacency,
  dir = tempdir(),
  iterations = 6000,
  show_plots = TRUE,
  verbose = TRUE,
  ignore_checks = FALSE,
  method = c("binomial", "poisson"),
  impute_lb = NULL,
  impute_ub = NULL,
  seed = NULL,
  initial_values = NULL,
  priors = NULL
) {
  method <- match.arg(method)
  pars <- c("lambda", "beta", "Z", "sig2", "tau2")
  RSTr_obj <- initialize_model(
    name = name,
    data = data,
    adjacency = adjacency,
    dir = dir,
    show_plots = show_plots,
    ignore_checks = ignore_checks,
    method = method,
    impute_lb = impute_lb,
    impute_ub = impute_ub,
    seed = seed,
    initial_values = initial_values,
    priors = priors,
    model = "ucar",
    pars = pars,
    restricted = FALSE
  )
  RSTr_obj <- update_model(RSTr_obj, iterations, show_plots, verbose)
  RSTr_obj
}

#' Initialize Restricted UCAR model
#' @rdname ucar
#' @export
rucar <- function(
  name,
  data,
  adjacency,
  dir = tempdir(),
  A = NULL,
  m0 = NULL,
  iterations = 6000,
  show_plots = TRUE,
  verbose = TRUE,
  ignore_checks = FALSE,
  method = c("binomial", "poisson"),
  impute_lb = NULL,
  impute_ub = NULL,
  seed = NULL,
  initial_values = NULL,
  priors = NULL
) {
  method <- match.arg(method)
  pars <- c("lambda", "beta", "Z", "sig2", "tau2")
  RSTr_obj <- initialize_model(
    name = name,
    data = data,
    adjacency = adjacency,
    dir = dir,
    show_plots = show_plots,
    ignore_checks = ignore_checks,
    method = method,
    impute_lb = impute_lb,
    impute_ub = impute_ub,
    seed = seed,
    initial_values = initial_values,
    priors = priors,
    model = "ucar",
    pars = pars,
    restricted = TRUE,
    A = A,
    m0 = m0
  )
  RSTr_obj <- update_model(RSTr_obj, iterations, show_plots, verbose)
  RSTr_obj
}



#' Initialize MCAR model
#' @rdname ucar
#' @export
mcar <- function(
  name,
  data,
  adjacency,
  dir = tempdir(),
  iterations = 6000,
  show_plots = TRUE,
  verbose = TRUE,
  ignore_checks = FALSE,
  method = c("binomial", "poisson"),
  impute_lb = NULL,
  impute_ub = NULL,
  seed = NULL,
  initial_values = NULL,
  priors = NULL
) {
  method <- match.arg(method)
  pars <- c("lambda", "beta", "Z", "G", "tau2")
  RSTr_obj <- initialize_model(
    name = name,
    data = data,
    adjacency = adjacency,
    dir = dir,
    show_plots = show_plots,
    ignore_checks = ignore_checks,
    method = method,
    impute_lb = impute_lb,
    impute_ub = impute_ub,
    seed = seed,
    initial_values = initial_values,
    priors = priors,
    model = "mcar",
    pars = pars
  )
  RSTr_obj <- update_model(RSTr_obj, iterations, show_plots, verbose)
  RSTr_obj
}

#' Initialize MSTCAR model
#' @rdname ucar
#' @export
mstcar <- function(
  name,
  data,
  adjacency,
  dir = tempdir(),
  iterations = 6000,
  show_plots = TRUE,
  verbose = TRUE,
  ignore_checks = FALSE,
  method = c("binomial", "poisson"),
  impute_lb = NULL,
  impute_ub = NULL,
  seed = NULL,
  initial_values = NULL,
  priors = NULL,
  update_rho = FALSE
) {
  method <- match.arg(method)
  pars <- c("lambda", "beta", "Z", "G", "Ag", "tau2")
  if (update_rho) pars <- c(pars, "rho")
  RSTr_obj <- initialize_model(
    name = name,
    data = data,
    adjacency = adjacency,
    dir = dir,
    show_plots = show_plots,
    ignore_checks = ignore_checks,
    method = method,
    impute_lb = impute_lb,
    impute_ub = impute_ub,
    seed = seed,
    initial_values = initial_values,
    priors = priors,
    model = "mstcar",
    pars = pars,
    update_rho = update_rho
  )
  RSTr_obj <- update_model(RSTr_obj, iterations, show_plots, verbose)
  RSTr_obj
}
