#' Run Gibbs sampler
#'
#' \code{run_sampler()} generates samples for model \code{name} in \code{dir}. The model used to generate samples (e.g., MSTCAR, MCAR, UCAR) along with the model's other parameters are specified in \code{initialize_*()}.
#' @param name Name of model and corresponding folder
#' @param dir Directory where model lives
#' @param iterations Specifies number of iterations to run
#' @param show_plots If set to \code{FALSE}, hides traceplots
#' @param show_progress If set to \code{FALSE}, hides progress bar
#' @param discard_burnin If set to \code{TRUE}, won't save burn-in samples
#' @returns No output, saves sampler output to \code{dir}
#' @examples
#' data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
#' adj_min <- list(2, 1)
#' initialize_mstcar("test", data_min, adj_min, tempdir())
#' run_sampler("test", show_plots = FALSE, show_progress = FALSE)
#' \dontshow{
#' unlink(paste0(tempdir(), "\\test"), recursive = TRUE)
#' }
#' @useDynLib RSTr, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom RcppDist bayeslm
#' @importFrom RcppArmadillo fastLm
#' @export
run_sampler <- function(name, dir = tempdir(), iterations = 6000, show_plots = TRUE, show_progress = TRUE, discard_burnin = FALSE) {
  if (show_plots) {
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
  }
  if (substr(dir, nchar(dir), nchar(dir)) != "/") {
    dir <- paste0(dir, "/")
  }
  iterations <- iterations - iterations %% 100

  sampler_start <- Sys.time()
  data <- readRDS(paste0(dir, name, "/data.Rds"))
  current_sample <- readRDS(paste0(dir, name, "/current_sample.Rds"))
  spatial_data <- readRDS(paste0(dir, name, "/spatial_data.Rds"))
  priors <- readRDS(paste0(dir, name, "/priors.Rds"))
  params <- readRDS(paste0(dir, name, "/params.Rds"))

  par_up <- names(current_sample)
  model <- params$model
  miss <- which(!is.finite(data$Y))
  total <- params$total
  start_batch <- params$batch
  t_accept <- priors$lambda_sd
  t_accept[] <- 0
  batches <- seq(start_batch + 1, start_batch + iterations / 100)

  if (model == "mstcar") {
    rho_up <- params$rho_up
    r_accept <- priors$rho_sd
    if (rho_up) r_accept[] <- 0
    if (!rho_up) par_up <- par_up[-which(par_up == "rho")]
  }
  if (show_plots) {
    plots <- vector("list", length(par_up))
    names(plots) <- par_up
    plot_its <- NULL
  }
  message("Starting sampler on Batch ", start_batch + 1, " at ", format(Sys.time(), "%a %b %d %X"))
  for (batch in batches) {
    T_inc <- 100
    if (show_progress) display_progress(batch, max(batches), total, 0, T_inc, sampler_start)
    output <- vector("list", length(par_up))
    names(output) <- par_up
    t_accept[] <- 0
    current_sample$lambda <- log_logit(current_sample$lambda, params$method)
    for (it in 1:T_inc) {
      if (length(miss)) data$Y <- impute_missing_events(data, current_sample, params, miss)
      current_sample$lambda <- update_lambda(current_sample, spatial_data, priors, params, data, t_accept)
      current_sample$Z <- update_Z(current_sample, spatial_data, params)
      current_sample$tau2 <- update_tau2(current_sample, spatial_data, priors, params)
      current_sample$beta <- update_beta(current_sample, spatial_data, params)
      if (model == "ucar") current_sample$sig2 <- update_sig2(current_sample, spatial_data, priors, params)
      if (model %in% c("mcar", "mstcar")) current_sample$G <- update_G(current_sample, spatial_data, priors, params)
      if (model == "mstcar") {
        current_sample$Ag <- update_Ag(current_sample, priors)
        if (rho_up) current_sample$rho <- update_rho(current_sample, spatial_data, priors, r_accept)
      }
      if (it %% 10 == 0) {
        output <- append_to_output(output, current_sample)
        if (show_plots) plots <- append_to_plots(plots, current_sample, params$method)
      }
      if (show_progress) display_progress(batch, max(batches), total, it, T_inc, sampler_start)
    }
    current_sample$lambda <- exp_expit(current_sample$lambda, params$method)
    output$lambda <- exp_expit(output$lambda, params$method)
    if (model == "mstcar") if (rho_up) priors$rho_sd <- tune_metropolis_sd(priors$rho_sd, r_accept / T_inc)
    priors$lambda_sd <- tune_metropolis_sd(priors$lambda_sd, t_accept / T_inc)
    total <- total + T_inc
    params$total <- total
    params$batch <- batch
    saveRDS(params, paste0(dir, name, "/params.Rds"))
    saveRDS(priors, paste0(dir, name, "/priors.Rds"))
    saveRDS(current_sample, paste0(dir, name, "/current_sample.Rds"))
    save_output(output, batch, dir, name, discard_burnin)

    if (show_plots) {
      output_its <- seq((batch - 1) * 100 + 10, batch * 100, 10)
      plot_its <- c(plot_its, output_its)
      grid <- c(2, 3)
      if (model == "mstcar") if (rho_up) grid <- c(2, 4)
      graphics::par(mfrow = grid)
      # Gradually remove values in burn-in, then plot
      if (plot_its[1] < 2000) {
        plots <- lapply(plots, \(par) par[-(1:5)])
        plot_its <- plot_its[-(1:5)]
      }
      lapply(
        par_up,
        \(par) plot(plot_its, plots[[par]], type = "l", main = par, xlab = "Iteration", ylab = "Value")
      )
    }
  }
  message("Model finished at ", format(Sys.time(), "%a %b %d %X"))
}

#' Impute event values
#' @noRd
impute_missing_events <- function(data, current_sample, params, miss) {
  if (params$method == "binomial") {
    rate <- expit(current_sample$lambda[miss])
    rp <- stats::runif(
      length(miss),
      stats::pbinom(params$impute_lb - 0.1, round(data$n[miss]), rate),
      stats::pbinom(params$impute_ub + 0.1, round(data$n[miss]), rate)
    )
    data$Y[miss] <- stats::qbinom(rp, round(data$n[miss]), rate)
  }
  if (params$method == "poisson") {
    rate <- data$n[miss] * exp(current_sample$lambda[miss])
    rp <- stats::runif(
      length(miss),
      stats::ppois(params$impute_lb - 0.1, rate),
      stats::ppois(params$impute_ub + 0.1, rate)
    )
    data$Y[miss] <- stats::qpois(rp, rate)
  }
  data$Y
}

#' Tune metropolis standard deviation
#' @noRd
tune_metropolis_sd <- function(sd, accept) {
  accept <- pmin(pmax(accept, 1 / 6), 0.75)
  sd <- ifelse(accept > 0.5, sd * accept / 0.5, sd)
  sd <- ifelse(accept < 0.35, sd * accept / 0.35, sd)
  sd
}

#' Append new values to output
#' @noRd
append_to_output <- function(output, current_sample) {
  onames <- names(output)
  output_mar <- sapply(current_sample, \(par) length(dim(par)) + 1) 
  output <- lapply(
    names(output),
    \(par) abind::abind(output[[par]], current_sample[[par]], along = output_mar[par])
  )
  names(output) <- onames
  output
}

#' Append new values to plots
#' @noRd
append_to_plots <- function(plots, current_sample, method) {
  current_sample$lambda <- exp_expit(current_sample$lambda, method)
  pnames <- names(plots)
  plots <- lapply(names(plots), \(par) c(plots[[par]], current_sample[[par]][1]))
  names(plots) <- pnames
  plots
}
