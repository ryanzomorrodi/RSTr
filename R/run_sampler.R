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
  inits <- readRDS(paste0(dir, name, "/inits.Rds"))
  spatial_data <- readRDS(paste0(dir, name, "/spatial_data.Rds"))
  priors <- readRDS(paste0(dir, name, "/priors.Rds"))
  params <- readRDS(paste0(dir, name, "/params.Rds"))

  par_up <- names(inits)
  model <- params$model
  miss <- which(!is.finite(data$Y))
  total <- params$total
  start_batch <- params$batch
  t_accept <- priors$theta_sd
  t_accept[] <- 0
  batches <- seq(start_batch + 1, start_batch + iterations / 100)

  if (model == "mstcar") {
    rho_up <- params$rho_up
    r_accept <- priors$rho_sd
    if (rho_up) r_accept[] <- 0
    if (!rho_up) par_up <- par_up[-which(par_up == "rho")]
  }
  output_mar <- list(
    "ucar" = c("theta" = 4, "beta" = 4, "sig2" = 3, "tau2" = 3, "Z" = 4),
    "mcar" = c("theta" = 4, "beta" = 4, "G" = 4, "tau2" = 3, "Z" = 4),
    "mstcar" = c("theta" = 4, "beta" = 4, "G" = 4, "tau2" = 2, "Ag" = 3, "Z" = 4, "rho" = 2)
  )[[model]]

  if (show_plots) {
    plots <- vector("list", length(par_up))
    names(plots) <- par_up
    plot_its <- NULL
  }

  message("Starting sampler on Batch ", start_batch + 1, " at ", format(Sys.time(), "%a %b %d %X"))
  for (batch in batches) {
    T_inc <- 100
    if (show_progress) {
      display_progress(batch, max(batches), total, 0, T_inc, sampler_start)
    }
    output <- vector("list", length(par_up))
    names(output) <- par_up
    if (model == "mstcar") if (rho_up) r_accept[] <- 0
    t_accept[] <- 0

    for (it in 1:T_inc) {
      if (length(miss)) {
        data$Y <- impute_missing_events(data, inits, params, miss)
      }
      if (model == "ucar") {
        inits$Z <- update_Z_ucar(inits, spatial_data)
        inits$theta <- update_theta_ucar(inits, data, priors, spatial_data, params, t_accept)
        if (params$restricted) {
          inits$sig2 <- update_sig2_ucar_restricted(inits, spatial_data, params, priors)
          inits$tau2 <- update_tau2_ucar_restricted(inits, spatial_data, params, priors)
          inits$beta <- update_beta_ucar_restricted(inits, spatial_data, params)
        } else {
          inits$sig2 <- update_sig2_ucar(inits, spatial_data, priors)
          inits$tau2 <- update_tau2_ucar(inits, spatial_data, priors)
          inits$beta <- update_beta_ucar(inits, spatial_data)
        }

      } else if (model == "mcar") {
        inits$G <- update_G_mcar(inits, priors, spatial_data)
        inits$Z <- update_Z_mcar(inits, spatial_data)
        inits$beta <- update_beta_mcar(inits, spatial_data)
        inits$tau2 <- update_tau2_mcar(inits, priors, spatial_data)
        inits$theta <- update_theta_mcar(inits, data, priors, spatial_data, params, t_accept)
        
      } else if (model == "mstcar") {
        inits$beta <- update_beta_mst(inits, spatial_data)
        inits$Z <- update_Z_mst(inits, spatial_data)
        inits$G <- update_G_mst(inits, priors, spatial_data)
        inits$Ag <- update_Ag_mst(inits, priors)
        inits$tau2 <- update_tau2_mst(inits, priors, spatial_data)
        inits$theta <- update_theta_mst(inits, data, priors, spatial_data, params, t_accept)
        if (rho_up) inits$rho <- update_rho_mst(inits, priors, spatial_data, r_accept)
      }
      if (it %% 10 == 0) {
        output <- append_to_output(output, inits, output_mar)
        if (show_plots) {
          plots <- append_to_plots(plots, inits)
        }
      }
      if (show_progress) {
        display_progress(batch, max(batches), total, it, T_inc, sampler_start)
      }
    }
    if (model == "mstcar") if (rho_up) priors$rho_sd <- tune_metropolis_sd(priors$rho_sd, r_accept / T_inc)
    priors$theta_sd <- tune_metropolis_sd(priors$theta_sd, t_accept / T_inc)
    total <- total + T_inc
    params$total <- total
    params$batch <- batch
    saveRDS(params, paste0(dir, name, "/params.Rds"))
    saveRDS(priors, paste0(dir, name, "/priors.Rds"))
    saveRDS(inits, paste0(dir, name, "/inits.Rds"))
    save_output(output, batch, dir, name, discard_burnin)

    if (show_plots) {
      output_its <- seq((batch - 1) * 100 + 10, batch * 100, 10)
      plot_its <- c(plot_its, output_its)
      grid <- c(2, 3)
      if (model == "mstcar") if (rho_up) grid <- c(2, 4)
      graphics::par(mfrow = grid)
      # Gradually remove plots in burn-in, then plot
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
impute_missing_events <- function(data, inits, params, miss) {
  if (params$method == "binomial") {
    rate <- expit(inits$theta[miss])
    rp <- stats::runif(
      length(miss),
      stats::pbinom(params$impute_lb - 0.1, round(data$n[miss]), rate),
      stats::pbinom(params$impute_ub + 0.1, round(data$n[miss]), rate)
    )
    data$Y[miss] <- stats::qbinom(rp, round(data$n[miss]), rate)
  }
  if (params$method == "poisson") {
    rate <- data$n[miss] * exp(inits$theta[miss])
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
append_to_output <- function(output, inits, output_mar) {
  onames = names(output)
  output <- lapply(names(output), \(par) abind::abind(output[[par]], inits[[par]], along = output_mar[par]))
  names(output) <- onames
  output
}

#' Append new values to plots
#' @noRd
append_to_plots <- function(plots, inits) {
  pnames = names(plots)
  plots <- lapply(names(plots), \(par) c(plots[[par]], inits[[par]][1]))
  names(plots) <- pnames
  plots
}
