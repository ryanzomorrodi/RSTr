#' Gibbs sampler
#' @useDynLib RSTr, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom RcppDist bayeslm
#' @importFrom RcppArmadillo fastLm
#'
#' @noRd
run_sampler_m <- function(name, dir, iterations, .show_plots, .show_progress, .discard_burnin) {
  sampler_start <- Sys.time()
  data <- readRDS(paste0(dir, name, "/data.Rds"))
  inits <- readRDS(paste0(dir, name, "/inits.Rds"))
  spatial_data <- readRDS(paste0(dir, name, "/spatial_data.Rds"))
  priors <- readRDS(paste0(dir, name, "/priors.Rds"))
  params <- readRDS(paste0(dir, name, "/params.Rds"))

  miss <- which(!is.finite(data$Y))
  total <- params$total
  start_batch <- params$batch
  t_accept <- priors$theta_sd
  t_accept[] <- 0
  batches <- seq(start_batch + 1, start_batch + iterations / 100)
  message("Starting sampler on Batch ", start_batch + 1, " at ", format(Sys.time(), "%a %b %d %X"))
  par_up <- names(inits)
  output_mar = c("theta" = 3, "beta" = 3, "G" = 3, "tau2" = 2, "Z" = 3) # unique to MCAR
  
  if (.show_plots) {
    plots <- vector("list", length(par_up))
    names(plots) <- par_up
    plot_its <- NULL
  }
  for (batch in batches) {
    T_inc <- 100
    if (.show_progress) {
      display_progress(batch, max(batches), total, 0, T_inc, sampler_start)
    }
    output <- vector("list", length(par_up))
    names(output) <- par_up
    t_accept[] <- 0
    for (it in 1:T_inc) {
      if (length(miss)) {
        data$Y <- impute_missing_events(data, inits, params, miss)
      }
      ### Unique to MCAR ###
      inits$beta <- update_beta_m(inits, spatial_data)
      inits$Z <- update_Z_m(inits, spatial_data)
      inits$G <- update_G_m(inits, priors, spatial_data)
      inits$tau2 <- update_tau2_m(inits, priors, spatial_data)
      inits$theta <- update_theta_m(inits, data, priors, spatial_data, params, t_accept)
      ###
      if (it %% 10 == 0) {
        output <- append_to_output(output, inits, output_mar)
        if (.show_plots) {
          plots <- append_to_plots(plots, inits)
        }
      }
      if (.show_progress) {
        display_progress(batch, max(batches), total, it, T_inc, sampler_start)
      }
    }
    priors$theta_sd <- tune_metropolis_sd(priors$theta_sd, t_accept / T_inc)
    total <- total + T_inc
    params$total <- total
    params$batch <- batch
    saveRDS(params, paste0(dir, name, "/params.Rds"))
    saveRDS(priors, paste0(dir, name, "/priors.Rds"))
    saveRDS(inits, paste0(dir, name, "/inits.Rds"))
    save_output(output, batch, dir, name, .discard_burnin)
    if (.show_plots) {
      output_its <- seq((batch - 1) * 100 + 10, batch * 100, 10)
      plot_its <- c(plot_its, output_its)
      grid <- c(2, 3)
      par(mfrow = grid)
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
