#' Gibbs sampler
#' @useDynLib RSTr, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom RcppDist bayeslm
#' @importFrom RcppArmadillo fastLm
#'
#' @noRd
run_sampler_mst <- function(name, dir, iterations, .show_plots, .show_progress, .discard_burnin) {
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

  rho_up <- params$rho_up # unique to MSTCAR
  r_accept <- priors$rho_sd # unique to MSTCAR
  if (rho_up) r_accept[] <- 0 # unique to MSTCAR

  

  message("Starting sampler on Batch ", start_batch + 1, " at ", format(Sys.time(), "%a %b %d %X"))
  par_up <- names(inits)

  if (!rho_up) par_up <- par_up[-which(par_up == "rho")] # unique to MSTCAR
  output_mar = c("theta" = 4, "beta" = 4, "G" = 4, "tau2" = 2, "Ag" = 3, "Z" = 4, "rho" = 2) # unique to MSTCAR
  
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
    if (rho_up) r_accept[] <- 0 # unique to MSTCAR
    t_accept[] <- 0
    
    for (it in 1:T_inc) {
      if (length(miss)) {
        data$Y <- impute_missing_events(data, inits, params, miss)
      }
      ### unique to MSTCAR ###
      inits$beta <- update_beta_mst(inits, spatial_data)
      inits$Z <- update_Z_mst(inits, spatial_data)
      inits$G <- update_G_mst(inits, priors, spatial_data)
      inits$Ag <- update_Ag_mst(inits, priors)
      inits$tau2 <- update_tau2_mst(inits, priors, spatial_data)
      inits$theta <- update_theta_mst(inits, data, priors, spatial_data, params, t_accept)
      if (rho_up) inits$rho <- update_rho_mst(inits, priors, spatial_data, r_accept)
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
    if (rho_up) priors$rho_sd <- tune_metropolis_sd(priors$rho_sd, r_accept / T_inc) # unique to MSTCAR
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
      grid <- c(2, ifelse(rho_up, 4, 3)) # unique to MSTCAR
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
