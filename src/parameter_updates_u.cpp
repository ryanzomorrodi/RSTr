#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::vec update_Z_u(List inits, List spatial_data) {
  vec Z = inits["Z"];
  double sig2 = inits["sig2"];
  vec theta = inits["theta"];
  vec beta = inits["beta"];
  double tau2 = inits["tau2"];
  field<uvec> adjacency = spatial_data["adjacency"];
  vec num_adj = spatial_data["num_adj"];
  field<uvec> island_region = spatial_data["island_region"];
  uvec island_id = spatial_data["island_id"];
  uword num_region = Z.n_elem;
  uword num_island = island_region.n_elem;
  for (uword reg = 0; reg < num_region; reg++) {
    double var_Z  = 1 / (1 / tau2 + num_adj[reg] / sig2);
    double mean_Z = var_Z * ((theta[reg] - beta[island_id[reg]]) / tau2 + sum(Z.elem(adjacency[reg])) / sig2);
    Z[reg] = R::rnorm(mean_Z, sqrt(var_Z));
  }
  vec Zi(num_island);
  for (uword isl = 0; isl < num_island; isl++) {
    Zi(isl) = mean(Z.elem(island_region[isl]));
  }
  Z -= Zi.elem(island_id);
  return Z;
}

//[[Rcpp::export]]
double update_sig2_u(List inits, List spatial_data, List params, List priors) {
  double sig2 = inits["sig2"];
  vec Z = inits["Z"];
  vec beta = inits["beta"];
  double tau2 = inits["tau2"];
  field<uvec> adjacency = spatial_data["adjacency"];
  vec num_adj = spatial_data["num_adj"];
  field<uvec> island_region = spatial_data["island_region"];
  uvec num_island_region = spatial_data["num_island_region"];
  double A = params["A"];
  double m0 = params["m0"];
  String method = params["method"];
  double sig_a = priors["sig_a"];
  double sig_b = priors["sig_b"];
  uword num_region = Z.n_elem;
  uword num_island = island_region.n_elem;
  double sum_adj   = 0;
  for (uword reg = 0; reg < num_region; reg++) {
    sum_adj += Z[reg] * sum(Z.elem(adjacency[reg]));
  }
  double a_sig = (num_region - num_island) / 2 + sig_a;
  double b_sig = 1 / ((sum(pow(Z, 2) % num_adj) - sum_adj) / 2 + sig_b);
  if (A >= 1e2) {
    sig2 = 1 / R::rgamma(a_sig, b_sig);
  } else if (A < 1e2) {
    double sig_thres = 0;
    if (method == "binomial") {
      double pi = sum(beta % num_island_region / num_region);
      pi        = exp(pi) / (1 + exp(pi));
      sig_thres = (1 / ((A + pi) * (1 - pi)) - tau2 * (1 + 1 / m0)) * m0;
    } else if (method == "poisson") {
      sig_thres = (log(1 / A + 1) - tau2 * (1 + 1 / m0)) * m0;
    }
    sig_thres = (sig_thres < 0) ? 0 : sig_thres;
    double u  = R::runif(0, R::pgamma(1 / sig_thres, a_sig, b_sig, true, false));
    sig2      = 1 / R::qgamma(u, a_sig, b_sig, true, false);
  } 
  return sig2;
}

//[[Rcpp::export]]
double update_tau2_u(List inits, List spatial_data, List params, List priors) {
  double tau2 = inits["tau2"];
  vec theta = inits["theta"];
  vec beta = inits["beta"];
  vec Z = inits["Z"];
  double sig2 = inits["sig2"];
  uvec num_island_region = spatial_data["num_island_region"];
  uvec island_id = spatial_data["island_id"];
  double A = params["A"];
  double m0 = params["m0"];
  String method = params["method"];
  double tau_a = priors["tau_a"];
  double tau_b = priors["tau_b"];
  uword num_region = Z.n_elem;
  double a_tau     = num_region / 2 + tau_a;
  double b_tau     = 1 / (sum(pow(theta - beta.elem(island_id) - Z, 2)) / 2 + tau_b);
  if (A >= 1e2) {
    tau2 = 1 / R::rgamma(a_tau, b_tau);
  } else if (A < 1e2) {
    double tau_thres = 0;
    if (method == "binomial") {
      double pi = sum(beta % num_island_region / num_region);
      pi        = exp(pi) / (1 + exp(pi));
      tau_thres = (1 / ((A + pi) * (1 - pi)) - sig2 / m0) / (1 + 1 / m0);
    } else if (method == "poisson") {
      tau_thres = (log(1 / A + 1) - sig2 / m0) / (1 + 1 / m0);
    }
    tau_thres = (tau_thres < 0) ? 0 : tau_thres;
    double u  = R::runif(0, R::pgamma(1 / tau_thres, a_tau, b_tau, true, false));
    tau2      = 1 / R::qgamma(u, a_tau, b_tau, true, false);
  }
  return tau2;
}

//[[Rcpp::export]]
arma::vec update_theta_u(List inits, List data, List priors, List spatial_data, List params, arma::vec& t_accept) {
  vec theta = inits["theta"];
  vec Z = inits["Z"];
  vec beta = inits["beta"];
  double tau2 = inits["tau2"];
  vec Y = data["Y"];
  vec n = data["n"];
  vec theta_sd = priors["theta_sd"];
  uvec island_id = spatial_data["island_id"];
  String method = params["method"];
  uword num_region = Z.n_elem;
  for (uword reg = 0; reg < num_region; reg++) {
    double theta_star = R::rnorm(theta(reg), theta_sd(reg));
    double rk1 = Y(reg) * (theta_star - theta(reg));
    double rk2 = 0;
    if (method == "binomial") {
      rk2 = n(reg) * (log(1 + exp(theta_star)) - log(1 + exp(theta(reg))));
    } else if (method == "poisson") {
      rk2 = n(reg) * (exp(theta_star) - exp(theta(reg)));
    }
    double rk3a = pow(theta_star - beta(island_id[reg]) - Z(reg), 2);
    double rk3b = pow(theta(reg) - beta(island_id[reg]) - Z(reg), 2);
    double rk   = exp(rk1 - rk2 - 1 / (2 * tau2) * (rk3a - rk3b));
    if (rk >= R::runif(0, 1)) {
      t_accept(reg)++;
      theta(reg) = theta_star;
    }
  }
  return theta;
}

//[[Rcpp::export]]
arma::vec update_beta_u(List inits, List spatial_data, List params) {
  vec beta = inits["beta"];
  vec theta = inits["theta"];
  vec Z = inits["Z"];
  double tau2 = inits["tau2"];
  double sig2 = inits["sig2"];
  field<uvec> island_region = spatial_data["island_region"];
  double A = params["A"];
  double m0 = params["m0"];
  String method = params["method"];
  uword num_island = island_region.n_elem;
  for (uword isl = 0; isl < num_island; isl++) {
    uword num_island_region = island_region[isl].n_elem;
    double sd_beta   = sqrt(tau2 / num_island_region);
    double mean_beta = mean(theta.elem(island_region[isl]) - Z.elem(island_region[isl]));
    if (A >= 1e2) {
      beta[isl] = R::rnorm(mean_beta, sd_beta);
    } else if (A < 1e2) {
      if (method == "binomial") {
        double var_t      = tau2 + (tau2 + sig2) / m0;
        double pi_beta    = pow(A - 1, 2) + 4 * (A - 1 / var_t);
        double beta_thres = ((1 - A) + sqrt(pi_beta)) / 2;
        beta_thres        = log(beta_thres / (1 - beta_thres));
        beta_thres        = (beta_thres < 0) ? 0 : beta_thres;
        double beta_max   = R::pnorm(beta_thres, mean_beta, sd_beta, true, false);
        if (beta_max > 0) {
          double u  = R::runif(0, beta_max);
          beta[isl] = R::qnorm(u, mean_beta, sd_beta, true, false);
        }
      } else if (method == "poisson") {
        beta[isl] = R::rnorm(mean_beta, sd_beta);
      }
    }
  }
  return beta;
}