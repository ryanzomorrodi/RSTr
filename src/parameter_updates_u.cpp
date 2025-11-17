#include <RcppArmadillo.h>
#include "cpp_helpers.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::cube update_Z_ucar(List inits, List spatial_data) {
  cube Z = inits["Z"];
  mat sig2 = inits["sig2"];
  cube theta = inits["theta"];
  cube beta = inits["beta"];
  mat tau2 = inits["tau2"];
  field<uvec> adjacency = spatial_data["adjacency"];
  vec num_adj = spatial_data["num_adj"];
  field<uvec> island_region = spatial_data["island_region"];
  uvec island_id = spatial_data["island_id"];
  uword num_region = Z.n_rows;
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  uword num_island = island_region.n_elem;
  for (uword reg = 0; reg < num_region; reg++) {
    for (uword grp = 0; grp < num_group; grp++) {
      for (uword time = 0; time < num_time; time++) {
        double sum_adj = sum(get_subregs(Z, adjacency(reg), grp, time));
        double tmb = (theta(reg, grp, time) - beta(island_id(reg), grp, time));
        double var_Z  = 1 / (1 / tau2(grp, time) + num_adj(reg) / sig2(grp, time));
        double mean_Z = var_Z * (tmb / tau2(grp, time) + sum_adj / sig2(grp, time));
        Z(reg, grp, time) = R::rnorm(mean_Z, sqrt(var_Z));
      }
    }
  }
  cube Zkt(num_island, num_group, num_time);
  for (uword isl = 0; isl < num_island; isl++) {
    Zkt.row(isl) = mean(get_regs(Z, island_region(isl)), 0);
  }
  Z -= get_regs(Zkt, island_id);
  return Z;
}

//[[Rcpp::export]]
arma::mat update_sig2_ucar(List inits, List spatial_data, List priors) {
  mat sig2 = inits["sig2"];
  cube Z = inits["Z"];
  field<uvec> adjacency = spatial_data["adjacency"];
  vec num_adj = spatial_data["num_adj"];
  field<uvec> island_region = spatial_data["island_region"];
  double sig_a = priors["sig_a"];
  double sig_b = priors["sig_b"];
  uword num_region = Z.n_rows;
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  uword num_island = island_region.n_elem;
  for (uword grp = 0; grp < num_group; grp++) {
    for (uword time = 0; time < num_time; time++) {
      double sum_adj = 0;
      for (uword reg = 0; reg < num_region; reg++) {
        sum_adj += Z(reg, grp, time) * sum(get_subregs(Z, adjacency(reg), grp, time));
      }
      double sig_shape = (num_region - num_island) / 2 + sig_a;
      double sig_scale = 1 / ((sum(pow(get_row(Z, grp, time), 2) % num_adj) - sum_adj) / 2 + sig_b);
      sig2(grp, time) = 1 / R::rgamma(sig_shape, sig_scale);
    }
  }
  return sig2;
}

//[[Rcpp::export]]
arma::mat update_sig2_ucar_restricted(List inits, List spatial_data, List params, List priors) {
  mat sig2 = inits["sig2"];
  cube Z = inits["Z"];
  cube beta = inits["beta"];
  mat tau2 = inits["tau2"];
  field<uvec> adjacency = spatial_data["adjacency"];
  vec num_adj = spatial_data["num_adj"];
  field<uvec> island_region = spatial_data["island_region"];
  uvec num_island_region = spatial_data["num_island_region"];
  mat A = params["A"];
  double m0 = params["m0"];
  String method = params["method"];
  double sig_a = priors["sig_a"];
  double sig_b = priors["sig_b"];
  uword num_region = Z.n_rows;
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  uword num_island = island_region.n_elem;
  for (uword grp = 0; grp < num_group; grp++) {
    for (uword time = 0; time < num_time; time++) {
      double sum_adj = 0;
      for (uword reg = 0; reg < num_region; reg++) {
        sum_adj += Z(reg, grp, time) * sum(get_subregs(Z, adjacency(reg), grp, time));
      }
      double sig_shape = (num_region - num_island) / 2 + sig_a;
      double sig_scale = 1 / ((sum(pow(get_row(Z, grp, time), 2) % num_adj) - sum_adj) / 2 + sig_b);
      double sig_thres = 0;
      if (method == "binomial") {
        double pi = sum(get_row(beta, grp, time) % num_island_region / num_region);
        pi = exp(pi) / (1 + exp(pi));
        sig_thres = (1 / ((A(grp, time) + pi) * (1 - pi)) - tau2(grp, time) * (1 + 1 / m0)) * m0;
      } else if (method == "poisson") {
        sig_thres = (log(1 / A(grp, time) + 1) - tau2(grp, time) * (1 + 1 / m0)) * m0;
      }
      sig_thres = (sig_thres < 0) ? 0 : sig_thres;
      double u = R::runif(0, R::pgamma(1 / sig_thres, sig_shape, sig_scale, true, false));
      sig2(grp, time) = 1 / R::qgamma(u, sig_shape, sig_scale, true, false);
    }
  }
  return sig2;
}

//[[Rcpp::export]]
arma::mat update_tau2_ucar(List inits, List spatial_data, List priors) {
  mat tau2 = inits["tau2"];
  cube theta = inits["theta"];
  cube beta = inits["beta"];
  cube Z = inits["Z"];
  uvec island_id = spatial_data["island_id"];
  double tau_a = priors["tau_a"];
  double tau_b = priors["tau_b"];
  uword num_region = Z.n_rows;
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  double tau_shape = num_region / 2 + tau_a;
  cube square_resid = pow(theta - get_regs(beta, island_id) - Z, 2);
  for (uword grp = 0; grp < num_group; grp++) {
    for (uword time = 0; time < num_time; time++) {
      double tau_scale = 1 / (sum(get_row(square_resid, grp, time)) / 2 + tau_b);
      tau2(grp, time) = 1 / R::rgamma(tau_shape, tau_scale);
    }
  }
  return tau2;
}

//[[Rcpp::export]]
arma::mat update_tau2_ucar_restricted(List inits, List spatial_data, List params, List priors) {
  mat tau2 = inits["tau2"];
  cube theta = inits["theta"];
  cube beta = inits["beta"];
  cube Z = inits["Z"];
  mat sig2 = inits["sig2"];
  uvec num_island_region = spatial_data["num_island_region"];
  uvec island_id = spatial_data["island_id"];
  mat A = params["A"];
  double m0 = params["m0"];
  String method = params["method"];
  double tau_a = priors["tau_a"];
  double tau_b = priors["tau_b"];
  uword num_region = Z.n_rows;
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  double tau_shape = num_region / 2 + tau_a;
  cube square_resid = pow(theta - get_regs(beta, island_id) - Z, 2);
  for (uword grp = 0; grp < num_group; grp++) {
    for (uword time = 0; time < num_time; time++) {
      double tau_scale = 1 / (sum(get_row(square_resid, grp, time)) / 2 + tau_b);
      double tau_thres = 0;
      if (method == "binomial") {
        double pi = sum(get_row(beta, grp, time) % num_island_region / num_region);
        pi = exp(pi) / (1 + exp(pi));
        tau_thres = (1 / ((A(grp, time) + pi) * (1 - pi)) - sig2(grp, time) / m0) / (1 + 1 / m0);
      } else if (method == "poisson") {
        tau_thres = (log(1 / A(grp, time) + 1) - sig2(grp, time) / m0) / (1 + 1 / m0);
      }
      tau_thres = (tau_thres < 0) ? 0 : tau_thres;
      double u = R::runif(0, R::pgamma(1 / tau_thres, tau_shape, tau_scale, true, false));
      tau2(grp, time) = 1 / R::qgamma(u, tau_shape, tau_scale, true, false);
    }
  }
  return tau2;
}

//[[Rcpp::export]]
arma::cube update_theta_ucar(List inits, List data, List priors, List spatial_data, List params, arma::cube& t_accept) {
  cube theta = inits["theta"];
  cube Z = inits["Z"];
  cube beta = inits["beta"];
  mat tau2 = inits["tau2"];
  cube Y = data["Y"];
  cube n = data["n"];
  cube theta_sd = priors["theta_sd"];
  uvec island_id = spatial_data["island_id"];
  String method = params["method"];
  uword num_region = Z.n_rows;
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  for (uword time = 0; time < num_time; time++) {
    for (uword reg = 0; reg < num_region; reg++) {
      for (uword grp = 0; grp < num_group; grp++) {
        double theta_star = R::rnorm(theta(reg, grp, time), theta_sd(reg, grp, time));
        double rk1 = Y(reg, grp, time) * (theta_star - theta(reg, grp, time));
        double rk2 = 0;
        if (method == "binomial") {
          rk2 = n(reg, grp, time) * (log(1 + exp(theta_star)) - log(1 + exp(theta(reg, grp, time))));
        } 
        if (method == "poisson") {
          rk2 = n(reg, grp, time) * (exp(theta_star) - exp(theta(reg, grp, time)));
        }
        double rk3a = pow(theta_star            - beta(island_id(reg), grp, time) - Z(reg, grp, time), 2);
        double rk3b = pow(theta(reg, grp, time) - beta(island_id(reg), grp, time) - Z(reg, grp, time), 2);
        double rk   = exp(rk1 - rk2 - 1 / (2 * tau2[grp]) * (rk3a - rk3b));
        if (rk >= R::runif(0, 1)) {
          t_accept(reg, grp, time)++;
          theta(reg, grp, time) = theta_star;
        }
      }
    }
  }
  return theta;
}

//[[Rcpp::export]]
arma::cube update_beta_ucar(List inits, List spatial_data) {
  cube beta = inits["beta"];
  cube theta = inits["theta"];
  cube Z = inits["Z"];
  mat tau2 = inits["tau2"];
  field<uvec> island_region = spatial_data["island_region"];
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  uword num_island = island_region.n_elem;
  cube tmZ = theta - Z;
  for (uword isl = 0; isl < num_island; isl++) {
    uword num_island_region = island_region(isl).n_elem;
    for (uword grp = 0; grp < num_group; grp++) {
      for (uword time = 0; time < num_time; time++) {
        double sd_beta = sqrt(tau2(grp, time) / num_island_region);
        double mean_beta = mean(get_subregs(tmZ, island_region(isl), grp, time));
        beta(isl, grp, time) = R::rnorm(mean_beta, sd_beta);
      }
    }
  }
  return beta;
}

//[[Rcpp::export]]
arma::cube update_beta_ucar_restricted(List inits, List spatial_data, List params) {
  cube beta = inits["beta"];
  cube theta = inits["theta"];
  cube Z = inits["Z"];
  mat tau2 = inits["tau2"];
  mat sig2 = inits["sig2"];
  field<uvec> island_region = spatial_data["island_region"];
  mat A = params["A"];
  double m0 = params["m0"];
  String method = params["method"];
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  uword num_island = island_region.n_elem;
  cube tmZ = theta - Z;
  for (uword isl = 0; isl < num_island; isl++) {
    uword num_island_region = island_region(isl).n_elem;
    mat var_t = tau2 + (tau2 + sig2) / m0;
    for (uword grp = 0; grp < num_group; grp++) {
      for (uword time = 0; time < num_time; time++) {
        double sd_beta = sqrt(tau2(grp, time) / num_island_region);
        double mean_beta = mean(get_subregs(tmZ, island_region(isl), grp, time));
        if (method == "binomial") {
          double pi_beta = pow(A(grp, time) - 1, 2) + 4 * (A(grp, time) - 1 / var_t(grp, time));
          double beta_thres = ((1 - A(grp, time)) + sqrt(pi_beta)) / 2;
          beta_thres = log(beta_thres / (1 - beta_thres));
          beta_thres = (beta_thres < 0) ? 0 : beta_thres;
          double beta_max = R::pnorm(beta_thres, mean_beta, sd_beta, true, false);
          if (beta_max > 0) {
            double u = R::runif(0, beta_max);
            beta(isl, grp, time) = R::qnorm(u, mean_beta, sd_beta, true, false);
          }
        } else if (method == "poisson") {
          beta(isl, grp, time) = R::rnorm(mean_beta, sd_beta);
        }
      }
    }
    
  }
  return beta;
}

