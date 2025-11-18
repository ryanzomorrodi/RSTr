//[[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "cpp_helpers.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::cube update_beta_mcar(List inits, List spatial_data) {
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
arma::cube update_Z_mcar(List inits, List spatial_data) {
  cube Z = inits["Z"];
  cube G = inits["G"];
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
  field<mat> Z_cov(max(num_adj) + 1, num_time);
  field<mat> Z_coveig(max(num_adj) + 1, num_time);
  cube rate_diff = theta - get_regs(beta, island_id);
  vec unique_num_adj = unique(num_adj);
  for (uword time = 0; time < num_time; time++) {
    vec taut = tau2.col(time);
    mat Gt = G.slice(time);
    for (uword count : unique_num_adj) {
      Z_cov(count, time) = inv(diagmat(1 / taut) + count * Gt);
      Z_coveig(count, time) = geteig(Z_cov(count, time));
    }
    for (uword reg = 0; reg < num_region; reg++) {
      vec sum_adj = sum(get_subgrp(Z, adjacency(reg), time), 0).t();
      vec Z_mean = Z_cov(num_adj(reg), time) * (get_grp(rate_diff, reg, time) / taut + Gt * sum_adj);
      vec Z_new  = cpp_rmvnorm(Z_mean, Z_coveig(num_adj(reg), time));
      for (uword grp = 0; grp < num_group; grp++) {
        Z(reg, grp, time) = Z_new(grp);
      }
    }
  }
  cube Zkt(num_island, num_group, num_time);
  for (uword isl = 0; isl < num_island; isl++) {
    Zkt.row(isl) = mean(get_regs(Z, island_region[isl]), 0);
  }
  Z -= get_regs(Zkt, island_id);
  return Z;
}

//[[Rcpp::export]]
arma::mat update_Z_m(List inits, List spatial_data) {
  mat Z = inits["Z"];
  mat G = inits["G"];
  mat theta = inits["theta"];
  mat beta = inits["beta"];
  vec tau2 = inits["tau2"];
  field<uvec> adjacency = spatial_data["adjacency"];
  vec num_adj = spatial_data["num_adj"];
  field<uvec> island_region = spatial_data["island_region"];
  uvec island_id = spatial_data["island_id"];
  uword num_region = Z.n_rows;
  uword num_group  = Z.n_cols;
  uword num_island = island_region.n_elem;
  field<mat> Z_cov(max(num_adj) + 1);
  field<mat> Z_coveig(max(num_adj) + 1);
  vec unique_num_adj = unique(num_adj);
  for (uword count : unique_num_adj) {
    Z_cov   (count) = inv(diagmat(1 / tau2) + count * G);
    Z_coveig(count) = geteig(Z_cov(count));
  }
  mat rate_diff = theta - beta.rows(island_id);
  for (uword reg = 0; reg < num_region; reg++) {
    vec sum_adj  = sum(Z.rows(adjacency(reg)), 0).t();
    vec Z_mean   = Z_cov(num_adj(reg)) * (rate_diff.row(reg).t() / tau2 + G * sum_adj);
    Z.row(reg)   = cpp_rmvnorm(Z_mean, Z_coveig(num_adj(reg))).t();
  }
  mat Zkt(num_island, num_group);
  for (uword isl = 0; isl < num_island; isl++) {
    Zkt.row(isl) = mean(Z.rows(island_region(isl)), 0);
  }
  Z -= Zkt.rows(island_id);
  return Z;
}

//[[Rcpp::export]]
arma::cube update_G_mcar(List inits, List priors, List spatial_data) {
  cube G = inits["G"];
  cube Z = inits["Z"];
  double G_df = priors["G_df"];
  mat G_scale = priors["G_scale"];
  field<uvec> adjacency = spatial_data["adjacency"];
  uword num_island = spatial_data["num_island"];
  uword num_region = Z.n_rows;
  uword num_time = Z.n_slices;
  double df_G = num_region - num_island + G_df;
  for (uword time = 0; time < num_time; time++) {
    mat scale_G = G_scale;
    for (uword reg = 0; reg < num_region; reg++) {
      vec Zit = get_grp(Z, reg, time);
      vec sum_adj = sum(get_subgrp(Z, adjacency(reg), time), 0).t();
      scale_G += adjacency(reg).n_elem * Zit * Zit.t() - sum_adj * Zit.t();
    }
    G.slice(time) = riwish(df_G, scale_G);
  }
  return G;
}

//[[Rcpp::export]]
arma::mat update_tau2_mcar(List inits, List priors, List spatial_data) {
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
arma::cube update_theta_mcar(List inits, List data, List priors, List spatial_data, List params, arma::cube& t_accept) {
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
