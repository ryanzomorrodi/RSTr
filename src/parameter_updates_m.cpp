//[[Rcpp::depends(RcppArmadillo, RcppDist)]]
#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "cpp_helpers.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::mat update_beta_m(List inits, List spatial_data) {
  mat beta = inits["beta"];
  mat theta = inits["theta"];
  mat Z = inits["Z"];
  rowvec tau2 = inits["tau2"];
  field<uvec> island_region = spatial_data["island_region"];
  uword num_group  = Z.n_cols;
  uword num_island = island_region.n_elem;
  for (uword isl = 0; isl < num_island; isl++) {
    rowvec sd_beta   = sqrt(tau2 / island_region[isl].n_elem);
    rowvec mean_beta = mean(theta.rows(island_region[isl]) - Z.rows(island_region[isl]), 0);
    beta.row(isl)    = rowvec(num_group, fill::randn) % sd_beta + mean_beta;
  }
  return beta;
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
    vec sum_adj  = sum(Z.rows(adjacency[reg]), 0).t();
    vec Z_mean   = Z_cov(num_adj[reg]) * (rate_diff.row(reg).t() / tau2 + G * sum_adj);
    Z.row(reg)   = cpp_rmvnorm(Z_mean, Z_coveig(num_adj[reg])).t();
  }
  mat Zkt(num_island, num_group);
  for (uword isl = 0; isl < num_island; isl++) {
    Zkt.row(isl) = mean(Z.rows(island_region[isl]), 0);
  }
  Z -= Zkt.rows(island_id);
  return Z;
}

//[[Rcpp::export]]
arma::mat update_G_m(List inits, List priors, List spatial_data) {
  mat G = inits["G"];
  mat Z = inits["Z"];
  double G_df = priors["G_df"];
  mat G_scale = priors["G_scale"];
  field<uvec> adjacency = spatial_data["adjacency"];
  uword num_island = spatial_data["num_island"];
  uword num_region = Z.n_rows;
  double df_G = num_region - num_island + G_df;
  mat scale_G = G_scale;
  for (uword reg = 0; reg < num_region; reg++) {
    vec sum_adj = sum(Z.rows(adjacency[reg]), 0).t();
    scale_G += adjacency[reg].n_elem * Z.row(reg).t() * Z.row(reg) - sum_adj * Z.row(reg);
  }
  G = riwish(df_G, scale_G);
  return G;
}

//[[Rcpp::export]]
arma::rowvec update_tau2_m(List inits, List priors, List spatial_data) {
  rowvec tau2 = inits["tau2"];
  mat theta = inits["theta"];
  mat beta = inits["beta"];
  mat Z = inits["Z"];
  double tau_a = priors["tau_a"];
  double tau_b = priors["tau_b"];
  uvec island_id = spatial_data["island_id"];
  uword num_region = Z.n_rows;
  uword num_group  = Z.n_cols;
  double tau_shape = num_region / 2 + tau_a;
  mat square_resid = pow(theta - beta.rows(island_id) - Z, 2) / 2;
  rowvec tau_scale = 1 / (sum(square_resid, 0) + tau_b);
  for (uword grp = 0; grp < num_group; grp++) {
    tau2[grp] = 1 / R::rgamma(tau_shape, tau_scale[grp]);
  }
  return tau2;
}

//[[Rcpp::export]]
arma::mat update_theta_m(List inits, List data, List priors, List spatial_data, List params, arma::mat& t_accept) {
  mat theta = inits["theta"];
  mat Z = inits["Z"];
  mat beta = inits["beta"];
  vec tau2 = inits["tau2"];
  mat Y = data["Y"];
  mat n = data["n"];
  mat theta_sd = priors["theta_sd"];
  uvec island_id = spatial_data["island_id"];
  String method = params["method"];
  uword num_region = Z.n_rows;
  uword num_group  = Z.n_cols;
  for (uword reg = 0; reg < num_region; reg++) {
    for (uword grp = 0; grp < num_group; grp++) {
      double theta_star = R::rnorm(theta(reg, grp), theta_sd(reg, grp));
      double rk1 = Y(reg, grp) * (theta_star - theta(reg, grp));
      double rk2 = 0;
      if (method == "binomial") {
        rk2 = n(reg, grp) * (log(1 + exp(theta_star)) - log(1 + exp(theta(reg, grp))));
      } 
      if (method == "poisson") {
        rk2 = n(reg, grp) * (exp(theta_star) - exp(theta(reg, grp)));
      }
      double rk3a = pow(theta_star      - beta(island_id[reg], grp) - Z(reg, grp), 2);
      double rk3b = pow(theta(reg, grp) - beta(island_id[reg], grp) - Z(reg, grp), 2);
      double rk   = exp(rk1 - rk2 - 1 / (2 * tau2[grp]) * (rk3a - rk3b));
      if (rk >= R::runif(0, 1)) {
        t_accept(reg, grp)++;
        theta(reg, grp) = theta_star;
      }
    }
  }
  return theta;
}
