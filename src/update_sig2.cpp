#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "cpp_helpers.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
void update_sig2_default(List& RSTr_obj) {
  List current_sample = RSTr_obj["current_sample"];
  mat sig2 = current_sample["sig2"];
  cube Z = current_sample["Z"];
  List spatial_data = RSTr_obj["spatial_data"];
  field<uvec> adjacency = spatial_data["adjacency"];
  vec num_adj = spatial_data["num_adj"];
  field<uvec> island_region = spatial_data["island_region"];
  List priors = RSTr_obj["priors"];
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
  current_sample["sig2"] = sig2;
  RSTr_obj["current_sample"] = current_sample;
}

//[[Rcpp::export]]
void update_sig2_ucar_restricted(List& RSTr_obj) {
  List current_sample = RSTr_obj["current_sample"];
  mat sig2 = current_sample["sig2"];
  cube Z = current_sample["Z"];
  cube beta = current_sample["beta"];
  mat tau2 = current_sample["tau2"];
  List spatial_data = RSTr_obj["spatial_data"];
  field<uvec> adjacency = spatial_data["adjacency"];
  vec num_adj = spatial_data["num_adj"];
  field<uvec> island_region = spatial_data["island_region"];
  uvec num_island_region = spatial_data["num_island_region"];
  List params = RSTr_obj["params"];
  mat A = params["A"];
  double m0 = params["m0"];
  String method = params["method"];
  List priors = RSTr_obj["priors"];
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
  current_sample["sig2"] = sig2;
  RSTr_obj["current_sample"] = current_sample;
}