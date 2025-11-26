#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "cpp_helpers.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
void update_tau2_default(List& RSTr_obj) {
  List current_sample = RSTr_obj["current_sample"];
  mat tau2 = current_sample["tau2"];
  cube lambda = current_sample["lambda"];
  cube beta = current_sample["beta"];
  cube Z = current_sample["Z"];
  List spatial_data = RSTr_obj["spatial_data"];
  uvec island_id = spatial_data["island_id"];
  List priors = RSTr_obj["priors"];
  double tau_a = priors["tau_a"];
  double tau_b = priors["tau_b"];
  uword num_region = Z.n_rows;
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  double tau_shape = num_region / 2 + tau_a;
  cube square_resid = pow(lambda - get_regs(beta, island_id) - Z, 2);
  for (uword grp = 0; grp < num_group; grp++) {
    for (uword time = 0; time < num_time; time++) {
      double tau_scale = 1 / (sum(get_row(square_resid, grp, time)) / 2 + tau_b);
      tau2(grp, time) = 1 / R::rgamma(tau_shape, tau_scale);
    }
  }
  current_sample["tau2"] = tau2;
  RSTr_obj["current_sample"] = current_sample;
}

//[[Rcpp::export]]
void update_tau2_ucar_restricted(List& RSTr_obj) {
  List current_sample = RSTr_obj["current_sample"];
  mat tau2 = current_sample["tau2"];
  cube lambda = current_sample["lambda"];
  cube beta = current_sample["beta"];
  cube Z = current_sample["Z"];
  mat sig2 = current_sample["sig2"];
  List spatial_data = RSTr_obj["spatial_data"];
  uvec num_island_region = spatial_data["num_island_region"];
  uvec island_id = spatial_data["island_id"];
  List params = RSTr_obj["params"];
  mat A = params["A"];
  double m0 = params["m0"];
  String method = params["method"];
  List priors = RSTr_obj["priors"];
  double tau_a = priors["tau_a"];
  double tau_b = priors["tau_b"];
  uword num_region = Z.n_rows;
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  double tau_shape = num_region / 2 + tau_a;
  cube square_resid = pow(lambda - get_regs(beta, island_id) - Z, 2);
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
  current_sample["tau2"] = tau2;
  RSTr_obj["current_sample"] = current_sample;
}

//[[Rcpp::export]]
void update_tau2_mstcar(List& RSTr_obj) {
  List current_sample = RSTr_obj["current_sample"];
  mat tau2 = current_sample["tau2"];
  cube lambda = current_sample["lambda"];
  cube beta = current_sample["beta"];
  cube Z = current_sample["Z"];
  List priors = RSTr_obj["priors"];
  double tau_a = priors["tau_a"];
  double tau_b = priors["tau_b"];
  List spatial_data = RSTr_obj["spatial_data"];
  uvec island_id = spatial_data["island_id"];
  uword num_region = Z.n_rows;
  uword num_group  = Z.n_cols;
  uword num_time   = Z.n_slices;
  double tau_shape = num_region * num_time / 2 + tau_a;
  cube square_resid = pow(lambda - get_regs(beta, island_id) - Z, 2) / 2;
  for (uword grp = 0; grp < num_group; grp++) {
    double tau_scale = 1 / (accu(square_resid.col(grp)) + tau_b);
    tau2[grp] = 1 / R::rgamma(tau_shape, tau_scale);
  }
  current_sample["tau2"] = tau2;
  RSTr_obj["current_sample"] = current_sample;
}