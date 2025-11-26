#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "cpp_helpers.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
void update_beta_default(List& RSTr_obj) {
  List current_sample = RSTr_obj["current_sample"];
  cube beta = current_sample["beta"];
  cube lambda = current_sample["lambda"];
  cube Z = current_sample["Z"];
  mat tau2 = current_sample["tau2"];
  List spatial_data = RSTr_obj["spatial_data"];
  field<uvec> island_region = spatial_data["island_region"];
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  uword num_island = island_region.n_elem;
  for (uword isl = 0; isl < num_island; isl++) {
    uword num_island_region = island_region[isl].n_elem;
    mat var_beta = sqrt(tau2 / num_island_region);
    mat mean_beta = mean(get_regs(lambda, island_region[isl]) - get_regs(Z, island_region[isl]), 0);
    if (num_time == 1) mean_beta = mean_beta.t();
    beta.row(isl) = mat(num_group, num_time, fill::randn) % var_beta + mean_beta;
  }
  current_sample["beta"] = beta;
  RSTr_obj["current_sample"] = current_sample;
}

//[[Rcpp::export]]
void update_beta_ucar_restricted(List& RSTr_obj) {
  List current_sample = RSTr_obj["current_sample"];
  cube beta = current_sample["beta"];
  cube lambda = current_sample["lambda"];
  cube Z = current_sample["Z"];
  mat tau2 = current_sample["tau2"];
  mat sig2 = current_sample["sig2"];
  List spatial_data = RSTr_obj["spatial_data"];
  field<uvec> island_region = spatial_data["island_region"];
  List params = RSTr_obj["params"];
  mat A = params["A"];
  double m0 = params["m0"];
  String method = params["method"];
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  uword num_island = island_region.n_elem;
  cube tmZ = lambda - Z;
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
  current_sample["beta"] = beta;
  RSTr_obj["current_sample"] = current_sample;
}

//[[Rcpp::export]]
void update_beta_mstcar(List& RSTr_obj) {
  List current_sample = RSTr_obj["current_sample"];
  cube beta = current_sample["beta"];
  cube lambda = current_sample["lambda"];
  cube Z = current_sample["Z"];
  vec tau2 = current_sample["tau2"];
  List spatial_data = RSTr_obj["spatial_data"];
  field<uvec> island_region = spatial_data["island_region"];
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  uword num_island = island_region.n_elem;
  for (uword isl = 0; isl < num_island; isl++) {
    uword num_island_region = island_region[isl].n_elem;
    mat var_beta = repmat(sqrt(tau2 / num_island_region), 1, num_time);
    mat mean_beta = mean(get_regs(lambda, island_region[isl]) - get_regs(Z, island_region[isl]), 0);
    if (num_time == 1) mean_beta = mean_beta.t();
    beta.row(isl) = mat(num_group, num_time, fill::randn) % var_beta + mean_beta;
  }
  current_sample["beta"] = beta;
  RSTr_obj["current_sample"] = current_sample;
}