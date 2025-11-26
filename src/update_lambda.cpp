#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "cpp_helpers.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
void update_lambda(List& RSTr_obj) {
  List current_sample = RSTr_obj["current_sample"];
  cube lambda = current_sample["lambda"];
  cube Z = current_sample["Z"];
  cube beta = current_sample["beta"];
  mat tau2 = current_sample["tau2"];
  List data = RSTr_obj["data"];
  cube Y = data["Y"];
  cube n = data["n"];
  List priors = RSTr_obj["priors"];
  cube lambda_sd = priors["lambda_sd"];
  cube lambda_accept = priors["lambda_accept"];
  List spatial_data = RSTr_obj["spatial_data"];
  uvec island_id = spatial_data["island_id"];
  List params = RSTr_obj["params"];
  String method = params["method"];
  uword num_region = Z.n_rows;
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  for (uword time = 0; time < num_time; time++) {
    for (uword reg = 0; reg < num_region; reg++) {
      for (uword grp = 0; grp < num_group; grp++) {
        double lambda_star = R::rnorm(lambda(reg, grp, time), lambda_sd(reg, grp, time));
        double rk1 = Y(reg, grp, time) * (lambda_star - lambda(reg, grp, time));
        double rk2 = 0;
        if (method == "binomial") {
          rk2 = n(reg, grp, time) * (log(1 + exp(lambda_star)) - log(1 + exp(lambda(reg, grp, time))));
        } 
        if (method == "poisson") {
          rk2 = n(reg, grp, time) * (exp(lambda_star) - exp(lambda(reg, grp, time)));
        }
        double rk3a = pow(lambda_star            - beta(island_id(reg), grp, time) - Z(reg, grp, time), 2);
        double rk3b = pow(lambda(reg, grp, time) - beta(island_id(reg), grp, time) - Z(reg, grp, time), 2);
        double rk = exp(rk1 - rk2 - 1 / (2 * tau2[grp]) * (rk3a - rk3b));
        if (rk >= R::runif(0, 1)) {
          lambda(reg, grp, time) = lambda_star;
          lambda_accept(reg, grp, time) ++;
        }
      }
    }
  }
  priors["lambda_accept"] = lambda_accept;
  RSTr_obj["priors"] = priors;
  current_sample["lambda"] = lambda;
  RSTr_obj["current_sample"] = current_sample;
}