#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "cpp_helpers.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
void update_Ag(List& RSTr_obj) {
  List current_sample = RSTr_obj["current_sample"];
  mat Ag = current_sample["Ag"];
  cube G = current_sample["G"];
  List priors = RSTr_obj["priors"];
  mat Ag_scale = priors["Ag_scale"];
  double G_df = priors["G_df"];
  double Ag_df = priors["Ag_df"];
  uword num_time  = G.n_slices;
  uword num_group = G.n_rows;
  mat Ag_covar(num_group, num_group, fill::zeros);
  Ag_covar += inv(Ag_scale);
  for (uword time = 0; time < num_time; time++) {
    Ag_covar += inv(G.slice(time));
  }
  Ag = rwish(num_time * G_df + Ag_df, inv(Ag_covar));
  current_sample["Ag"] = Ag;
  RSTr_obj["current_sample"] = current_sample;
}