#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "cpp_helpers.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
void update_rho(List& RSTr_obj) {
  List current_sample = RSTr_obj["current_sample"];
  vec rho = current_sample["rho"];
  cube G = current_sample["G"];
  cube Z = current_sample["Z"];
  List priors = RSTr_obj["priors"];
  double rho_a = priors["rho_a"];
  double rho_b = priors["rho_b"];
  vec rho_sd = priors["rho_sd"];
  vec rho_accept = priors["rho_accept"];
  List spatial_data = RSTr_obj["spatial_data"];
  field<uvec> adjacency = spatial_data["adjacency"];
  uword num_island = spatial_data["num_island"];
  uword num_region = Z.n_rows;
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  vec logit_rho = log(rho / (1 - rho));
  vec rand = Rcpp::rnorm(num_group, 0, 1);
  vec expit_rho = rand % rho_sd + logit_rho;
  vec rho_star_0 = 1 / (1 + exp(-expit_rho));
  vec r(num_group, fill::zeros);
  vec ra(num_group, fill::zeros);
  vec rb(num_group, fill::zeros);
  vec rc(num_group, fill::zeros);
  cube Zm(num_region, num_group, num_time);
  for (uword reg = 0; reg < num_region; reg++) {
    Zm.row(reg) = Z.row(reg) - mean(get_regs(Z, adjacency[reg]), 0);
  }
  for (uword grp = 0; grp < num_group; grp++) {
    vec rho_star = rho;
    rho_star[grp] = rho_star_0[grp];
    ra[grp] = (1 - pow(rho[grp], 2)) / (1 - pow(rho_star[grp], 2));
    field<mat> Sein_rho = Sig_eta_i(G, rho);
    field<mat> Sein_rho_star = Sig_eta_i(G, rho_star);
    field<mat> Sein_diff(num_time, num_time);    
    for (uword time1 = 0; time1 < num_time; time1++) {
      uword time1_l = (time1 == 0) ? 0 : time1 - 1;
      uword time2_u = (time1 == num_time - 1) ? num_time : time1 + 2;
      for (uword time2 = time1_l; time2 < time2_u; time2++) {
        Sein_diff(time2, time1) = Sein_rho_star(time2, time1) - Sein_rho(time2, time1);
      }
    }
    for (uword reg = 0; reg < num_region; reg++) {
      for (uword time1 = 0; time1 < num_time; time1++) {
        uword time2_l = (time1 == 0) ? 0 : time1 - 1;
        uword time2_u = (time1 == num_time - 1) ? num_time : time1 + 2;
        for (uword time2 = time2_l; time2 < time2_u; time2++) {
          mat rb_ik = get_grp(Z, reg, time2).t() * Sein_diff(time2, time1) * get_grp(Zm, reg, time1);
          rb[grp] += adjacency[reg].n_elem * rb_ik[0] / 2;
        }
      }
    }
    rc[grp] = pow(rho_star[grp] / rho[grp], rho_a) * pow((1 - rho_star[grp]) / (1 - rho[grp]), rho_b);
    r[grp] = exp((num_region - num_island) * (num_time - 1) / 2 * log(ra[grp]) - rb[grp]) * rc[grp];
    if (r[grp] >= R::runif(0, 1)) {
      rho[grp] = rho_star[grp];
      rho_accept[grp]++;
    }
  }
  priors["rho_accept"] = rho_accept;
  RSTr_obj["priors"] = priors;
  current_sample["rho"] = rho;
  RSTr_obj["current_sample"] = current_sample;
}