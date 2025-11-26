#include <RcppArmadillo.h>
#include <RcppDist.h>
#include "cpp_helpers.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
void update_G_default(List& RSTr_obj) {
  List current_sample = RSTr_obj["current_sample"];
  cube G = current_sample["G"];
  cube Z = current_sample["Z"];
  List priors = RSTr_obj["priors"];
  double G_df = priors["G_df"];
  mat G_scale = priors["G_scale"];
  List spatial_data = RSTr_obj["spatial_data"];
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
  current_sample["G"] = G;
  RSTr_obj["current_sample"] = current_sample;
}

//[[Rcpp::export]]
void update_G_mstcar(List& RSTr_obj) {
  List current_sample = RSTr_obj["current_sample"];
  cube G = current_sample["G"];
  cube Z = current_sample["Z"];
  mat Ag = current_sample["Ag"];
  vec rho = current_sample["rho"];
  List priors = RSTr_obj["priors"];
  double G_df = priors["G_df"];
  List spatial_data = RSTr_obj["spatial_data"];
  field<uvec> adjacency = spatial_data["adjacency"];
  uword num_island = spatial_data["num_island"];
  uword num_region = Z.n_rows;
  uword num_group = Z.n_cols;
  uword num_time = Z.n_slices;
  cube Ags(num_group, num_group, num_time, fill::zeros);
  Ags.each_slice() += Ag;
  vec r  = rho;
  vec sr = sqrt(1 - pow(rho, 2));

  for (uword reg = 0; reg < num_region; reg++) {
    double num_adj = adjacency[reg].n_elem;
    mat Zmikt = Z.row(reg) - mean(get_regs(Z, adjacency[reg]), 0);
    if (num_time == 1) Zmikt = Zmikt.t();
    mat Zt = get_grp(Z, reg, 0).t();
    Ags.slice(0) += num_adj * Zmikt.col(0) * Zt;
    for (uword time = 1; time < num_time; time++) {
      vec Zt  = 1 / sr % get_grp(Z, reg, time);
      vec Ztl = r / sr % get_grp(Z, reg, time - 1);
      vec Zm  = 1 / sr % Zmikt.col(time);
      vec Zml = r / sr % Zmikt.col(time - 1);
      Ags.slice(time) += num_adj * ((Zm - Zml) * (Zt - Ztl).t());
    }
  }
  for (uword time = 0; time < num_time; time++) {
    G.slice(time) = riwish((num_region - num_island) + G_df, Ags.slice(time));
  }
  current_sample["G"] = G;
  RSTr_obj["current_sample"] = current_sample;
}